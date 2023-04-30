import music21 as m21

from music21.chord import Chord
from music21.note import Note
from music21.key import Key
from music21.stream import Stream

import numpy as np
from pathlib import Path


def read_progression(
    cname: Path,
    dname: Path,
    scale_duration=False,
    max_chords=None,
    force_octave=None,
):
    durations = np.loadtxt(dname)
    if scale_duration:
        durations = [1 for _ in durations]
    chords = np.loadtxt(cname, dtype=int, delimiter=" ", max_rows=max_chords)
    for c, d in zip(chords, durations):
        if force_octave is not None:
            for i in range(c.shape[0]):
                c[i] = np.mod(c[i], 12)
                c[i] += (force_octave - (c[i] > 9)) * 12
        if scale_duration:
            try:
                dur_type = m21.duration.quarterLengthToClosestType(d)
                m21.musicxml.m21ToXml.typeToMusicXMLType(dur_type[0])
                duration = m21.duration.Duration(type=dur_type[0])
            except Exception:
                # skip if duration is too short
                continue
        else:
            duration = m21.duration.Duration(d)
        try:
            chord = Chord(c.tolist(), duration=duration)
            if force_octave is not None:
                chord.inversion(0)
            yield chord
        except:
            yield Note(c, duration=duration)


if __name__ == "__main__":
    progression = Stream()
    # for n in read_voice("chords.txt", "durations.txt", False):
    #     progression.append(n)
    # progression.show()
    # output = Stream()
    for c in read_progression(
        "output/chords.txt", "output/durations.txt", True, force_octave=5
    ):
        symbol = m21.harmony.chordSymbolFromChord(c)
        progression.append(symbol)
        progression.append(c)

    progression = m21.harmony.realizeChordSymbolDurations(progression)
    # print(progression)
    progression.show()
