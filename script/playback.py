import music21 as m21

from music21.chord import Chord
from music21.key import Key
from music21.stream import Stream

import numpy as np
from pathlib import Path


def read_progression(cname: Path, dname: Path, scale_duration=False, max_chords=None):
    durations = np.loadtxt(dname)
    if scale_duration:
        durations = [
            max(d, 1 / 8) for d in 4 * durations / np.max(durations)
        ]  # make min length 64th note
    chords = np.loadtxt(cname, dtype=int, delimiter=" ", max_rows=max_chords)
    for c, d in zip(chords, durations):
        try:
            dur_type = m21.duration.quarterLengthToClosestType(d)
            m21.musicxml.m21ToXml.typeToMusicXMLType(dur_type[0])
        except Exception:
            # skip if duration is too short
            continue
        try:
            clist = c.tolist()
            yield m21.chord.Chord(
                clist, duration=m21.duration.Duration(type=dur_type[0])
            )
        except:
            yield m21.note.Note(c, duration=m21.duration.Duration(type=dur_type[0]))


if __name__ == "__main__":
    progression = m21.stream.Stream()
    # for n in read_voice("chords.txt", "durations.txt", False):
    #     progression.append(n)
    # progression.show()
    for c in read_progression(
        "output/chords-section.txt", "output/durations-section.txt", False
    ):
        progression.append(c)
    # print(progression)
    progression.show()
