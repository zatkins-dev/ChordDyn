import music21 as m21
import numpy as np
from pathlib import Path


def read_progression(cname: Path, dname: Path, scale_duration=False):
    durations = np.loadtxt(dname)
    if scale_duration:
        durations = [
            max(d, 1 / 8) for d in 4 * durations / np.max(durations)
        ]  # make min length 64th note
    print(durations)
    chords = np.loadtxt(cname, dtype=int, delimiter=" ")
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


def diff(c1, c2):
    return list(set(c1).difference(set(c2)))


def shared(c1, c2):
    return list(set(c1).intersection(set(c2)))


def read_voice(cname: Path, dname: Path, ret_changes: bool = True):
    durations = np.loadtxt(dname)
    durations = 4 * durations / np.max(durations)  # make max length 4
    chords = np.loadtxt(cname, dtype=int, delimiter=" ")
    changed_notes_fwd_idx = [
        c1.tolist().index(diff(c1, c2)[0])
        for c1, c2 in zip(chords[1:, :], chords[:-1, :])
    ]
    changed_notes_bwd_idx = [
        c2.tolist().index(diff(c2, c1)[0])
        for c1, c2 in zip(chords[1:, :], chords[:-1, :])
    ]
    print(changed_notes_fwd_idx)
    changed_notes_fwd = [c[i] for c, i in zip(chords[:-1, :], changed_notes_fwd_idx)]
    changed_notes_bwd = [c[i] for c, i in zip(chords[1:, :], changed_notes_bwd_idx)]
    same_notes_fwd = [
        c[[j != i for j in range(len(c))]]
        for c, i in zip(chords[:-1, :], changed_notes_fwd_idx)
    ]
    same_notes_bwd = [
        c[[j != i for j in range(len(c))]]
        for c, i in zip(chords[1:, :], changed_notes_bwd_idx)
    ]
    if ret_changes:
        for n1, n2, d in zip(
            changed_notes_fwd[1:], changed_notes_bwd[:-1], durations[1:-1]
        ):
            try:
                dur_type = m21.duration.quarterLengthToClosestType(d)
                m21.musicxml.m21ToXml.typeToMusicXMLType(dur_type[0])
            except Exception:
                # skip if duration is too short
                continue
            print(n1, n2, dur_type[0])
            yield m21.chord.Chord(
                [int(n1), int(n2)],
                duration=m21.duration.Duration(type=dur_type[0]),
            )
    else:
        for cf, cb, d in zip(same_notes_fwd[1:], same_notes_bwd[:-1], durations[1:-1]):
            try:
                dur_type = m21.duration.quarterLengthToClosestType(d)
                m21.musicxml.m21ToXml.typeToMusicXMLType(dur_type[0])
            except Exception:
                # skip if duration is too short
                continue
            yield m21.chord.Chord(
                cf.tolist(),
                duration=m21.duration.Duration(type=dur_type[0]),
            )


if __name__ == "__main__":
    progression = m21.stream.Stream()
    # for n in read_voice("chords.txt", "durations.txt", False):
    #     progression.append(n)
    # progression.show()
    for c in read_progression("chords.txt", "durations.txt", False):
        progression.append(c)
    # print(progression)
    progression.show()
