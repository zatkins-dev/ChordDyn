import music21 as m21
from music21.chord import Chord
from music21.key import Key
from music21.stream import Stream
from music21.note import Note, Pitch
from itertools import permutations
from playback import read_progression
import numpy as np
from collections.abc import Iterable
from tonal_interval_space import *

np.set_printoptions(precision=4, linewidth=90)


HARMONIC_FUNCTION_D_MAX = 1.6863383443715223


def complex_array_str(x):
    return (
        "["
        + " ".join(
            f"{np.real(xi):0.3g}"
            + (f"+{np.imag(xi):0.3g}j" if np.abs(np.imag(xi)) > 1e-10 else "")
            for xi in x
        )
        + "]"
    )


def guess_key(chord_stream):
    s1 = Stream()
    for chord in chord_stream:
        c = Chord(chord.normalOrder).simplifyEnharmonics().closedPosition(forceOctave=4)
        c.inversion(0)
        s1.append(c)
    return s1.analyze("key")


def dist_mod_12(p, q):
    match (p, q):
        case (Iterable(), Iterable()):
            p = np.asarray(p)
            q = np.asarray(q)
            return np.min([np.mod(p - q, 12), np.mod(q - p, 12)], axis=0)
        case _:
            return np.min([np.mod(p - q, 12), np.mod(q - p, 12)])


def get_voices(c1: Chord, c2: Chord):
    pc1 = c1.orderedPitchClasses
    pc2 = c2.orderedPitchClasses
    num_pitches = len(pc1)
    choices = [p for p in permutations(range(num_pitches), num_pitches)]
    semitone_dist = [
        [dist_mod_12(pc2[i], pc1[j]) for (i, j) in zip(range(num_pitches), p)]
        for p in choices
    ]
    total_dist = [sum(st_dist) for st_dist in semitone_dist]
    min_ind = np.argmin(total_dist)
    min_permutation = choices[min_ind]
    nl1 = np.asarray([pc1[j] for j in min_permutation])
    nl2 = np.asarray(pc2)
    return nl1, nl2


def perceptual_unrelatedness(T, Tm1):
    if T == Tm1:
        return 0.0
    return (TIV.euclidean(T, Tm1) - D_MIN_TRIAD_TRIAD) / (
        D_MAX_TRIAD_TRIAD - D_MIN_TRIAD_TRIAD
    )


def dissonance(T):
    return T.dissonance()


def key_unrelatedness(T, k):
    # (1 - np.cos(TIV.cosine(T, tiv(k)))) / 2
    return scaled_cosine_dist_triad_key(T, tiv(k))


def harmonic_function_unrelatedness(T: TIV, k):
    # return (1 - np.cos(TIV.cosine(T, Tf(T, k)))) / 2
    Tkey = tiv(k)
    return (
        1 - np.cos(TIV.cosine(T - Tkey, Tf(T, k) - Tkey))
    ) / 2  # / HARMONIC_FUNCTION_D_MAX


def anti_voice_leading(c: Chord, cm1: Chord, k: Key):
    v1, v2 = get_voices(cm1, c)
    s = dist_mod_12(v1, v2)
    num = 0.0
    if sum(s) == 0:
        return 0.0
    for pc, dist in zip(v2, s):
        cos_dist = scaled_cosine_dist_note_key(tiv([pc]), tiv(k))
        num += cos_dist * np.exp(0.05 * (dist - 5)) / 3
    return num * perceptual_unrelatedness(tiv(c), tiv(cm1))


def tonal_fitness(ci, cim1, k):
    Ti = tiv(ci)
    Tim1 = tiv(cim1) if cim1 else None
    c = dissonance(Ti)
    m = anti_voice_leading(ci, cim1, k)
    h = 0
    return np.array([4.22, 2.14, 2.06, 3.78]).dot(np.array([delta, c, m, h]))


def get_time_signature(window):
    match window:
        case 3:
            return m21.meter.TimeSignature("3/4")
        case 6:
            return m21.meter.TimeSignature("slow 6/8")
        case _:
            return m21.meter.TimeSignature("4/4")


def get_whole_rest(window, ts):
    return m21.note.Rest(window * ts.beatLengthToQuarterLengthRatio, fullMeasure="auto")


def pitch_in_key(pitch: Pitch, key: Key):
    degree, accidental = key.getScaleDegreeAndAccidentalFromPitch(
        Pitch(pitch.pitchClass),
        comparisonAttribute="pitchClass",
    )

    # if accidental is not None and np.abs(accidental.alter) > 1:
    #     degree, accidental = key.getScaleDegreeAndAccidentalFromPitch(
    #         Pitch(pitch.pitchClass + np.sign(accidental.alter) * 1),
    #         comparisonAttribute="pitchClass",
    #     )
    #     return key.pitchFromDegree(degree)
    new_pitch = key.pitchFromDegree(degree)
    new_pitch.octave = pitch.octave
    if accidental is None:
        return new_pitch
    if accidental.alter < 0:
        if new_pitch.accidental is None:
            if key.sharps > 0:  # sharps
                print(f"changing flat to sharp")
                new_pitch.transpose(-2, inPlace=1)
                accidental.alter += 2
            new_pitch.accidental = accidental
        elif new_pitch.accidental.alter < 0:
            new_pitch.transpose(-1, inPlace=1)
        else:
            new_pitch.accidental = None
        return new_pitch
    if accidental.alter > 0:
        if new_pitch.accidental is None:
            if key.sharps <= 0:  # flats
                print(f"changing sharp to flat")
                new_pitch.transpose(2, inPlace=1)
                accidental.alter -= 2
            new_pitch.accidental = accidental
        elif new_pitch.accidental.alter > 0:
            new_pitch.transpose(1, inPlace=1)
        else:
            new_pitch.accidental = None
        return new_pitch


def annotate(progression, window, closedPosition=False):
    aprog = Stream()
    part = m21.stream.Part()
    aprog.append(part)

    ts = get_time_signature(window)
    part.timeSignature = ts
    if closedPosition:
        part.append(m21.clef.TrebleClef())
    else:
        part.append(m21.clef.BassClef())
    # part.append(ts)
    one_beat = m21.duration.Duration(1 * ts.beatLengthToQuarterLengthRatio)
    for i in range(0, len(progression), window):
        s = Stream(progression[i : i + window])
        k: Key = guess_key(s)
        print(k)
        part.append(k)

        s1 = Stream()
        for c in s:
            new_c = Chord([pitch_in_key(n.pitch, k) for n in c.notes])
            print(f"{c}, {new_c}")
            assert set(c.pitchClasses).intersection(new_c.pitchClasses) == set(
                c.pitchClasses
            )
            if closedPosition:
                s1.append(new_c.closedPosition(forceOctave=4))
            else:
                s1.append(new_c)
        s = s1

        for chord in s:
            rn = m21.roman.romanNumeralFromChord(chord, k)
            rn.simplifyEnharmonics(inPlace=True)
            # rn.lyric = rn.figure
            # rn.duration = one_beat
            # part.append(rn)
            c = Chord(chord.pitches, duration=one_beat)
            c.lyric = rn.romanNumeral
            print(f"{c} : {rn.figureAndKey}")
            part.append(c)
        if i < len(progression) - window:
            part.append(get_whole_rest(window, ts))
    aprog.makeMeasures(inPlace=True)
    aprog.show("text")

    return aprog


def analyze_progressions(full_progression, window):
    P_tiv = [tiv(c) for c in full_progression]
    for i in range(len(full_progression) - window):
        s = Stream(full_progression[i : (i + window)])
        k = guess_key(s)
        s1 = Stream()
        for c in s:
            new_c = Chord([pitch_in_key(n.pitch, k) for n in c.notes])
            s1.append(new_c.closedPosition(forceOctave=4))
        s = s1
        rns = [m21.roman.romanNumeralFromChord(c, k).romanNumeral for c in s]
        print(f"key {k}: {' '.join(rns)}")
        print(f"   Tkey: {complex_array_str(tiv(k).vector)}")
        for j in range(0, window):
            Tj = P_tiv[i + j]
            Tjm1 = P_tiv[i + (j - 1) % window]
            delta = perceptual_unrelatedness(Tj, Tjm1)
            xi = dissonance(Tj)
            ell = key_unrelatedness(Tj, k)
            phi = harmonic_function_unrelatedness(Tj, k)
            vl = anti_voice_leading(
                full_progression[i + j], full_progression[i + (j - 1) % window], k
            )
            print(f"  {rns[j]} {s[j].orderedPitchClasses} ({s[j].pitchedCommonName})")
            print(f"    T(c) = {complex_array_str(Tj.vector)}")
            if j > 0:
                print(f"    perceptual unrelatedness        = {delta}")
                print(f"    melodic unattraction            = {vl}")
            print(f"    dissonance                      = {xi}")
            print(f"    key unrelatedness               = {ell}")
            print(f"    harmonic function unrelatedness = {phi}")

        print(f"  Progression:")

        print()


if __name__ == "__main__":
    s = Stream()
    for c in read_progression(
        "output/chords.txt", "output/durations.txt", max_chords=15
    ):
        s.append(c)
    analyze_progressions(s, 3)
    # print(f"max dist: {max_distance()}")
    # print(f"max dist note from major: {max_distance_from_key(1)}")
    # print(f"max dist note from minor: {max_distance_from_key(1, mode='minor')}")

    # print(f"max dist triad from major: {max_distance_from_key(3)}")
    # print(f"max dist triad from minor: {max_distance_from_key(3, mode='minor')}")
    # print(f"min dist note from major: {min_distance_from_key(1)}")
    # print(f"min dist note from minor: {min_distance_from_key(1, mode='minor')}")

    # print(f"min dist triad from major: {min_distance_from_key(3)}")
    # print(f"min dist triad from minor: {min_distance_from_key(3, mode='minor')}")
    # annotated = annotate(s, 3, closedPosition=True)
    # annotated.show()
