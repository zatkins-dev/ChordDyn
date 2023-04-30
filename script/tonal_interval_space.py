from music21.chord import Chord
from music21.key import Key
from music21.note import Note
from music21.roman import RomanNumeral
import numpy as np
from itertools import combinations
from collections.abc import Iterable
from TIVlib import TIV, TIVCollection

TIV_WEIGHTS = np.array([2, 11, 17, 16, 19, 7])
D_MAX_TRIAD_TRIAD = 36.427610264768354
D_MIN_TRIAD_TRIAD = 12.566546832245992
D_MAX_NOTE_KEY = 39.62532649435763
D_MIN_NOTE_KEY = 29.495270674859842
D_MAX_TRIAD_KEY = 29.35600628754387
D_MIN_TRIAD_KEY = 9.49705496319371
C_MAX = 32.863353450309965
COS_NOTE_KEY_MIN = 1.0929416422891878
COS_NOTE_KEY_MAX = 2.12247555540512
COS_TRIAD_KEY_MIN = 0.6428199666592707
COS_TRIAD_KEY_MAX = 2.653664640179884


def scaled_cosine_dist_note_key(Tnote, Tkey):
    # return (TIV.cosine(Tnote, Tkey) - COS_NOTE_KEY_MIN) / (
    #     COS_NOTE_KEY_MAX - COS_NOTE_KEY_MIN
    # )
    return (1 - np.cos(TIV.cosine(Tnote, Tkey))) / 2


def scaled_cosine_dist_triad_key(Ttriad, Tkey):
    return (TIV.cosine(Ttriad, Tkey) - COS_TRIAD_KEY_MIN) / (
        COS_TRIAD_KEY_MAX - COS_TRIAD_KEY_MIN
    )


def chroma_vector(value):
    v = np.zeros((12,))
    match value:
        case Chord() as c:
            v[c.orderedPitchClasses] = 1
        case Note() as n:
            v[n.pitch.pitchClass] = 1
        case Key() as k:
            v[np.array([p.pitchClass for p in k.pitches])] = 1
        case Iterable() as l:
            v[np.asarray(l, dtype=int)] = 1
    return v


def tiv(o):
    return TIV.from_chroma(chroma_vector(o))


def tiv_norm(T):
    return T.dissonance()


def triad_dist(T1, T2):
    return TIV.euclidean(T1, T2)


def note_key_dist(pc, k):
    Tnote = tiv(Note(pitchClass=pc))
    return (TIV.euclidean(Tnote, tiv(k)) - D_MIN_NOTE_KEY) / (
        D_MAX_NOTE_KEY - D_MIN_NOTE_KEY
    )


def triad_key_dist(T, k):
    return (TIV.euclidean(T, tiv(k)) - D_MIN_TRIAD_KEY) / (
        D_MAX_TRIAD_KEY - D_MIN_TRIAD_KEY
    )


def max_distance():
    all_tivs = [tiv(c) for c in combinations(range(12), 3)]
    return max([TIV.euclidean(tiv1, tiv2) for tiv1, tiv2 in combinations(all_tivs, 2)])


def max_distance_from_key(num_notes, mode="major"):
    if mode == "major":
        key = Key("C")
    elif mode == "minor":
        key = Key("a")
    else:
        raise ValueError("mode must be 'major' or 'minor'")
    chords = [c for c in combinations(range(12), num_notes)]
    dists = [TIV.euclidean(tiv(c), tiv(key)) for c in chords]
    print(chords[np.argmax(dists)])
    return np.max(dists)


def min_distance():
    all_tivs = [tiv(c) for c in combinations(range(12), 3)]
    return min([TIV.euclidean(tiv1, tiv2) for tiv1, tiv2 in combinations(all_tivs, 2)])


def min_distance_from_key(num_notes, mode="major"):
    if mode == "major":
        key = Key("C")
    elif mode == "minor":
        key = Key("a")
    else:
        raise ValueError("mode must be 'major' or 'minor'")
    chords = [c for c in combinations(range(12), num_notes)]
    dists = [TIV.euclidean(tiv(c), tiv(key)) for c in chords]
    print(chords[np.argmin(dists)])
    return np.min(dists)


def angle(v1, v2):
    inner = cos_angle(v1, v2)
    if np.abs(inner) - 1 < 1e-12:
        return 0 if inner > 0 else np.pi
    return np.arccos(inner)


def tonic(k: Key):
    if k.mode == "minor":
        return RomanNumeral("i", k)
    return RomanNumeral("I", k)


def dominant(k: Key):
    if k.mode == "minor":
        return RomanNumeral("v", k)
    return RomanNumeral("V", k)


def subdominant(k: Key):
    if k.mode == "minor":
        return RomanNumeral("iv", k)
    return RomanNumeral("IV", k)


def Tf(T, k, return_rn=False):
    Tkey = tiv(k)
    Tfs = [tiv(tonic(k)), tiv(dominant(k)), tiv(subdominant(k))]
    dist = [1 - np.cos(TIV.cosine(T - Tkey, tf - Tkey)) for tf in Tfs]
    # print(dist)
    ind = np.argmin(dist)
    labels = ["I", "V", "IV"] if k.mode == "major" else ["i", "v", "iv"]
    if return_rn:
        return Tfs[ind], labels[ind]
    else:
        return Tfs[ind]


def build_tis_tree(progression):
    pass
