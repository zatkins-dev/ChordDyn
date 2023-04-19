from music21.chord import Chord
from music21.key import Key
from music21.note import Note
from music21.roman import RomanNumeral
import numpy as np

TIV_WEIGHTS = np.array([2, 11, 17, 16, 19, 7])
TIV_WEIGHTS_PLOT = np.array([0.049, 0.132, 0.203, 0.221, 0.246, 0.121])
C_MAX = 32.863353450309965


def chroma_vector(value):
    v = np.zeros((12,))
    match value:
        case Chord() as c:
            v[c.orderedPitchClasses] = 1
        case Note() as n:
            v[n.pitch.pitchClass] = 1
        case Key() as k:
            v[np.array([p.pitchClass for p in k.pitches])] = 1
        case list() as l:
            v[np.asarray(l)] = 1
    return v


def tiv(o):
    c = chroma_vector(o)
    T = np.fft.rfft(c)
    return T[1:7] / T[0] * TIV_WEIGHTS / C_MAX


def cos_angle(v1, v2):
    return np.real(np.vdot(v2, v1)) / (np.linalg.norm(v1) * np.linalg.norm(v2))


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
    Tfs = [tiv(tonic(k)) - Tkey, tiv(dominant(k)) - Tkey, tiv(subdominant(k)) - Tkey]
    T_ref = T - tiv(k)
    dist = [angle(T_ref, tf) for tf in Tfs]
    ind = np.argmin(dist)
    if return_rn:
        labels = ["I", "V", "IV"] if k.mode == "major" else ["i", "v", "iv"]
        return Tfs[ind], labels[ind]
    else:
        return Tfs[ind]


def build_tis_tree(progression):
    pass
