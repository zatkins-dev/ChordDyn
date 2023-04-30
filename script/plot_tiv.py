from TIVlib import TIV
from tonal_interval_space import *


def my_cos(tiv1, tiv2):
    return np.arccos(
        np.clip(
            np.real(np.vdot(tiv1.vector, tiv2.vector))
            / (np.linalg.norm(tiv1.vector) * np.linalg.norm(tiv2.vector)),
            -1,
            1,
        )
    )


if __name__ == "__main__":
    T = tiv([0, 4, 7])  # C major
    Tm = tiv([0, 4, 9])  # a minor
    G = tiv([7, 11, 2])  # G
    F = tiv([5, 9, 0])  # F
    Tkey = tiv(Key("C"))
    print(TIV.cosine(G - Tkey, F - Tkey) / np.pi * 180)
