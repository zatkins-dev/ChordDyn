import music21 as m21
from playback import read_progression


def guess_key(chord_stream):
    return chord_stream.analyze("key")


def analyze_keys(full_progression, window):
    for i in range(len(full_progression) - window):
        s = m21.stream.Stream(full_progression[i : (i + window)])
        k = guess_key(s)
        common_names = ", ".join(c.pitchedCommonName for c in s)
        rns = " ".join(str(m21.roman.romanNumeralFromChord(c, k).figure) for c in s)
        print(common_names)
        print(f"key {k}: {rns}")


if __name__ == "__main__":
    s = m21.stream.Stream()
    for c in read_progression("chords.txt", "durations.txt"):
        s.append(c)
    analyze_keys(s, 4)
