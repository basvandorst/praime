"""Frozen predictions for prime triples (p, p+a, p+b): the distribution of the number of
even clocks (0 to 3) for each clock, from the same model as predict.py (classes mod 5 * 2^K spread evenly, coins independent
given the classes).  Written and committed before count3.py was run.
"""
from fractions import Fraction as F
from predict import single

PATTERNS = [(2, 6), (4, 6)]


def triple(a, b, K):
    M = 1 << K
    off = (0, a, b)
    cls5 = [c for c in range(1, 5) if all((c + o) % 5 for o in off)]
    acc = {"%s%d" % (f, j): F(0) for f in "HL" for j in range(4)}
    n = 0
    for r in range(1, M, 2):
        for c in cls5:
            hs, ls = zip(*(single((r + o) % M, (c + o) % 5, K) for o in off))
            for f, ps in (("H", hs), ("L", ls)):
                dist = [F(1)]
                for q in ps:
                    dist = [(dist[j] if j < len(dist) else 0) * (1 - q) + (dist[j - 1] * q if j else 0)
                            for j in range(len(dist) + 1)]
                for j in range(4):
                    acc["%s%d" % (f, j)] += dist[j]
            n += 1
    return {k: v / n for k, v in acc.items()}


if __name__ == "__main__":
    for a, b in PATTERNS:
        x, y = triple(a, b, 10), triple(a, b, 14)
        assert x == y
        print("(p, p+%d, p+%d): " % (a, b) + ", ".join("%s = %s = %.6f" % (k, v, float(v)) for k, v in x.items()))
