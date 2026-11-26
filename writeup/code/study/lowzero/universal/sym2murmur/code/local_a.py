"""alpha(l, e, v) = (1/phi(l^e)) sum_{b unit mod l^e} e(-2 l^v b / l^e) S(l^{2v} b^2, 1; l^e)   (complete sums, floats)."""
import numpy as np
from math import gcd

def kloost(m, c):
    x = np.array([x for x in range(c) if gcd(x, c) == 1])
    xi = np.array([pow(int(v), -1, c) for v in x])
    return np.exp(2j * np.pi * ((m * x + xi) % c) / c).sum()

def alpha(l, e, v):
    c = l ** e
    bs = [b for b in range(c) if gcd(b, c) == 1]
    tot = 0
    for b in bs:
        tot += np.exp(-2j * np.pi * ((2 * l ** v * b) % c) / c) * kloost((l ** (2 * v) * b * b) % c, c)
    return tot / len(bs)

if __name__ == "__main__":
    for l in (2, 3, 5, 7, 11):
        for e in (1, 2, 3):
            if l ** e > 400: continue
            row = []
            for v in range(0, 4):
                a = alpha(l, e, v)
                row.append("%+.6f%+.6fi" % (a.real, a.imag))
            print("l=%d e=%d: " % (l, e) + "  ".join(row), " [1/phi=%.6f]" % (1 / (l ** e - l ** (e - 1))))
