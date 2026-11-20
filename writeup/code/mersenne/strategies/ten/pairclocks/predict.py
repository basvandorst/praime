"""Predictions for prime pairs (p, p+d): how often do both members have an even clock?

Two clocks of a prime q > 5:
  H(q): the order of 2 mod q is even      <=> q divides some 2^n + 1   (Hasse: density 17/24)
  L(q): the rank z(q) of q in Fibonacci is even <=> q divides some Lucas number (Lagarias: 2/3)

Single-prime model (derivation in MERSENNE.md section 53), e = v2(q - 1):
  H: q = 3, 5 mod 8 -> 1;  q = 7 mod 8 -> 0;  q = 1 mod 8 -> 1 - 2^-(e-1)
  L: q = 3 mod 4 -> 1;  q = 1 mod 4 and q = 2, 3 mod 5 -> 0;  q = 1 mod 4 and q = 1, 4 mod 5 -> 1 - 2^(1-e)
Pair model (the conjecture): p is spread evenly over the classes mod 5 * 2^K that a prime pair can
occupy, and given the classes of p and p+d the two clocks are independent.

Exact rationals: enumerate p mod 2^K; the one class where q = 1 mod 2^K (e unknown) gets the exact
average over e >= K (geometric). K = 12 and K = 16 must give the same fraction.
"""
from fractions import Fraction as F
import sys


def single(q2, q5, K):
    """P(H), P(L) for a prime q with q mod 2^K = q2 (odd) and q mod 5 = q5 (nonzero)."""
    r = (q2 - 1) % (1 << K)
    if r == 0:
        deep = F(4, 3) / (1 << K)
        eH = 1 - deep
        eL = 1 - deep
        e = None
    else:
        e = (r & -r).bit_length() - 1
    m8 = q2 % 8
    if m8 in (3, 5):
        h = F(1)
    elif m8 == 7:
        h = F(0)
    else:
        h = eH if e is None else 1 - F(1, 1 << (e - 1))
    if q2 % 4 == 3:
        lu = F(1)
    elif q5 in (2, 3):
        lu = F(0)
    else:
        lu = eL if e is None else 1 - F(2, 1 << e)
    return h, lu


def pair(d, K):
    M = 1 << K
    cls5 = [c for c in range(1, 5) if (c + d) % 5 != 0]
    tot = 0
    acc = {k: F(0) for k in ("H2", "H1", "H0", "L2", "L1", "L0")}
    for a in range(1, M, 2):
        b = (a + d) % M
        for c in cls5:
            h1, l1 = single(a, c, K)
            h2, l2 = single(b, (c + d) % 5, K)
            tot += 1
            acc["H2"] += h1 * h2
            acc["H0"] += (1 - h1) * (1 - h2)
            acc["L2"] += l1 * l2
            acc["L0"] += (1 - l1) * (1 - l2)
    for f in "HL":
        acc[f + "2"] /= tot
        acc[f + "0"] /= tot
        acc[f + "1"] = 1 - acc[f + "2"] - acc[f + "0"]
    return acc


def singles(K):
    M = 1 << K
    h = lu = F(0)
    n = 0
    for a in range(1, M, 2):
        for c in range(1, 5):
            x, y = single(a, c, K)
            h += x
            lu += y
            n += 1
    return h / n, lu / n


if __name__ == "__main__":
    ds = list(range(2, 32, 2))
    h12, l12 = singles(12)
    h16, l16 = singles(16)
    assert h12 == h16 == F(17, 24) and l12 == l16 == F(2, 3), (h12, l12)
    print("single primes: H =", h12, " L =", l12, "(Hasse 17/24, Lagarias 2/3)")
    print("naive independent pair: H both = %s = %.6f, L both = %s = %.6f" % (
        h12 ** 2, float(h12 ** 2), l12 ** 2, float(l12 ** 2)))
    print("d  | H both      H one       H none     | L both      L one       L none")
    for d in ds:
        a = pair(d, 12)
        b = pair(d, 16)
        assert a == b, d
        print("%2d | %-10s %-10s %-10s | %-10s %-10s %-10s" % (
            d, a["H2"], a["H1"], a["H0"], a["L2"], a["L1"], a["L0"]))
        print("   | %.6f   %.6f   %.6f   | %.6f   %.6f   %.6f" % tuple(
            float(a[k]) for k in ("H2", "H1", "H0", "L2", "L1", "L0")))
    sys.stdout.flush()
