"""Check 7: n^2 + 1. (a) count even n in F1, F2, F3 with n^2+1 (probable) prime, one by one with gmpy2.is_prime.
(b) S_f({0}) by my own product to 1e7; omega_B(p) by brute force vs the rule
    omega = |B|(1 + chi4(p)) minus coincidences at p | (b - b')^2 + 4 (and p | b - b'), for p < 2000."""
import math, itertools
import gmpy2
import numpy as np
from myS import PRI


def count(lo, hi):
    n0 = lo + (lo % 2)
    return sum(1 for n in range(n0, hi, 2) if gmpy2.is_prime(n * n + 1))


def omega_brute(B, p):
    return sum(1 for x in range(p) if any((((x + b) ** 2 + 1) % p) == 0 for b in B))


def omega_rule(B, p):
    if p == 2:
        return len({(1 - b) % 2 for b in B})
    roots = 2 if p % 4 == 1 else 0
    if roots == 0:
        return 0
    # union of shifted root pairs {r - b}: subtract coincidences
    pts = set()
    i = int(gmpy2.isqrt(p))
    r = next(x for x in range(p) if (x * x + 1) % p == 0)
    for b in B:
        pts |= {(r - b) % p, (-r - b) % p}
    return len(pts)


def omega_text_rule(B, p):
    """the write-up's large-p rule: generic |B|(1+chi4(p)) minus one per pair with p | (b-b')^2 + 4."""
    if p % 4 == 3:
        return 0
    om = 2 * len(B)
    for b, bb in itertools.combinations(B, 2):
        if ((b - bb) ** 2 + 4) % p == 0:
            om -= 1
    return om


if __name__ == '__main__':
    Ps = PRI[PRI < 10**7]
    chi = np.where(Ps % 4 == 1, 1, np.where(Ps % 4 == 3, -1, 0))
    w = 1 + chi; w[0] = 1
    Sf = math.exp(float(np.sum(np.log1p(-w / Ps) - np.log1p(-1.0 / Ps))))
    print('(b) S_f({0}) product to 1e7: %.6f   (literature Bateman-Horn constant of n^2+1: 1.3728134628)' % Sf)
    bad = bad_text = 0; tested = 0
    for B in ([0, 2], [0, 4], [0, 2, 6], [0, 10], [0, 2, 4, 8, 12], [0, 6, 20]):
        span = max(B)
        for p in [int(q) for q in PRI[PRI < 2000]]:
            ob = omega_brute(B, p)
            bad += ob != omega_rule(B, p)
            if p > span and p > 2:
                tested += 1
                bad_text += ob != omega_text_rule(B, p)
    print('(b) omega brute vs union-of-shifted-roots: %d mismatches; vs text rule for p > span: %d of %d' % (bad, bad_text, tested))
    for name, lo, hi in (('F1', 10**6, 3 * 10**6), ('F2', 10**8, 10**8 + 10**6), ('F3', 10**10, 10**10 + 2 * 10**5)):
        print('(a) %s [%d, %d): even n with n^2+1 prime: %d' % (name, lo, hi, count(lo, hi)), flush=True)
