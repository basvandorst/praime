"""Insertion-response counts per class x mod 210p (x coprime to 210) for single primes p, to measure how the
explained response variance V grows as residue information mod p is added (reviewer's own code).
usage: mywin2.py OUTDIR X0 [X0 ...]"""
import sys, os, math
import numpy as np
def small_primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.flatnonzero(s)

def sieve(lo, hi, base):
    """primes in [lo, hi), odd-only segmented sieve"""
    out = []
    SEG = 2 * 10 ** 7
    for a in range(lo, hi, SEG):
        b = min(hi, a + SEG)
        a0 = a | 1                                   # first odd >= a
        n = (b - a0 + 1) // 2                        # odd numbers a0, a0+2, ..., < b
        if n <= 0: continue
        seg = np.ones(n, bool)
        for p in base[1:].tolist():                  # skip 2
            if p * p >= b: break
            # first odd multiple of p that is >= max(a0, p*p)
            st = max(p * p, ((a0 + p - 1) // p) * p)
            if st % 2 == 0: st += p
            seg[(st - a0) // 2::p] = False
        out.append(a0 + 2 * np.flatnonzero(seg))
    P = np.concatenate(out)
    return P[P >= max(lo, 3)]

def red(q): return [a for a in range(q) if math.gcd(a, q) == 1]

OUT = sys.argv[1]; WIN = [int(float(a)) for a in sys.argv[2:]]
SPAN = 10 ** 8; PAD = 4000
PL = [11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 151, 199, 251, 307, 401, 503]
base = small_primes(int(math.isqrt(max(WIN) + SPAN + PAD)) + 2)
ok210 = np.ones(210, bool)
for p in (2, 3, 5, 7): ok210[::p] = False
QS = (3, 5, 7)
for x0 in WIN:
    P = sieve(x0 - PAD, x0 + SPAN + PAD, base)
    res = {}
    for q in QS:
        C = len(red(q))
        for p in PL: res['c_%d_%d' % (q, p)] = np.zeros(210 * p * C * C, np.int64)
        res['c_%d_1' % q] = np.zeros(210 * C * C, np.int64)
    CH = 10 ** 7
    for a in range(x0, x0 + SPAN, CH):
        xs = a + np.arange(CH, dtype=np.int64); xs = xs[ok210[xs % 210]]
        j = np.searchsorted(P, xs); isp = P[j] == xs
        prevp = P[j - 1]; nextp = np.where(isp, P[j + 1], P[j])
        for q in QS:
            r = red(q); C = len(r); lut = -np.ones(q, np.int64); lut[r] = np.arange(C)
            lr = lut[prevp % q] * C + lut[nextp % q]
            res['c_%d_1' % q] += np.bincount((xs % 210) * C * C + lr, minlength=210 * C * C)
            for p in PL:
                M = 210 * p
                res['c_%d_%d' % (q, p)] += np.bincount((xs % M) * C * C + lr, minlength=M * C * C)
    np.savez_compressed(os.path.join(OUT, 'v_%d.npz' % x0), **res)
    print('window %d done' % x0, flush=True)
