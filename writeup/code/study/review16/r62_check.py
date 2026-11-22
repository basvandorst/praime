import math, itertools, sys, json
import numpy as np
from sieve import base_primes, primes_in
L0, L1, SEG = 10 ** 10, 11 * 10 ** 9, 10 ** 8
KMAX = 4
def Hb(p):
    p = np.asarray(p, float).ravel(); p = p[p > 0]; p = p / p.sum(); return float(-(p * np.log2(p)).sum())
def mi(J):
    return Hb(J.sum(axis=-1)) + Hb(J.reshape(-1, 15).sum(axis=0)) - Hb(J)
bp = base_primes(int(math.isqrt(L1)) + 2)
cnt = {k: np.zeros(15 ** (k + 1), np.int64) for k in range(1, KMAX + 1)}
tail = None; ngaps = 0
for s in range(L0, L1, SEG):
    P = primes_in(s, s + SEG, bp)
    if tail is not None: P = np.concatenate([tail, P])
    c = ((np.diff(P) // 2) % 15).astype(np.int64)
    # windows: only those whose last gap is new in this segment
    new0 = 0 if tail is None else len(tail) - 1   # index of first new gap in c
    ngaps += len(c) - new0
    for k in range(1, KMAX + 1):
        ends = np.arange(max(k, new0), len(c))
        code = np.zeros(len(ends), np.int64)
        for j in range(k, -1, -1):
            code = code * 15 + c[ends - j]
        cnt[k] += np.bincount(code, minlength=15 ** (k + 1))
    tail = P[-(KMAX + 1):]
    print(s, len(P), flush=True)
lam = float(np.mean(np.log(np.linspace(L0, L1, 1001))))
out = {'lambda': lam, 'gaps': ngaps}
for k in range(1, KMAX + 1):
    J = cnt[k].reshape((15,) * (k + 1)).astype(float); N = J.sum()
    I = mi(J)
    nz = np.count_nonzero(J); nzp = np.count_nonzero(J.sum(-1)); nzn = np.count_nonzero(J.reshape(-1, 15).sum(0))
    mm = ((nzp - 1) + (nzn - 1) - (nz - 1)) / (2 * N * math.log(2))
    out[k] = {'N': int(N), 'plugin': I, 'MM': I + mm, 'mm_corr': mm, 'cells': int(nz)}
    print(k, out[k], flush=True)
np.savez_compressed('r62_counts.npz', **{'k%d' % k: cnt[k] for k in cnt})
json.dump(out, open('r62_check.json', 'w'), indent=1)
