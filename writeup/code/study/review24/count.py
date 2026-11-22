"""Independent recount: last-digit transitions (p_n, p_{n+k}), k=1..3, p_{n+k}-p_n <= DMAX, p_n in [L0,L1),
p_{n+k} < L1 (as r85 does). Also per-block counts (block by p_n) for batch-means standard errors.
usage: count.py L0 L1 DMAX NB out.npz"""
import sys, math, time
import numpy as np
L0 = int(float(sys.argv[1])); L1 = int(float(sys.argv[2])); DMAX = int(sys.argv[3]); NB = int(sys.argv[4]); OUT = sys.argv[5]
t0 = time.time()
R = math.isqrt(L1) + 1
s = np.ones(R + 1, bool); s[:2] = False
for i in range(2, math.isqrt(R) + 1):
    if s[i]: s[i*i::i] = False
base = np.nonzero(s)[0][1:].astype(np.int64)   # odd primes
del s
dig = np.full(10, -1); dig[[1, 3, 7, 9]] = [0, 1, 2, 3]
cnt = np.zeros((3, NB, 16), np.int64)
SEG = 10**8
tail = np.zeros(0, np.int64)
lo = L0 | 1
while lo < L1:
    hi = min(lo + SEG, L1)
    n = (hi - lo + 1) // 2
    m = np.ones(n, bool)                          # m[i] <-> lo + 2i
    first = np.maximum(base * base, -(-lo // base) * base)
    first = np.where(first % 2 == 0, first + base, first)
    idx = (first - lo) // 2
    for p, i in zip(base.tolist(), idx.tolist()):
        if i < n: m[i::p] = False
    P = lo + 2 * np.nonzero(m)[0].astype(np.int64)
    P = P[P < L1]
    allp = np.concatenate([tail, P])
    nt = len(tail)
    for k in (1, 2, 3):
        # pairs whose second element is new (index >= nt)
        j = np.arange(max(nt, k), len(allp))
        a = allp[j - k]; b = allp[j]
        ok = (b - a) <= DMAX
        a = a[ok]; b = b[ok]
        blk = ((a - L0) * NB) // (L1 - L0)
        code = dig[a % 10] * 4 + dig[b % 10]
        np.add.at(cnt[k-1], (blk, code), 1)
    tail = allp[-3:]
    lo = hi if hi % 2 == 1 else hi + 1
    print(lo, round(time.time() - t0), flush=True)
np.savez(OUT, cnt=cnt, L0=L0, L1=L1, DMAX=DMAX)
print('done', round(time.time() - t0), 's')
