"""Independent check of charlog.py (sympy factorisation, mpmath; no FFT, no L-values).
For sample n, the class sum C(n) = sum_{m = n (mod q)} (Lambda(m)/ln m) m^-s is computed from its definition
(the first J members factorised by sympy, the rest bounded), times n^s, and compared with D(n) from the run.
Usage: charlog_check.py charlog_X.npy q s J"""
import sys, math
import numpy as np
import mpmath as mp
from sympy import factorint
mp.mp.dps = 30
A = np.load(sys.argv[1]); q = int(sys.argv[2]); s = mp.mpf(sys.argv[3]); J = int(sys.argv[4])
ns, D = A[0].astype(np.int64), A[1]
rng = np.random.default_rng(3)
pick = sorted(set(rng.choice(len(ns), 40, replace=False).tolist()) | {0, 1, len(ns) - 1, int(np.argmax(A[2]))})
worst = 0
for i in pick:
    n = int(ns[i]); acc = mp.mpf(0)
    for j in range(J):
        m = n + j * q
        f = factorint(m)
        if len(f) == 1:
            acc += mp.mpf(1) / next(iter(f.values())) * mp.mpf(m) ** (-s)
    m0 = n + J * q
    rest = mp.mpf(m0) ** (-s) + mp.mpf(m0) ** (1 - s) / (q * (s - 1))   # bound for j >= J
    exact_lo = acc * mp.mpf(n) ** s; exact_hi = (acc + rest) * mp.mpf(n) ** s
    d = float(D[i])
    off = 0.0 if exact_lo <= d <= exact_hi else float(min(abs(d - exact_lo), abs(d - exact_hi)))
    worst = max(worst, off)
    if i in pick[:6] or off > 1e-3:
        print(f'n = {n}: D from the run {d:.12f}; exact class sum x n^s in [{mp.nstr(exact_lo, 13)}, {mp.nstr(exact_hi, 13)}]; outside by {off:.2e}')
print(f'{len(pick)} numbers checked; largest distance of D from the exact interval: {worst:.2e}')
