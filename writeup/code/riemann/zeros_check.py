"""Checks zeros written by zeros.js against mpmath.zetazero at sample indices.

usage: python3 zeros_check.py zeros.f64 [samples]
"""
import random
import sys
import time

import mpmath
import numpy as np

mpmath.mp.dps = 20
Z = np.fromfile(sys.argv[1], dtype=np.float64)
ns = int(sys.argv[2]) if len(sys.argv) > 2 else 60
random.seed(5)
fixed = [1, 10, 31, 32, 50, 100, 150, 200, 300, 500, 700, 1000, 1500, 2000, 3000, 5000, 7000, 10000,
         20000, 50000, 100000, 300000, 1000000, 1999999, 2000000]
idx = sorted({n for n in fixed if n <= len(Z)} | {random.randint(31, len(Z)) for _ in range(ns)})
worst, at, t0 = 0.0, 0, time.time()
for n in idx:
    e = abs(float(mpmath.zetazero(n).imag) - Z[n - 1])
    if e > worst:
        worst, at = e, n
print(f'{len(Z)} zeros, {len(idx)} checked against mpmath.zetazero: '
      f'max error {worst:.2e} at zero #{at} [{time.time() - t0:.0f} s]')
