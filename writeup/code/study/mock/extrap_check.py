#!/usr/bin/env python3
"""extrap_check.py -- (1) finite-N U(N) cumulants vs N = infinity; (2) linear-in-1/log q extrapolation of the model-R
deficit R_m (a quantity that tends to 0) across the fresh moduli: the method's own systematic for reading T3."""
import json, math, sys
import numpy as np
sys.path.insert(0, '.')
from mock_test import toeplitz_cumulants
print('finite N (Toeplitz exact) vs N=300:')
for s in (1.0, 1.4):
    ref = toeplitz_cumulants(300, s)
    for N in (9, 12, 14):
        k = toeplitz_cumulants(N, s)
        print(f'  s={s} N={N}: k2={k[0]:.5f} k3={k[1]:+.5f} k4={k[2]:+.5f}   (N=300: {ref[0]:.5f} {ref[1]:+.5f} {ref[2]:+.5f})')
cells = json.load(open('reg/predictions_frozen.json'))['cells']
print('model-R deficit: weighted linear fit R_m = alpha + beta/log q over the fresh moduli (alpha should be 0)')
for m in (2, 3, 4):
    for s in (0.8, 1.0, 1.1, 1.2, 1.3, 1.4):
        cl = [c for c in cells.values() if abs(c['s'] - s) < 1e-9 and f'R{m}' in c]
        x = np.array([1 / c['logq'] for c in cl]); y = np.array([c[f'R{m}'] for c in cl])
        w = np.array([1 / c[f'se{m}_exp'] ** 2 for c in cl])
        A = np.vstack([np.ones_like(x), x]).T
        coef = np.linalg.solve(A.T @ (A * w[:, None]), A.T @ (w * y))
        U = cl[0][f'U{m}']
        print(f'  m={m} s={s}: R range {y.min():+.5f}..{y.max():+.5f}  alpha={coef[0]:+.5f}  beta={coef[1]:+.4f}  (U={U:+.5f}, n={len(cl)})')
print('expected se of alpha (extrapolated structured part S_4 at 1/log q -> 0) from the expected per-modulus se:')
for s in (1.1, 1.2, 1.3, 1.4):
    cl = [c for c in cells.values() if abs(c['s'] - s) < 1e-9]
    x = np.array([1 / c['logq'] for c in cl]); w = np.array([1 / c['se4_exp'] ** 2 for c in cl])
    A = np.vstack([np.ones_like(x), x]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    print(f'  s={s}: se(alpha)={math.sqrt(cov[0,0]):.5f}  U4={cl[0]["U4"]:+.5f}  Dpilot4={cl[0]["Dpilot4"]:+.5f}')
