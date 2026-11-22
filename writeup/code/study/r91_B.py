"""Round 91: the second-order constant B, derived. Var(psi(x+H) - psi(x) - H) = H(ln x - 1) (diagonal: E Lambda^2 - 1)
+ 2 sum_{h<H} (H - h)(S(h) - 1); the pair sum is evaluated numerically and fitted to -1/2 H ln H + c H.
Result: c = (2 - gamma - ln 2pi)/2, so Var = H(ln(x/H) + 1 - gamma - ln 2pi), B = 1 - gamma - ln 2pi. usage: r91_B.py"""
import json, math, os
import numpy as np
import common as C

C2 = 0.6601618158468696; G = 0.5772156649015329
out = {}
for H in (10 ** 4, 10 ** 5, 10 ** 6):
    h = np.arange(1, H); S = np.where(h % 2 == 0, 2 * C2, 0.0)
    for p in C.small_primes(H).tolist():
        if p > 2:
            S[p - 1::p] *= (p - 1) / (p - 2)
    F = float(np.sum((H - h) * (S - 1.0)))
    out[str(H)] = (F + 0.5 * H * math.log(H)) / H
out['(2 - gamma - ln 2pi)/2'] = (2 - G - math.log(2 * math.pi)) / 2
out['B = 1 - gamma - ln 2pi'] = 1 - G - math.log(2 * math.pi)
print(out)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study', 'r91_B.json'), 'w'), indent=1)
