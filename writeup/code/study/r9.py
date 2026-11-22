"""Round 9 (preregister/round9.txt): the variance law along quadratics on fresh data (n in [1e6, 1.05e6)),
Ulam's 195 diagonals and Euler's family n^2 + n + c (c odd, 1..199). Primality by gmpy2.is_prime (no sieve)."""
import os, json, math, time
import numpy as np
import gmpy2
from scipy.signal import fftconvolve
import probes as PR
import theory_geo as TG

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
N0, N = 1_000_000, 50_000
FAM = {'U': PR.curves(), 'E': [(1, 1, c) for c in range(1, 200, 2)]}
# 4n^2 + 4n + c = f_{-4,c}(n + 1): the b = 4 curves are shifted copies of the b = -4 curves (review, point 10)
FAM['U159'] = [f for f in PR.curves() if f[1] != 4]
if os.environ.get('R9_ONLY'):
    FAM = {k: FAM[k] for k in os.environ['R9_ONLY'].split(',')}


def one(f):
    a, b, c = f
    n = np.arange(N0, N0 + N, dtype=np.int64)
    obs = sum(1 for k in range(N0, N0 + N) if gmpy2.is_prime(a * k * k + b * k + c))
    v = (a * n * n + b * n + c).astype(np.float64)
    w = 1.0 / np.log(v)
    C = PR.bh(f)
    E = C * w.sum()
    q = C * w
    R = TG.ratio_Cd(f, N)
    cor = fftconvolve(w, w[::-1])[N:]                  # lags 1..N-1
    var = float(np.sum(q * (1 - q))) + 2 * float(np.sum((R[1:N] - 1) * C * C * cor))
    return obs, E, var / E, C


res = {}
for fam, curves in FAM.items():
    t0 = time.time()
    rows = [one(f) for f in curves]
    obs = np.array([r[0] for r in rows], float); E = np.array([r[1] for r in rows]); VE = np.array([r[2] for r in rows])
    chi = (obs - E) ** 2 / E
    K = len(curves)
    dev = (obs - E) / E
    res[fam] = {'curves': K, 'primes': int(obs.sum()), 'chi2': float(chi.mean()), 'theory': float(VE.mean()),
                'z': float((chi.mean() - VE.mean()) / (chi.std(ddof=1) / math.sqrt(K))),
                'mean_rel_dev_from_BH': float(dev.mean()), 'z_BH': float(dev.mean() / (dev.std(ddof=1) / math.sqrt(K))),
                'C_range': [float(min(r[3] for r in rows)), float(max(r[3] for r in rows))],
                'secs': round(time.time() - t0)}
    print(fam, res[fam], flush=True)
json.dump(res, open(os.path.join(OUT, 'r9%s.json' % ('_' + os.environ['R9_ONLY'] if os.environ.get('R9_ONLY') else '')), 'w'), indent=1)
