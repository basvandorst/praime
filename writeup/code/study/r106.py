"""Round 106 (preregister/round106.txt): the variance of psi(n + h) - psi(n) - h over n in [N/2, N], N = 1e7, divided by
Montgomery-Soundararajan's mean of h (ln(n/h) + B): (a) from a sieve, (b) from the explicit formula with 10,000 zeros
(stratified random x, chunked) plus the diagonal tail above the last zero. usage: r106.py"""
import os, json, math, time
import numpy as np
from scipy import integrate
import common as C
HERE = os.path.dirname(os.path.abspath(__file__))
N = 10 ** 7; HS = (1000, 3000, 10000); B = 1 - 0.5772156649015329 - math.log(2 * math.pi)
t0 = time.time()
# (a) psi on 0 .. N + max(h)
M = N + max(HS) + 2
lam = np.zeros(M)
P = C.small_primes(M).astype(np.int64)
lp = np.log(P.astype(float)); lam[P] = lp
for p, l in zip(P[P * P < M].tolist(), lp[P * P < M].tolist()):
    q = p * p
    while q < M: lam[q] += l; q *= p
psi = np.cumsum(lam); del lam
n = np.arange(N // 2, N + 1)
res = {'N': N, 'rows': []}
Z = np.array([float(l) for l in open(os.path.join(HERE, '..', '..', 'data', 'riemann', 'zeros_10000_flint.txt')) if l.strip()])
T = float(Z[-1])
rng = np.random.default_rng(106); K = 40000
xs = N / 2 + (np.arange(K) + rng.random(K)) * (N / 2) / K
for h in HS:
    d = psi[n + h] - psi[n] - h
    Vm = float(np.mean(d * d)); F = float(np.mean(h * (np.log(n / h) + B)))
    # (b) zero sum S(x) = sum over rho of ((x+h)^rho - x^rho)/rho, rho = 1/2 +- i gamma  ->  2 Re sum over gamma > 0
    S = np.zeros(K)
    lx, lxh = np.log(xs), np.log(xs + h)
    for i0 in range(0, len(Z), 400):
        g = Z[i0:i0 + 400]; rho = 0.5 + 1j * g
        A = (np.exp(np.outer(lxh, rho)) - np.exp(np.outer(lx, rho))) / rho[None, :]
        S += 2 * np.real(A.sum(axis=1))
    Vz_low = float(np.mean(S * S))
    # diagonal tail above T: 2 * x * int_T^inf |(1+d)^rho - 1|^2/|rho|^2 * N0'(t) dt, averaged over x (d = h/x)
    tails = []
    for x in np.linspace(N / 2, N, 9):
        dd = h / x; u = math.log1p(dd); s = math.sqrt(1 + dd)
        f = lambda t: ((1 + dd) + 1 - 2 * s * math.cos(t * u)) / (0.25 + t * t) * math.log(t / (2 * math.pi)) / (2 * math.pi)
        per = 2 * math.pi / u; lo = T; tot = 0.0
        while lo < T + 60 * per:
            tot += integrate.quad(f, lo, lo + per, limit=100)[0]; lo += per
        tot += (2 + dd) * (math.log(lo / (2 * math.pi)) + 1) / (2 * math.pi * lo)
        tails.append(2 * x * tot)
    Vz = Vz_low + float(np.mean(tails))
    row = {'h': h, 'measured_ratio': Vm / F, 'zeros_ratio': Vz / F, 'zeros_low_only_ratio': Vz_low / F, 'tail_ratio': float(np.mean(tails)) / F,
           'V_measured': Vm, 'V_zeros': Vz, 'F': F}
    res['rows'].append(row)
    print('h %5d  measured %.4f  zeros %.4f (low zeros %.4f + tail %.4f)' % (h, Vm / F, Vz / F, Vz_low / F, float(np.mean(tails)) / F), flush=True)
res['T'] = T; res['seconds'] = time.time() - t0
json.dump(res, open(os.path.join(HERE, '..', '..', 'data', 'study', 'r106.json'), 'w'), indent=1)
