"""Independent k-th neighbour law, y = 7 (P = 210). Writes law.json.
Variants: full (order 2, exponentiated), noR (drop pair term = independent
Bernoulli with HL p_i), lin2 (order 2 Bahadur, not exponentiated),
S1 (S_{>7} == 1: sieve by q<=7 only), cramer (all odd n prime w.p. 2/lambda)."""
import json, math, sys
Y = int(sys.argv[1]) if len(sys.argv) > 1 else 7
SMALLP = [q for q in (2, 3, 5, 7, 11, 13, 17) if q <= Y]
import numpy as np
from math import gcd

LAM = math.log(1.05e9)
P = 1
for _q in SMALLP: P *= _q
DMAX = 240
K = 4                      # coefficients z^0..z^3
QMAX = 10**5

def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]

PR = primes_upto(QMAX)
QS = PR[(PR > Y) & (PR <= DMAX)].astype(float)          # explicit residue counting
QL = PR[PR > DMAX].astype(float)                          # all residues distinct
def Lbig(m):  # sum over DMAX<q<=QMAX plus tail -C(m,2)/(Q ln Q)
    return float(np.sum(np.log1p(-m / QL) - m * np.log1p(-1 / QL))) - m*(m-1)/2/(QMAX*math.log(QMAX))
LB = {m: Lbig(m) for m in (2, 3, 4)}
PHI = 1.0
for q in SMALLP: PHI *= 1 - 1/q
RHO0 = 1 / (LAM * PHI)
lq = np.log1p(-1 / QS)

def logS_small(nu, m):
    return np.sum(np.log1p(-nu / QS) - m * lq, axis=-1)

def series_mul(a, b):
    c = np.zeros(K)
    for i in range(K):
        for j in range(K - i): c[i+j] += a[i]*b[j]
    return c

def series_exp(f):
    g = np.zeros(K); g[0] = 1.0
    h = f.copy(); h[0] = 0
    for n in range(1, K):
        g[n] = sum(k*h[k]*g[n-k] for k in range(1, n+1)) / n
    return g * math.exp(f[0])

def coeffs(p, R, mode):
    """Coefficients z^0..z^3 of the generating function for interior chances p, pair matrix R."""
    n = len(p)
    t = p / (1 - p)
    # ln(1-p+pz) series
    L = np.zeros((n, K)); L[:, 0] = np.log1p(-p)
    for m in range(1, K): L[:, m] = (-1)**(m+1) * t**m / m
    f = L.sum(0)
    if mode == 'noR' or n == 0:
        return series_exp(f)
    s = np.stack([(-t)**m for m in range(K)], 1)            # 1/(1+tz)
    w = np.zeros((n, K)); w[:, 0] = -t * s[:, 0]
    for m in range(1, K): w[:, m] = t * (s[:, m-1] - s[:, m])
    pair = np.zeros(K)
    for a in range(K):
        Rw = R @ w[:, a]
        for b in range(K - a): pair[a+b] += 0.5 * w[:, b] @ Rw
    if mode == 'full':
        return series_exp(f + pair)
    if mode == 'lin2':
        one = np.zeros(K); one[0] = 1
        return series_mul(series_exp(f), one + pair)
    raise ValueError(mode)

classes = [a for a in range(P) if gcd(a, P) == 1]
modes = ['full', 'noR', 'lin2', 'S1', 'cramer']
out = {m: np.zeros((DMAX + 1, K)) for m in modes}   # P(D is k-th gap), averaged over classes
for a in classes:
    for D in range(2, DMAX + 1, 2):
        if gcd(a + D, P) != 1: continue
        I = np.array([i for i in range(2, D, 2) if gcd(a + i, P) == 1], dtype=float)
        rD = D % QS
        nu2 = 1 + (rD != 0)
        l2 = logS_small(nu2, 2) + LB[2]
        if len(I):
            ri = I[:, None] % QS[None, :]
            ni = (ri != 0) & (ri != rD)
            l3 = logS_small(nu2 + ni, 3) + LB[3]
            nu4 = nu2 + ni[:, None, :] + ni[None, :, :] - (ni[:, None, :] & ni[None, :, :] & (ri[:, None, :] == ri[None, :, :]))
            l4 = logS_small(nu4, 4) + LB[4]
            p = RHO0 * np.exp(l3 - l2)
            R = np.exp(l4 + l2 - l3[:, None] - l3[None, :]) - 1
            np.fill_diagonal(R, 0.0)
        else:
            p = np.zeros(0); R = np.zeros((0, 0))
        base = RHO0 * math.exp(l2)
        for m in ('full', 'noR', 'lin2'):
            out[m][D] += base * coeffs(p, R, m)
        pS = np.full(len(I), RHO0)
        out['S1'][D] += RHO0 * coeffs(pS, None, 'noR')
for m in ('full', 'noR', 'lin2', 'S1'): out[m] /= len(classes)
# Cramer with parity: every odd number prime w.p. 2/lambda independently (no class sieve)
r = 2 / LAM
for D in range(2, DMAX + 1, 2):
    nint = D // 2 - 1
    for k in range(K):
        out['cramer'][D, k] = r * math.comb(nint, k) * r**k * (1 - r)**(nint - k)

res = {'lambda': LAM, 'rho0': RHO0, 'LB': LB}
for m in modes:
    res[m] = {}
    for Dm in (200, 240):
        A = out[m][:Dm + 1]
        tot = A.sum(0)
        res[m][Dm] = {'mass': tot.tolist(),
                      'rep': {k+1: {q: float(A[::q][:, k].sum() / tot[k]) if True else 0
                                    for q in (3, 5, 7, 11, 13)} for k in range(K)}}
    res[m]['perD'] = out[m].tolist()
# Note A[::q] picks D = 0, q, 2q, ... (D=0 row is zero)
json.dump(res, open(f'/home/bas/scratch/rebuild-C/law_y{Y}.json', 'w'), indent=1)
for m in modes:
    print(m, 'mass<=240', np.round(res[m][240]['mass'], 6))
    for k in range(1, K+1):
        print('  k', k, {q: round(v, 5) for q, v in res[m][240]['rep'][k].items()})
