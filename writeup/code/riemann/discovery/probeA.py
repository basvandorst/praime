"""Probe A: how much of each L-function's low zeros is fixed by the character's values at small primes?
Data: all zeros with |gamma| <= 15.41 of the 395,136 Dirichlet L-functions mod 395,137.
Weil's explicit formula says a smooth zero statistic sum_gamma f(gamma) is exactly linear in Re chi(p^k):
  sum_gamma f(gamma) = (f-dependent constant, by parity) - 2 sum_{p^k} Re chi(p^k) ln p fhat(k ln p) / p^(k/2) / (2 pi)
so its regression on Re chi(p^k) must approach R^2 = 1 with predicted coefficients. The lowest zero is not smooth;
how much of it the first primes fix is the measured quantity."""
import math, sys
import numpy as np
D = np.load('lfam395137f.npz'); q = int(D['q']); pw = D['pw']; kz = D['k'].astype(np.int64); gz = D['gamma']
N = q - 1
ind = np.empty(q, dtype=np.int64); ind[pw] = np.arange(N)
ks = np.arange(1, N)                                   # non-principal characters
kbar = (-ks) % N
# lowest |gamma| per character: min over its zeros with gamma > 0 and its conjugate's
low = np.full(N, np.inf); np.minimum.at(low, kz, gz)
lowest = np.minimum(low[ks], low[kbar])
parity = ks % 2
# smooth statistic: S_f = sum over all zeros with |gamma| <= T of f(gamma), f(g) = exp(-(g/2)^2) (tail beyond T is negligible)
f = np.exp(-(gz / 2.0) ** 2)
Sp = np.bincount(kz, weights=f, minlength=N)
Sf = Sp[ks] + Sp[kbar]
primes = [p for p in range(2, 200) if all(p % r for r in range(2, int(p ** 0.5) + 1))]
def feats(P, kmax=1):
    cols = [np.ones(len(ks)), parity.astype(float)]
    names = ['1', 'odd']
    for p in primes:
        if p > P: break
        for e in range(1, kmax + 1):
            th = 2 * math.pi * ((ks * (e * ind[p] % N)) % N) / N
            cols.append(np.cos(th)); names.append(f'Re chi({p}^{e})')
    return np.stack(cols, 1), names
def r2(Xm, y):
    c, *_ = np.linalg.lstsq(Xm, y, rcond=None); res = y - Xm @ c
    return 1 - res.var() / y.var(), c
out = []
for target, y in (('smooth sum of exp(-(gamma/2)^2) over the zeros', Sf), ('lowest zero |gamma_1|', lowest), ('log of the lowest zero', np.log(lowest))):
    row = []
    for P in (2, 3, 5, 7, 13, 29, 59, 113, 199):
        Xm, names = feats(P, kmax=2)
        R, c = r2(Xm, y); row.append(f'p <= {P}: {R:.3f}')
    out.append(f'{target}: R^2 from Re chi(p), Re chi(p^2) and parity: ' + '; '.join(row))
# predicted coefficients for the smooth statistic: d S_f / d Re chi(p) = -2 ln p fhat(ln p)/sqrt p / (2 pi), fhat(x) = int f(g) e^{-igx} dg = 2 sqrt(pi) e^{-x^2}
Xm, names = feats(29, kmax=1)
R, c = r2(Xm, Sf)
pred = {p: -2 * math.log(p) * 2 * math.sqrt(math.pi) * math.exp(-math.log(p) ** 2) / math.sqrt(p) / (2 * math.pi) for p in primes if p <= 29}
out.append('smooth statistic, fitted against predicted coefficient of Re chi(p): ' + '; '.join(f'{p}: {c[2 + i]:+.5f} vs {pred[p]:+.5f}' for i, p in enumerate([p for p in primes if p <= 29])))
# the extremes: characters with chi(p) = 1 for many small p
agree = np.zeros(len(ks), dtype=int)
for p in primes[:10]:
    agree += ((ks * ind[p]) % N == 0)
for a in range(0, 4):
    s = agree == a
    if s.sum(): out.append(f'characters with chi(p) = 1 for exactly {a} of the first 10 primes: {int(s.sum())}, median lowest zero {np.median(lowest[s]):.4f}')
s = agree >= 4
out.append(f'characters with chi(p) = 1 for 4 or more of the first 10 primes: {int(s.sum())}, median lowest zero {np.median(lowest[s]) if s.sum() else float("nan"):.4f}')
# which single prime matters most for the lowest zero
Xm, names = feats(29, kmax=1); R, c = r2(Xm, lowest)
out.append('lowest zero, coefficient per Re chi(p) (p <= 29): ' + ', '.join(f'{names[i]} {c[i]:+.4f}' for i in range(2, len(names))))
out.append(f'lowest zero over all characters: min {lowest.min():.6f} (character {int(ks[np.argmin(lowest)])}), median {np.median(lowest):.4f}, max {lowest.max():.4f} (character {int(ks[np.argmax(lowest)])})')
print('\n'.join(out))
open('probeA.txt', 'w').write('\n'.join(out) + '\n')
