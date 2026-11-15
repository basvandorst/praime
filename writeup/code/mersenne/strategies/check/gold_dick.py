"""Independent checks of the goldbach and dickman notes (numpy, direct counting; no FFT)."""
import numpy as np, math, sympy
N = 500000
s = np.ones(N + 1, bool); s[:2] = False
for i in range(2, int(N ** .5) + 1):
    if s[i]: s[i * i::i] = False
P = np.flatnonzero(s)
g = np.zeros(N + 1, np.int64)
for p in P:                                  # direct: for each p, count q >= p with p + q = n
    q = P[P >= p]; n = p + q; n = n[n <= N]; np.add.at(g, n, 1)
ev = np.arange(4, N + 1, 2)
print('even n with g = 0:', int((g[ev] == 0).sum()))
big = ev[ev >= 1000]; m = g[big].min(); print('min g above 1000:', int(m), 'at', big[g[big] == m].tolist()[:3], '; g(1112) by sympy:', sum(1 for p in sympy.primerange(2, 557) if sympy.isprime(1112 - p)))
v = g[ev[ev >= 6]] * np.log(ev[ev >= 6]) ** 2 / ev[ev >= 6]; e6 = ev[ev >= 6]
print('band ratio (mult of 3 vs not):', round(float(v[e6 % 3 == 0].mean() / v[e6 % 3 != 0].mean()), 3), '; mult of 5 vs not (among non-mult of 3):', round(float(v[(e6 % 5 == 0) & (e6 % 3 != 0)].mean() / v[(e6 % 5 != 0) & (e6 % 3 != 0)].mean()), 3), '(HL 4/3)')
for M in (10**4, 10**5, 5 * 10**5, 10**6, 10**7):
    lp = np.zeros(M + 1, np.int64)
    ss = np.ones(M + 1, bool); ss[:2] = False
    for i in range(2, int(M ** .5) + 1):
        if ss[i]: ss[i * i::i] = False
    for p in np.flatnonzero(ss): lp[p::p] = p
    n = np.arange(2, M + 1)
    print(f'N={M:>9,}: share P(n) <= sqrt n {100 * np.mean(lp[2:] ** 2 <= n):.2f}%, P(n) <= n^(1/3) {100 * np.mean(lp[2:].astype(float) ** 3 <= n):.2f}%')
print('rho(2) = 1 - ln 2 =', round(1 - math.log(2), 4))
