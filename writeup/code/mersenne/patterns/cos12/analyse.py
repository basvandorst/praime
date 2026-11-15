"""Cosine rows with period 12: the lanes of the primes, consecutive primes, twin pairs, the exponents."""
import numpy as np, math
from scipy.stats import binom
def sieve(N):
    s = np.ones(N + 1, bool); s[:2] = False
    for i in range(2, int(N ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return np.flatnonzero(s)
UP = {1, 11}
for N in (500000, 10 ** 7):
    P = sieve(N); P = P[P >= 5]; r = P % 12
    up = np.isin(r, [1, 11]); cnt = {a: int((r == a).sum()) for a in (1, 5, 7, 11)}
    print(f'\n== primes 5..{N:,}: {len(P):,}; by remainder mod 12: {cnt}; up lane (1, 11) {up.sum():,} = {up.mean():.4f}, down (5, 7) {(~up).sum():,}')
    lead = np.cumsum(np.where(up, -1, 1))
    print(f'   down minus up along the way: final {lead[-1]}, minimum {lead.min()}, below zero at {int((lead < 0).sum())} of {len(lead):,} primes (first at p = {int(P[np.argmax(lead < 0)]) if (lead < 0).any() else "-"})')
    print(f'   class 1 against the mean of 5, 7, 11: {cnt[1] - (cnt[5] + cnt[7] + cnt[11]) / 3:+.0f}; sqrt(x)/ln x = {math.sqrt(N) / math.log(N):.0f}')
    g = np.diff(P); same = up[1:] == up[:-1]
    print(f'   consecutive primes in the same lane: {same.mean():.4f} (switch {1 - same.mean():.4f}); a random lane order would give {up.mean() ** 2 + (1 - up.mean()) ** 2:.4f}')
    for G in (2, 4, 6, 8, 10, 12):
        m = g % 12 == G % 12
        print(f'     gap = {G} mod 12: {m.mean():.4f} of all gaps, same lane {same[m].mean():.3f}')
    T = np.zeros((4, 4), int); idx = {1: 0, 5: 1, 7: 2, 11: 3}
    for a, b in zip(r[:-1], r[1:]): T[idx[a], idx[b]] += 1
    exp = np.outer(T.sum(1), T.sum(0)) / T.sum()
    print('   consecutive remainders (rows: this prime 1,5,7,11; columns: next), observed / independent:')
    for i, a in enumerate((1, 5, 7, 11)): print('     ', a, ' '.join(f'{T[i, j]:>7}/{exp[i, j]:>8.0f}' for j in range(4)))
    tw = P[:-1][g == 2]; tu = int((tw % 12 == 11).sum()); td = int((tw % 12 == 5).sum())
    print(f'   twin pairs (p, p + 2): {len(tw):,}; both up (p = 11 mod 12) {tu:,}, both down (p = 5 mod 12) {td:,}; other lanes: {len(tw) - tu - td}')
E = [5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933]
c = {a: sum(1 for k in E if k % 12 == a) for a in (1, 5, 7, 11)}
Ecl = {1: 11.81, 5: 13.39, 7: 11.85, 11: 11.27}
pu = (Ecl[1] + Ecl[11]) / sum(Ecl.values()); nu = c[1] + c[11]; n = len(E)
print(f'\n== exponents 5 .. 82,589,933 (complete): {n}; by remainder mod 12 {c}; expected (heuristic, a = 6/2) {Ecl}')
print(f'   up lane {nu} of {n} ({nu / n:.2f}), expected share {pu:.3f}; P(up <= {nu}) = {binom.cdf(nu, n, pu):.4f}; two-sided {min(1, 2 * binom.cdf(nu, n, pu)):.4f}')
for lo, hi, nm in ((5, 10 ** 5, 'below 10^5'), (10 ** 5, 10 ** 9, 'from 10^5')):
    S = [k for k in E if lo <= k < hi]; u = sum(1 for k in S if k % 12 in UP)
    print(f'   {nm}: up {u} of {len(S)}; P(up <= {u}) = {binom.cdf(u, len(S), pu):.3f}')
print(f'   no. 52, 136,279,841 = {136279841 % 12} mod 12 ({"up" if 136279841 % 12 in UP else "down"})')
