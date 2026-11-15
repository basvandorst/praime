"""Exponent-set check with the independent LL of llind.py.
(1) the 28 listed exponents < 10^5 all pass LL (proved prime M_p);
(2) every prime p < 2*10^4 is either removed by a proved factor (TF kmin file: 2kp+1 | M_p, or a P-1 factor g | M_p)
    or LL-tested here; the exponents found must be exactly the listed ones below 2*10^4.
Above 2*10^4 the completeness of the list is the published GIMPS verification, not checked here."""
import json, numpy as np, time
from llind import ll, primes_below
MEXP = [2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243]
t = time.time()
bad = [p for p in MEXP if not ll(p)]
print(f'(1) LL on the 28 listed exponents: failures {bad} ({time.time() - t:.0f} s)')
P, km = np.load('kmin_100000_1000000.npy'); kmin = {int(p): int(k) for p, k in zip(P, km)}
pm = {}
for f in ['pm1_1e4.json', 'pm1_1e5a.json']:
    for p, r in json.load(open(f))['res'].items():
        if r['g'] is not None: pm[int(p)] = int(r['g'])
X = 20000; t = time.time(); found = []; nll = 0
for p in primes_below(X):
    M = (1 << p) - 1
    if p > 2 and kmin[p]:
        q = 2 * kmin[p] * p + 1; assert pow(2, p, q) == 1 and q < M; continue
    if p in pm:
        g = pm[p]; assert 1 < g < M and M % g == 0; continue
    nll += 1
    if ll(p): found.append(p)
print(f'(2) p < {X}: {nll} LL tests after the proved removals; exponents {found}; '
      f'match the list: {found == [m for m in MEXP if m < X]} ({time.time() - t:.0f} s)')
