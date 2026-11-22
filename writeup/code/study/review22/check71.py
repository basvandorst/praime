"""Round 71 check, own code: (1) the law value at d = 6 from the Hardy-Littlewood series (no interior twin is possible),
(2) the HL expectation of the twin count in [1e15, 1e15+1e10), (3) recount of twins and twin gaps in 1e8 sub-windows,
checked with gmpy2, against the stored law. usage: check71.py"""
import json, math, sys
import numpy as np
import gmpy2
import mpmath
sys.path.insert(0, '.')
from sieve import primes_in, base_primes

r = json.load(open('/home/bas/praime/writeup/data/study/r71.json'))
out = {}
# (1) d = 6: P(next twin start at 6 | twin at 0) = S({0,2,6,8}) / (S({0,2}) lam^2)
ps = base_primes(10 ** 7)
lr = 0.0
for q in ps.tolist():
    n4 = len({0, 2 % q, 6 % q, 8 % q}); n2 = len({0, 2 % q})
    lr += math.log(1 - n4 / q) - math.log(1 - n2 / q) - 2 * math.log(1 - 1 / q)
# tail beyond 1e7: each q > 8 has n4 = 4, n2 = 2: log(1-4/q) - log(1-2/q) - 2 log(1-1/q) ~ -(1/q^2)(8-2-1)... tiny; estimate
tail = -sum((16 - 4 - 2) / 2 / q ** 2 for q in [1e7]) * 0  # negligible (< 1e-7)
for i, lam in enumerate([math.log(1e15), math.log(1e15 + 5e9), math.log(1e15 + 1e10)]):
    out['law6_own_lam%d' % i] = math.exp(lr) / lam ** 2
out['law6_stored'] = r['law']['6']
# (2) HL expectation of twin pairs: 2 C2 * integral dx / ln^2 x
C2 = 0.6601618158468696
I = mpmath.quad(lambda x: 1 / mpmath.log(x) ** 2, [1e15, 1e15 + 1e10])
out['twins_HL_expect'] = float(2 * C2 * I); out['twins_stored'] = r['twins']
# (3) recount in sub-windows
law = {int(k): v[0] for k, v in r['law'].items()}
sub = {}
for L in (10 ** 15, 10 ** 15 + 5 * 10 ** 9):
    P = primes_in(L, 10 ** 8)
    # gmpy2 cross-check on a stretch
    x = int(P[5000]); ch = [x]
    for _ in range(3000):
        x = int(gmpy2.next_prime(x)); ch.append(x)
    assert np.array_equal(P[5000:5000 + 3001], np.array(ch))
    tw = P[:-1][np.diff(P) == 2]
    dd = np.diff(tw)
    assert np.all(dd % 6 == 0)
    bands = [(lo, lo + 24) for lo in range(6, 600, 30)] + [(lo, lo + 54) for lo in range(606, 1200, 60)]
    N = len(dd); chi = 0.0
    for lo, hi in bands:
        c = int(np.sum((dd >= lo) & (dd <= hi))); e = N * sum(law.get(d, 0.0) for d in range(lo, hi + 1, 6))
        chi += (c - e) ** 2 / e
    sub[str(L)] = {'primes': int(len(P)), 'twins': int(len(tw)), 'gaps': int(N), 'chi2_30_bands': chi,
                   'share_gt_1200': float(np.mean(dd > 1200)), 'twins_HL': float(2 * C2 * mpmath.quad(lambda x: 1 / mpmath.log(x) ** 2, [L, L + 10 ** 8]))}
out['subwindows'] = sub
out['law_mass_to_1200_stored'] = r['compare']['law_mass_to_1200']
print(json.dumps(out, indent=1))
json.dump(out, open('check71.json', 'w'), indent=1)
