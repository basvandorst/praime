"""Spot the primes: in a window [N, N + L) a detective sees the true primes and one fake set and must say which is
which, using one statistic at a time and one rule fixed in advance: the primes are the more regular set, so it picks the
set with the smaller statistic (fewer elements divisible by a prime <= 31; smaller gap2, rep, bask, bask1k).  (A first
rule, 'closer to the primes' value in [N0/2, 2 N0]', was dropped: that value drifts with n, and the rule scored below
chance.)  The reference values are still printed.  Fakes: numbers without prime factors in W (W = {2}, {2,3,5}, {2..13}, primes <= 100), chance
prod p/(p-1) / ln n.  Statistics: div (elements with a prime factor <= 31), gap2 (E[g^2]/E[g]^2 of the gaps), rep (share
of repeated last digits), bask, bask1k (variance / mean of the counts in non-overlapping baskets of 100 and of 1,000, detrended by h / ln n).  Accuracy over
R windows with N uniform in [N0/2, N0].  Prediction written before the run: div unmasks W = {2} and {2,3,5} at once;
gap2 and bask unmask the deeper wheels, bask best; nothing is perfect against the primes <= 100 at L = 10^4.
usage: spot.py N0 L R out.txt"""
import sys, math
import numpy as np
N0, L, RR = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); out = open(sys.argv[4], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
M = 2 * N0 + L + 10
sv = np.ones(M + 1, bool); sv[:2] = False
for i in range(2, int(M ** .5) + 1):
    if sv[i]: sv[i * i::i] = False
P = np.nonzero(sv)[0]; small = P[P <= 31]
def stats(S, lo, hi):
    S = S[(S >= lo) & (S < hi)]
    div = int(np.sum(np.any(S[:, None] % small[None, :] == 0, axis=1) & (S > 31)))
    g = np.diff(S).astype(float); gap2 = float((g * g).mean() / g.mean() ** 2) if len(g) > 2 else float('nan')
    d = S % 10; ok = np.isin(d, (1, 3, 7, 9)); dd = d[ok]; rep = float(np.mean(dd[1:] == dd[:-1])) if len(dd) > 2 else float('nan')
    edges = np.arange(lo, hi + 1, 100); c = np.histogram(S, bins=edges)[0].astype(float); mid = (edges[:-1] + edges[1:]) / 2
    e = 100 / np.log(mid); bask = float(np.var(c - e) / e.mean())
    edges = np.arange(lo, hi + 1, 1000); c = np.histogram(S, bins=edges)[0].astype(float); mid = (edges[:-1] + edges[1:]) / 2
    e = 1000 / np.log(mid); bask1k = float(np.var(c - e) / e.mean())
    return dict(div=div, gap2=gap2, rep=rep, bask=bask, bask1k=bask1k)
ref = stats(P, N0 // 2, 2 * N0)
pr('N0 = %d, L = %d, %d windows; reference values of the primes in [N0/2, 2 N0]: gap2 %.4f, rep %.4f, bask %.4f, bask1k %.4f' % (N0, L, RR, ref['gap2'], ref['rep'], ref['bask'], ref['bask1k']))
rng = np.random.default_rng(1818)
for name, W in (('W={2}', [2]), ('W={2,3,5}', [2, 3, 5]), ('W={2..13}', [2, 3, 5, 7, 11, 13]), ('W=primes<=100', list(P[P <= 100]))):
    W = np.array(W); c = float(np.prod(W / (W - 1.0)))
    acc = {k: 0.0 for k in ref}; acc['div'] = 0.0; vals = {k: [] for k in ('div', 'gap2', 'rep', 'bask', 'bask1k')}
    for r in range(RR):
        lo = int(rng.integers(N0 // 2, N0)); hi = lo + L; n = np.arange(lo, hi)
        ok = np.ones(L, bool)
        for p in W: ok &= (n % p != 0)
        F = n[ok & (rng.random(L) < np.minimum(1, c / np.log(n)))]
        sp = stats(P, lo, hi); sf = stats(F, lo, hi)
        for k in ('gap2', 'rep', 'bask', 'bask1k'):
            a, b = sp[k], sf[k]; acc[k] += 1.0 if a < b else (0.5 if a == b else 0.0); vals[k].append(sf[k])
        acc['div'] += 1.0 if sp['div'] < sf['div'] else (0.5 if sp['div'] == sf['div'] else 0.0); vals['div'].append(sf['div'])
    pr('%-15s accuracy: div %.2f  gap2 %.2f  rep %.2f  bask %.2f  bask1k %.2f  | fake means: div %.1f gap2 %.4f rep %.4f bask %.4f bask1k %.4f' % (name, acc['div'] / RR, acc['gap2'] / RR, acc['rep'] / RR, acc['bask'] / RR, acc['bask1k'] / RR, np.mean(vals['div']), np.mean(vals['gap2']), np.mean(vals['rep']), np.mean(vals['bask']), np.mean(vals['bask1k'])))
