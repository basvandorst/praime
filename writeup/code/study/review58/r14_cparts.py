"""Where does the eps-c correlation come from?  Per genus-2 curve: log of the local factor F_p(1/p)/(1 - 1/p) at p <= 13
(own code: local polynomial from lfuneuler, lambda(p^j) by recursion, F_p(X) = sum_j lambda(p^j) mu(p^j) X^j), the bad-prime
factors separately, the sign of the linear (toric) piece at conductor-exponent-1 primes, and an own c for 3 curves
(exact factors p <= 3000, first order to 5e4, plain partial products and the registered log-Cesaro rule)."""
import json, math, sys
import numpy as np, cypari2
sys.path.insert(0, '.')
from r1_offsets import load, jk_mean
pari = cypari2.Pari(); pari.allocatemem(500 * 10 ** 6)
rows = [r for r in load() if r['cls'] == 'g2']
def cc(x):
    try: return complex(x)
    except Exception: return complex(float(x))
def Fp(L, p):
    P = 1 / pari.lfuneuler(L, p)
    co = [cc(c).real for c in pari.Vecrev(P)] if str(pari.type(P)) == 't_POL' else [cc(P).real]
    an = [c * p ** (-j / 2) for j, c in enumerate(co)]          # analytic normalisation (k = 2): coefficient of X^j times p^{-j/2}
    J = 60; lam = [1.0] + [0.0] * J
    for j in range(1, J + 1):
        lam[j] = -sum(an[i] * lam[j - i] for i in range(1, min(j, len(an) - 1) + 1))
    F = sum(lam[j] * an[j] * p ** (-j) for j in range(len(an)))
    return F, co
out = []
for r in rows:
    L = pari('lfungenus2([%s,%s])' % (r['spec'][1], r['spec'][2])); N = int(r['Q'])
    d = {}
    for p in (2, 3, 5, 7, 11, 13):
        F, co = Fp(L, p); d[p] = (math.log(abs(F) / (1 - 1 / p)), F < 0, N % p == 0, co)
    bad = [int(q) for q in pari.factor(N)[0]]
    bd = {}
    for q in bad:
        F, co = Fp(L, q); bd[q] = (F / (1 - 1 / q), co)
    out.append(dict(key=r['key'], eps=r['eps'], c=r['c'], small=d, bad=bd))
eps = np.array([o['eps'] for o in out]); c = np.array([o['c'] for o in out])
print('mean log c (|c|) by eps: +1 %.3f, -1 %.3f' % (np.mean(np.log(np.abs(c[eps > 0]))), np.mean(np.log(np.abs(c[eps < 0])))))
tot = np.zeros(len(out))
for p in (2, 3, 5, 7, 11, 13):
    v = np.array([o['small'][p][0] for o in out]); tot += v
    neg = np.array([o['small'][p][1] for o in out])
    print('p = %2d: mean log F_p/(1-1/p): eps+1 %+.3f, eps-1 %+.3f ; F_p < 0 for %d curves (eps of those: %s)' % (
        p, v[eps > 0].mean(), v[eps < 0].mean(), neg.sum(), list(eps[neg].astype(int))))
print('p <= 13 together: eps+1 %+.3f, eps-1 %+.3f  (difference %.3f; log c difference %.3f)' % (tot[eps > 0].mean(), tot[eps < 0].mean(),
      tot[eps > 0].mean() - tot[eps < 0].mean(), np.mean(np.log(np.abs(c[eps > 0]))) - np.mean(np.log(np.abs(c[eps < 0])))))
# bad primes: exponent-1 primes have a linear toric factor (1 -+ X); its sign vs eps
print('bad primes (q, local factor arithmetic, F_q/(1-1/q)) for the first 12 curves:')
for o in out[:12]:
    print('  %-14s eps %+d c %+.3f  %s' % (o['key'], o['eps'], o['c'], {q: (round(v[0], 4), [round(x, 3) for x in v[1]]) for q, v in o['bad'].items()}))
# own c for three curves
pr = [int(p) for p in pari.primes([2, 50000])]
for r in rows[:2] + rows[-1:]:
    L = pari('lfungenus2([%s,%s])' % (r['spec'][1], r['spec'][2]))
    an = pari.lfunan(L, 50000); lam = np.array([float(an[n - 1]) for n in range(1, 50001)]) / np.sqrt(np.arange(1, 50001))
    lg = 0.0; sgn = 1; cum = []
    for p in pr:
        if p <= 3000:
            F, _ = Fp(L, p); sgn *= (1 if F > 0 else -1); lg += math.log(abs(F) / (1 - 1 / p))
        else:
            lg += math.log((1 - lam[p - 1] ** 2 / p) / (1 - 1 / p))
        cum.append(lg)
    cum = np.array(cum); P = np.array(pr); sel = P >= 5000; w = 1.0 / P[sel]; corr = 2 / (np.sqrt(P[sel]) * np.log(P[sel]))
    ces = sgn * math.exp(np.sum((cum[sel] + corr) * w) / w.sum())
    print('own c %-14s: partial products at 1e4 %+.4f, 5e4 %+.4f, log-Cesaro %+.4f ; registered %+.4f' % (
        r['key'], sgn * math.exp(cum[P <= 10000][-1]), sgn * math.exp(cum[-1]), ces, r['c']))
json.dump([dict(key=o['key'], eps=o['eps'], c=o['c'], small={str(p): v[:3] for p, v in o['small'].items()},
                bad={str(q): [v[0], v[1]] for q, v in o['bad'].items()}) for o in out], open('r14_cparts.json', 'w'), indent=1)
