"""Per genus-2 curve (resumable, one json line per curve in r4_arith.jsonl):
(1) moments of lambda_p over good p <= 5e4: E[l^2], E[l^4], mean of sum alpha^2 = 2 lambda(p^2) - lambda(p)^2 (p <= 223), mean l;
(2) bad primes: conductor exponent and PARI's local factor (arithmetic normalisation);
(3) for every 11th curve: functional-equation check of my theta: e^{2 i theta(t)} L(1/2+it)/conj L(1/2+it) = eps at t = 7.3, 23.1, 51.7."""
import json, math, sys, os
import numpy as np, cypari2
sys.path.insert(0, '.')
from r1_offsets import theta, load
pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6); pari.default('realprecision', 38)
rows = [r for r in load() if r['cls'] == 'g2']
def cc(x):
    try: return complex(x)
    except Exception: return complex(float(x))
done = set()
if os.path.exists('r4_arith.jsonl'):
    for line in open('r4_arith.jsonl'): done.add(json.loads(line)['key'])
fo = open('r4_arith.jsonl', 'a')
PM = 50000
primes = [int(p) for p in pari.primes([2, PM])]
for idx, r in enumerate(rows):
    if r['key'] in done: continue
    f, h = r['spec'][1], r['spec'][2]
    L = pari('lfungenus2([%s,%s])' % (f, h))
    N = int(pari.lfunparams(L)[0]); eps = cc(pari.lfunrootres(L)[2]).real
    an = pari.lfunan(L, PM)
    a = np.array([float(an[n - 1]) for n in range(1, PM + 1)]); lam = a / np.sqrt(np.arange(1, PM + 1))
    fa = pari.factor(N); bad = [(int(fa[0][i]), int(fa[1][i])) for i in range(len(fa[0]))]
    good = np.array([p for p in primes if N % p]); l = lam[good - 1]
    g2 = good[good * good <= PM]; sa2 = 2 * lam[g2 ** 2 - 1] - lam[g2 - 1] ** 2
    loc = {str(p): [e, str(1 / pari.lfuneuler(L, p))] for p, e in bad}
    o = dict(key=r['key'], rnd=r['rnd'], N=N, Q_file=r['Q'], eps=eps, eps_file=r['eps'], c_file=r['c'], bad=bad, loc=loc,
             m2=float(np.mean(l ** 2)), m4=float(np.mean(l ** 4)), m6=float(np.mean(l ** 6)), m1=float(np.mean(l)), msa2=float(np.mean(sa2)),
             lam_small={str(p): float(lam[p - 1]) for p in (2, 3, 5, 7, 11, 13)},
             frac_zero=float(np.mean(np.abs(a[good - 1]) < 0.5)))
    if idx % 11 == 0:
        Li = pari.lfuninit(L, [1, 55]); fe = []
        for t in (7.3, 23.1, 51.7):
            v = cc(pari.lfun(Li, pari('1 + %s*I' % t)))      # arithmetic centre s = 1 (k = 2)
            z = complex(np.exp(2j * theta(t, N, [0.5, 0.5, 1.5, 1.5])) * v / v.conjugate()); fe.append([z.real, z.imag])
        o['fe'] = fe
    fo.write(json.dumps(o) + '\n'); fo.flush()
    print('%-14s N %5d eps %+d (file %+d) bad %s m2 %.3f m4 %.3f m6 %.2f msa2 %+.3f m1 %+.4f zero-frac %.3f c %+.3f %s' % (
        o['key'], N, eps, o['eps_file'], bad, o['m2'], o['m4'], o['m6'], o['msa2'], o['m1'], o['frac_zero'], o['c_file'],
        ('FE ' + str(np.round(np.array(o['fe']), 10).tolist())) if 'fe' in o else ''), flush=True)
