"""z and chi^2 of each fake against the r81 law and against the fake's own exact prediction (cramer, sieve13)."""
import json, math, glob
import numpy as np
from math import comb, gcd

R = json.load(open('/home/bas/praime/writeup/data/study/r81.json'))
LAW = {int(k): {int(q): v for q, v in d.items()} for k, d in R['theory'].items()}
KS = (1, 2, 3, 4); QS = (3, 5, 7, 11, 13); DMAX = 240
Ds = np.arange(2, DMAX + 1, 2)
LAMS = np.log(np.linspace(1e9, 1.1e9, 9))         # average over the window, weighted by density

def own_cramer():
    out = {k: {q: 0.0 for q in QS} for k in KS}
    for k in KS:
        w = np.zeros(len(Ds))
        for lam in LAMS:
            r = 2 / lam
            w += np.array([comb(D // 2 - 1, k - 1) * r ** k * (1 - r) ** (D // 2 - k) for D in Ds]) / lam
        for q in QS:
            out[k][q] = w[Ds % q == 0].sum() / w.sum()
    return out

def own_sieve13():
    PP = 30030; cop = np.array([gcd(x, PP) == 1 for x in range(PP + DMAX + 1)])
    A = np.nonzero(cop[:PP])[0]
    # m[a, D] = # coprime positions strictly between a and a+D; adm[a, D]
    cs = np.cumsum(cop)
    M = np.array([[cs[a + D - 1] - cs[a] for D in Ds] for a in A])
    ADM = np.array([[cop[a + D] for D in Ds] for a in A])
    out = {k: {} for k in KS}
    for k in KS:
        w = np.zeros(len(Ds))
        for lam in LAMS:
            r = 1 / (lam * 5760 / 30030)
            binom = np.vectorize(lambda m: comb(int(m), k - 1))(M).astype(float)
            w += (ADM * binom * r ** k * (1 - r) ** (M - k + 1)).sum(0) / lam
        for q in QS:
            out[k][q] = w[Ds % q == 0].sum() / w.sum()
    return out

OWN = {'cramer': own_cramer(), 'sieve13': own_sieve13()}
# power: law vs the sieve13 exact prediction, in sigma units at n = 4.8M-like counts
def zs(f, pred):
    Z = {}
    for k in KS:
        for q in QS:
            n = f['nin'][str(k)]; m = f['cnt'][str(k)][str(q)] / n; se = math.sqrt(m * (1 - m) / n)
            Z[(k, q)] = ((pred[k][q] - m) / se, m)
    return Z

def show(name, Z):
    z = np.array([v[0] for v in Z.values()])
    print('  vs %-8s max|z| %6.1f  chi2 %9.1f (20 dof)  #|z|>2: %d' % (name, np.abs(z).max(), (z ** 2).sum(), (np.abs(z) > 2).sum()))
    return z

for fn in sorted(glob.glob('/home/bas/scratch/null-D/out_*.json')):
    f = json.load(open(fn)); print(f['model'], 'seed', f['seed'], 'N', f['nprimes'])
    Zl = zs(f, LAW); show('law', Zl)
    print('   z(law) by k (q=3,5,7,11,13): ' + ' | '.join(' '.join('%+.1f' % Zl[(k, q)][0] for q in QS) for k in KS))
    if f['model'] in OWN:
        Zo = zs(f, OWN[f['model']]); show('own', Zo)
# the primes themselves (r81 data)
zp = np.array([R['data'][str(k)][str(q)]['z'] for k in KS for q in QS])
print('primes (r81.json)  vs law  max|z| %.1f chi2 %.1f' % (np.abs(zp).max(), (zp ** 2).sum()))
print('\nlaw minus sieve13-own prediction, in units of the primes\' se (r81.json):')
for k in KS:
    print(' k=%d ' % k + ' '.join('q%d %+.4f (%+.1f se)' % (q, LAW[k][q] - OWN['sieve13'][k][q], (LAW[k][q] - OWN['sieve13'][k][q]) / R['data'][str(k)][str(q)]['se']) for q in QS))
print('\nprimes vs sieve13-own: ' + ' '.join('%+.1f' % ((OWN['sieve13'][k][q] - R['data'][str(k)][str(q)]['p']) / R['data'][str(k)][str(q)]['se']) for k in KS for q in QS))
json.dump({m: {str(k): {str(q): v for q, v in d.items()} for k, d in o.items()} for m, o in OWN.items()}, open('/home/bas/scratch/null-D/own_pred.json', 'w'), indent=1)
