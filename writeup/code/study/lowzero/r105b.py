"""Round 105 follow-up (not pre-registered): the race between chi_D = -1 and +1 primes predicted from the zeros alone.
Normalised race E(x) = 1 + 2 sum_{0 < gamma <= T} Re(x^(i gamma)/rho), rho = 1/2 + i gamma (Rubinstein-Sarnak); the
-1 primes lead where E > 0. Predicted share of primes p <= 1e5 with E(p) > 0 against the measured share (round 105),
with the zeros up to T = 30 and with the lowest zero alone. usage: r105b.py"""
import json, math, os
import numpy as np, cypari2
from sympy import primerange
from scipy.stats import spearmanr, pearsonr
pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6)
H = os.path.dirname(os.path.abspath(__file__))
rows = json.load(open(os.path.join(H, 'r105.json')))['rows']
P = np.array(list(primerange(2, 10 ** 5)), float); lnP = np.log(P)
for r in rows:
    z = np.array([float(t) for t in pari.lfunzeros(pari.lfuncreate(r['D']), 30)])
    z = z[z > 0]
    def share(zz):
        E = np.ones_like(lnP)
        for g in zz:
            rho = complex(0.5, g)
            E += 2 * np.real(np.exp(1j * g * lnP) / rho)
        return float(np.mean(E > 0))
    r['pred_all'] = share(z); r['pred_one'] = share(z[:1]); r['nzeros'] = int(len(z))
out = {}
for name, sel in (('negative', lambda r: r['D'] < 0), ('positive', lambda r: r['D'] > 0), ('all', lambda r: True)):
    R = [r for r in rows if sel(r)]
    b = np.array([r['bias'] for r in R]); pa = np.array([r['pred_all'] for r in R]); p1 = np.array([r['pred_one'] for r in R])
    out[name] = {'n': len(R), 'pearson_all': pearsonr(pa, b)[0], 'spearman_all': spearmanr(pa, b)[0], 'mean_abs_err_all': float(np.mean(np.abs(pa - b))),
                 'pearson_one': pearsonr(p1, b)[0], 'mean_abs_err_one': float(np.mean(np.abs(p1 - b)))}
    print(name, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in out[name].items()})
json.dump({'tests': out, 'rows': rows}, open(os.path.join(H, 'r105b.json'), 'w'), indent=0)
