"""DRAFT. The blind test: refuses unless sha256(predictions_g2b.json) equals the registered hash.
Primary D1: b = sum Re(G conj law)/sum law^2 over COMPLETE L (sharp (0,150]), jackknife se over L:
  LAW iff |b - 1| <= 2.5 se and b - 0.929 >= 2.5 se; DEFICIT iff |b - 0.929| <= 2.5 se and 1 - b >= 2.5 se; else UNDECIDED.
Reported: Tukey(0.1), band (10,150], lines pooled and lines minus r = 1, the proportional and additive models on (10,150].
--dry: runs the same readings on the 100 genus-2 L-functions of rounds 151/152 and review 58 (no hash check)."""
import sys, os, json, math, hashlib, glob
import numpy as np
from scipy.special import loggamma
HERE = os.path.dirname(os.path.abspath(__file__)); REG = 'REPLACE_WITH_REGISTERED_SHA256'
def theta(t, Q, mu):
    t = np.asarray(t, float); out = t / 2 * np.log(Q / np.pi ** len(mu))
    for m in mu: out = out + np.imag(loggamma((0.5 + m + 1j * t) / 2))
    return out
def tukey(x, f):
    w = np.ones_like(x); e = x < f; w[e] = np.sin(math.pi * x[e] / (2 * f)) ** 2; e = x > 1 - f; w[e] = np.sin(math.pi * (1 - x[e]) / (2 * f)) ** 2; return w
def slope(G, P):
    G = np.asarray(G, complex); P = np.asarray(P, complex); x = (G * P.conj()).real; y = np.abs(P) ** 2; n = len(x)
    loo = (x.sum() - x) / (y.sum() - y); return x.sum() / y.sum(), math.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2))
LINES = {'2': 2.0, '1/2': 0.5, '3': 3.0, '1/3': 1 / 3}
items = []
if '--dry' in sys.argv:
    sys.path.insert(0, os.path.dirname(HERE)); from g2data import load
    sys.path.insert(0, os.path.dirname(HERE))
    for r in load():
        items.append(dict(key=r['key'], g=r['gg'], Q=r['Q'], mu=r['mu'], eps=r['eps'], c=r['c'], lines=None))
else:
    raw = open(os.path.join(HERE, 'predictions_g2b.json'), 'rb').read(); h = hashlib.sha256(raw).hexdigest()
    if h != REG: sys.exit('REFUSED: sha256 %s != registered %s' % (h, REG))
    P = json.loads(raw)['L']; excl = []
    for key, e in P.items():
        fn = os.path.join(HERE, 'zeros_g2b', key + '.npz')
        if not os.path.exists(fn): excl.append((key, 'MISSING')); continue
        z = np.load(fn); info = json.loads(str(z['info']))
        if info['status'] != 'COMPLETE': excl.append((key, info['status'])); continue
        items.append(dict(key=key, g=z['g'][(z['g'] > 0) & (z['g'] <= 150)], Q=e['Q'], mu=e['mu'], eps=e['eps'], c=e['c'], lines=e['lines']))
    print('excluded:', excl)
T = 150.0
for it in items:
    ph = np.exp(2j * theta(it['g'], it['Q'], it['mu'])); it['G'] = complex(ph.sum()); it['law'] = -it['eps'] * it['c'] * T / (2 * math.pi)
    it['Gt'] = complex((tukey(it['g'] / T, 0.1) * ph).sum()); it['Gb'] = complex(ph[it['g'] > 10].sum())
    it['Gl'] = {k: complex((ph * np.exp(-1j * it['g'] * math.log(r))).sum()) for k, r in LINES.items()}
b, se = slope([i['G'] for i in items], [i['law'] for i in items])
dec = 'LAW' if (abs(b - 1) <= 2.5 * se and b - 0.929 >= 2.5 * se) else ('DEFICIT' if (abs(b - 0.929) <= 2.5 * se and 1 - b >= 2.5 * se) else 'UNDECIDED')
print('D1 sharp (0,150]: n = %d  b = %.4f +- %.4f  -> %s' % (len(items), b, se, dec))
print('   Tukey(0.1): %.4f +- %.4f ; band (10,150]: %.4f +- %.4f' % (*slope([i['Gt'] for i in items], [0.9 * i['law'] for i in items]),
      *slope([i['Gb'] for i in items], [i['law'] * 140 / 150 for i in items])))
if items[0]['lines'] is not None:
    x = np.array([sum((i['Gl'][k] * np.conj(i['lines'][k])).real for k in LINES) for i in items]); y = np.array([sum(abs(i['lines'][k]) ** 2 for k in LINES) for i in items])
    n = len(x); bl = x.sum() / y.sum(); loo = (x.sum() - x) / (y.sum() - y); sl = math.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2))
    print('   lines pooled: %.4f +- %.4f ; lines - r=1: %+.4f +- %.4f' % (bl, sl, bl - b, math.hypot(sl, se)))
y = np.array([(i['Gb'] - i['law'] * 140 / 150).real for i in items]); Pn = np.array([-i['law'] * 140 / 150 for i in items])
for nm, X in (('additive', np.ones((len(y), 1))), ('proportional', Pn[:, None]), ('both', np.vstack([np.ones_like(y), Pn]).T)):
    beta = np.linalg.lstsq(X, y, rcond=None)[0]; print('   model %-12s on (10,150]: coef %s rss %.0f' % (nm, np.array2string(beta, precision=3), ((y - X @ beta) ** 2).sum()))
