"""H4: own check of stored zeros for one heavy curve: recompute zeros on (a, b] at realprecision 57 (stored: default
precision of the registered routine), compare positions; independent count of sign changes of the Hardy function on a grid;
window sums from both lists. Root number set from lfunrootres before lfuninit (lfungenus2 leaves it unset).
usage: h4_zerocheck.py key a b"""
import sys, os, math, time; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, cypari2
from g2data import load, theta
key, a, b = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
r = [x for x in load() if x['key'] == key][0]
pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6); pari.default('realprecision', 57)
t0 = time.time()
L = pari('my(L = lfungenus2([%s,%s])); L[6] = lfunrootres(L)[3]; L' % (r['spec'][1], r['spec'][2]))
Li = pari.lfuninit(L, [1, int(b) + 3]); t1 = time.time() - t0
z = np.sort(np.array([float(x) for x in pari.lfunzeros(Li, [a, b])]))
gs = r['g'][(r['g'] > a) & (r['g'] <= b)]
print('%s (%g, %g]: stored %d zeros, recomputed at prec 57: %d zeros (init %.0f s, total %.0f s)' % (key, a, b, len(gs), len(z), t1, time.time() - t0))
if len(z) == len(gs): print('  max |diff| = %.2e' % np.max(np.abs(z - gs)))
else: print('  stored only:', [x for x in gs if np.min(np.abs(z - x)) > 1e-6], ' new only:', [x for x in z if np.min(np.abs(gs - x)) > 1e-6])
tt = np.arange(a + 0.005, b, 0.02); hv = np.array([float(pari.lfunhardy(Li, t)) for t in tt])
sc = int(np.sum(np.sign(hv[1:]) != np.sign(hv[:-1])))
print('  sign changes of the Hardy function on a 0.02 grid: %d (zeros %d); %.0f s' % (sc, len(gs), time.time() - t0))
for lab, zz in (('stored', gs), ('recomputed', z)):
    s = np.exp(2j * theta(zz, r['Q'], r['mu'])).sum(); print('  window sum %s %+.6f %+.6fi ; law %+.3f' % (lab, s.real, s.imag, -r['eps'] * r['c'] * (b - a) / (2 * math.pi)))
