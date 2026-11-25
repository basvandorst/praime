"""Independent check of stored genus-2 zeros: recompute with lfuninit at higher precision (realprecision 57) on (a, b],
compare positions; independent count of sign changes of the Hardy function lfunhardy on a fine grid; window phasor sums from
both lists.  usage: r8_zerocheck.py key a b"""
import json, math, sys, time
import numpy as np, cypari2
sys.path.insert(0, '.')
from r1_offsets import theta, load
key, a, b = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
prec = int(sys.argv[4]) if len(sys.argv) > 4 else 57
r = [x for x in load() if x['key'] == key][0]
pari = cypari2.Pari(); pari.allocatemem(900 * 10 ** 6); pari.default('realprecision', prec)
t0 = time.time()
L = pari('lfungenus2([%s,%s])' % (r['spec'][1], r['spec'][2]))
Li = pari.lfuninit(L, [1, int(b) + 3])
z = np.sort(np.array([float(x) for x in pari.lfunzeros(Li, [a, b])]))
g = r['g']; gs = g[(g > a) & (g <= b)]
print('%s window (%g,%g]: stored %d zeros, recomputed (prec %d) %d zeros, %.0fs' % (key, a, b, len(gs), prec, len(z), time.time() - t0))
if len(z) == len(gs):
    print('  max |diff| = %.3e' % np.max(np.abs(z - gs)))
else:
    print('  stored only:', [x for x in gs if np.min(np.abs(z - x)) > 1e-6]); print('  new only:', [x for x in z if np.min(np.abs(gs - x)) > 1e-6])
# sign changes of the Hardy function on a grid (step 0.01) as an independent count
tt = np.arange(a, b + 1e-9, 0.01)
hv = np.array([float(pari.lfunhardy(Li, t)) for t in tt])
sc = int(np.sum(np.sign(hv[1:]) != np.sign(hv[:-1])))
print('  sign changes of lfunhardy on step-0.01 grid: %d ; min |Z| at grid %.3e ; %.0fs' % (sc, np.min(np.abs(hv)), time.time() - t0))
for lab, zz in (('stored', gs), ('recomputed', z)):
    s = np.exp(2j * theta(zz, r['Q'], r['mu'])).sum()
    print('  window sum %s: %+.4f %+.4fi   pred %+.4f' % (lab, s.real, s.imag, -r['eps'] * r['c'] * (b - a) / (2 * math.pi)))
