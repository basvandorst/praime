# Round 170: weight 2, squarefree level, a_p whole numbers (Hasse), plus Steinberg columns.  Usage: integral.py <X,X,..> [center]
import sys, time, numpy as np
from maass import Problem
center = 'center' in sys.argv
for X in [float(v) for v in sys.argv[1].split(',')]:
    t0 = time.time(); P0 = Problem(2.0, 0, X, hol=True, J=3)
    angles = sorted({round(float(np.arccos(k/(2*np.sqrt(p)))), 12) for p in P0.primes for k in range(-int(2*np.sqrt(p)), int(2*np.sqrt(p)) + 1) if abs(k) < 2*np.sqrt(p)})
    th = np.array(angles); P = Problem(2.0, 0, X, hol=True, th=th, bad=True, center=center)
    mask = np.zeros((len(P.primes), len(th) + 2)); mask[:, -2:] = 1
    for i, p in enumerate(P.primes):
        for j, a in enumerate(th):
            k = 2*np.sqrt(p)*np.cos(a)
            if abs(k - round(k)) < 1e-8: mask[i, j] = 1
    lam, t, rd = P.solve(mask=mask)
    out = []
    for i, p in enumerate(P.primes):
        vals = np.concatenate([2*np.sqrt(p)*np.cos(th), [1.0, -1.0]]); w = P.w[i]
        out.append('%d:%s' % (p, '/'.join('%+d%s(%.2f)' % (round(vals[j]), 'S' if j >= len(th) else '', w[j]) for j in np.nonzero(w > 0.01)[0])))
    print('weight 2 X %g%s: smallest level %.5f (t* = %.6f, %d rounds, %.0fs); forms %s; a_p(weight): %s' % (X, ' central zero' if center else '', np.exp(-t), t, rd, time.time() - t0,
          ' '.join('%.1e' % v for v in P.mineigs(P.w)), ' '.join(out)), flush=True)
