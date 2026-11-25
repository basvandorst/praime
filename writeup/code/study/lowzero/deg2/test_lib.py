import sys, time, pickle, os, math
import numpy as np, cypari2
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deg2lib as D
pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6)
H = os.path.dirname(os.path.abspath(__file__))
zf = os.path.join(H, 'zgrid.pkl')
t = time.time()
if os.path.exists(zf):
    Z = pickle.load(open(zf, 'rb'))
else:
    Z = D.ZetaGrid(); pickle.dump(Z, open(zf, 'wb'))
print('zeta grid %.0fs' % (time.time() - t), flush=True)
for ainv in ([1, -1, 0, -17, 29],):
    t = time.time(); c = D.Curve(pari, ainv); print(ainv, 'N', c.N, 'bad', c.bad, 'LS(1)', c.LS1, '%.0fs' % (time.time() - t), flush=True)
    for L in (10.0, 13.0):
        r = np.array([1e-4, 1e-3, 0.01, 0.1, 0.5, 1.0])
        Cv, dg, of = D.C_full(Z, c, L, r=r, parts=True)
        print(' L', L, 'C(r) at', r, ':', ' '.join('%.2e' % v for v in Cv))
        print('    diag part', ' '.join('%.3e' % v for v in dg), ' off part', ' '.join('%.3e' % v for v in of))
    s = D.sigma2_curve(Z, c, 20, 400); g = D.sigma2_gue(); d1 = D.sigma2_deg1(Z, c.N, c.bad, 20, 400)
    k1 = D.sigma2_curve(Z, c, 20, 400, kmode='k1'); oz = D.sigma2_curve(Z, c, 20, 400, kmode='offzeta')
    idx = [int(round(l / D.W)) - 1 for l in (0.5, 1, 2, 4, 6, 8, 10)]
    print(' l      ', ' '.join('%7.1f' % D.LGRID[i] for i in idx))
    for nm, v in (('GUE', g), ('full', s), ('deg1', d1), ('k1', k1), ('offzeta', oz)):
        print(' %-7s' % nm, ' '.join('%7.4f' % v[i] for i in idx))
