"""Level predictions (no bad prime: zeta; and the effect of the exact two-point change per prime) by the second-order
deformation formula with the structure function read off the exact Bogomolny-Keating number variance:
   v(l) = (Sigma2_{L,B}(l) - Sigma2_GUE(l)) / (1 - sinc^2 l)      [closure: Delta Sigma2 = v R2_base, base = sine]
   dVar = int v p,   d<r~> = (second-order formula, lib/gapcalc.py).
This v is formal (negative at small l: the zeros are more rigid than GUE)."""
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/lib')
import numpy as np
from scipy.interpolate import CubicSpline
import loadlibs, gapcalc
bk = loadlibs.load_bklib(); tp = loadlibs.load_tplib(); G = tp.G
B0 = gapcalc.Base(); m1, V0, R0, _ = gapcalc.stats_with(B0, None)
LG = bk.LGRID; GUE = bk.sigma2_gue()

def v_from(L, S):
    d = bk.sigma2_full(G, L, tuple(S)) - GUE
    keep = LG >= 0.1; vv = d[keep] / (1 - np.sinc(LG[keep]) ** 2)
    c = np.polyfit(LG[keep][:5] ** 2, vv[:5], 1)            # v ~ a l^2 + b near 0 -> force through 0 below 0.1
    a = vv[0] / LG[keep][0] ** 2
    x = np.concatenate([[0.0, 0.05], LG[keep]]); y = np.concatenate([[0.0, a * 0.05 ** 2], vv])
    cs = CubicSpline(x, y)
    def v(u):
        u = np.abs(np.asarray(u, float)); out = cs(np.minimum(u, LG[-1])); return out
    return v

def level(L, S=()):
    return gapcalc.second_order(B0, v_from(L, S))

if __name__ == '__main__':
    out = {'V0': V0, 'R0': R0, 'rows': []}
    print('zeta (no bad prime): level of Var(s) and <r~> from the exact two-point function (second order, base sine)')
    print('   L     Var(s)     <r~>     dVar      d<r~>    | BBLM 1/L^2 and with Q: Var')
    lam = 1.573151; Q = 2.315846
    for L in (5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.3, 9.4, 10.2, 10.7, 11.0, 11.2, 24.47, 44.58, 46.83):
        dv, dr = level(L)
        bb = V0 - (2 * lam / L ** 2) * (1 + V0); bq = bb + 6 * Q * (1 + V0) / L ** 3
        out['rows'].append(dict(L=L, var=V0 + dv, r=R0 + dr, dvar=dv, dr=dr, bblm=bb, bblmq=bq))
        print('  %5.2f  %.5f  %.5f  %+.5f  %+.5f  | %.5f %.5f' % (L, V0 + dv, R0 + dr, dv, dr, bb, bq))
    print('\nper-prime changes from the exact two-point change (v_B - v_0), second-order formula, base sine (p missing minus present)')
    for L in (9.4, 10.7, 11.0):
        base = level(L)
        for p in (2, 3, 5, 7, 11, 13):
            d = level(L, (p,))
            print('  L=%.1f p=%2d  dVar %+.5f  d<r~> %+.5f' % (L, p, d[0] - base[0], d[1] - base[1]))
            out['rows'].append(dict(L=L, p=p, dvar_exact2pt=d[0] - base[0], dr_exact2pt=d[1] - base[1]))
    json.dump(out, open('level.json', 'w'), indent=1)
