"""oned.py: the one-dimensional pieces of the per-prime change of Var log|L'(rho)| (derivation.md sect. 3-4), at fixed L,
bad prime p added to the bad set S = {} (q = 1 arithmetic otherwise), with Gaussian damping exp(-(eps/E0)^2) of the
long-range tails (Abel limit; E0 = 15, 30, 45 raw).
  dmu  = int ln|r| dT2            (closed diagonal part + numerical off-diagonal part)
  dA   = int ln^2|r| dT2          (closed diagonal part + numerical off-diagonal part)
  dB   = -Li2(1/p)/2 + dB_off,  dB_off = -(pi^2/2) int |r| dT2_off   (sum rule int dT2_off = 0 checked)
  mu   = int ln|r| T2 (base S = {}) = mu_GUE + int ln|r| C(r)
The closed diagonal forms are checked against the damped numerical integrals of dT2_diag.
usage: oned.py OUT.json
"""
import sys, math, json
import numpy as np
from mpmath import polylog
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vlib

EG = vlib.EG; TWO_PI = vlib.TWO_PI


def rnodes(Rmax, h=0.125, ng=8):
    x, w = np.polynomial.legendre.leggauss(ng)
    edges = [0.0] + [h * 4.0 ** (-k) for k in range(8, 0, -1)] + list(np.arange(h, Rmax + 1e-9, h))
    edges = np.unique(np.array(edges)); a, b = edges[:-1], edges[1:]
    return ((b - a)[:, None] / 2 * (x[None, :] + 1) + a[:, None]).ravel(), ((b - a)[:, None] / 2 * w[None, :]).ravel()


def closed(p, L):
    lp = math.log(p); l1 = math.log(p / (p - 1))
    dA = 0.0; m = 1
    while p ** m < 1e18:
        dA += -(2 / L) * lp / (m * p ** m) * (EG + math.log(TWO_PI * m * lp / L)); m += 1
    return dict(dmu_diag=lp * l1 / L, dA_diag=dA, dB_diag=-0.5 * float(polylog(2, 1.0 / p)), dV3_Q=lp * l1 ** 2 / L,
                first_order=-0.5 * float(polylog(2, 1.0 / p)) + lp * l1 * (2 + l1) / L)


def main(out):
    E0s = (15.0, 30.0, 45.0); Ls = (9.5, 10.0, 10.5, 11.0, 11.1); Ps = (2, 3, 5, 7, 11, 13)
    res = {'E0': E0s, 'rows': []}
    # mu_GUE check (no damping needed; tail of -sinc^2 ln r beyond R added analytically: -(ln R + 1)/(2 pi^2 R))
    r, w = rnodes(400.0)
    mug = -2 * np.sum(w * np.log(r) * np.sinc(r) ** 2) - 2 * (math.log(400) + 1) / (2 * math.pi ** 2 * 400)
    print('mu_GUE numerical %.6f   closed gamma + ln 2pi - 1 = %.6f' % (mug, EG + math.log(TWO_PI) - 1))
    res['mu_GUE'] = EG + math.log(TWO_PI) - 1; res['mu_GUE_num'] = mug
    for L in Ls:
        Rmax = 150.0 * L / TWO_PI; r, w = rnodes(Rmax); e = TWO_PI * r / L; lr = np.log(r)
        D = [np.exp(-(e / E0) ** 2) for E0 in E0s]
        C = vlib.C_of_r(L, (), r)
        mu_ar = [2 * np.sum(w * lr * C * d) for d in D]
        print('L=%.2f  mu = mu_GUE + arith: arith(E0=15/30/45) %s' % (L, ' '.join('%+.5f' % x for x in mu_ar)), flush=True)
        for p in Ps:
            c = closed(p, L)
            ho = vlib.dT2_off(L, p, r); hd = vlib.dT2_diag(L, p, r)
            row = dict(L=L, p=p, **c, mu_arith=mu_ar)
            row['dmu_off'] = [2 * np.sum(w * lr * ho * d) for d in D]
            row['dA_off'] = [2 * np.sum(w * lr ** 2 * ho * d) for d in D]
            row['dB_off'] = [-(math.pi ** 2 / 2) * 2 * np.sum(w * r * ho * d) for d in D]
            row['sumrule_off'] = [2 * np.sum(w * ho * d) for d in D]
            row['dmu_diag_num'] = [2 * np.sum(w * lr * hd * d) for d in D]
            row['dA_diag_num'] = [2 * np.sum(w * lr ** 2 * hd * d) for d in D]
            row['dB_diag_num'] = [-(math.pi ** 2 / 2) * 2 * np.sum(w * r * hd * d) for d in D]
            res['rows'].append(row)
            f = lambda k: ' '.join('%+.5f' % x for x in row[k])
            print('  p=%2d dmu diag %.5f (num %s) off %s | dA diag %+.5f (num %s) off %s | dB diag %.5f (num %s) off %s | sum rule %s' % (
                p, c['dmu_diag'], f('dmu_diag_num'), f('dmu_off'), c['dA_diag'], f('dA_diag_num'), f('dA_off'), c['dB_diag'],
                f('dB_diag_num'), f('dB_off'), f('sumrule_off')), flush=True)
    json.dump(res, open(out, 'w'), indent=1)


if __name__ == '__main__':
    main(sys.argv[1])
