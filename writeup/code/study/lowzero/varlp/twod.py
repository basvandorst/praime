"""twod.py: the three-point piece of the per-prime change of Var log|L'(rho)|:
    dV3 = int int ln|u| ln|v| dT3(0, u, v) du dv,   dT3 = dT3_Q + dT3_P   (bad prime p added to S = {}),
dT3_Q in closed form (ln p/L) ln^2(p/(p-1)) (checked here numerically), dT3_P by 2-D quadrature over the half plane
u >= 0 (T3(0,u,v) = T3(0,-u,-v)), Gauss panels of width 1 (NG nodes) graded at u = 0 and v = 0, |u|, |v| <= V, with the
product damping d(u) d(v), d(x) = exp(-(2 pi x/(L E0))^2), and with the first-order far-field model subtracted:
    M(u, v) = -dC(v) D(u) - dC(u) D(v) - dC(u) D(v - u),   D(x) = sin(2 pi x)/(pi x),   dC = dT2_diag,
whose undamped integral is known exactly (derivation.md sect. 4):  2 (gamma + ln 2pi) dmu_diag - dA_diag + J,
J = int ln|u| Ci(2 pi |u|) dC(u) du.   Result(E0) = [int int (dT3_P - M) d d]_grid + M_exact; E0 -> infinity by
Richardson in 1/E0^2.
usage: twod.py L V NG PB OUT.json
"""
import sys, math, json, time
import numpy as np
from scipy.special import sici
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vlib

EG = vlib.EG; TWO_PI = vlib.TWO_PI
PS_ = (2, 3, 5, 7, 11, 13)
SIGF = (5.0, 4.0, 3.3, 2.7)          # damping widths sigma = V/SIGF (unfolded)


def axis(V, ng, half):
    x, w = np.polynomial.legendre.leggauss(ng)
    edges = [0.0] + [4.0 ** (-k) for k in range(7, 0, -1)] + list(np.arange(1.0, V + 1e-9, 1.0))
    edges = np.unique(np.array(edges)); a, b = edges[:-1], edges[1:]
    X = ((b - a)[:, None] / 2 * (x[None, :] + 1) + a[:, None]).ravel(); W = ((b - a)[:, None] / 2 * w[None, :]).ravel()
    if half: return X, W
    return np.concatenate([-X[::-1], X]), np.concatenate([W[::-1], W])


def dQ(L, p, u, v, w):
    """change of the Q-part of T3 when p is added to the bad set: -(2/L^3) Re sum_3 q_p(e1, e2)"""
    k = TWO_PI / L; pos = {'u': u * k, 'v': v * k, 'w': w * k}; tot = 0
    lp = math.log(p)
    for (A, B1n, B2n) in vlib.TERMS:
        al = -pos[A]; e1 = al + pos[B1n]; e2 = al + pos[B2n]
        X = p * np.exp(1j * e1 * lp); Y = p * np.exp(1j * e2 * lp)
        tot = tot + (-(lp ** 3) / ((X - 1) * (Y - 1)))
    return -(2 / L ** 3) * tot.real


def main(L, V, ng, PB, out):
    t0 = time.time()
    E0S = tuple(TWO_PI * V / (L * c) for c in SIGF)
    U, WU = axis(V, ng, True); Vv, WV = axis(V, ng, False)
    PS = vlib.PrimeSums(PB); Slist = [()] + [(p,) for p in PS_]
    dU = [np.exp(-(TWO_PI * U / (L * E0)) ** 2) for E0 in E0S]; dV = [np.exp(-(TWO_PI * Vv / (L * E0)) ** 2) for E0 in E0S]
    lnU = np.log(U); lnV = np.log(np.abs(Vv))
    D = lambda x: 2 * np.sinc(2 * x)
    nE = len(E0S); nP = len(PS_)
    accP = np.zeros((nP, nE)); accM = np.zeros((nP, nE)); accQ = np.zeros((nP, nE)); accPM_V = np.zeros((nP, 3))
    dCv = [vlib.dT2_diag(L, p, Vv) for p in PS_]
    for i, u in enumerate(U):
        z = np.zeros_like(Vv)
        Ve = Vv.copy(); Ve[np.abs(Ve - u) < 1e-9] += 1e-6               # T3 is continuous; avoid the 0/0 of single terms
        args = vlib.b1_args(z, z + u, Ve, L)                      # positions (0, u, v)
        B1all = [PS.B1_all(a, b) for (a, b) in args]
        T = vlib.T3P(z, z + u, Ve, L, Slist, B1all)
        dCu = [vlib.dT2_diag(L, p, np.array([u]))[0] for p in PS_]
        for j, p in enumerate(PS_):
            dP = T[j + 1] - T[0]
            M = -dCv[j] * D(u) - dCu[j] * D(Vv) - dCu[j] * D(Vv - u)
            q = dQ(L, p, z, z + u, Vv)
            base = WU[i] * lnU[i] * WV * lnV
            for k in range(nE):
                wk = 2 * base * dU[k][i] * dV[k]                   # factor 2: the half plane u < 0
                accP[j, k] += np.sum(wk * dP); accM[j, k] += np.sum(wk * M); accQ[j, k] += np.sum(wk * q)
            for k, Vc in enumerate((V / 4, V / 2, V)):            # undamped, truncated at |u|,|v| <= Vc (diagnostic)
                if u <= Vc: accPM_V[j, k] += np.sum(2 * base * (dP - M) * (np.abs(Vv) <= Vc))
        if i % 100 == 0: print('row %d/%d  %.0f s' % (i, len(U), time.time() - t0), flush=True)
    # exact (undamped) integral of the model M
    r, w = axis(4000.0 * 0 + 400.0, 8, True)
    rows = []
    for j, p in enumerate(PS_):
        lp = math.log(p); l1 = math.log(p / (p - 1)); dmu = lp * l1 / L
        dA = sum(-(2 / L) * lp / (m * p ** m) * (EG + math.log(TWO_PI * m * lp / L)) for m in range(1, 80) if p ** m < 1e30)
        Jv = 2 * np.sum(w * np.log(r) * sici(TWO_PI * r)[1] * vlib.dT2_diag(L, p, r))
        Mex = 2 * (EG + math.log(TWO_PI)) * dmu - dA + Jv
        resk = [accP[j, k] - accM[j, k] + Mex for k in range(nE)]
        x = np.array([1 / E0 ** 2 for E0 in E0S]); A = np.vstack([np.ones_like(x), x]).T  # Richardson in 1/E0^2
        ext = float(np.linalg.lstsq(A, np.array(resk), rcond=None)[0][0])
        extP = float(np.linalg.lstsq(A, accP[j], rcond=None)[0][0])
        extQ = float(np.linalg.lstsq(A, accQ[j], rcond=None)[0][0])
        row = dict(L=L, p=p, V=V, ng=ng, PB=PB, E0=list(E0S), dV3P_damped=list(accP[j]), M_damped=list(accM[j]), M_exact=Mex,
                   J=float(Jv), dV3P_sub=resk, dV3P_extrap=ext, dV3P_raw_extrap=extP, dV3Q_num=list(accQ[j]), dV3Q_num_extrap=extQ,
                   dV3Q_closed=lp * l1 ** 2 / L, PminusM_undamped_trunc=list(accPM_V[j]))
        rows.append(row)
        print('L=%.2f p=%2d  dV3_P: subtracted %s -> %.5f | raw damped %s -> %.5f | M exact %.5f | dV3_Q num -> %.5f closed %.5f' % (
            L, p, ' '.join('%.5f' % x for x in resk), ext, ' '.join('%.5f' % x for x in accP[j]), extP, Mex, extQ, lp * l1 ** 2 / L), flush=True)
    json.dump(rows, open(out, 'w'), indent=1)
    print('time %.0f s' % (time.time() - t0))


if __name__ == '__main__':
    main(float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5])
