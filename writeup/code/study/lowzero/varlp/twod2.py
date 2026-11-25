"""twod2.py: as twod.py (the three-point piece dV3 of the per-prime change of Var log|L'(rho)|), for many bad sets at once:
single primes (with the first-order far-field model M subtracted, Richardson in 1/E0^2) and pairs (the non-additive
part T(p,q) - T(p) - T(q) + T(), whose far field has no first-order term; damped and extrapolated).
Same quadrature as twod.py (Gauss panels of width 1, NG nodes, graded at u = 0 and v = 0; half plane u >= 0; |u|,|v| <= V;
product damping with widths V/5, V/4, V/3.3, V/2.7). The S-independent parts of the six P-terms are computed once per row.
usage: twod2.py L V NG PB OUT.json
"""
import sys, os, math, json, time
import numpy as np
from scipy.special import sici
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vlib
from twod import axis, dQ, SIGF

EG = vlib.EG; TWO_PI = vlib.TWO_PI
SINGLES = (2, 3, 5, 7, 11, 13, 17, 23, 31, 47, 67, 97, 151, 251, 397, 631, 1009, 1297)
PAIRS = ((2, 3), (2, 5), (2, 7), (2, 11), (2, 13), (3, 5), (3, 7), (3, 11), (5, 7))


def t3p_multi(u, v, w, L, Slist, B1all):
    k = TWO_PI / L; pos = {'u': u * k, 'v': v * k, 'w': w * k}
    pre = []; ib = 0
    for (A, B1n, B2n) in vlib.TERMS:
        al = -pos[A]; e1 = al + pos[B1n]; e2 = al + pos[B2n]; ed = pos[B2n] - pos[B1n]
        for (f1, f2, fd) in ((e1, e2, ed), (e2, e1, -ed)):
            a = np.abs(f1)
            pre.append(dict(f1=f1, f2=f2, fd=fd, a=a, B=B1all[ib], ph=np.exp(-1j * L * f1), g2=vlib.g_reg(f2), gd=vlib.g_reg(fd),
                            Z2m1=vlib.X('Z2m1', a), Phi=vlib.X('Phi', a), Phim1=vlib.X('Phim1', a),
                            P3G=-1 / (1j * f2) + 1 / (1j * fd), pos1=f1 > 0))
            ib += 1
    res = []
    for S in Slist:
        tot = 0
        for d in pre:
            B = d['B']; g2 = d['g2']; gd = d['gd']
            Psi = 1.0; Psim1 = 0.0
            for p in S:
                lp = math.log(p)
                B = B - vlib.b1_p(p, d['f1'], d['f2'])
                g2 = g2 + lp / (p ** (1 + 1j * d['f2']) - 1); gd = gd + lp / (p ** (1 + 1j * d['fd']) - 1)
                fp = vlib.bklib.f_p(p, d['a']); F = 1 + d['a'] ** 2 * fp; Psim1 = Psim1 * F + fp; Psi = Psi * F
            rego = d['Z2m1'] * d['Phi'] * Psi + d['Phim1'] * Psi + Psim1
            Zc = np.where(d['pos1'], np.conj(rego), rego)
            h = B + g2 - gd
            tot = tot - d['ph'] * (Zc * (d['P3G'] + h) + h / d['f1'] ** 2)
        res.append((2 / L ** 3 * tot).real)
    return res


def main(L, V, ng, PB, out):
    t0 = time.time()
    E0S = tuple(TWO_PI * V / (L * c) for c in SIGF)
    U, WU = axis(V, ng, True); Vv, WV = axis(V, ng, False)
    PS = vlib.PrimeSums(PB)
    Slist = [()] + [(p,) for p in SINGLES] + list(PAIRS); idx = {S: i for i, S in enumerate(Slist)}
    dU = [np.exp(-(TWO_PI * U / (L * E0)) ** 2) for E0 in E0S]; dV = [np.exp(-(TWO_PI * Vv / (L * E0)) ** 2) for E0 in E0S]
    lnU = np.log(U); lnV = np.log(np.abs(Vv)); D = lambda x: 2 * np.sinc(2 * x)
    nE = len(E0S); nS = len(SINGLES); nQ = len(PAIRS)
    accP = np.zeros((nS, nE)); accM = np.zeros((nS, nE)); accQ = np.zeros((nS, nE)); accNA = np.zeros((nQ, nE))
    dCv = [vlib.dT2_diag(L, p, Vv) for p in SINGLES]
    for i, u in enumerate(U):
        z = np.zeros_like(Vv); Ve = Vv.copy(); Ve[np.abs(Ve - u) < 1e-9] += 1e-6
        B1all = [PS.B1_all(a, b) for (a, b) in vlib.b1_args(z, z + u, Ve, L)]
        T = t3p_multi(z, z + u, Ve, L, Slist, B1all)
        base = WU[i] * lnU[i] * WV * lnV; wk = [2 * base * dU[k][i] * dV[k] for k in range(nE)]
        for j, p in enumerate(SINGLES):
            dP = T[idx[(p,)]] - T[0]
            dCu = vlib.dT2_diag(L, p, np.array([u]))[0]
            M = -dCv[j] * D(u) - dCu * D(Vv) - dCu * D(Vv - u)
            q = dQ(L, p, z, z + u, Vv)
            for k in range(nE):
                accP[j, k] += np.sum(wk[k] * dP); accM[j, k] += np.sum(wk[k] * M); accQ[j, k] += np.sum(wk[k] * q)
        for j, (p, q_) in enumerate(PAIRS):
            na = T[idx[(p, q_)]] - T[idx[(p,)]] - T[idx[(q_,)]] + T[0]
            for k in range(nE): accNA[j, k] += np.sum(wk[k] * na)
        if i % 100 == 0: print('row %d/%d  %.0f s' % (i, len(U), time.time() - t0), flush=True)
    x = np.array([1 / E0 ** 2 for E0 in E0S]); A = np.vstack([np.ones_like(x), x]).T
    ext = lambda y: float(np.linalg.lstsq(A, np.asarray(y), rcond=None)[0][0])
    r, w = axis(400.0, 8, True)
    res = {'L': L, 'V': V, 'ng': ng, 'PB': PB, 'E0': list(E0S), 'singles': [], 'pairs': []}
    for j, p in enumerate(SINGLES):
        lp = math.log(p); l1 = math.log(p / (p - 1)); dmu = lp * l1 / L
        dA = sum(-(2 / L) * lp / (m * p ** m) * (EG + math.log(TWO_PI * m * lp / L)) for m in range(1, 80) if p ** m < 1e30)
        Jv = 2 * np.sum(w * np.log(r) * sici(TWO_PI * r)[1] * vlib.dT2_diag(L, p, r))
        Mex = 2 * (EG + math.log(TWO_PI)) * dmu - dA + Jv
        sub = [accP[j, k] - accM[j, k] + Mex for k in range(nE)]
        row = dict(p=p, dV3P_sub=sub, dV3P=ext(sub), dV3P_raw=list(accP[j]), dV3P_raw_ext=ext(accP[j]), M_exact=Mex, J=float(Jv),
                   dV3Q_num=ext(accQ[j]), dV3Q_closed=lp * l1 ** 2 / L)
        res['singles'].append(row)
        print('L=%.2f p=%4d dV3_P %s -> %.5f (raw ext %.5f) | M exact %.5f J %+.5f | dV3_Q num %.5f closed %.5f' % (
            L, p, ' '.join('%.5f' % x for x in sub), row['dV3P'], row['dV3P_raw_ext'], Mex, Jv, row['dV3Q_num'], row['dV3Q_closed']), flush=True)
    for j, S in enumerate(PAIRS):
        res['pairs'].append(dict(S=list(S), na=list(accNA[j]), na_ext=ext(accNA[j])))
        print('L=%.2f pair %s non-additive dV3_P %s -> %+.5f' % (L, S, ' '.join('%+.5f' % x for x in accNA[j]), ext(accNA[j])), flush=True)
    json.dump(res, open(out, 'w'), indent=1)
    print('time %.0f s' % (time.time() - t0))


if __name__ == '__main__':
    main(float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5])
