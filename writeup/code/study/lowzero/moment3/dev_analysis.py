"""dev_analysis.py: development set (8 characters, 1000 <= |D| <= 1300; these count as used).
Measures the integrals and compares with the recipe and the comparators; estimates the noise.
usage: dev_analysis.py VALDIR
"""
import sys, os, math, itertools, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import tmlib as T

DEV = [1001, -1003, 1005, -1048, 1068, -1195, 1253, -1267]
T0, H, NPT = 100.0, 0.05, 18001
t = T0 + H * np.arange(NPT)


def w_taper(tt, ramp=100.0):
    tt = np.asarray(tt, dtype=float)
    a = np.clip((tt - 100.0) / ramp, 0, 1); b = np.clip((1000.0 - tt) / ramp, 0, 1)
    return np.sin(0.5 * math.pi * a) ** 2 * np.sin(0.5 * math.pi * b) ** 2


def trap(y, h=H):
    return h * (y.sum() - 0.5 * (y[0] + y[-1]))


def main(valdir):
    V = {D: np.load(os.path.join(valdir, 'absL_%d.npy' % D)) ** 2 for D in DEV + [1]}
    wt = w_taper(t)
    out = {'single': [], 'zetaL': [], 'pair': [], 'triple': []}
    print('quadrature check: trapezoid h = 0.05 vs h = 0.1 (every other point), triple (1001, -1003, 1005):')
    y = V[1001] * V[-1003] * V[1005]
    a1 = trap(y); a2 = 2 * H * (y[::2].sum() - 0.5 * (y[0] + y[-1]))
    print('   %.10e %.10e rel diff %.1e' % (a1, a2, (a1 - a2) / a1))
    yw = y * wt; a1 = trap(yw); a2 = 2 * H * (yw[::2].sum() - 0.5 * (yw[0] + yw[-1]))
    print('   tapered: rel diff %.1e' % ((a1 - a2) / a1))

    print('\nSINGLE L (theorem: (phi(q)/q)(ln(qt/2pi) + 2 gamma + 2 sum_{p|q} ln p/(p-1)))')
    for D in DEV:
        tr = T.Triple((D,), pmax=10**5)
        m = trap(V[D]); p = T.integrate(tr.P, 100, 1000)
        mw = trap(V[D] * wt); pw = T.integrate_w(lambda x: tr.P(x) * w_taper(x), 100, 1000)
        out['single'].append({'D': D, 'meas': m, 'pred': p, 'meas_w': mw, 'pred_w': pw})
        print('  %6d  ratio sharp %.5f  tapered %.5f' % (D, m / p, mw / pw))
    print('\nZETA * L(chi_D) (Topacogullari Thm 1.5 polynomial = recipe, checked to 2e-12)')
    for D in DEV:
        tr = T.Triple((1, D), pmax=10**5)
        y = V[1] * V[D]
        m = trap(y); p = T.integrate(tr.P, 100, 1000)
        mw = trap(y * wt); pw = T.integrate_w(lambda x: tr.P(x) * w_taper(x), 100, 1000)
        out['zetaL'].append({'D': D, 'meas': m, 'pred': p, 'meas_w': mw, 'pred_w': pw, 'diag': T.integrate(tr.diag, 100, 1000)})
        print('  %6d  ratio sharp %.4f  tapered %.4f   (AFE-diagonal/recipe %.4f)' % (D, m / p, mw / pw, out['zetaL'][-1]['diag'] / p))
    print('\nPAIRS (Topacogullari Thm 1.3 polynomial = recipe)')
    for Dp in itertools.combinations(DEV, 2):
        tr = T.Triple(Dp, pmax=10**5)
        y = V[Dp[0]] * V[Dp[1]]
        m = trap(y); p = T.integrate(tr.P, 100, 1000)
        mw = trap(y * wt); pw = T.integrate_w(lambda x: tr.P(x) * w_taper(x), 100, 1000)
        d = T.integrate(tr.diag, 100, 1000); ld = T.integrate(tr.lead, 100, 1000)
        win = [trap(y[k * 2000:(k + 1) * 2000 + 1]) / T.integrate(tr.P, 100 + 100 * k, 200 + 100 * k) for k in range(9)]
        out['pair'].append({'Ds': Dp, 'meas': m, 'pred': p, 'meas_w': mw, 'pred_w': pw, 'diag': d, 'lead': ld, 'win': win})
        print('  %6d %6d  ratio sharp %.4f  tapered %.4f   diag/full %.4f  lead/full %.4f  h0 %.4f' % (Dp + (m / p, mw / pw, d / p, ld / p, tr.hm[0])))
    print('\nTRIPLES')
    for Dt in itertools.combinations(DEV, 3):
        F = T.Family(Dt)
        y = V[Dt[0]] * V[Dt[1]] * V[Dt[2]]
        m = trap(y); mw = trap(y * wt)
        pr = F.predictions(); prw = F.predictions(w=w_taper)
        win = [trap(y[k * 2000:(k + 1) * 2000 + 1]) / T.integrate(F.full.P, 100 + 100 * k, 200 + 100 * k) for k in range(9)]
        out['triple'].append({'Ds': Dt, 'meas': m, 'meas_w': mw, 'pred': pr, 'pred_w': prw, 'win': win, 'h0': float(F.full.hm[0]),
                              'A3_over_Aij': float(F.full.hm[0] / F.hm_no3[0])})
        print('  %6d %6d %6d  r_full %.3f (tap %.3f)  r_diag %.3f  r_kirk %.3f  r_no3 %.3f  r_ind %.3f  r_lead %.3f' % (
            Dt + (m / pr['full'], mw / prw['full'], m / pr['diag'], m / pr['kirk'], m / pr['no3'], m / pr['ind'], m / pr['lead'])), flush=True)
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dev_out.json'), 'w'), indent=0)


if __name__ == '__main__':
    main(sys.argv[1])
