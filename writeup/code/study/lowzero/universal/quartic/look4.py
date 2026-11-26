# residual pulls, X-trend (halves), comparators with the amplitude fitted for every model; cubic check against m_table
import sys, os, json, time; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pred4 as P
H = os.path.dirname(os.path.abspath(__file__))
ue = np.array([0.0] + [(2 * j + 1) / 24 for j in range(24)]); ye = ue**3; fi = np.round(ue * 96).astype(int)
def coarse(x): return np.add.reduceat(x, fi[:-1], axis=-1)[..., :len(ue) - 1]
def jack(num, den):
    n = len(num); N = num.sum(0); Dn = den.sum(0); m = N / Dn; mj = (N - num) / (Dn - den)
    return m, np.sqrt((n - 1) / n * ((mj - mj.mean(0))**2).sum(0))
def fit(m, s, pr, sl):
    w = 1 / s[sl]**2; a = np.sum(pr[sl] * m[sl] * w) / np.sum(pr[sl]**2 * w); return a, 1 / np.sqrt(np.sum(pr[sl]**2 * w)), np.sum((m[sl] - a * pr[sl])**2 * w)
t0 = time.time(); base = P.quartic_model(T=1000.0); cub = P.cubic_model(T=1000.0); print('models %.0f s' % (time.time() - t0))
# cubic check: bins of width 0.05 on (0, 3] against m_table.txt (round 161)
tab = np.loadtxt('/home/bas/praime/writeup/code/study/lowzero/universal/cubic/m_table.txt')
e5 = np.round(np.linspace(0, 3, 61), 10)
for sg, col in ((-1, 2), (1, 3)):
    mine = P.bin_means(cub, P.CUBIC_MU[sg], e5); ref = tab[:, col].reshape(60, 5).mean(1)
    print('cubic D%s0 vs m_table (60 bins of 0.05): max |diff| %.4f, rms %.4f; y^-1/3 coeff %.6f' % ('<' if sg < 0 else '>', np.abs(mine - ref).max(), np.sqrt(np.mean((mine - ref)**2)), P.tail_coeffs(cub, P.CUBIC_MU[sg])['0.3333']))
sl = slice(4, 24); res = {}
for r2 in (1, 2, 0):
    z = np.load(os.path.join(H, 'data', 'disc_r%d.npz' % r2)); num, den, D = coarse(z['num']), coarse(z['den']), np.abs(z['D'])
    m, s = jack(num, den); mus = P.SIG_MU[r2]; sg = -1 if r2 == 1 else 1
    comps = {'pred': P.bin_means(base, mus, ye), 'one_line': P.bin_means(P.single_line_model(P.A1_exact(), T=1000.0), mus, ye, squares=False),
             'deg2_kernel': P.bin_means(base, P.CUBIC_MU[sg], ye), 'cubic_transplant': P.bin_means(cub, P.CUBIC_MU[sg], ye),
             'no_square_term': P.bin_means(base, mus, ye, squares=False)}
    a, ae, c2 = fit(m, s, comps['pred'], sl)
    print('r2=%d: limit amplitude %.3f +- %.3f, shape chi2 %.1f/20; chi2(zero) %.1f' % (r2, a, ae, c2, np.sum((m[sl] / s[sl])**2)))
    print('   pulls (data - a*pred)/err, bins 4..23:', np.round((m[sl] - a * comps['pred'][sl]) / s[sl], 1).tolist())
    for nm, pr in comps.items():
        if nm == 'pred': continue
        b, be, cc = fit(m, s, pr, sl); print('   %-17s fitted amp %6.3f, shape chi2 %7.1f, delta %7.1f' % (nm, b, cc, cc - c2))
    for lo, hi in ((5000, 10000), (10000, 20000)):
        k = (D >= lo) & (D < hi); mm, ss = jack(num[k], den[k]); b, be, cc = fit(mm, ss, comps['pred'], sl)
        print('   half [%d,%d): %d fields, amplitude %.3f +- %.3f, shape chi2 %.1f' % (lo, hi, k.sum(), b, be, cc))
