# shape fits on discovery: one amplitude; oscillating part and y^(-1/3) tail with separate amplitudes; full jackknife covariance
import sys, os, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pred4 as P
H = os.path.dirname(os.path.abspath(__file__))
ue = np.array([0.0] + [(2 * j + 1) / 24 for j in range(24)]); ye = ue**3; fi = np.round(ue * 96).astype(int)
def coarse(x): return np.add.reduceat(x, fi[:-1], axis=-1)[..., :len(ue) - 1]
base = P.quartic_model(T=1000.0); sl = slice(4, 24)
for r2 in (1, 2, 0):
    z = np.load(os.path.join(H, 'data', 'disc_r%d.npz' % r2)); num, den = coarse(z['num'])[:, sl], coarse(z['den'])[:, sl]
    n = len(num); N = num.sum(0); Dn = den.sum(0); m = N / Dn; mj = (N - num) / (Dn - den); dm = mj - mj.mean(0); C = (n - 1) / n * dm.T @ dm
    Ci = np.linalg.inv(C); mus = P.SIG_MU[r2]
    Pf = P.bin_means(base, mus, ye)[sl]; c3 = P.tail_coeffs(base, mus)['0.3333']
    Pt = c3 * np.diff(1.5 * ye**(2 / 3))[sl] / np.diff(ye)[sl]; Po = Pf - Pt
    def gls(cols):
        A = np.array(cols).T; G = A.T @ Ci @ A; b = np.linalg.solve(G, A.T @ Ci @ m); r = m - A @ b
        return b, np.sqrt(np.diag(np.linalg.inv(G))), float(r @ Ci @ r)
    print('r2=%d (%d fields), 20 bins u in (7/24, 47/24], full covariance:' % (r2, n))
    print('   chi2 zero %.1f; chi2 pred (no fit) %.1f' % (m @ Ci @ m, (m - Pf) @ Ci @ (m - Pf)))
    for name, cols in (('a*pred', [Pf]), ('a*osc + b*tail', [Po, Pt]), ('osc only, a', [Po])):
        b, be, c2 = gls(cols); print('   %-15s coeffs %s +- %s  chi2 %.1f' % (name, np.round(b, 3), np.round(be, 3), c2))
    corr = np.corrcoef(m, Pf)[0, 1]; print('   correlation data-pred over the 20 bins: %.3f' % corr)
