import sys, os, time, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pred4 as P
H = os.path.dirname(os.path.abspath(__file__))
ue = np.array([0.0] + [(2 * j + 1) / 24 for j in range(24)]); ye = ue**3       # coarse edges in u
fine_idx = np.round(ue * 96).astype(int)
def coarse(x): return np.add.reduceat(x, fine_idx[:-1], axis=-1)[..., :len(ue) - 1]
def jack(num, den):
    n = len(num); N = num.sum(0); Dn = den.sum(0); m = N / Dn; mj = (N - num) / (Dn - den)
    return m, np.sqrt((n - 1) / n * ((mj - mj.mean(0))**2).sum(0))
t0 = time.time(); mods = {T: P.quartic_model(T=T) for T in (750.0, 1500.0)}; print('models %.0f s' % (time.time() - t0))
out = {}
for r2 in (0, 1, 2):
    z = np.load(os.path.join(H, 'data', 'disc_r%d.npz' % r2)); num, den = coarse(z['num']), coarse(z['den'])
    m, s = jack(num, den); mus = P.SIG_MU[r2]
    p15 = P.bin_means(mods[1500.0], mus, ye); p7 = P.bin_means(mods[750.0], mus, ye)
    out[r2] = dict(m=m.tolist(), s=s.tolist(), pred=p15.tolist(), numerr=np.abs(p15 - p7).tolist())
    print('r2=%d: %d fields; u_lo u_hi  y_hi  data +- err   pred (T=1500, |T=750 diff|)' % (r2, len(z['D'])))
    for j in range(len(m)): print('  %.3f %.3f %6.3f  %8.3f +- %6.3f   %8.3f (%.1e)' % (ue[j], ue[j + 1], ye[j + 1], m[j], s[j], p15[j], abs(p15[j] - p7[j])))
    sel = slice(1, None); w = 1 / s[sel]**2
    amp = np.sum(p15[sel] * m[sel] * w) / np.sum(p15[sel]**2 * w); print('  amplitude (bins 2..): %.3f; chi2 pred %.1f, chi2 zero %.1f, %d bins' % (amp, np.sum((m[sel] - p15[sel])**2 * w), np.sum(m[sel]**2 * w), len(w)))
json.dump(out, open(os.path.join(H, 'look1.json'), 'w'))
