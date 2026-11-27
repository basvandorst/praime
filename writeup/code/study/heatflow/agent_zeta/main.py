import numpy as np, sys
from analyze import *
np.seterr(divide='ignore', invalid='ignore')

run = sys.argv[1]; f0, f1 = float(sys.argv[2]), float(sys.argv[3])
tlist = [float(v) for v in sys.argv[4].split(',')]
R = np.load(run); ts = R['t']; X = R['x']; x0 = X[0]
lo, hi = window(x0, f0, f1); Delta = hi - lo
hw = np.sqrt(1 + (BETA / np.pi) ** 2) * 2 * np.pi / Delta
idx = (x0 >= lo) & (x0 <= hi)
print('%s M=%d window %.1f..%.1f (zeros %d) L %.3f..%.3f hw %.2e' % (run, len(x0), lo, hi, idx.sum(),
      np.log(lo / 2 / np.pi), np.log(hi / 2 / np.pi), hw))
ns = np.arange(2, 201); kind = np.array([classify(n) for n in ns])
om = np.linspace(np.log(1.6), np.log(205), 6000)
dmin = np.min(np.abs(om[:, None] - np.log(np.arange(1, 3000.))[None, :]), 1)
xs = np.exp(om[dmin > 1.5 * hw])
show = [2, 3, 5, 7, 11, 13, 4, 9, 6, 10, 12, 14, 15, 21]
for t in tlist:
    it = int(np.argmin(np.abs(ts - t))); t = ts[it]; g = X[it]
    s = np.mean((g - x0)[idx]); sd = np.std((g - x0)[idx])
    w = kaiser(g, lo, hi)
    Fn = Fsum(g, w, ns.astype(float)); Fx = Fsum(g, w, xs)
    dn = dirlog(t, 400); P = template(t, ns, lo, hi, dn)
    Rn = Fn * np.exp(-1j * s * np.log(ns)) / P
    # noise: rms |F| at non-integer x, in 6 omega bins, interpolated to ln n
    lx = np.log(xs); e = np.linspace(lx.min(), lx.max(), 7); bc = 0.5 * (e[1:] + e[:-1])
    rms = np.array([np.sqrt(np.mean(np.abs(Fx[(lx >= e[k]) & (lx <= e[k + 1])]) ** 2)) for k in range(6)])
    sig = np.interp(np.log(ns), bc, rms) / np.sqrt(2)          # per real component
    print('\nt=%.2f  mean shift s=%.6f (-pi t/8=%.6f, s/t=%.6f) spread of displacement %.4f' % (t, s, -np.pi * t / 8, s / t if t else 0, sd))
    print('  noise rms|F| non-integer x by omega bin:', ' '.join('%.2e' % v for v in rms), ' (|P| range %.2f..%.2f)' % (np.abs(P[P != 0]).min() if t else 0, np.abs(P).max()))
    for n in show:
        i = n - 2
        print('   n=%3d %-5s P=%9.4f  ReR-1=%+.2e ImR=%+.2e  |R|-1=%+.2e  err=%.1e  SNR=%.0f' % (n, kind[i], P[i], Rn[i].real - 1,
              Rn[i].imag, abs(Rn[i]) - 1, sig[i] / abs(P[i]), abs(P[i]) / sig[i]))
    if t == 0:
        c = kind == 'comp'
        print('  t=0 composites: max |F| %.2e, rms %.2e' % (np.abs(Fn[c]).max(), np.sqrt(np.mean(np.abs(Fn[c]) ** 2))))
        continue
    a = (Fn * np.exp(-1j * s * np.log(ns))).real; q = (Fn * np.exp(-1j * s * np.log(ns))).imag
    for k in ['prime', 'ppow', 'comp']:
        m = kind == k
        W = P[m] ** 2 / sig[m] ** 2
        b = np.sum(a[m] * P[m] / sig[m] ** 2) / W.sum(); eb = 1 / np.sqrt(W.sum())
        chi = np.sum((a[m] - b * P[m]) ** 2 / sig[m] ** 2) / (m.sum() - 1)
        r = Rn[m]
        print('  pooled %-5s n<=200 (%3d): slope %.6f +- %.1e  (%.0f err from 1, %.0f from 0), chi2/dof %.1f; '
              'ReR-1 rms %.1e max %.1e; ImR rms %.1e; unweighted mean ReR %.5f'
              % (k, m.sum(), b, eb, abs(b - 1) / eb, b / eb, chi, np.sqrt(np.mean((r.real - 1) ** 2)),
                 np.max(np.abs(r.real - 1)), np.sqrt(np.mean(r.imag ** 2)), r.real.mean()))
    worst = np.argsort(-np.abs(Rn - 1))[:5]
    print('  worst n:', ', '.join('%d(%s) %+.1e%+.1ei' % (ns[j], kind[j], Rn[j].real - 1, Rn[j].imag) for j in worst))
