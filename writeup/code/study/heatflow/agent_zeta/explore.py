import numpy as np, sys
from analyze import *

R = np.load(sys.argv[1]); ts = R['t']; X = R['x']; x0 = X[0]
lo, hi = window(x0, 0.2, 0.8); Delta = hi - lo
hw = np.sqrt(1 + (BETA / np.pi) ** 2) * 2 * np.pi / Delta
idx = (x0 >= lo) & (x0 <= hi)
print('window %.1f..%.1f  Delta %.0f  zeros %d  L %.2f..%.2f  main-lobe half width %.2e' % (lo, hi, Delta, idx.sum(), np.log(lo / 2 / np.pi), np.log(hi / 2 / np.pi), hw))
om = np.linspace(np.log(1.6), np.log(205), 6000); xs = np.exp(om)
mint = np.arange(1, 3000)
dmin = np.min(np.abs(om[:, None] - np.log(mint)[None, :]), 1)
ok = dmin > 1.5 * hw
print('noise points kept', ok.sum())
ns = np.arange(2, 201)
for it in [0, 4, 8, 12]:
    t = ts[it]; g = X[it]
    s = np.mean((g - x0)[idx])
    w = kaiser(g, lo, hi)
    Fn = Fsum(g, w, ns.astype(float)); Fx = Fsum(g, w, xs[ok])
    dn = dirlog(t, 400); P = template(t, ns, lo, hi, dn)
    Rn = Fn * np.exp(-1j * s * np.log(ns)) / P
    print('t=%.2f shift s=%.5f (-pi t/8 = %.5f)' % (t, s, -np.pi * t / 8))
    # noise spectrum in bins of omega
    edges = np.linspace(om[0], om[-1], 13)
    b = np.digitize(om[ok], edges)
    print('  noise rms |F| by omega bin:', ' '.join('%.1f' % np.sqrt(np.mean(np.abs(Fx[b == k]) ** 2)) for k in range(1, 13)))
    print('  sqrt(sum w^2)=%.1f' % np.sqrt((w ** 2).sum()))
    for n in [2, 3, 5, 7, 11, 13, 4, 9, 6, 10, 12, 14, 15, 21]:
        print('   n=%3d %-5s P=%9.3f  F=%9.3f%+9.3fi  R=%7.3f%+7.3fi' % (n, classify(n), P[n - 2], Fn[n - 2].real, Fn[n - 2].imag, Rn[n - 2].real, Rn[n - 2].imag))
