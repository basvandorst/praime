"""Validation on the USED set r138_Fn (allowed): own X12 and own predictor, against the review's retrospective numbers
(s8_out_0.5.txt / s8_out_1.0.txt: slope 1.0037 +- 0.0067 (sigma 0.5), 1.0090 +- 0.0083 (sigma 1.0), R2 0.911; one-point
slope 1.008 +- 0.007). usage: fn_validate.py P"""
import sys, math, time, os
import numpy as np
from scipy.special import loggamma
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hlib import *

P = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
HERE = os.path.dirname(os.path.abspath(__file__))
ZF = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'data', 'riemann', 'lzeros')) + '/r138_Fn.npz'
t0 = time.time()


def theta(t, q, a):
    t = np.asarray(t, float)
    return t / 2 * math.log(q / math.pi) + np.imag(loggamma(0.25 + a / 2 + 0.5j * t))


z = np.load(ZF); Z = {int(k[1:]): np.sort(z[k]) for k in z.files}
Ds = sorted(Z)
# completeness (round 138 rule)
keep = []
for D in Ds:
    g = Z[D]; q = abs(D); a = 0 if D > 0 else 1; t = np.linspace(100, 600, 200)
    s = np.searchsorted(g, t) - theta(t, q, a) / math.pi
    b = g[(g >= 100) & (g <= 600)]; u = b / (2 * math.pi) * np.log(q * b / (2 * math.pi * math.e))
    if np.abs(s).max() > 1.5 or np.diff(u).max() > 3.0:
        print('flagged', D, round(float(np.abs(s).max()), 2), round(float(np.diff(u).max()), 2))
    else:
        keep.append(D)
Ds = keep; n = len(Ds); I, J = np.triu_indices(n, 1)
print('sets', n, 'pairs', len(I))
kap = np.array([phi_ratio(abs(D)) for D in Ds])

# one-point sums
G1 = np.array([np.exp(2j * theta(Z[D][(Z[D] > 100) & (Z[D] <= 600)], abs(D), 0 if D > 0 else 1)).sum() for D in Ds])
pred1 = -kap * 500 / (2 * math.pi)
sl = (G1.real * pred1).sum() / (pred1 ** 2).sum(); res = G1.real - sl * pred1
print('one-point: slope Re G on -kappa T/2pi %.4f +- %.4f (OLS through 0); Im G mean %.2f sd %.2f' % (
    sl, math.sqrt((res ** 2).sum() / (n - 1) / (pred1 ** 2).sum()), G1.imag.mean(), G1.imag.std()))
jk = []
for f in range(n):
    m = np.arange(n) != f; jk.append((G1.real[m] * pred1[m]).sum() / (pred1[m] ** 2).sum())
jk = np.array(jk); print('   jackknife se %.4f' % math.sqrt((n - 1) / n * ((jk - jk.mean()) ** 2).sum()))

# X12 (review's definition: t-grid [100, 600], dt 0.02, Gaussians of width s = sigma/sqrt 2)
pr, X = real_table(pari_instance(), Ds, P)
Jm, H0 = J_real_all(X, pr, (0.5, 1.0))
print('J computed %.0fs' % (time.time() - t0))
for sig in (0.5, 1.0):
    s = sig / math.sqrt(2); dt = 0.02; tg = np.arange(100.0, 600.0 + dt / 2, dt); T = 500.0
    Gm = np.zeros((n, len(tg)), complex); half = int(6 * s / dt) + 1
    for i, D in enumerate(Ds):
        g = Z[D]; b = g[(g > 100 - 6 * s) & (g < 600 + 6 * s)]
        ph = np.exp(2j * theta(b, abs(D), 0 if D > 0 else 1)); c = np.round((b - 100) / dt).astype(int)
        for off in range(-half, half + 1):
            j = c + off; ok = (j >= 0) & (j < len(tg))
            np.add.at(Gm[i], j[ok], ph[ok] * np.exp(-(tg[j[ok]] - b[ok]) ** 2 / (2 * s * s)))
    Xm = (Gm @ Gm.conj().T) * dt / (math.sqrt(math.pi) * s * T)
    x = Xm[I, J]; y = x.real
    predH = kap[I] * kap[J] * Jm[sig][I, J] / (4 * math.pi ** 2)
    pred1x = kap[I] * kap[J] * sig * math.sqrt(2 * math.pi) / (4 * math.pi ** 2)

    def jfit(cols, names):
        Zm = np.column_stack(cols); c = np.linalg.lstsq(Zm, y, rcond=None)[0]; reps = []
        for f in range(n):
            k = ~((I == f) | (J == f)); reps.append(np.linalg.lstsq(Zm[k], y[k], rcond=None)[0])
        reps = np.array(reps); se = np.sqrt((n - 1) / n * ((reps - reps.mean(0)) ** 2).sum(0)); r = y - Zm @ c
        print('   ' + ', '.join('%s %.4f +- %.4f' % (nm, v, e) for nm, v, e in zip(names, c, se)) + ';  residual sd %.4f, R2 %.3f' % (r.std(), 1 - r.var() / y.var()))
    one = np.ones(len(y))
    print('sigma %.1f (P = %d): mean Re X %.4f, mean Im X %.4f, sd Re %.4f, sd Im %.4f' % (sig, P, y.mean(), x.imag.mean(), y.std(), x.imag.std()))
    jfit([one, predH], ['icpt', 'c(H)'])
    jfit([one, pred1x], ['icpt', 'c(H=1)'])
    jfit([one, pred1x, predH - pred1x], ['icpt', 'c(H=1)', 'c(H-1)'])
    jfit([predH], ['c(H), no icpt'])
    np.savez(os.path.join(HERE, 'fn_val_%.1f_P%d.npz' % (sig, P)), Ds=np.array(Ds), x=x, predH=predH, pred1=pred1x, H0=H0[I, J])
print('done %.0fs' % (time.time() - t0))
