"""(B) the Gram-weighted cross statistic (retrospective on NN; proposed as the clean re-test).
X_12 = (1/T) sum_{gamma1, gamma2} e^{2i theta_1(gamma1) - 2i theta_2(gamma2)} K(gamma2 - gamma1),  K(u) = exp(-u^2/2 sigma^2),
heights in [100, 600]. The weight cancels the phase of the swapped term exactly (for every pair, any conductors), so
E X_12 = (c/4 pi^2) kappa_1 kappa_2 int H(u) K(u) du  (+ leakage of d1 d2 ~ K^(L) ~ e^{-sigma^2 L^2/2}, + diagonal
orbits n ~ e^L, both mean zero over pairs). Registered constant c = 1. Compared: H from the recipe (Euler product to
PH), H == 1 (product of the one-point Gram terms), H == H(0)."""
import sys, math, time
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__file__))
from mylib import *
from swapmine import H_and_dH, hp_factors

sig = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
PH = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
t0 = time.time()
Z = load_zeros(); S1 = np.load(OUT + '/s1.npz'); S2 = np.load(OUT + '/s2_pairs.npz'); Ds = list(S1['Ds'])
I, J = S2['I'], S2['J']; Xk = np.load(OUT + '/kron.npy'); P = np.load(OUT + '/primes.npy')
s = sig / math.sqrt(2); dt = 0.02; tg = np.arange(100.0, 600.0 + dt / 2, dt); T = 500.0
G = np.zeros((len(Ds), len(tg)), complex); half = int(6 * s / dt) + 1
for i, D in enumerate(Ds):
    g = Z[D]; b = g[(g > 100 - 6 * s) & (g < 600 + 6 * s)]; ph = np.exp(2j * theta(b, abs(D), 0 if D > 0 else 1))
    c = np.round((b - 100) / dt).astype(int)
    for off in range(-half, half + 1):
        j = c + off; ok = (j >= 0) & (j < len(tg))
        np.add.at(G[i], j[ok], ph[ok] * np.exp(-(tg[j[ok]] - b[ok]) ** 2 / (2 * s * s)))
X = (G @ G.conj().T) * dt / (math.sqrt(math.pi) * s * T)
x = X[I, J]
print('sigma %.2f: X_12 computed (%.0fs); mean Re %.4f, mean Im %.4f, sd Re %.4f' % (sig, time.time() - t0, x.real.mean(), x.imag.mean(), x.real.std()), flush=True)
# predictors
kap = np.array([phi_ratio(abs(D)) for D in Ds]); kk = kap[I] * kap[J]
ug = np.arange(0, 6 * sig + 1e-9, 0.05); wq = np.full(len(ug), 0.05); wq[0] = 0.025; Kw = np.exp(-ug ** 2 / (2 * sig ** 2)) * wq
m = P <= PH; Pm = P[m]; Xm = Xk[:, m].astype(int)
intHK = np.zeros(len(I)); H0 = np.zeros(len(I))
for n, (a, b) in enumerate(zip(I, J)):
    Pf, sg = hp_factors(Xm[a] * Xm[b], Pm)
    Hv, _ = H_and_dH(ug, Pf, sg)
    intHK[n] = 2 * (Hv * Kw).sum().real; H0[n] = Hv[0].real
    if n % 4000 == 0: print(n, '%.0fs' % (time.time() - t0), flush=True)
c0 = 1 / (4 * math.pi ** 2); Kint = sig * math.sqrt(2 * math.pi)
pred_H = c0 * kk * intHK; pred_1 = c0 * kk * Kint; pred_0 = c0 * kk * H0 * Kint
np.savez(OUT + '/s8_gramx_%.2f.npz' % sig, X=X, x=x, pred_H=pred_H, pred_1=pred_1, pred_0=pred_0, H0=H0)
y = x.real
def jfit(cols, names):
    Zm = np.column_stack(cols); c = np.linalg.lstsq(Zm, y, rcond=None)[0]; reps = []
    for f in range(len(Ds)):
        k = ~((I == f) | (J == f)); reps.append(np.linalg.lstsq(Zm[k], y[k], rcond=None)[0])
    reps = np.array(reps); nj = len(reps); se = np.sqrt((nj - 1) / nj * ((reps - reps.mean(0)) ** 2).sum(0))
    r = y - Zm @ c
    print('   ' + ', '.join('%s %.4f +- %.4f' % (nm, v, e) for nm, v, e in zip(names, c, se)) + ';  residual sd %.4f, R2 %.3f' % (r.std(), 1 - r.var() / y.var()))
print('fits of Re X_12 over all %d pairs (jackknife over L-functions):' % len(y))
one = np.ones(len(y))
jfit([one, pred_H], ['icpt', 'c(recipe H)'])
jfit([one, pred_1], ['icpt', 'c(H=1)'])
jfit([one, pred_0], ['icpt', 'c(H=H(0))'])
jfit([one, pred_1, pred_H - pred_1], ['icpt', 'c(H=1)', 'c(H-1)'])
jfit([pred_H], ['c(recipe H), no icpt'])
print('mean Re X / mean pred_H: %.4f; / mean pred_1: %.4f' % (y.mean() / pred_H.mean(), y.mean() / pred_1.mean()))
# by H(0) bins
qs = np.quantile(H0, [0, 0.2, 0.4, 0.6, 0.8, 1])
for lo, hi in zip(qs[:-1], qs[1:]):
    b = (H0 >= lo) & (H0 <= hi)
    print('   H(0) in [%.2f, %.2f]: n %5d, mean Re X %.4f, mean pred_H %.4f, mean pred_1 %.4f' % (lo, hi, b.sum(), y[b].mean(), pred_H[b].mean(), pred_1[b].mean()))
print('done %.0fs' % (time.time() - t0))
