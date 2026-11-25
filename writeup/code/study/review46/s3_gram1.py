"""(B) one-point Gram term per L-function: G_D = sum_{100 < gamma <= 600} e^{2 i theta_D(gamma)}.
Recipe / stationary phase: <e^{2i theta}(L'/L)(1/2 + it)> = -eps k/phi(k), so the zero density carries
-(1/pi)(k/phi(k)) Re[eps e^{-2i theta(t)}] and G_D = -(k/phi(k)) (600 - 100)/(2 pi) + fluctuation (eps = 1 for real chi).
Also: the same with a smooth taper, and per sign of D, and the local amplitude in 5 height bands."""
import sys, math
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__file__))
from mylib import *

Z = load_zeros(); S = np.load(OUT + '/s1.npz'); Ds = list(S['Ds'])
kap = np.array([phi_ratio(abs(D)) for D in Ds])
G = []; Gt = []; bands = []
edges = np.linspace(100, 600, 6)
for D in Ds:
    g = Z[D]; q = abs(D); a = 0 if D > 0 else 1
    b = g[(g > 100) & (g <= 600)]; e = np.exp(2j * theta(b, q, a))
    G.append(e.sum())
    x = (b - 100) / 500; tap = np.sin(math.pi * x) ** 2          # Hann taper, mean 1/2
    Gt.append((e * tap).sum() / 0.5)
    bands.append([e[(b > lo) & (b <= hi)].sum() for lo, hi in zip(edges[:-1], edges[1:])])
G = np.array(G); Gt = np.array(Gt); bands = np.array(bands)
pred = -kap * 500 / (2 * math.pi)
print('sets %d; kappa range %.3f-%.3f' % (len(Ds), kap.min(), kap.max()))
print('Re G: mean %.2f, mean prediction %.2f; ratio of means %.4f' % (G.real.mean(), pred.mean(), G.real.mean() / pred.mean()))
sl = (G.real * pred).sum() / (pred * pred).sum(); res = G.real - sl * pred
se = math.sqrt((res ** 2).sum() / (len(G) - 1) / (pred * pred).sum())
print('slope of Re G on prediction (through 0): %.4f +- %.4f;  Im G mean %.2f (sd %.2f)' % (sl, se, G.imag.mean(), G.imag.std()))
slt = (Gt.real * pred).sum() / (pred * pred).sum(); rest = Gt.real - slt * pred
print('tapered: slope %.4f +- %.4f' % (slt, math.sqrt((rest ** 2).sum() / (len(G) - 1) / (pred * pred).sum())))
# OLS with intercept
Xd = np.column_stack([np.ones(len(G)), pred]); c = np.linalg.lstsq(Xd, G.real, rcond=None)[0]
print('OLS Re G = %.2f + %.4f * prediction' % (c[0], c[1]))
for sgn in (1, -1):
    m = np.array([np.sign(D) == sgn for D in Ds])
    print('  D %s 0: n=%d slope %.4f' % ('>' if sgn > 0 else '<', m.sum(), (G.real[m] * pred[m]).sum() / (pred[m] ** 2).sum()))
for k in np.unique(np.round(kap, 3)):
    m = np.abs(kap - k) < 5e-4
    if m.sum() >= 5: print('  kappa %.3f: n=%3d  mean Re G / prediction %.3f' % (k, m.sum(), G.real[m].mean() / pred[m].mean()))
pb = -kap[:, None] * 100 / (2 * math.pi) * np.ones((1, 5))
print('per height band (100-200,...,500-600) slope:', ' '.join('%.3f' % ((bands[:, j].real * pb[:, j]).sum() / (pb[:, j] ** 2).sum()) for j in range(5)))
print('residual sd of Re G after the law: %.1f (sqrt N = %.1f)' % (res.std(), math.sqrt(np.mean([((Z[D] > 100) & (Z[D] <= 600)).sum() for D in Ds]))))
np.savez(OUT + '/s3_gram1.npz', G=G, Gt=Gt, kap=kap, bands=bands)
