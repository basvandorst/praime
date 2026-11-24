"""RETROSPECTIVE. Where does gamma < 1 in test B come from?
 (1) gamma - 1 is linear in the residuals e_i = S_i - full_i: gamma - 1 = sum_i c_i . e_i. Share of each L-function.
 (2) Completeness of every zero set (B and the rounds 117/120 sample): count against the exact smooth count
     N0(t) = theta_chi(t)/pi + 1 (true gamma factor, parity a) averaged over the band -> number of missing (or extra) zeros;
     smallest unfolded gaps; and gamma after dropping suspicious sets.
 (3) Null distribution of gamma - pooled (full law true, residual rows resampled), B.
usage: influence.py -> influence_out.txt"""
import os, sys, json, math
import numpy as np
from scipy.special import loggamma
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'bk'))
from common import *

Ds, S, F, Dg = load_B(); X = design(Ds); n = len(Ds)
A = np.linalg.pinv(X.T @ X) @ X.T
off = (A @ F)[1:] - (A @ Dg)[1:]
den = (off ** 2).sum()
# gamma - 1 = sum over cells off * (A e)[1:] / den = sum_i e_i . (A[1:, i]^T off) / den
E = S - F
contrib = np.array([(A[1:, i][:, None] * off * E[i][None, :]).sum() / den for i in range(n)])
g = gamma(S, F, Dg, X)
print('[1] B: gamma %.4f; sum of contributions %+.4f (gamma - 1 = %+.4f)' % (g, contrib.sum(), g - 1))
o = np.argsort(contrib)
print('    most negative contributions: ' + ' '.join('%d:%+.4f' % (Ds[i], contrib[i]) for i in o[:8]))
print('    most positive contributions: ' + ' '.join('%d:%+.4f' % (Ds[i], contrib[i]) for i in o[-8:]))
cs = np.sort(contrib)
print('    gamma - 1 from the 10 most negative: %+.4f; from all others: %+.4f; sd of a contribution %.4f, expected sd of the sum %.4f' % (
    cs[:10].sum(), cs[10:].sum(), contrib.std(), contrib.std() * math.sqrt(n)))

def theta_count(t, k, a):
    s = 0.5 + 1j * np.asarray(t, float)
    return ((np.asarray(t) / 2) * np.log(k / np.pi) + np.imag(loggamma((s + a) / 2))) / np.pi + 1.0

def completeness(g, D, lo=100, hi=600):
    g = np.sort(g); a = 0 if D > 0 else 1; k = abs(D)
    ts = np.linspace(lo, hi, 2001)
    cnt = np.searchsorted(g, ts, side='right')
    dev = cnt - theta_count(ts, k, a)          # S(t) + corrections; mean of S(t) over a long range is ~ 0
    gg = g[(g >= lo) & (g <= hi)]; u = gg / (2 * math.pi) * np.log(k * gg / (2 * math.pi * math.e))
    return float(dev.mean()), float(dev[-200:].mean() - dev[:200].mean()), float(np.diff(u).min())

zB = np.load(os.path.join(LZ, 'r134_B.npz'))
cB = np.array([completeness(zB['D%d' % D], D) for D in Ds])
Z2 = load_r120_zeros(); D2 = sorted([int(k[1:]) for k in Z2], key=lambda D: (abs(D), D))
c2 = np.array([completeness(Z2['D%d' % D], D) for D in D2])
for nm, c in (('B', cB), ('rounds 117/120', c2)):
    print('[2] %s: mean of N - N0 over 100-600: mean %+.3f, sd %.3f, min %+.3f, max %+.3f; |mean| > 0.75 in %d sets; drift (end - start) sd %.3f, max |.| %.2f; smallest unfolded gap %.4f' % (
        nm, c[:, 0].mean(), c[:, 0].std(), c[:, 0].min(), c[:, 0].max(), int((np.abs(c[:, 0]) > 0.75).sum()), c[:, 1].std(), np.abs(c[:, 1]).max(), c[:, 2].min()))
print('    B: correlation of contribution with |mean(N - N0)|: %.3f, with drift %.3f, with the smallest gap %.3f' % (
    np.corrcoef(contrib, np.abs(cB[:, 0]))[0, 1], np.corrcoef(contrib, cB[:, 1])[0, 1], np.corrcoef(contrib, cB[:, 2])[0, 1]))
for thr in (0.75, 0.6, 0.5):
    keep = np.nonzero(np.abs(cB[:, 0]) <= thr)[0]
    print('    B without sets with |mean(N - N0)| > %.2f (%d dropped): gamma %.4f' % (thr, n - len(keep), gamma(S[keep], F[keep], Dg[keep], X[keep])))
keep = np.nonzero(np.abs(cB[:, 1]) <= 0.8)[0]
print('    B without sets with |drift| > 0.8 (%d dropped): gamma %.4f' % (n - len(keep), gamma(S[keep], F[keep], Dg[keep], X[keep])))
# robust versions: drop the k most influential sets of either sign symmetrically
for k in (5, 10):
    keep = np.setdiff1d(np.arange(n), np.concatenate([o[:k], o[-k:]]))
    print('    B dropping the %d most negative and %d most positive contributors: gamma %.4f' % (k, k, gamma(S[keep], F[keep], Dg[keep], X[keep])))

# (3) null distribution of gamma and of gamma - pooled
rng = np.random.default_rng(3)
res = (S - F) - X @ (A @ (S - F))
gn, dn = [], []
for _ in range(1000):
    Sn = F + res[rng.integers(0, n, n)]
    a = gamma(Sn, F, Dg, X); b = pooled(Sn, F, Dg); gn.append(a); dn.append(a - b)
gn = np.array(gn); dn = np.array(dn)
d_obs = g - pooled(S, F, Dg)
print('[3] B null (full law true): gamma %.4f +- %.4f; gamma - pooled %.4f +- %.4f; observed gamma - pooled %+.4f (z %.2f); corr(gamma, pooled) under null %.2f' % (
    gn.mean(), gn.std(), dn.mean(), dn.std(), d_obs, (d_obs - dn.mean()) / dn.std(), np.corrcoef(gn, gn - dn)[0, 1]))
np.savez(os.path.join(HERE, 'influence.npz'), contrib=contrib, cB=cB, c2=c2, D=np.array(Ds), D2=np.array(D2))
