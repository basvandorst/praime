"""Null distribution of the window over-dispersion statistic (chi2 of 14 window means) under curve-level random signs
(keeps every within-curve correlation between windows, removes any common structure); genus 2 (55) and Sym^3 (75, 8 windows)."""
import math, sys
import numpy as np
sys.path.insert(0, '.')
from r1_offsets import theta, load
rows = load()
def W(S, edges):
    out = np.zeros((len(S), len(edges) - 1), complex)
    for i, r in enumerate(S):
        g = r['g']; ph = np.exp(2j * theta(g, r['Q'], r['mu']))
        for j in range(len(edges) - 1):
            a, b = edges[j], edges[j + 1]; m = (g > a) & (g <= b)
            out[i, j] = ph[m].sum() + r['eps'] * r['c'] * (b - a) / (2 * math.pi)
    return out
rng = np.random.default_rng(7)
for cls, edges in (('g2', np.arange(10, 151, 10)), ('sym3', np.arange(10, 51, 5))):
    S = [r for r in rows if r['cls'] == cls]; X = W(S, edges); n = len(S)
    for part in ('real', 'imag'):
        Y = X.real if part == 'real' else X.imag
        def stat(Z):
            m = Z.mean(0); s = Z.std(0, ddof=1) / math.sqrt(n); return np.sum((m / s) ** 2)
        def stat_c(Z):   # remove a common constant per window length (allows a linear drift) before the chi2
            m = Z.mean(0); s = Z.std(0, ddof=1) / math.sqrt(n); w = 1 / s ** 2; mu = np.sum(w * m) / np.sum(w); return np.sum(((m - mu) / s) ** 2)
        obs = stat(Y); obs_c = stat_c(Y)
        null = np.array([stat(Y * rng.choice([-1, 1], size=(n, 1))) for _ in range(4000)])
        null_c = np.array([stat_c(Y * rng.choice([-1, 1], size=(n, 1))) for _ in range(4000)])
        print('%-5s %s: chi2 obs %.1f, null mean %.1f, 95%% %.1f, 99.9%% %.1f, p = %.4f ;  drift-removed obs %.1f, p = %.4f' % (
            cls, part, obs, null.mean(), np.quantile(null, 0.95), np.quantile(null, 0.999), np.mean(null >= obs), obs_c, np.mean(null_c >= obs_c)))
