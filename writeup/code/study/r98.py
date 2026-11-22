"""Round 98 (preregister/round98.txt): RMS over ln X of the swing of R(X, d), (2/d) sqrt(2Q),
Q = sum_{j != k} |c_j|^2 |c_k|^2/(4 + (g_j - g_k)^2) over both signs, for the actual zeros, Poisson sets and GUE sets
with the same mean count. usage: r98.py"""
import os, json, math, time
import numpy as np
from scipy.linalg import eigh_tridiagonal

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', 'data', 'study')
Z = np.array([float(l) for l in open(os.path.join(HERE, '..', '..', 'data', 'riemann', 'zeros_10000_flint.txt')) if l.strip()])
DS = (0.1, 0.03, 0.01)
t = np.linspace(2 * math.pi, Z[-1] + 3000, 2_400_000)
Nmean = np.maximum.accumulate((t / (2 * math.pi)) * np.log(t / (2 * math.pi * math.e)) + 7 / 8)


def rms_all(g):
    G = np.concatenate([g, -g]); rho = 0.5 + 1j * G; out = {}
    for d in DS:
        c2 = np.abs(((1 + d) ** rho - 1) / rho) ** 2; Q = 0.0
        for i0 in range(0, len(G), 500):
            dg = G[i0:i0 + 500, None] - G[None, :]
            M = c2[i0:i0 + 500, None] * c2[None, :] / (4 + dg * dg)
            ii = np.arange(i0, min(i0 + 500, len(G))); M[ii - i0, ii] = 0.0
            Q += float(M.sum())
        out[str(d)] = 2 / d * math.sqrt(2 * Q)
    return out


def poisson(seed):
    rng = np.random.default_rng(seed); lo, hi = max(Nmean[0], 0.0), float(np.interp(Z[-1], t, Nmean))
    n = rng.poisson(hi - lo)
    return np.sort(np.interp(rng.uniform(lo, hi, n), Nmean, t))


def gue(seed, n=30_000, keep=10_000):
    rng = np.random.default_rng(seed)
    dgl = rng.normal(0, math.sqrt(2), n)
    off = np.sqrt(rng.chisquare(2 * np.arange(n - 1, 0, -1)))          # beta = 2
    ev = np.sort(eigh_tridiagonal(dgl, off, eigvals_only=True)) / math.sqrt(2 * n)   # semicircle on [-sqrt2, sqrt2]... scale below
    R = float(np.max(np.abs(ev))) * (1 + 1e-9)                       # edge of the spectrum (semicircle radius), from the sample
    x = np.clip(ev / R, -1, 1)
    F = 0.5 + (x * np.sqrt(1 - x * x) + np.arcsin(x)) / math.pi       # semicircle CDF
    u = F * n
    mid = (n - keep) // 2
    u = u[mid:mid + keep] - u[mid]                                   # unfolded, mean spacing 1, starting at 0
    lo = float(np.interp(Z[0], t, Nmean))
    return np.interp(lo + u, Nmean, t)                                # mapped through the zeros' mean count


if __name__ == '__main__':
    import sys
    PART = os.path.join(OUT, 'r98_parts.json')
    mode = sys.argv[1]
    if mode in ('actual', 'poisson', 'gue'):
        s0, s1 = (int(sys.argv[2]), int(sys.argv[3])) if mode != 'actual' else (0, 1)
        out = {}
        for s in range(s0, s1):
            g = Z if mode == 'actual' else (poisson(s) if mode == 'poisson' else gue(s))
            out['%s_%d' % (mode, s)] = rms_all(g)
        tmp = PART + '.%s_%d' % (mode, s0)
        json.dump(out, open(tmp, 'w')); print(mode, s0, s1, 'done', flush=True)
    else:
        res = {}
        for f in sorted(os.listdir(OUT)):
            if f.startswith('r98_parts.json.'):
                res.update(json.load(open(os.path.join(OUT, f))))
        summ = {}
        for d in map(str, DS):
            p = np.array([v[d] for k, v in res.items() if k.startswith('poisson')]); g = np.array([v[d] for k, v in res.items() if k.startswith('gue')])
            a = res['actual_0'][d]
            summ[d] = {'actual': a, 'n_poisson': len(p), 'poisson_mean': float(p.mean()), 'poisson_sd': float(p.std(ddof=1)), 'n_gue': len(g),
                       'gue_mean': float(g.mean()), 'gue_sd': float(g.std(ddof=1)),
                       'z_vs_poisson': float((a - p.mean()) / p.std(ddof=1)), 'z_vs_gue': float((a - g.mean()) / g.std(ddof=1))}
            print(d, {k: round(v, 4) for k, v in summ[d].items()}, flush=True)
        json.dump({'sets': res, 'summary': summ, 'note': 'GUE edge R taken from the sample (fixed before the run; the first draft used sqrt 2, wrong scale)'},
                  open(os.path.join(OUT, 'r98.json'), 'w'), indent=1)
