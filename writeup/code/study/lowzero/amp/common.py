"""Shared helpers: BK predictions per L-function (as the repo's predict.py, from the cached grid), the number variance
of a zero set exactly as r134.py, and the amplitude gamma of the bad-prime off-diagonal part (as r134.py / r134b.py)."""
import os, sys, json, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'bk'))
import bklib
from bklib import sigma2_full, dS_diag, LGRID, W, L_nodes, TWO_PI

REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..', '..'))
LZ = os.path.join(REPO, 'writeup/data/riemann/lzeros')
LOW = os.path.join(REPO, 'writeup/code/study/lowzero')
LS = [0.5 * i for i in range(1, 17)]
IDX = [int(round(l / W)) - 1 for l in LS]
PR6 = (2, 3, 5, 7, 11, 13)


def prime_factors(n):
    n = abs(n); ps = []; p = 2
    while p * p <= n:
        if n % p == 0:
            ps.append(p)
            while n % p == 0: n //= p
        p += 1
    if n > 1: ps.append(n)
    return ps


def predict(G, Ds, t1, t2, n=8, off_scale=1.0, off_dL=0.0, extra=None):
    """per-L-function predicted Sigma2 on the 16 primary lengths, averaged over [t1, t2] with the density weight.
    Returns (full, diag). off_scale multiplies the bad-prime off-diagonal change; off_dL evaluates that change at L + off_dL.
    extra(G, L, B) -> array on LGRID added to full (for alternative models)."""
    F, Dg = [], []
    for D in Ds:
        q = abs(D); B = tuple(prime_factors(D)); Ln, wn = L_nodes(q, t1, t2, n=n)
        f = 0; d = 0
        for L, w in zip(Ln, wn):
            base = sigma2_full(G, L, ())
            dg = base + sum(np.array([dS_diag(p, L, l) for l in LGRID]) for p in B)
            if off_dL == 0.0:
                fu = sigma2_full(G, L, B)
                off = fu - dg
            else:
                L2 = L + off_dL
                off = sigma2_full(G, L2, B) - (sigma2_full(G, L2, ()) + sum(np.array([dS_diag(p, L2, l) for l in LGRID]) for p in B))
            fu = dg + off_scale * off
            if extra is not None:
                fu = fu + extra(G, L, B)
            f = f + w * fu; d = d + w * dg
        F.append(f[IDX]); Dg.append(d[IDX])
    return np.array(F), np.array(Dg)


def sigma2_set(g, q, lo, hi):
    """r134.py's statistic: g sorted positive ordinates; windows from the first to the last zero in [lo, hi]."""
    g = np.sort(g); g = g[(g >= lo) & (g <= hi)]
    u = g / (2 * math.pi) * np.log(q * g / (2 * math.pi * math.e)); v = []
    for l in LS:
        x = np.arange(u[0], u[-1] - l, 0.01); c = np.searchsorted(u, x + l) - np.searchsorted(u, x); v.append(np.var(c))
    return np.array(v)


def design(Ds, primes=PR6):
    return np.array([[1] + [int(abs(D) % p == 0) for p in primes] for D in Ds], float)


def gamma(S, F, Dg, X, ix=None, rows=slice(None), cols=slice(None)):
    ix = np.arange(len(S)) if ix is None else ix
    Ai = np.linalg.pinv(X[ix].T @ X[ix]) @ X[ix].T
    fit = (Ai @ S[ix])[1:]; pf = (Ai @ F[ix])[1:]; pd = (Ai @ Dg[ix])[1:]
    off = (pf - pd)[rows, cols]; res = (fit - pd)[rows, cols]
    return float((res * off).sum() / (off ** 2).sum())


def gamma_boot(S, F, Dg, X, nb=1000, seed=7, **kw):
    rng = np.random.default_rng(seed); n = len(S); b = []
    for _ in range(nb):
        try: b.append(gamma(S, F, Dg, X, ix=rng.integers(0, n, n), **kw))
        except Exception: pass
    return gamma(S, F, Dg, X, **kw), float(np.std(b, ddof=1))


def pooled(S, F, Dg, ix=None):
    ix = np.arange(len(S)) if ix is None else ix
    y = (S - Dg)[ix]; x = (F - Dg)[ix]; y = y - y.mean(0); x = x - x.mean(0)
    return float((x * y).sum() / (x * x).sum())


def load_B():
    P = json.load(open(os.path.join(LOW, 'bk', 'predictions_per_lfunction.json')))
    R = json.load(open(os.path.join(LOW, 'r134_B.json')))
    idx = [P['l'].index(l) for l in LS]
    Ds = P['B']['D']; S = np.array(R['sigma2']); F = np.array(P['B']['full'])[:, idx]; Dg = np.array(P['B']['diag'])[:, idx]
    return Ds, S, F, Dg


def load_r120_zeros():
    Z = {}
    for f in ('r117_part0.npz', 'r117_part1.npz', 'r120_part0.npz', 'r120_part1.npz'):
        z = np.load(os.path.join(LZ, f)); Z.update({k: z[k] for k in z.files})
    return Z
