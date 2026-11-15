"""Models with the S-null log-odds as an offset, and the leave-families-out pipeline.

S null (round 5): w_i = min(c_f * lnQ_i/lnR_i, 0.95) for sieve survivors, c_f = O_f / sum_f lnQ/lnR
(fitted per family on the labels at hand; a family with O_f = 0 has w = 0 and adds nothing).
Models: logit p = logit w + b0 + x.beta
  L2: NLL + lam/2 |beta|^2       L1: NLL + alpha |beta|_1      (b0 unpenalised; L-BFGS-B)
  GB: gradient-boosted depth-2 trees on the binary columns, Newton leaves -G/(H + lam_leaf), shrinkage lr.
Score: out-of-family gain in bits = sum_i [y log2(p/w) + (1-y) log2((1-p)/(1-w))] over held-out rows,
and the hit-weighted mean within-family AUC of logit p (model) against logit w (null).
"""
import numpy as np
import scipy.sparse as sp
from scipy.optimize import minimize

L2_GRID = [30.0, 100.0, 300.0, 1000.0]
L1_GRID = [4.0, 6.0, 8.0, 12.0]
GB_ROUNDS = [10, 25, 50, 100]
GB_LR, GB_LAM, GB_MINH = 0.1, 10.0, 2.0
NFOLD = 5
FOLD_SEED = 6


def folds(fam):
    ids = np.unique(fam)
    rng = np.random.default_rng(FOLD_SEED)
    perm = rng.permutation(ids)
    return {int(f): i % NFOLD for i, f in enumerate(perm)}


def s_null(y, wraw, fam):
    w = np.zeros(len(y))
    for f in np.unique(fam):
        m = fam == f
        O = y[m].sum()
        w[m] = np.minimum(O / wraw[m].sum() * wraw[m], 0.95)
    return w


def sigmoid(z):
    return 0.5 * (1 + np.tanh(0.5 * z))


def logit(p):
    return np.log(p) - np.log1p(-p)


def nll_grad(theta, X, y, off):
    b0, beta = theta[0], theta[1:]
    z = off + b0 + X @ beta
    # log(1 + e^z) - y z, stable
    nll = np.sum(np.logaddexp(0, z) - y * z)
    r = sigmoid(z) - y
    return nll, np.concatenate(([r.sum()], X.T @ r))


def fit_l2(X, y, off, lam, theta0=None):
    p = X.shape[1]
    th = np.zeros(p + 1) if theta0 is None else theta0.copy()

    def f(t):
        v, g = nll_grad(t, X, y, off)
        v += 0.5 * lam * t[1:] @ t[1:]
        g[1:] += lam * t[1:]
        return v, g
    res = minimize(f, th, jac=True, method="L-BFGS-B", options=dict(maxiter=500, gtol=1e-6))
    return res.x


def fit_l1(X, y, off, alpha, theta0=None):
    p = X.shape[1]
    # variables: b0 (free), u >= 0, v >= 0; beta = u - v
    if theta0 is None:
        x0 = np.zeros(1 + 2 * p)
    else:
        x0 = np.concatenate(([theta0[0]], np.maximum(theta0[1:], 0), np.maximum(-theta0[1:], 0)))

    def f(x):
        b0, u, v = x[0], x[1:p + 1], x[p + 1:]
        t = np.concatenate(([b0], u - v))
        val, g = nll_grad(t, X, y, off)
        val += alpha * (u.sum() + v.sum())
        return val, np.concatenate(([g[0]], g[1:] + alpha, -g[1:] + alpha))
    bounds = [(None, None)] + [(0, None)] * (2 * p)
    res = minimize(f, x0, jac=True, method="L-BFGS-B", bounds=bounds, options=dict(maxiter=1000, gtol=1e-6))
    x = res.x
    return np.concatenate(([x[0]], x[1:p + 1] - x[p + 1:]))


def fit_gb(Xb, y, off, rounds, Xtest_list=()):
    """Depth-2 boosted trees on binary Xb (CSR). Returns trees and, for each matrix in Xtest_list,
    the cumulative score at every round in GB_ROUNDS (dict round -> array)."""
    n, p = Xb.shape
    XT = Xb.T.tocsr()
    XC = Xb.tocsc()
    F = np.zeros(n)
    trees = []
    out = [dict() for _ in Xtest_list]
    acc = [np.zeros(T.shape[0]) for T in Xtest_list]
    TC = [T.tocsc() for T in Xtest_list]

    def best(G1, H1, G, H):
        G0, H0 = G - G1, H - H1
        ok = (H1 >= GB_MINH) & (H0 >= GB_MINH)
        gain = np.where(ok, G1 ** 2 / (H1 + GB_LAM) + G0 ** 2 / (H0 + GB_LAM) - G ** 2 / (H + GB_LAM), -np.inf)
        j = int(np.argmax(gain))
        return (j, gain[j]) if np.isfinite(gain[j]) and gain[j] > 0 else (None, 0.0)

    def col(M, j):
        return M[:, j].toarray().ravel() > 0.5

    for t in range(1, rounds + 1):
        pr = sigmoid(off + F)
        g = pr - y
        h = pr * (1 - pr)
        GH = XT @ np.column_stack([g, h])
        j1, _ = best(GH[:, 0], GH[:, 1], g.sum(), h.sum())
        if j1 is None:
            leafv = -g.sum() / (h.sum() + GB_LAM) * GB_LR
            tree = (None, None, None, leafv)
            F += leafv
            for i in range(len(acc)):
                acc[i] += leafv
        else:
            m1 = col(XC, j1)
            V = np.column_stack([g * ~m1, h * ~m1, g * m1, h * m1])
            GH2 = XT @ V
            subs = []
            for side, mask in ((0, ~m1), (1, m1)):
                Gs, Hs = g[mask].sum(), h[mask].sum()
                j2, _ = best(GH2[:, 2 * side], GH2[:, 2 * side + 1], Gs, Hs)
                if j2 is None:
                    subs.append((None, -Gs / (Hs + GB_LAM) * GB_LR, None))
                else:
                    m2 = col(XC, j2)
                    a = mask & m2; b_ = mask & ~m2
                    va = -g[a].sum() / (h[a].sum() + GB_LAM) * GB_LR
                    vb = -g[b_].sum() / (h[b_].sum() + GB_LAM) * GB_LR
                    subs.append((j2, vb, va))      # (feature, value if 0, value if 1)
            tree = (j1, subs[0], subs[1], None)

            def apply(M, Mc):
                mm = col(Mc, j1)
                v = np.zeros(M.shape[0])
                for side, mask in ((0, ~mm), (1, mm)):
                    j2, v0, v1 = subs[side]
                    if j2 is None:
                        v[mask] = v0
                    else:
                        m2 = col(Mc, j2)
                        v[mask & ~m2] = v0
                        v[mask & m2] = v1
                return v
            F += apply(Xb, XC)
            for i, (T, Tc) in enumerate(zip(Xtest_list, TC)):
                acc[i] += apply(T, Tc)
        trees.append(tree)
        if t in GB_ROUNDS:
            for i in range(len(acc)):
                out[i][t] = acc[i].copy()
    return trees, out


def gb_apply(trees, Xb, upto):
    Xc = Xb.tocsc()
    v = np.zeros(Xb.shape[0])
    for tree in trees[:upto]:
        j1, s0, s1, leaf = tree
        if j1 is None:
            v += leaf
            continue
        mm = Xc[:, j1].toarray().ravel() > 0.5
        for mask, (j2, v0, v1) in ((~mm, s0), (mm, s1)):
            if j2 is None:
                v[mask] += v0
            else:
                m2 = Xc[:, j2].toarray().ravel() > 0.5
                v[mask & ~m2] += v0
                v[mask & m2] += v1
    return v


def gain_bits(y, p, w):
    return float(np.sum(y * (np.log2(p) - np.log2(w)) + (1 - y) * (np.log2(1 - p) - np.log2(1 - w))))


def auc(score, y):
    pos, neg = score[y == 1], score[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    allv = np.concatenate([pos, neg])
    order = np.argsort(allv, kind="mergesort")
    ranks = np.empty(len(allv)); ranks[order] = np.arange(1, len(allv) + 1)
    # average ranks for ties
    s = allv[order]
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (i + j + 2) / 2
        i = j + 1
    return (ranks[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def fam_auc(score, y, fam):
    num = den = 0.0
    for f in np.unique(fam):
        m = fam == f
        a = auc(score[m], y[m])
        if np.isfinite(a):
            num += a * y[m].sum(); den += y[m].sum()
    return num / den


def pipeline(X, Xb, y, wraw, fam, fold_of, return_preds=False):
    """Leave-families-out CV of every configuration. Returns dict config -> (gain bits, dAUC)."""
    w = s_null(y, wraw, fam)
    act = w > 0
    off = np.full(len(y), -np.inf); off[act] = logit(w[act])
    fo = np.array([fold_of[int(f)] for f in fam])
    configs = [f"L2 lam={l:g}" for l in L2_GRID] + [f"L1 alpha={a:g}" for a in L1_GRID] + \
              [f"GB rounds={r}" for r in GB_ROUNDS]
    pred = {c: np.zeros(len(y)) for c in configs}
    for k in range(NFOLD):
        tr = act & (fo != k); te = act & (fo == k)
        Xtr, Xte = X[tr], X[te]
        th = None
        for l in sorted(L2_GRID, reverse=True):
            th = fit_l2(Xtr, y[tr], off[tr], l, th)
            pred[f"L2 lam={l:g}"][te] = th[0] + Xte @ th[1:]
        th = None
        for a in sorted(L1_GRID, reverse=True):
            th = fit_l1(Xtr, y[tr], off[tr], a, th)
            pred[f"L1 alpha={a:g}"][te] = th[0] + Xte @ th[1:]
        _, outs = fit_gb(Xb[tr], y[tr], off[tr], max(GB_ROUNDS), [Xb[te]])
        for r in GB_ROUNDS:
            pred[f"GB rounds={r}"][te] = outs[0][r]
    res = {}
    base_auc = fam_auc(off[act], y[act], fam[act])
    for c in configs:
        p = sigmoid(off[act] + pred[c][act])
        res[c] = (gain_bits(y[act], p, w[act]), fam_auc(off[act] + pred[c][act], y[act], fam[act]) - base_auc)
    res["_base_auc"] = (0.0, base_auc)
    if return_preds:
        return res, pred, w
    return res
