"""The probes: one statistic family per research line. Each takes the member mask M of the window [L, L+W)
and returns {name: value}. theory_values(L, W) returns the predicted value where a formula exists.

Lines: N numerical, G geometric, S statistical, A algebraic, D dynamic, I information-theoretic.
(The computational line is the engine itself: nulls, holdouts and the second implementation.)
"""
import math
import numpy as np
import theory as T

GAPS = [2, 4, 6, 8, 10, 12, 18, 30]
VAR_H = [50, 500, 5000]
SAME_Q = [3, 4, 5, 7, 8]
MI_LAGS = [1, 2, 3, 5, 10]


# ---------- geometric: the quadratics 4n^2 + bn + c that Ulam's spiral draws as diagonals ----------
def curves():
    out = []
    for b in (-4, -2, 0, 2, 4):
        for c in range(1, 81, 2):
            D = b * b - 16 * c
            if D >= 0 and math.isqrt(D) ** 2 == D:
                continue                       # reducible
            out.append((4, b, c))
    return out


_BH = {}


def bh(f):
    if f not in _BH:
        _BH[f] = T.bh_constant(*f)
    return _BH[f]


def _curve_ns(f, L, W):
    a, b, c = f
    # n >= 1 with L <= f(n) < L+W (f increasing for n >= 1 since b >= -4)
    lo = max(1, int(math.isqrt(max(L // a, 0))) - 2)
    n = np.arange(lo, int(math.isqrt((L + W) // a)) + 3, dtype=np.int64)
    v = a * n * n + b * n + c
    keep = (v >= L) & (v < L + W)
    return n[keep], v[keep]


def geometric(M, L, W):
    obs, exp_bh, exp_flat = [], [], []
    for f in curves():
        n, v = _curve_ns(f, L, W)
        if len(v) == 0:
            continue
        o = int(M[v - L].sum())
        e = float(np.sum(1.0 / np.log(v.astype(np.float64))))
        obs.append(o); exp_flat.append(e); exp_bh.append(bh(f) * e)
    obs, eb, ef = np.array(obs, float), np.array(exp_bh), np.array(exp_flat)
    chi2 = float(np.mean((obs - eb) ** 2 / eb))
    slope = float(np.sum(obs * eb) / np.sum(eb * eb))          # least squares obs = s * BH
    corr = float(np.corrcoef(obs / ef, eb / ef)[0, 1])
    return {'G.chi2_BH': chi2, 'G.slope_BH': slope, 'G.corr_BH': corr}


# ---------- numerical, statistical, algebraic, dynamic, information ----------
def _mi(x, y, kx, ky):
    """plug-in mutual information in bits with the Miller-Madow correction."""
    n = len(x)
    J = np.bincount(x * ky + y, minlength=kx * ky).reshape(kx, ky).astype(float)
    px, py = J.sum(1), J.sum(0)
    nz = J > 0
    mi = np.sum(J[nz] * np.log(J[nz] * n / np.outer(px, py)[nz])) / n
    dof = (np.count_nonzero(px) - 1) * (np.count_nonzero(py) - 1)
    return float((mi - dof / (2 * n)) / math.log(2))


def _markov_gain(s, k, K):
    """bits/symbol: log-loss of order-(k-1) minus order-k Markov model, fitted on the first half (add-1/2),
    scored on the second half."""
    n = len(s)
    h = n // 2

    def ctx(order):
        c = np.zeros(n, dtype=np.int64)
        for j in range(1, order + 1):
            c[order:] = c[order:] * K + s[order - j:n - j]
        return c

    def loss(order):
        c = ctx(order)
        idx = np.arange(order, n)
        tr, te = idx[idx < h], idx[idx >= h]
        nc = K ** order
        tab = np.bincount(c[tr] * K + s[tr], minlength=nc * K).reshape(nc, K) + 0.5
        tab /= tab.sum(1, keepdims=True)
        te = te[te >= max(order, 3)]
        return float(-np.mean(np.log2(tab[c[te], s[te]])))

    return loss(k - 1) - loss(k)


def sequential(M, L, W):
    P = np.nonzero(M)[0].astype(np.int64) + L
    g = np.diff(P)
    out = {}
    # numerical
    out['N.eqgap'] = float(np.mean(g[1:] == g[:-1]))
    out['N.updown'] = float(np.mean(g[1:] > g[:-1]) - np.mean(g[1:] < g[:-1]))
    # statistical: gap shares
    for h in GAPS:
        out['S.gap%d' % h] = float(np.mean(g == h))
    # statistical: counts in short intervals, detrended by li (expected = sum 1/ln n)
    for H in VAR_H:
        nb = W // H
        cnt = np.add.reduceat(M[:nb * H].astype(np.int64), np.arange(0, nb * H, H))
        mid = L + np.arange(nb, dtype=np.float64) * H + H / 2
        e = H / np.log(mid)
        out['S.var%d' % H] = float(np.var(cnt - e) / np.mean(e))
    # algebraic: next prime in the same class mod q
    for q in SAME_Q:
        r = P % q
        out['A.same%d' % q] = float(np.mean(r[1:] == r[:-1]))
    # algebraic: Chebyshev race, increments from L, regression on sqrt(x)/ln x - sqrt(L)/ln L
    xs = np.unique(np.geomspace(L + W / 1000, L + W - 1, 400).astype(np.int64))
    pos = np.searchsorted(P, xs, side='right')
    E = np.sqrt(xs) / np.log(xs) - math.sqrt(L) / math.log(L)
    for q, N, R in ((4, (3,), 1), (3, (2,), 1), (8, (3, 5, 7), 1)):
        cR = np.concatenate(([0], np.cumsum(P % q == R)))[pos]
        d = np.zeros(len(xs))
        for a in N:
            cN = np.concatenate(([0], np.cumsum(P % q == a)))[pos]
            d += (cN - cR) / len(N)
        out['A.race%d' % q] = float(np.sum(d * E) / np.sum(E * E))
    # dynamic: Markov predictability of the gap class (g/2 mod 3) and of the capped gap
    r3 = ((g // 2) % 3).astype(np.int64)
    gc = np.minimum(g // 2, 20).astype(np.int64)
    out['D.gain1_mod6'] = _markov_gain(r3, 1, 3)
    out['D.gain2_mod6'] = _markov_gain(r3, 2, 3)
    out['D.gain1_gap'] = _markov_gain(gc, 1, 21)
    out['D.gain2_gap'] = _markov_gain(gc, 2, 21)
    # information: mutual information between gap classes k steps apart
    for k in MI_LAGS:
        out['I.mi6_lag%d' % k] = _mi(r3[:-k], r3[k:], 3, 3)
    out['I.migap_lag1'] = _mi(gc[:-1], gc[1:], 21, 21)
    out['I.migap_lag2'] = _mi(gc[:-2], gc[2:], 21, 21)
    return out


def all_stats(M, L, W):
    d = sequential(M, L, W)
    d.update(geometric(M, L, W))
    return d


def theory_values(L, W):
    t = {}
    t['N.eqgap'] = T.eqgap_share(L, W)
    t['N.updown'] = 0.0
    for h in GAPS:
        t['S.gap%d' % h] = T.gap_share(h, L, W)
    for H in VAR_H:
        t['S.var%d' % H] = T.ms_ratio(H, L, W)
    for q in SAME_Q:
        t['A.same%d' % q] = T.same_class_share(q, L, W)
    for q in (4, 3, 8):
        t['A.race%d' % q] = 1.0
    t['G.chi2_BH'] = 1.0
    t['G.slope_BH'] = 1.0
    return t
