"""FROZEN validation harness for 'exotic arithmetic features' of Mersenne exponents.

Protocol (fixed before any real feature was scored):
  primes 3 <= p < 1e7 (odd primes); null w(p) = e^g log2(a p)/p, a = 6 (p = 1 mod 4), 2 (p = 3 mod 4);
  w = 1 for p <= 7 (deterministic, so 3, 5, 7 carry no information and are left out of the scoring).
  Discovery  D: 11 <= p < 1e5   (24 exponents 13 ... 86243)
  Validation V: 1e5 <= p < 1e7   (10 exponents 110503 ... 6972593)
  Final exam F: 1e7 <= p < B     (13 exponents) -- NOT touched here; see final_exam.py.
For each feature v (float array over the primes, NaN = undefined):
  candidate classes come from the discovery VALUES only (label free): v <= t, v > t for discovery
  quantiles t (5% .. 95%), or every level/threshold when there are <= 12 levels; a class needs
  expected E_d >= 0.5.  The class with the smallest discovery Poisson-binomial p (ties: larger O/E)
  is chosen on discovery labels; it is scored on validation: O_v, E_v, exact Poisson-binomial
  one-sided p_v = P(X >= O_v), and the log-likelihood gain (bits) of the discovery-fitted two-rate
  model (rate (O+.5)/(E+.5) in and out of the class) over the null.
The identical pipeline runs on NS synthetic label sets drawn from the null (labels for p < 1e7).
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GAMMA = 0.5772156649015329
EXPS = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281,
        3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243,
        110503, 132049, 216091, 756839, 859433, 1257787, 1398269, 2976221, 3021377,
        6972593, 13466917, 20996011, 24036583, 25964951, 30402457, 32582657,
        37156667, 42643801, 43112609, 57885161, 74207281, 77232917, 82589933]
B = 83880077
S1, S2 = 10**5, 10**7
NS = 400            # synthetic null label sets
K = 70              # PMF truncation (E <= 25 everywhere; tail beyond 70 < 1e-20)
J = 8               # power sums used for the small-w part
WBIG = 0.02         # discovery primes with w > WBIG go through an exact DP


def odd_primes_below(n):
    s = np.ones(n // 2, dtype=bool)          # s[i] <-> 2i+1
    s[0] = False
    for i in range(1, int(n ** 0.5) // 2 + 1):
        if s[i]:
            q = 2 * i + 1
            s[q * q // 2::q] = False
    return 2 * np.flatnonzero(s).astype(np.int64) + 1


def null_w(p):
    p = np.asarray(p, dtype=np.float64)
    a = np.where(np.mod(p, 4) == 1, 6.0, 2.0)
    w = np.exp(GAMMA) * np.log2(a * p) / p
    w = np.minimum(w, 1.0)
    w[p <= 7] = 1.0
    return w


P = odd_primes_below(S2)                       # 664,578 odd primes < 1e7
W = null_w(P)
DMASK = (P >= 11) & (P < S1)
VMASK = (P >= S1) & (P < S2)
REAL = np.flatnonzero(np.isin(P, EXPS))       # odd exponents < 1e7 (positions in P)


def label_sets(ns=NS, seed=20261115):
    path = os.path.join(HERE, f"labels_{ns}_{seed}.npz")
    if os.path.exists(path):
        z = np.load(path)
        return [REAL] + [z[f"s{i}"] for i in range(ns)]
    rng = np.random.default_rng(seed)
    sets = [np.flatnonzero(rng.random(len(P)) < W) for _ in range(ns)]
    np.savez(path, **{f"s{i}": s for i, s in enumerate(sets)})
    return [REAL] + sets


# ---------------------------------------------------------------- Poisson-binomial machinery
def _exp_series(c):
    """c: (nc, K+1) coefficients of log G; returns coefficients of G = exp(log G)."""
    nc = c.shape[0]
    g = np.zeros((nc, K + 1))
    g[:, 0] = np.exp(c[:, 0])
    kc = c * np.arange(K + 1)[None, :]
    for n in range(1, K + 1):
        g[:, n] = (kc[:, 1:n + 1] * g[:, n - 1::-1][:, :n]).sum(axis=1) / n
    return g


def pmf_from_sums(logprod, S):
    """logprod: (nc,) sum log(1-w); S: (nc, J) power sums of r = w/(1-w)."""
    nc = len(logprod)
    c = np.zeros((nc, K + 1))
    c[:, 0] = logprod
    for j in range(1, J + 1):
        c[:, j] = (-1) ** (j + 1) * S[:, j - 1] / j
    return _exp_series(c)


def in_class(v, cls):
    typ, t = cls
    with np.errstate(invalid="ignore"):
        if typ == 0:
            return v <= t
        if typ == 1:
            return v > t
        return v == t


def candidate_classes(vd):
    vd = vd[~np.isnan(vd)]
    u = np.unique(vd)
    if len(u) <= 1:
        return []
    if len(u) <= 12:
        th = u[:-1]
        cl = [(0, t) for t in th] + [(1, t) for t in th]
        if len(u) > 2:
            cl += [(2, c) for c in u]
    else:
        th = np.unique(np.quantile(vd, np.arange(0.05, 0.951, 0.05)))
        th = th[th < u[-1]]
        cl = [(0, t) for t in th] + [(1, t) for t in th]
    return cl


def class_stats(v, w, classes, scale=1.0):
    """Per class: E = sum w, logprod = sum log(1-w), S_j = sum r^j (j=1..J); all times scale."""
    r = w / (1.0 - w)
    cols = np.empty((len(v), J + 2))
    cols[:, 0] = w
    cols[:, 1] = np.log1p(-w)
    rj = np.ones_like(r)
    for j in range(J):
        rj = rj * r
        cols[:, 2 + j] = rj
    ok = ~np.isnan(v)
    vs, cs = v[ok], cols[ok]
    o = np.argsort(vs, kind="stable")
    vs, cs = vs[o], cs[o]
    pre = np.vstack([np.zeros((1, J + 2)), np.cumsum(cs, axis=0)])
    tot = pre[-1]
    out = np.empty((len(classes), J + 2))
    for i, (typ, t) in enumerate(classes):
        hi = np.searchsorted(vs, t, "right")
        if typ == 0:
            out[i] = pre[hi]
        elif typ == 1:
            out[i] = tot - pre[hi]
        else:
            out[i] = pre[hi] - pre[np.searchsorted(vs, t, "left")]
    return out * scale


def disc_pmf(vD, wD, classes):
    big = wD > WBIG
    st = class_stats(np.where(big, np.nan, vD), wD, classes)
    pmf = pmf_from_sums(st[:, 1], st[:, 2:])
    vb, wb = vD[big], wD[big]
    for x, wi in zip(vb, wb):
        m = np.array([bool(in_class(np.array([x]), c)[0]) for c in classes], dtype=float) * wi
        sh = np.zeros_like(pmf)
        sh[:, 1:] = pmf[:, :-1]
        pmf = pmf * (1 - m)[:, None] + sh * m[:, None]
    E = st[:, 0] + np.array([wb[in_class(vb, c)].sum() for c in classes])
    return pmf, E


def survival(pmf):
    sf = np.cumsum(pmf[:, ::-1], axis=1)[:, ::-1]
    return np.minimum(sf, 1.0)


def evaluate(v, sets=None, sample_mask=None, sample_frac=1.0):
    """Score one feature on every label set. v: float array over P (NaN = undefined).
    sample_mask: if given, validation expectations come from these validation primes only,
    weighted by 1/sample_frac (v must be computed there and at every validation positive)."""
    if sets is None:
        sets = label_sets()
    v = np.asarray(v, dtype=np.float64)
    vD, wD = v[DMASK], W[DMASK]
    classes = candidate_classes(vD)
    if not classes:
        return None
    pmfD, ED = disc_pmf(vD, wD, classes)
    keep = ED >= 0.5
    classes = [c for c, k in zip(classes, keep) if k]
    if not classes:
        return None
    pmfD, ED = pmfD[keep], ED[keep]
    sfD = survival(pmfD)
    EDtot = wD.sum()
    vm = VMASK if sample_mask is None else (VMASK & sample_mask)
    scale = 1.0 if sample_mask is None else 1.0 / sample_frac
    stV = class_stats(v[vm], W[vm], classes, scale)
    EV = stV[:, 0]
    EVtot = W[vm].sum() * scale
    sfV = survival(pmf_from_sums(stV[:, 1], stV[:, 2:]))
    nc, ns = len(classes), len(sets)
    res = {k: np.zeros(ns) for k in ("cls", "Od", "Ed", "pd", "Ov", "Ev", "pv", "bits")}
    dpos_all = [s[DMASK[s]] for s in sets]
    vpos_all = [s[VMASK[s]] for s in sets]
    for si in range(ns):
        xd, xv = v[dpos_all[si]], v[vpos_all[si]]
        Od = np.array([in_class(xd, c).sum() for c in classes])
        Ov = np.array([in_class(xv, c).sum() for c in classes])
        pd = sfD[np.arange(nc), np.minimum(Od, K)]
        ratio = Od / ED
        key = np.lexsort((-ratio, np.round(np.log(pd), 9)))
        b = key[0]
        r_in = (Od[b] + 0.5) / (ED[b] + 0.5)
        r_out = (len(xd) - Od[b] + 0.5) / (EDtot - ED[b] + 0.5)
        Ein, Eout = EV[b], EVtot - EV[b]
        Oin, Oout = Ov[b], len(xv) - Ov[b]
        bits = (Oin * np.log(r_in) - (r_in - 1) * Ein + Oout * np.log(r_out) - (r_out - 1) * Eout) / np.log(2)
        res["cls"][si] = b
        res["Od"][si], res["Ed"][si], res["pd"][si] = Od[b], ED[b], pd[b]
        res["Ov"][si], res["Ev"][si], res["pv"][si] = Ov[b], EV[b], sfV[b, min(Ov[b], K)]
        res["bits"][si] = bits
    res["classes"] = classes
    return res


def describe(cls):
    typ, t = cls
    return ("<=", ">", "==")[typ] + " %.6g" % t


def final_exam_score(vF, wF, posF, cls):
    """Score ONE frozen class on the final-exam range (called only by final_exam.py)."""
    vF = np.asarray(vF, dtype=np.float64)
    st = class_stats(vF, wF, [cls])
    sf = survival(pmf_from_sums(st[:, 1], st[:, 2:]))
    O = int(in_class(vF[posF], cls).sum())
    return O, st[0, 0], sf[0, min(O, K)]
