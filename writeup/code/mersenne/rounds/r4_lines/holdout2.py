"""Experiment 4b: holdout tests with power.  Training = the 28 exponents
below 10^5, holdout = the 23 in [10^5, B).  Null for the holdout: training
fixed, holdout replaced by 2000 Wagstaff sets of exactly 23 primes in
[10^5, B) (same exact Poisson-splitting sampler as gen_null.py).

(a) Columns: for each of the 10 most striking w on the training data
    (columns_h5.json 'top', ranked by max z; chosen without the holdout),
    C_w = number of pairs (h, x) in the same column with h a holdout exponent
    and x any other exponent.  Tests: w = top-1 alone, and sum over the 10.
(b) Lines: at the training-selected w (smallest per-w p in lines_h5) and
    at the page default w = 708: T_new = collinear triples with at least one
    holdout point = T(51) - T(28).
(c) Lines: among all w = 2..3000 and all lines through >= 3 training
    exponents that have cells in [10^5, B), the one with the smallest
    training Poisson tail; count holdout exponents on it, against the
    Poisson-binomial over the primes on its holdout cells.
"""
import json, os, math
from itertools import combinations
import numpy as np
from common import EXPS, B, HERE, load_primes
from holdout import line_points, prob_of, pb_tail, lines_through

p_all, pr_all = load_primes()
LO = 10**5
train = np.array([e for e in EXPS if e < LO], dtype=np.int64)
hold = np.array([e for e in EXPS if e >= LO], dtype=np.int64)
NH = 2000
rng = np.random.default_rng(777)
sel = p_all >= LO
rp, rq = p_all[sel], pr_all[sel]
mu = -np.log1p(-rq); cmu = np.cumsum(mu); tot = cmu[-1]
H = np.empty((NH, 23), dtype=np.int64); got = 0
while got < NH:
    idx = np.unique(np.searchsorted(cmu, rng.random(rng.poisson(tot)) * tot, side="right"))
    idx = idx[idx < len(rp)]
    if len(idx) == 23:
        H[got] = rp[idx]; got += 1
HS = np.vstack([hold[None, :], H])               # row 0 = real holdout
S = HS.shape[0]
out = {"NH": NH}


def exch(s):
    s = np.asarray(s, dtype=np.float64)
    return float((s >= s[0] - 1e-12).sum() / len(s))


# (a) columns
cj = json.load(open(os.path.join(HERE, "columns_h5.json")))
wtop = [r["w"] for r in cj["top"][:10]]
Csum = np.zeros(S)
res = {}
for w in wtop:
    rt = train % w
    rh = HS % w                                  # (S, 23)
    c_tr = (rh[:, :, None] == rt[None, None, :]).sum(axis=(1, 2))
    srt = np.sort(rh, axis=1)
    # pairs inside the holdout: sum over runs C(m,2)
    eq = (rh[:, :, None] == rh[:, None, :]).sum(axis=(1, 2)) - 23
    C = c_tr + eq // 2
    res[w] = {"real": int(C[0]), "null_mean": float(C[1:].mean()), "p": exch(C)}
    Csum += C
out["columns"] = {"w_list": wtop, "per_w": res,
                  "sum": {"real": int(Csum[0]), "null_mean": float(Csum[1:].mean()),
                          "p": exch(Csum)}}
print("columns:", out["columns"])


# (b) lines, T_new at fixed w
def T_sets(X, w):
    Sx, K = X.shape
    I, J = np.triu_indices(K, 1)
    r = (X - 1) // w; c = (X - 1) % w
    dr = (r[:, J] - r[:, I]).ravel(); dc = (c[:, J] - c[:, I]).ravel()
    # order points by n so dr >= 0: X rows are sorted
    g = np.gcd(dr, np.abs(dc))
    a = dr // g; b = dc // g + w
    SI = (np.arange(Sx)[:, None] * K + I[None, :]).ravel()
    key = (SI << 39) | (a << 13) | b
    key.sort()
    st = np.empty(len(key), dtype=bool); st[0] = True; st[1:] = key[1:] != key[:-1]
    cnt = np.bincount(np.cumsum(st) - 1).astype(np.float64)
    gset = (key[st] >> 39) // K
    return np.bincount(gset, cnt * (cnt - 1) / 2, minlength=Sx)


lj = json.load(open(os.path.join(HERE, "lines_h5.json")))
w_sel = lj["summary"]["ALL"]["w"]
out["lines_Tnew"] = {}
T28 = None
for w in sorted({w_sel, 708}):
    T28 = T_sets(train[None, :], w)[0]
    tn = []
    for c0 in range(0, S, 500):
        blk = HS[c0:c0 + 500]
        X = np.sort(np.hstack([np.repeat(train[None, :], len(blk), 0), blk]), axis=1)
        tn.append(T_sets(X, w) - T28)
    tn = np.concatenate(tn)
    out["lines_Tnew"][w] = {"real": float(tn[0]), "null_mean": float(tn[1:].mean()),
                            "p": exch(tn)}
    print("T_new w=", w, out["lines_Tnew"][w])

# (c) most surprising extendable line over all w
best = None
nlines = 0
p_tr = p_all[p_all < LO]; q_tr = pr_all[p_all < LO]
p_h = p_all[p_all >= LO]; q_h = pr_all[p_all >= LO]


def holdout_cells(key, w):
    """(n values, probs) of the primes on the line in [LO, B)."""
    dr, dc, inv = key
    if dc == 0:                                  # column c = -inv: n = c+1 mod w
        m = (p_h % w) == ((-inv + 1) % w)
        return p_h[m], q_h[m]
    ext = line_points(key, w, LO, B)
    pr = prob_of(ext)
    return ext[pr > 0], pr[pr > 0]


for w in range(2, 3001):
    colmass = None
    for key, mem in lines_through(train.tolist(), w).items():
        dr, dc, inv = key
        if dc == 0:
            if colmass is None:
                colmass = np.bincount(p_h % w, weights=q_h, minlength=w)
            if colmass[(-inv + 1) % w] == 0:
                continue
            m = (p_tr % w) == ((-inv + 1) % w)
            tail = pb_tail(q_tr[m], len(mem))
        else:
            ext = line_points(key, w, LO, B)
            if len(ext) == 0 or prob_of(ext).sum() == 0:
                continue
            tail = pb_tail(prob_of(line_points(key, w, 1, LO)), len(mem))
        nlines += 1
        if best is None or tail < best[0]:
            best = (tail, w, key, sorted(mem))
tail, w, key, mem = best
ext, pe = holdout_cells(key, w)
hits = [int(e) for e in hold if e in set(ext.tolist())]
out["best_extendable_line"] = {
    "n_candidate_lines": nlines, "w": w, "line(dr,dc,inv)": key, "train_members": mem,
    "train_tail": tail, "holdout_cells": int(len(ext)),
    "expected_holdout": float(pe.sum()), "holdout_hits": hits,
    "p": pb_tail(pe, len(hits))}
print("best extendable line:", out["best_extendable_line"])
with open(os.path.join(HERE, "holdout2.json"), "w") as f:
    json.dump(out, f, indent=1, default=str)
