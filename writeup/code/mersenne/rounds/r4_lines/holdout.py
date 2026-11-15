"""Experiment 4: holdout.  Choose on the 28 exponents below 10^5 (null h5),
test on the 23 exponents in [10^5, B).

Columns: the w with the smallest per-w p (ties: largest z) in columns_h5.json;
  the most populated column(s) at that w, and all occupied columns.
Lines: the w with the smallest per-w p in lines_h5_perw.json (ties: largest z),
  plus the page default w = 708; at that w, among all lines through >= 3
  training exponents, the one with the smallest Poisson tail
  P(Poisson(lambda) >= m), lambda = sum of Wagstaff P over primes on the line
  below 10^5.  Extend the line over the grid to [10^5, B) and count exponents.
Null for the test: (a) unconditional Poisson-binomial over the primes on the
extension; (b) conditional on 23 holdout exponents: Binomial(23, q), q = the
line's share of the Wagstaff mass in [10^5, B).
"""
import json, os, math
from itertools import combinations
from fractions import Fraction
import numpy as np
from scipy.stats import binom
from common import EXPS, B, HERE, load_primes

p_all, pr_all = load_primes()
LO = 10**5
train = [e for e in EXPS if e < LO]
hold = [e for e in EXPS if LO <= e < B]
assert len(train) == 28 and len(hold) == 23
mass_hold = pr_all[(p_all >= LO)].sum()


def prob_of(ns):
    ns = np.asarray(ns, dtype=np.int64)
    if len(ns) == 0:
        return np.zeros(0)
    i = np.searchsorted(p_all, ns)
    i = np.minimum(i, len(p_all) - 1)
    ok = p_all[i] == ns
    return np.where(ok, pr_all[i], 0.0)


def pb_tail(probs, k):
    """P(sum Bernoulli(probs) >= k), exact DP truncated at k."""
    if k <= 0:
        return 1.0
    dist = np.zeros(k + 1); dist[0] = 1.0
    for q in probs:
        if q <= 0:
            continue
        new = dist * (1 - q)
        new[1:] += dist[:-1] * q
        new[k] += dist[k] * q          # absorbing state ">= k"
        dist = new
    return float(dist[k])


def test_set(label, members_hold, probs_ext):
    k = len(members_hold)
    lam = float(np.sum(probs_ext))
    q = lam / mass_hold
    res = {"label": label, "holdout_hits": k, "hits": members_hold,
           "expected_uncond": lam, "p_uncond": pb_tail(probs_ext, k),
           "q_share": q, "expected_cond23": 23 * q,
           "p_cond23": float(binom.sf(k - 1, 23, q))}
    print(res)
    return res


# ---------------- lines ----------------
def lines_through(pts_n, w):
    """All lines (as member sets) through >= 3 of the points, square layout w."""
    pos = {n: ((n - 1) // w, (n - 1) % w) for n in pts_n}
    lines = {}
    for a, b in combinations(pts_n, 2):
        (r1, c1), (r2, c2) = pos[a], pos[b]
        dr, dc = r2 - r1, c2 - c1
        g = math.gcd(dr, abs(dc))
        dr, dc = dr // g, dc // g
        if dr < 0 or (dr == 0 and dc < 0):
            dr, dc = -dr, -dc
        # line invariant: dc*r - dr*c
        key = (dr, dc, dc * r1 - dr * c1)
        lines.setdefault(key, set()).update((a, b))
    return {k: v for k, v in lines.items() if len(v) >= 3}


def line_points(key, w, n_lo, n_hi):
    """Integers n in [n_lo, n_hi) whose cell (r, c) lies on the line."""
    dr, dc, inv = key
    if dr == 0:                     # row r = -inv/dc ... dc = 1: inv = r
        r = inv // dc
        lo, hi = max(n_lo, r * w + 1), min(n_hi, r * w + w + 1)
        return np.arange(lo, hi, dtype=np.int64) if lo < hi else np.zeros(0, np.int64)
    if dc == 0:                     # column c = -inv/dr (dr = 1)
        c = -inv // dr
        n0 = c + 1
        start = n0 + max(0, -(-(n_lo - n0) // w)) * w
        return np.arange(start, n_hi, w, dtype=np.int64)
    out = []
    # points r = r0 + t*dr, c = c0 + t*dc; find one lattice point on the line
    for c0 in range(w):
        if (inv + dr * c0) % dc == 0:
            r0 = (inv + dr * c0) // dc
            break
    ts = []
    # c = c0 + t*dc in [0, w)
    tmin = math.ceil((0 - c0) / dc) if dc > 0 else math.ceil((w - 1 - c0) / dc)
    tmax = math.floor((w - 1 - c0) / dc) if dc > 0 else math.floor((0 - c0) / dc)
    t = np.arange(tmin, tmax + 1, dtype=np.int64)
    r = r0 + t * dr
    c = c0 + t * dc
    n = r * w + c + 1
    n = n[(r >= 0) & (n >= n_lo) & (n < n_hi)]
    return np.sort(n)



if __name__ == "__main__":
    out = {}
    # ---------------- columns ----------------
    cj = json.load(open(os.path.join(HERE, "columns_h5.json")))
    names = ["D_low", "M_high", "CHI_high", "G_high"]
    # rebuild per-w best from 'top' (top-15 by max z) + summary
    w_col = cj["summary"]["ALL"]["w"]
    w_colz = cj["summary"]["ALL"]["w_z"]
    print("column selection: min-p w =", w_col, " max-z w =", w_colz)
    for w in sorted({w_col, w_colz}):
        cnt = {}
        for e in train:
            cnt.setdefault(e % w, []).append(e)
        mx = max(len(v) for v in cnt.values())
        top = sorted(r for r, v in cnt.items() if len(v) == mx)
        occ = sorted(cnt)
        print(f"w={w}: top columns {[(r, cnt[r]) for r in top]}")
        ext_top = p_all[(p_all >= LO) & np.isin(p_all % w, top)]
        pr_top = prob_of(ext_top)
        out[f"col_w{w}_top"] = test_set(f"columns {top} mod {w}",
                                        [e for e in hold if e % w in top], pr_top)
        ext_occ = p_all[(p_all >= LO) & np.isin(p_all % w, occ)]
        out[f"col_w{w}_occ"] = test_set(f"{len(occ)} occupied columns mod {w}",
                                        [e for e in hold if e % w in occ], prob_of(ext_occ))
        # pairs: holdout exponents sharing a column with each other or the training set
        allc = {}
        for e in train + hold:
            allc.setdefault(e % w, []).append(e)
        out[f"col_w{w}_classes_with_2plus"] = {str(r): v for r, v in allc.items() if len(v) > 1}

    lj = json.load(open(os.path.join(HERE, "lines_h5_perw.json")))
    lnames = ["T", "Q", "Mx"]
    bestp = min(min(r["p_" + s] for s in lnames) for r in lj)
    cands = [(max(r["z_" + s] for s in lnames if r["p_" + s] == bestp), r["w"]) for r in lj
             if min(r["p_" + s] for s in lnames) == bestp]
    w_line = max(cands)[1]
    print("line selection: min p =", bestp, " w =", w_line, " (tied w:", len(cands), ")")
    for w in (w_line, 708):
        L = lines_through(train, w)
        best = None
        for key, mem in L.items():
            pts = line_points(key, w, 1, LO)
            assert all(m in set(pts.tolist()) for m in mem), (key, mem)
            pr = prob_of(pts)
            tail = pb_tail(pr, len(mem))
            if best is None or tail < best[0]:
                best = (tail, key, sorted(mem), float(pr.sum()), len(pts))
        tail, key, mem, lam, npts = best
        print(f"w={w}: {len(L)} lines with >=3 training exponents; most surprising "
              f"{key} members {mem} lambda_train={lam:.3f} tail={tail:.3g} ({npts} cells < 1e5)")
        ext = line_points(key, w, LO, B)
        hits = [e for e in hold if e in set(ext.tolist())]
        out[f"line_w{w}"] = {"train_line": {"key": key, "members": mem,
                                            "lambda_train": lam, "train_tail": tail,
                                            "n_lines_3plus": len(L)},
                             "cells_in_holdout_range": int(len(ext)),
                             "test": test_set(f"line {key} at w={w}", hits, prob_of(ext))}
    with open(os.path.join(HERE, "holdout.json"), "w") as f:
        json.dump(out, f, indent=1, default=str)
