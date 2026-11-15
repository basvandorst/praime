"""Experiment 5: independent recheck in plain Python integers.

- primes: sympy.isprime on the real exponents (all < 2^64, so a proof).
- columns: D, M, chi^2 (as an exact Fraction) via collections.Counter at the
  reported best w; the per-w p recomputed over the same null sets with the
  same exchangeable rule.
- lines: collinear triples by exact cross products over all triples,
  4-sets over all quadruples, most points on a line by grouping pairs on
  (reduced slope Fraction or 'vertical', intercept Fraction); per-w p for T
  recomputed on all null sets at the reported best w.
- spiral: coordinates by walking the spiral step by step (no formula),
  triples by cross products, A0 and Dmax by counting.
"""
import json, os, sys
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import gcd
import numpy as np
from sympy import isprime, totient
from common import EXPS, HERE

assert all(isprime(e) for e in EXPS)
LIM = {"full": 10**9, "page": 500000, "h5": 100000}
report = {}


def col_stats(xs, w):
    c = Counter(x % w for x in xs)
    D = len(c)
    M = max(c.values())
    co = [x for x in xs if gcd(x, w) == 1]
    nc = len(co)
    cc = Counter(x % w for x in co)
    phi = int(totient(w))
    chi = Fraction(phi, nc) * sum(v * v for v in cc.values()) - nc
    return D, M, chi


def exch_p(vals, larger_is_striking=True):
    s = vals if larger_is_striking else [-v for v in vals]
    return sum(1 for v in s if v >= s[0]) / len(s)


def collinear(p, q, r):
    return (q[0] - p[0]) * (r[1] - p[1]) - (r[0] - p[0]) * (q[1] - p[1]) == 0


def line_stats_bruteforce(pts):
    T = sum(1 for a, b, c in combinations(pts, 3) if collinear(a, b, c))
    Q = sum(1 for a, b, c, d in combinations(pts, 4)
            if collinear(a, b, c) and collinear(a, b, d))
    lines = {}
    for a, b in combinations(pts, 2):
        if a[0] == b[0]:
            key = ("x", a[0])
        else:
            s = Fraction(b[1] - a[1], b[0] - a[0])
            key = (s, a[1] - s * a[0])
        lines.setdefault(key, set()).update((a, b))
    Mx = max(len(v) for v in lines.values())
    return T, Q, Mx


def square(n, w):
    return ((n - 1) // w, (n - 1) % w)


# ---- columns
for view in ("full", "page", "h5"):
    cj = json.load(open(os.path.join(HERE, f"columns_{view}.json")))
    null = np.load(os.path.join(HERE, f"null_{view}.npy")).tolist()
    real = [e for e in EXPS if e < LIM[view]]
    sets = [real] + null
    ws = sorted({cj["summary"][k]["w"] for k in cj["summary"]} |
                {cj["summary"][k]["w_z"] for k in cj["summary"]})
    for w in ws:
        st = [col_stats(s, w) for s in sets]
        D = [x[0] for x in st]; M = [x[1] for x in st]; C = [x[2] for x in st]
        rec = {"D": D[0], "M": M[0], "chi2": str(C[0]),
               "p_D": exch_p(D, False), "p_M": exch_p(M), "p_chi": exch_p(C)}
        # compare with the vectorised run
        r0 = [r for r in json.load(open(os.path.join(HERE, f"columns_{view}_perw.json"))) if r["w"] == w]
        if r0:
            r0 = r0[0]
            rec["match_vectorised"] = (r0["D"] == D[0] and r0["M"] == M[0]
                                       and abs(r0["CHI"] - float(C[0])) < 1e-6
                                       and abs(r0["p_D_low"] - rec["p_D"]) < 1e-12
                                       and abs(r0["p_M_high"] - rec["p_M"]) < 1e-12
                                       and abs(r0["p_CHI_high"] - rec["p_chi"]) < 1e-12)
        report[f"col_{view}_w{w}"] = rec
        print(view, w, rec, flush=True)

# ---- lines (square layout)
for view in ("page", "h5", "full"):
    f = os.path.join(HERE, f"lines_{view}_perw.json")
    if not os.path.exists(f):
        continue
    perw = {r["w"]: r for r in json.load(open(f))}
    lj = json.load(open(os.path.join(HERE, f"lines_{view}.json")))
    real = [e for e in EXPS if e < LIM[view]]
    null = np.load(os.path.join(HERE, f"null_{view}.npy"))[:lj["N"]].tolist()
    ws = sorted({lj["summary"][k]["w"] for k in lj["summary"]} | {708})
    for w in ws:
        T, Q, Mx = line_stats_bruteforce([square(n, w) for n in real])
        rec = {"T": T, "Q": Q, "Mx": Mx}
        v = perw[w]
        rec["match_vectorised"] = (v["T"] == T and v["Mx"] == Mx and v["Q"] == Q)
        if True:                 # recompute the per-w p of T over all null sets
            Ts = [T] + [sum(1 for a, b, c in combinations([square(n, w) for n in s], 3)
                            if collinear(a, b, c)) for s in null]
            rec["p_T"] = exch_p(Ts)
            rec["p_T_match"] = abs(rec["p_T"] - v["p_T"]) < 1e-12
        report[f"lines_{view}_w{w}"] = rec
        print("lines", view, w, rec, flush=True)

# ---- spiral by walking
def walk(nmax):
    pos = {1: (0, 0)}
    x = y = 0
    n = 1
    step = 1
    dirs = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    d = 0
    while n < nmax:
        for _ in range(2):
            dx, dy = dirs[d % 4]
            for _ in range(step):
                if n >= nmax:
                    break
                x += dx; y += dy; n += 1
                pos[n] = (x, y)
            d += 1
        step += 1
    return pos

pos = walk(600000)
for view in ("page", "full"):
    sj = json.load(open(os.path.join(HERE, f"spiral_{view}.json")))
    real = [e for e in EXPS if e < LIM[view]]
    if view == "page":
        # formula coordinates agree with the walk up to a rotation/reflection:
        # compare the invariant |x|,|y| sets and pairwise squared distances
        ok = all(sorted(map(abs, pos[e])) == sorted(map(abs, sj["coords"][str(e)]))
                 for e in real)
        dist_ok = all(
            (pos[a][0]-pos[b][0])**2 + (pos[a][1]-pos[b][1])**2 ==
            (sj["coords"][str(a)][0]-sj["coords"][str(b)][0])**2 +
            (sj["coords"][str(a)][1]-sj["coords"][str(b)][1])**2
            for a, b in combinations(real, 2))
        pts = [pos[e] for e in real]
    else:
        ok = dist_ok = None
        pts = [tuple(sj["coords"][str(e)]) for e in real]
    T, Q, Mx = line_stats_bruteforce(pts)
    A0 = sum(1 for x, y in pts if x == 0 or y == 0 or x == y or x == -y)
    Dmax = max(max(Counter(x - y for x, y in pts).values()),
               max(Counter(x + y for x, y in pts).values()))
    rec = {"coords_abs_match": ok, "distances_match": dist_ok, "T": T, "Mx": Mx, "A0": A0,
           "Dmax": Dmax,
           "match": (T == sj["stats"]["T"]["real"] and Mx == sj["stats"]["Mx"]["real"]
                     and A0 == sj["stats"]["A0"]["real"] and Dmax == sj["stats"]["Dmax"]["real"])}
    report[f"spiral_{view}"] = rec
    print("spiral", view, rec, flush=True)

with open(os.path.join(HERE, "recheck.json"), "w") as f:
    json.dump(report, f, indent=1)
