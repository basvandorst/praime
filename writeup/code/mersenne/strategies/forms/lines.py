"""Question 5, method A: bright straight lines in the quadform (P = 1, 2, 3)
and eisenstein layouts.

Lattice points (a, b) of the drawn region R:
  quadform P=1 : 0 <= b <= a,   a^2 + b^2 <= N        (smallest-b rep has b <= a)
  quadform P>=2: a, b >= 0,     a^2 + P b^2 <= N
  eisenstein   : a >= b >= 0,   a^2 + a b + b^2 <= N
Every prime with a representation has exactly one point in R (checked), so the
primes on a line are the lattice points of R on it whose norm is prime.
Lines: direction (da, db) primitive, 0 <= db <= 3, |da| <= 3; key = db*a - da*b.
For each line: number of lattice points, number of primes, the quadratic
f(t) = Q(a0 + t da, b0 + t db) along it, the Bateman-Horn constant
C(f) = prod_{p <= 10000} (1 - w(p)/p)/(1 - 1/p) with w(p) = #roots of f mod p
(brute force with numpy), and the prediction C(f) * sum 1/ln f(t) over the
lattice points of R on the line with f(t) >= 2.
Background: the chance that a point of R with norm n is prime is ~ B/ln n with
B = (share of primes represented) / (lattice points of R per unit of n):
  P=1: (1/2)/(pi/8) = 4/pi; P>=2: (1/(2h)) * 4 sqrt(P)/pi; eis: (1/2)/(pi/(6 sqrt 3)).
contrast = C(f)/B: > 1 means the line is lit more densely than its surroundings.
Writes lines.json; checked by lines_check.py.
"""
import json
from math import gcd, isqrt, log, pi, sqrt

import numpy as np

N = 500_000
PMAX_BH = 10_000


def sieve(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = False
    return s


ISP = sieve(N)
SMALLP = np.nonzero(sieve(PMAX_BH))[0]


def layout(name):
    """Return (Q(a,b) function on ints, lattice points a, b arrays, background B)."""
    A = np.arange(0, isqrt(N) + 1, dtype=np.int64)
    aa, bb = np.meshgrid(A, A, indexing="ij")
    aa, bb = aa.ravel(), bb.ravel()
    if name == "eis":
        Q = lambda a, b: a * a + a * b + b * b
        m = (aa >= bb) & (Q(aa, bb) <= N)
        B = 0.5 / (pi / (6 * sqrt(3)))
    else:
        P = name
        Q = lambda a, b, P=P: a * a + P * b * b
        m = Q(aa, bb) <= N
        if P == 1:
            m &= bb <= aa
            B = 0.5 / (pi / 8)
        else:
            h = {2: 1, 3: 1}[P]
            B = (1 / (2 * h)) * 4 * sqrt(P) / pi
    return Q, aa[m], bb[m], B


def coeffs(name, a0, b0, da, db):
    """f(t) = A t^2 + Bc t + C for the point (a0 + t da, b0 + t db)."""
    if name == "eis":
        A = da * da + da * db + db * db
        Bc = 2 * a0 * da + a0 * db + b0 * da + 2 * b0 * db
        C = a0 * a0 + a0 * b0 + b0 * b0
    else:
        P = name
        A = da * da + P * db * db
        Bc = 2 * a0 * da + 2 * P * b0 * db
        C = a0 * a0 + P * b0 * b0
    return A, Bc, C


def bh_constant(A, Bc, C):
    prod = 1.0
    for p in SMALLP:
        p = int(p)
        t = np.arange(p, dtype=np.int64)
        v = ((A % p) * t % p * t + (Bc % p) * t + (C % p)) % p
        w = int((v == 0).sum())
        prod *= (1 - w / p) / (1 - 1 / p)
    return prod


def describe(da, db, key):
    if db == 0:
        return f"b = {-key}"
    if da == 0:
        return f"a = {key}"
    # key = db*a - da*b
    left = (f"{db}a" if db != 1 else "a")
    if da > 0:
        right = f" - {da if da != 1 else ''}b"
    else:
        right = f" + {-da if da != -1 else ''}b"
    return f"{left}{right} = {key}"


def base_point(da, db, key):
    """Lattice point on the line with 0 <= b0 < db (db > 0), or a0 = 0 (db = 0)."""
    if db == 0:
        return 0, -key
    for b0 in range(db):
        if (key + da * b0) % db == 0:
            return (key + da * b0) // db, b0
    raise ValueError


def poly_str(A, Bc, C, var):
    s = f"{A}{var}^2" if A != 1 else f"{var}^2"
    if Bc:
        s += f" {'+' if Bc > 0 else '-'} {abs(Bc) if abs(Bc) != 1 else ''}{var}"
    if C:
        s += f" {'+' if C > 0 else '-'} {abs(C)}"
    return s


DIRS = [(da, db) for db in range(0, 4) for da in range(-3, 4)
        if gcd(abs(da), db) == 1 and (db > 0 or da == 1)]


def analyse(name):
    Q, a, b, B = layout(name)
    n = Q(a, b)
    prime = ISP[n]
    # sanity: each represented prime appears exactly once in R
    pn = n[prime]
    assert len(np.unique(pn)) == len(pn)
    lines = []
    for (da, db) in DIRS:
        key = db * a - da * b
        uk, inv = np.unique(key, return_inverse=True)
        npts = np.bincount(inv)
        npr = np.bincount(inv, weights=prime.astype(np.float64)).astype(int)
        for i in np.nonzero(npr >= 5)[0]:
            lines.append((int(npr[i]), int(npts[i]), da, db, int(uk[i])))
    lines.sort(reverse=True)
    # empirical background: primes in R / sum over points of R of 1/ln n
    nn = n[n >= 2].astype(np.float64)
    B_emp = prime.sum() / np.sum(1 / np.log(nn))

    def detail(npr, npts, da, db, k):
        a0, b0 = base_point(da, db, k)
        A, Bc, C = coeffs(name, a0, b0, da, db)
        on = (db * a - da * b) == k
        ts = ((a[on] - a0) // da) if da != 0 else ((b[on] - b0) // db)
        if db == 0:
            ts = a[on]
        f = A * ts * ts + Bc * ts + C
        assert np.array_equal(np.sort(f), np.sort(n[on]))
        fv = f[f >= 2].astype(np.float64)
        S = float(np.sum(1 / np.log(fv)))
        Cbh = bh_constant(A, Bc, C)
        var = "a" if db == 0 else ("b" if db == 1 else "t")
        return {
            "line": describe(da, db, k), "dir": [da, db], "key": k,
            "points": npts, "primes": npr,
            "n(t)": poly_str(A, Bc, C, var),
            "t_range": [int(ts.min()), int(ts.max())],
            "C_BH": Cbh, "sum_1_over_ln": S, "predicted": Cbh * S,
            "ratio_obs_pred": npr / (Cbh * S) if Cbh * S > 0 else None,
            "contrast_C_over_B": Cbh / B,
            "local_rate_over_B_emp": (npr / S) / B_emp if S > 0 else None,
        }

    UNIT = {(1, 0), (0, 1), (1, 1), (-1, 1)}
    unit_lines = [L for L in lines if (L[2], L[3]) in UNIT]
    top_count = [detail(*L) for L in unit_lines[:5]]

    def top_contrast_of(cands, k=5):
        # brightest by contrast among long lines (>= 100 lattice points, >= 10 primes);
        # quick C with primes <= 200 to preselect 40, then the exact C (p <= 10000)
        pre = []
        for L in cands:
            npr, npts, da, db, kk = L
            if npts < 100 or npr < 10:
                continue
            a0, b0 = base_point(da, db, kk)
            A, Bc, C = coeffs(name, a0, b0, da, db)
            prod = 1.0
            for p in SMALLP[:46]:
                p = int(p)
                t = np.arange(p, dtype=np.int64)
                w = int((((A % p) * t % p * t + (Bc % p) * t + C % p) % p == 0).sum())
                prod *= (1 - w / p) / (1 - 1 / p)
            pre.append((prod, L))
        pre.sort(reverse=True)
        return sorted((detail(*L) for _, L in pre[:40]),
                      key=lambda d: (-round(d["contrast_C_over_B"], 9), -d["primes"]))[:k]

    top_contrast = top_contrast_of(unit_lines)
    top_contrast_any = top_contrast_of(lines)
    named = {}
    wanted = {1: [(1, 0, -1), (1, 1, 1), (1, 0, -2)],
              2: [(1, 0, -1), (0, 1, 1), (1, 1, 1)],
              3: [(1, 0, -1), (0, 1, 1), (1, 1, 1), (0, 1, 2)],
              "eis": [(1, 0, -1), (1, 1, 1)]}[name]
    for (da, db, k) in wanted:
        on = (db * a - da * b) == k
        npts, npr = int(on.sum()), int(prime[on].sum())
        named[describe(da, db, k)] = detail(npr, npts, da, db, k)
    return {"B_theory": B, "B_empirical": float(B_emp),
            "primes_in_R": int(prime.sum()), "points_in_R": int(len(a)),
            "top_by_count": top_count, "top_by_contrast": top_contrast,
            "top_by_contrast_any_direction": top_contrast_any,
            "named": named}


def main():
    out = {}
    for name in (1, 2, 3, "eis"):
        r = analyse(name)
        out[str(name)] = r
        print(f"=== {name}: B={r['B_theory']:.4f} B_emp={r['B_empirical']:.4f} "
              f"primes in R={r['primes_in_R']} points={r['points_in_R']}")
        for tag in ("top_by_count", "top_by_contrast", "top_by_contrast_any_direction"):
            print(" ", tag)
            for d in r[tag]:
                print(f"    {d['line']:>14s} n={d['n(t)']:<22s} pts={d['points']:4d} "
                      f"primes={d['primes']:4d} pred={d['predicted']:7.1f} "
                      f"C={d['C_BH']:.4f} C/B={d['contrast_C_over_B']:.3f}")
        print("  named")
        for k, d in r["named"].items():
            print(f"    {d['line']:>14s} n={d['n(t)']:<22s} pts={d['points']:4d} "
                  f"primes={d['primes']:4d} pred={d['predicted']:7.1f} "
                  f"C={d['C_BH']:.4f} C/B={d['contrast_C_over_B']:.3f}")
    json.dump(out, open("lines.json", "w"), indent=1)


if __name__ == "__main__":
    main()
