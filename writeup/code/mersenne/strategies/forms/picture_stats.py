"""Question 6, method A: numbers that describe how each quadform picture looks.

For P = 1..100, region R as on the page (P=1: 0 <= b <= a; else a, b >= 0;
a^2 + P b^2 <= N):
  points      = lattice points of R (numpy count),
  dots        = primes drawn (from qf_count.json),
  density     = dots / points (share of lattice points lit as primes),
  B           = 2 sqrt(P)/(pi h) (P >= 2), 4/pi (P = 1): background constant,
                chance that a point of R with norm n is prime ~ B / ln n,
  B_emp       = dots / sum_{points of R, n >= 2} 1/ln n,
  rows        : n = a^2 + P b^2 for fixed b = 1..20 (the horizontal lines);
                Bateman-Horn constant C (primes <= 10000, root counts by the
                Legendre-symbol formula), contrast = C/B,
  diagonals   : a - b = 1..20, columns a = 1..20, the same way.
The page draws (a, b sqrt P), so for large P the picture is a stack of
horizontal rows sqrt(P) apart; the row contrasts are what the eye sees.
Writes picture_stats.json; checked by picture_check.py.
"""
import json
from math import isqrt, log, pi, sqrt

import numpy as np

N = 500_000


def sieve(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = False
    return s


SMALLP = [int(p) for p in np.nonzero(sieve(10_000))[0]]


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def omega(A, B, C, p):
    if p == 2:
        return sum(1 for t in range(2) if (A * t * t + B * t + C) % 2 == 0)
    if A % p == 0:
        if B % p:
            return 1
        return p if C % p == 0 else 0
    return 1 + legendre(B * B - 4 * A * C, p)


def bh(A, B, C):
    prod = 1.0
    for p in SMALLP:
        prod *= (1 - omega(A, B, C, p) / p) / (1 - 1 / p)
    return prod


def main():
    ref = {r["P"]: r for r in json.load(open("qf_count.json"))["rows"]}
    out = []
    for P in range(1, 101):
        h = ref[P]["h"]
        dots = ref[P]["drawn_all"]
        A = np.arange(isqrt(N) + 1, dtype=np.int64)
        aa, bb = np.meshgrid(A, A[: isqrt(N // P) + 1], indexing="ij")
        n = aa * aa + P * bb * bb
        m = n <= N
        if P == 1:
            m &= bb <= aa
        points = int(m.sum())
        nn = n[m & (n >= 2)].astype(np.float64)
        B_emp = dots / float(np.sum(1 / np.log(nn)))
        B = 4 / pi if P == 1 else 2 * sqrt(P) / (pi * h)
        bmax = isqrt(N // P)

        def row(b):     # n = a^2 + P b^2, a free
            lo = b if P == 1 else 0
            length = isqrt(N - P * b * b) - lo + 1
            return bh(1, 0, P * b * b), length

        def col(a):     # n = P b^2 + a^2, b free
            if P == 1:
                return None, 0
            return bh(P, 0, a * a), isqrt((N - a * a) // P) + 1

        def diag(c):    # a = b + c: n = (1+P) b^2 + 2 c b + c^2
            cnt = sum(1 for b in range(0, isqrt(N) + 1) if (b + c) ** 2 + P * b * b <= N)
            return bh(1 + P, 2 * c, c * c), cnt

        rows = {b: row(b) for b in range(1, min(20, bmax) + 1)}
        cols = {a: col(a) for a in range(1, 21)}
        diags = {c: diag(c) for c in range(1, 21)}

        def best(d, sign):
            cand = [(v[0] / B, k, v[0], v[1]) for k, v in d.items()
                    if v[0] is not None and v[1] >= 100]
            if not cand:
                return None
            x = max(cand) if sign > 0 else min(cand)
            return {"index": x[1], "C": x[2], "contrast": x[0], "points": x[3]}

        rec = {
            "P": P, "h": h, "points": points, "dots": dots,
            "density": dots / points, "B": B, "B_emp": B_emp,
            "rows_in_picture": bmax + 1,
            "row1_C": rows[1][0], "row1_contrast": rows[1][0] / B,
            "row1_len": rows[1][1],
            "col1_C": cols[1][0], "diag1_C": diags[1][0],
            "best_row": best(rows, 1), "worst_row": best(rows, -1),
            "best_col": best(cols, 1), "best_diag": best(diags, 1),
            "row_contrasts": {b: rows[b][0] / B for b in rows},
        }
        out.append(rec)
        br = rec["best_row"]
        print(f"P={P:3d} h={h:2d} pts={points:6d} dots={dots:5d} dens={dots/points:.4f} "
              f"B={B:.3f} Bemp={B_emp:.3f} rows={bmax+1:3d} "
              f"row1 C={rows[1][0]:.3f} ({rows[1][0]/B:.2f}xB) "
              f"best row b={br['index'] if br else '-'} "
              f"{(br['contrast'] if br else 0):.2f}xB", flush=True)
    json.dump(out, open("picture_stats.json", "w"), indent=1)


if __name__ == "__main__":
    main()
