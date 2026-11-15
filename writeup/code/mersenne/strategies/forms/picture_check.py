"""Independent check of picture_stats.py.

  points  : sum over b of the number of a (closed formula with isqrt),
  dots    : from qf_check.json (Cornacchia counts),
  B       : recomputed from h of PARI (qfbclassno),
  B_emp   : plain-Python sum of 1/ln n over the lattice points,
  row1 C, best-row C, best-diagonal C: root counts by brute force over
            t = 0..p-1 in plain Python (not the Legendre formula), p <= 10000
            from sympy.primerange.
Also: the fraction of drawn primes (p not dividing 4P) is below 1/(2h) for how
many P (from qf_count.json, recomputed from the Cornacchia counts).
"""
import json
from math import isqrt, log, pi, sqrt

import cypari2
import numpy as np
from sympy import primerange

pari = cypari2.Pari()
N = 500_000
SMALLP = list(primerange(2, 10_001))


def bh_brute(A, B, C):
    prod = 1.0
    for p in SMALLP:
        t = np.arange(p, dtype=np.int64)
        v = (((A % p) * t % p) * t + (B % p) * t + (C % p)) % p
        w = int((v == 0).sum())
        prod *= (1 - w / p) / (1 - 1 / p)
    return prod


def main():
    ref = {r["P"]: r for r in json.load(open("picture_stats.json"))}
    cnt = {r["P"]: r for r in json.load(open("qf_count.json"))["rows"]}
    cor = {r["P"]: r for r in json.load(open("qf_check.json"))}
    bad = 0
    below = 0
    rel = []
    for P in range(1, 101):
        r = ref[P]
        pts = 0
        b = 0
        while P * b * b <= N:
            amax = isqrt(N - P * b * b)
            lo = b if P == 1 else 0
            pts += max(0, amax - lo + 1)
            b += 1
        S = 0.0
        b = 0
        while P * b * b <= N:
            lo = b if P == 1 else 0
            for a in range(lo, isqrt(N - P * b * b) + 1):
                n = a * a + P * b * b
                if n >= 2:
                    S += 1 / log(n)
            b += 1
        h = int(pari.qfbclassno(-4 * P))
        dots = cor[P]["drawn_all"]
        B = 4 / pi if P == 1 else 2 * sqrt(P) / (pi * h)
        c1 = bh_brute(1, 0, P)
        br = r["best_row"]
        cbr = bh_brute(1, 0, P * br["index"] ** 2)
        bd = r["best_diag"]
        bc = r["best_col"]       # columns a = c: n = P b^2 + c^2 (none for P = 1 or P >= 50)
        cbc = bh_brute(P, 0, bc["index"] ** 2) if bc else None
        ok_col = (bc is None) or abs(cbc - bc["C"]) < 1e-9
        # diagonals a - b = c have < 100 points for P >= 50: none reported there
        cbd = bh_brute(1 + P, 2 * bd["index"], bd["index"] ** 2) if bd else None
        if not bd:
            bd = {"C": None}
        ok = (pts == r["points"] and dots == r["dots"] and abs(B - r["B"]) < 1e-12
              and abs(dots / S - r["B_emp"]) < 1e-9 and abs(c1 - r["row1_C"]) < 1e-9
              and abs(cbr - br["C"]) < 1e-9 and ok_col
              and (cbd is None and bd["C"] is None or abs(cbd - bd["C"]) < 1e-9))
        bad += not ok
        frac = cnt[P]["drawn_excl"] / cnt[P]["n_excl"]
        below += frac < 1 / (2 * h)
        rel.append((1 - frac * 2 * h, P))
        print(f"P={P:3d} pts={pts} dots={dots} B={B:.4f} B_emp={dots/S:.4f} "
              f"row1 C={c1:.4f} best row C={cbr:.4f} best diag C={cbd} best col C={cbc} "
              f"{'OK' if ok else 'MISMATCH'}", flush=True)
    rel.sort()
    res = {"mismatches": bad, "P with fraction below 1/(2h)": below,
           "smallest relative deficit": rel[0], "largest relative deficit": rel[-1],
           "densest P": max(ref.values(), key=lambda r: r["density"])["P"]}
    print(res)
    json.dump(res, open("picture_check.json", "w"), indent=1)


if __name__ == "__main__":
    main()
