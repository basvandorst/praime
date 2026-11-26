"""Class E optima near u -> 0: S(u) = sum_k e^{-ku/2} g(ku) ~ I0/u - g(0)/2 + O(u), I0 = int_0^oo e^{-v/2} g(v) dv.
Prints I0 (Simpson on the grid u = i/400, i <= 2400, from grid.pkl), S(u) at the first grid points, and
min_{u>0 on grid} S/scale; also g at the touch points' neighbours (g itself may be negative in class E).
Usage: python classE_check.py"""
import json, glob
import mpmath as mp
import lp_opt as LO

mp.mp.dps = 50
files = sorted(glob.glob("polish_E_d*.json"), key=lambda s: int(s.split("_d")[1][:-5]))
for f in files:
    r = json.load(open(f))
    b = [mp.mpf(x) for x in r["b"]]
    n = len(b)
    if n - 1 > LO.JMAX:
        continue
    h = LO.DELTA
    gv = [mp.fsum(b[j] * LO.G[i][j] for j in range(n)) for i in range(LO.IMAX + 1)]
    fv = [mp.exp(-i * h / 2) * gv[i] for i in range(LO.IMAX + 1)]
    simpson = h / 3 * (fv[0] + fv[-1] + 4 * mp.fsum(fv[1:-1:2]) + 2 * mp.fsum(fv[2:-1:2]))
    S = [None] + [mp.fsum(b[j] * row[j] for j in range(n)) for row in (LO.S_grid_row(i, n - 1) for i in range(1, 6))]
    gneg = min((gv[i], i * h) for i in range(LO.IMAX + 1))
    print("%-18s I0 = %s   S(1/400..5/400) = %s   min g on grid = %s at u = %.4f" % (
        f, mp.nstr(simpson, 12), ", ".join(mp.nstr(x, 8) for x in S[1:]), mp.nstr(gneg[0], 4), float(gneg[1])))
