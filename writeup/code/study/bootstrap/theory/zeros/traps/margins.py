"""D1 supplement: interior local maxima of rho(u) = eta h_k(u)/gt(u) on the scan grid, away from the main
component (a, b) of S (grid points with u < a or u > b, local maxima of rho in the grid order, excluding the
monotone tail of the main bump next to (a, b)).  A second component of S would show up as such a local
maximum with rho > 1; the largest one measures how far the scan is from finding one.  Also the minimum of
gt/(scale) on the grid (sanity of the optimum).  Usage: python margins.py > margins.out"""
import json, pickle, math
from fam import *

GRID = pickle.load(open("grid.pkl", "rb"))
ctx.prec = 400
GUf = [float(Fraction(s)) for s in GRID["u"]]
GROWS = [[arb(x) for x in row] for row in GRID["rows"]]


def fl(x):
    return float(x.mid().str(20, radius=False)) if isinstance(x, arb) else float(x)


print("deg  k  n  x   largest interior local max of rho off (a,b)  at u    | #local maxima > 1e-3")
for d in range(2, 26, 2):
    try:
        T = json.load(open("traps_d%d.json" % d))
    except FileNotFoundError:
        continue
    R = json.load(open("hp_d%d.json" % d))
    n = d + 1
    bg = [arb(s) for s in R["b"]]
    g = [dot(bg, row[:n]) for row in GROWS]
    for nd in T["nodes"]:
        k = nd["k"]
        bh = [arb(s) for s in R["h"][k - 1]["b"]]
        h = [dot(bh, row[:n]) for row in GROWS]
        for E in nd["etas"]:
            if "a" not in E:
                print("%3d %2d %3d %3d   %s" % (2 * d, k, nd["n"], E["factor"], E["status"]))
                continue
            eta = arb(E["eta"])
            a, b = float(E["a"]), float(E["b"])
            rho = [fl(eta * hv / gv) if gv > 0 else float("nan") for hv, gv in zip(h, g)]
            best, at, cnt = -1e300, None, 0
            for i in range(1, len(rho) - 1):
                u = GUf[i]
                if a <= u <= b:
                    continue
                # exclude the tail of the main bump: rho increasing towards (a, b)
                if rho[i] >= rho[i - 1] and rho[i] >= rho[i + 1]:
                    if (u < a and GUf[i + 1] >= a) or (u > b and GUf[i - 1] <= b):
                        continue
                    if rho[i] > 1e-3:
                        cnt += 1
                    if rho[i] > best:
                        best, at = rho[i], u
            print("%3d %2d %3d %3d   %12.4e                                %.3f  | %d" % (2 * d, k, nd["n"], E["factor"], best, at if at else -1, cnt))
