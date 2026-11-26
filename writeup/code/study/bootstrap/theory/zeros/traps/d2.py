"""D2: how the trap shrinks with the degree (reads hp_d*.json, traps_d*.json, basevals.json).
For k = 1..6 and every degree with K >= k: trap radius (eta = 2 eta_0) around log n_k, gt''(u_k), u_k - log n_k,
the comparison value sqrt(2 lambda/(mu_k gt''(u_k))) and the ratio radius_S/that (radius_S = half-width of S
around u_k).  Then g_d(log 6), g_d(log 10), g_d(log 12), b_0 per degree.   Usage: python d2.py > d2.out"""
import json, math
from fam import *

BV = json.load(open("basevals.json"))
GL = {int(n): [arb(s) for s in v] for n, v in BV["gl"].items()}
D = list(range(2, 26, 2))
HP = {d: json.load(open("hp_d%d.json" % d)) for d in D}
TR = {d: json.load(open("traps_d%d.json" % d)) for d in D}


def fl(s):
    return float(arb(s).mid().str(20, radius=False)) if isinstance(s, str) else float(s)


print("D2a. Trap at eta = 2 eta_0 (and 10 eta_0) around u_k, as a function of the degree 4K = 2d'.")
print("radS = max(|a-u_k|,|b-u_k|); pred = sqrt(2 lambda/(mu_k g''(u_k))); R = trap radius around log n_k")
for k in range(1, 7):
    print("k = %d" % k)
    print("  deg  n   u_k-log n     g''(u_k)        lambda       radS(x2)     pred         radS/pred  R(x2)        R(x10)")
    for d in D:
        R = HP[d]
        if R["K"] < k:
            continue
        nd = TR[d]["nodes"][k - 1]
        E2, E10 = nd["etas"]
        lam = fl(R["lam"]); mu = fl(R["mu"][k - 1]); gpp = fl(R["gpp"][k - 1])
        pred = math.sqrt(2 * lam / (mu * gpp))
        if "a" in E2:
            radS = max(abs(fl(E2["a_minus_uk"])), abs(fl(E2["b_minus_uk"])))
            r2 = "%.4e" % E2["radius_logn"]
        else:
            radS = float("nan"); r2 = E2["status"]
        r10 = ("%.4e" % E10["radius_logn"]) if "a" in E10 else E10["status"]
        print("  %3d %3d  %11.4e  %.8e  %.4e  %.4e  %.4e  %.6f   %s   %s" % (
            2 * d, nd["n"], fl(nd["u_minus_logn"]), gpp, lam, radS, pred, radS / pred, r2, r10))
print()
print("D2b. Optimum g_d at non-prime-powers, and b_0 (g(0) = 1).")
print(" deg   g(log 6)        g(log 10)       g(log 12)       b_0")
for d in D:
    b = [arb(s) for s in HP[d]["b"]]
    vals = [dot(b, GL[n][:d + 1]) for n in (6, 10, 12)]
    print(" %3d  %.9e  %.9e  %.9e  %.15f" % (2 * d, *[float(v.mid().str(15, radius=False)) for v in vals], fl(HP[d]["b"][0])))
