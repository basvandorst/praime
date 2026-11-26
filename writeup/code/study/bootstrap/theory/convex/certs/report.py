"""Human-readable dump of every polished optimum (polish_<CLASS>_d<d>.json): B coefficients and roots, touch points,
g(log n)/sqrt(n) for n = 2..12, cost carriers, dual multipliers.  Usage: python report.py > report.txt"""
import json, glob, math


def key(f):
    return (f.split("_")[1], int(f.split("_d")[1][:-5]))


for f in sorted(glob.glob("polish_*_d*.json"), key=key):
    r = json.load(open(f))
    print("=" * 100)
    print("class %s   deg B = %d   cost/g(0) = %s   (KKT value lambda = %s, |stationarity| = %s)" %
          (r["cls"], r["deg"], r["cost"], r["lam"], r["stationarity_residual"]))
    print("  arch-side cost (24-node moments) = %s ; prime - arch = %s" % (r["cost_arch"], r["prime_minus_arch"]))
    print("  touch points u_k (double zeros of %s):" % ("g" if r["cls"] == "A" else "S"))
    for u, e, mu in zip(r["touch_points"], r["touch_points_exp"], r["mu"]):
        print("     u = %s   e^u = %s   dual weight mu = %s" % (u[:22], e, mu))
    print("  all mu >= 0: %s ; B has a real root: %s ; min of F/scale on the u-grid (i/400): %s" %
          (r["all_mu_nonneg"], r["B_has_real_root"], r["grid_min_rel"]))
    print("  grid local minima (u, F/scale):", r["grid_local_minima"])
    print("  B(t) = sum_j b_j t^(2j), normalised g(0) = 1:")
    for j, b in enumerate(r["b"]):
        print("     b_%-2d = %s" % (j, b[:28]))
    print("  roots of B in t (principal sqrt of each root of q(s) = B(sqrt s); all roots are +- these):", ", ".join(r["roots_t"]))
    print("  g(log n)/sqrt(n):", "  ".join("%s: %s" % (k, v) for k, v in r["g_logn_over_sqrtn"].items()))
    if "S_logp" in r:
        print("  S(log p):", "  ".join("%s: %s" % (k, v) for k, v in r["S_logp"].items()))
    if r["cls"] == "A":
        sh = sorted(((int(k), float(v)) for k, v in r["cost_share"].items()), key=lambda kv: -abs(kv[1]))[:6]
        print("  share of cost by n (2 Lambda(n) n^-1/2 g(log n) / cost):", "  ".join("%d: %.4g" % kv for kv in sh))
    else:   # class E: cost = 2 sum_p log p S(log p); per-n shares cancel (g(log 2) < 0 < g(log 4)), so use S
        c = float(r["cost"])
        shp = {int(p): 2 * math.log(int(p)) * float(v) / c for p, v in r["S_logp"].items()}
        for k, v in r["cost_share"].items():
            if int(k) >= 17 and all(int(k) % q for q in range(2, int(k))):
                shp[int(k)] = float(v)     # p^2 > 40 terms are below the stored precision
        sh = sorted(shp.items(), key=lambda kv: -abs(kv[1]))[:6]
        print("  share of cost by prime p (2 log p S(log p) / cost):", "  ".join("%d: %.4g" % kv for kv in sh))
