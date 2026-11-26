"""Decay table: optimal cost/g(0) versus deg B (class A and E), from polish_*.json (KKT optima) and b1_cost (deg 0).
deg 4k+2 equals deg 4k (the top coefficient is forced to 0: B >= 0 needs b_d >= 0, g >= 0 at u -> oo needs
(-1)^d b_d >= 0, since g(u) ~ B(2 pi i e^u) (Phi*Phi)(u)).
Model check: at deg 4K the K touch points sit just above the first K prime powers (class A) / primes (class E);
we compare log10(cost) with -(n_{K+1} - 1), n_{K+1} = first prime power (A) / prime (E) not touched.
Usage: python decay.py"""
import json, os
import numpy as np

C0 = 0.075197132697534236812          # b1_cost.py (deg 0, both formulas)


def is_prime(p):
    return p > 1 and all(p % q for q in range(2, int(p ** .5) + 1))


def is_pp(n):
    for p in range(2, n + 1):
        if n % p == 0:
            while n % p == 0:
                n //= p
            return n == 1
    return False


PP = [n for n in range(2, 300) if is_pp(n)]
PR = [n for n in range(2, 300) if is_prime(n)]
rows = {}
for cls in "AE":
    rows[cls] = {0: dict(cost=C0, costs="%.20g" % C0, touch=[], share={})}
    for d in range(2, 60, 2):
        f = "polish_%s_d%d.json" % (cls, d)
        if os.path.exists(f):
            r = json.load(open(f))
            if cls == "A":
                share = {int(k): float(v) for k, v in r["cost_share"].items()}
            else:   # by prime: 2 log p S(log p)/cost (per-n shares cancel in class E)
                share = {int(p): 2 * np.log(int(p)) * float(v) / float(r["cost"]) for p, v in r["S_logp"].items()}
                share.update({int(k): float(v) for k, v in r["cost_share"].items() if int(k) >= 17 and is_prime(int(k))})
            rows[cls][2 * d] = dict(cost=float(r["cost"]), costs=r["cost"], touch=[float(x) for x in r["touch_points_exp"]],
                                    share=share,
                                    mu_ok=r["all_mu_nonneg"], Bpos=not r["B_has_real_root"],
                                    res=r["stationarity_residual"], gridmin=r["grid_min_rel"],
                                    arch=r["prime_minus_arch"], lam=r["lam"],
                                    limited=float(r["stationarity_residual"]) > 1e-10 * float(r["cost"]))
for cls, seq in (("A", PP), ("E", PR)):
    print("class %s   (deg 4K+2 = deg 4K)" % cls)
    print("deg  K | cost/g(0) (KKT optimum)   | ratio  | log10   | next n | log10+(n-1) | touch points e^u")
    prev = None
    for deg in sorted(rows[cls]):
        a = rows[cls][deg]
        K = deg // 4
        nxt = seq[K]
        ratio = "" if prev is None else "%7.1f" % (prev / a["cost"])
        print("%3d %2d | %-25s | %6s | %7.3f | %6d | %11.2f | %s" % (
            deg, K, "%.15e" % a["cost"] + (" *" if a.get("limited") else ""), ratio, np.log10(a["cost"]), nxt, np.log10(a["cost"]) + nxt - 1,
            ", ".join("%.10g" % t for t in a["touch"])))
        prev = a["cost"]
    degs = np.array([d_ for d_ in sorted(rows[cls]) if not rows[cls][d_].get("limited")])
    y = np.log10([rows[cls][d]["cost"] for d in degs])
    if len(degs) >= 4:
        ks = degs / 4
        p1, r1, *_ = np.polyfit(ks, y, 1, full=True)
        p2, r2, *_ = np.polyfit(ks, y, 2, full=True)
        print("  fit log10 cost = a + b K       : b = %.3f, rms %.3f" % (p1[0], np.sqrt(r1[0] / len(ks))))
        print("  fit log10 cost = a + b K + c K^2: c = %.4f, b = %.3f, rms %.3f" % (p2[0], p2[1], np.sqrt(r2[0] / len(ks))))
        print("  decrements of log10 cost per +4 degrees:", " ".join("%.2f" % x for x in -np.diff(y)))
        nx = np.array([seq[int(k)] for k in ks])
        p3, r3, *_ = np.polyfit(nx, y, 1, full=True)
        print("  fit log10 cost = a + s n_{K+1} (n_{K+1} = first untouched %s): s = %.3f, a = %.2f, rms %.3f" % (
            "prime power" if cls == "A" else "prime", p3[0], p3[1], np.sqrt(r3[0] / len(ks))))
        h = len(ks) // 2
        print("  local slope (decades per unit n) first half %.2f, second half %.2f" % (
            np.polyfit(nx[:h + 1], y[:h + 1], 1)[0], np.polyfit(nx[h:], y[h:], 1)[0]))
    if any(rows[cls][d_].get("limited") for d_ in rows[cls]):
        print("  * precision-limited (stationarity residual > 1e-10 cost at 90 digits): a feasible point, upper bound only; not in the fits")
    print("  main cost carriers (class A: by n; class E: by prime p):")
    for deg in sorted(rows[cls]):
        sh = sorted(rows[cls][deg]["share"].items(), key=lambda kv: -kv[1])[:4]
        if sh:
            print("    deg %2d: %s" % (deg, "  ".join("n=%d %.1f%%" % (n, 100 * s) for n, s in sh)))
    print("  checks: mu>=0, B>0 (no real root), stationarity |z|, min of F/scale on the u-grid, prime-arch with 24-node moments (see archcheck.out for the 48-node check)")
    for deg in sorted(rows[cls]):
        a = rows[cls][deg]
        if deg:
            print("    deg %2d: %s %s %s %s %s" % (deg, a["mu_ok"], a["Bpos"], a["res"], a["gridmin"], a["arch"]))
    print()
