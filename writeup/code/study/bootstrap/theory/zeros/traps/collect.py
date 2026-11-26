"""Collect: traps.out (D1, all degrees), d3_summary.out (D3 certified radii) and the W(h) cross-check
(ball W(h) of trapcert + DELTA * W(Phi*Phi)  vs  lambda h_k(0) + mu_k from the re-polish)."""
import json, os, glob
from mpmath import mp, mpf
mp.dps = 50
with open("traps.out", "w") as f:
    f.write("# D1, one line per (degree, k). Columns: n = prime power nearest e^{u_k}; u-logn = u_k - log n; mu = KKT weight;\n"
            "# 2L/sqn = 2 Lambda(n)/sqrt(n); g'' = gt''(u_k); W = W(h_k); kkt = W(h_k) - lambda h_k(0) - mu_k; h0 = h_k(0);\n"
            "# eta0 = lambda/W(h_k); eta1/eta0 (inf: B_h < 0 on R); for eta = 2 eta0 (x2) and 10 eta0 (x10): R = trap radius around log n,\n"
            "# lb/mu = weight lower bound / mu_k, extra = number of extra components of S found on the grid, rhomax = max eta h/gt on\n"
            "# the grid away from the main component, rho_j = max_j eta h''(u_j)/gt''(u_j) over the other nodes, logs[m] = log m in S;\n"
            "# verdict: isolates = S one interval containing no log m (m != n); single-but-wide = one interval containing another log m.\n")
    for d in range(2, 26, 2):
        if os.path.exists("traps_d%d.out" % d):
            f.write(open("traps_d%d.out" % d).read())
WPP = mpf("0.0240444194343550508370662176054")   # W(Phi*Phi) (certs README, B = 1)
H = json.load(open("hp_d24.json"))
lines = ["# D3 (degree 48), rigorous: S c (a_k, b_k); R_k >= max |a_k - log n|, |b_k - log n|; weight lower bound;",
         "# Wcheck = ball W(h) + DELTA W(Phi*Phi) - (lambda h_k(0) + mu_k) (independent code paths; rounding of h to SIGH digits included)"]
for k in range(1, 13):
    fn = "trapcert_k%d.json" % k
    if not os.path.exists(fn):
        lines.append("k %2d: no certificate" % k)
        continue
    C = json.load(open(fn))
    hk = H["h"][k - 1]
    Wb = mpf(C["W_h"].strip("[").split("+/-")[0])
    chk = Wb + (lambda f: mpf(f.numerator) / f.denominator)(__import__("fractions").Fraction(C["DELTA"])) * WPP - (mpf(H["lam"]) * mpf(hk["h0"]) + mpf(H["mu"][k - 1]))
    wfn = "trapweight_k%d.json" % k
    wl = json.load(open(wfn))["weight_lower"] if os.path.exists(wfn) else "(crude) " + C["weight_lower"]
    lines.append("k %2d n %2d eta'/eta0 %.4f  a-log n %s  b-log n %s  R_k <= %s  weight >= %s  (mu_k %s)  intervals %d  Wcheck %.1e  time %.0fs" % (
        k, C["n_k"], C["eta_prime_over_eta0"], C["a_minus_logn"], C["b_minus_logn"], C["R_k_upper"], wl,
        C["mu_k"][:12], C["S_intervals"], chk, C["time_s"]))
open("d3_summary.out", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
