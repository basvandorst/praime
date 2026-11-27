"""Read the repo scans (after my own numbers were written): at each LMFDB r_j <= 15, the two neighbouring grid points,
t_low there, mean a_p vs truth; also check t*_odd >= t*_even pointwise (parity-dominance theorem) on the scans."""
import json, numpy as np
D = "/home/bas/praime/writeup/data/study/maass/"
key = json.load(open(D + "lmfdb_key.json"))
for X in (8, 16, 32):
    S = {e: sorted([json.loads(l) for l in open(D + "scan_eps%d_X%d_k4_J61.jsonl" % (e, X))], key=lambda d: d["r"]) for e in (0, 1)}
    print("X=%d  rows: even %d odd %d" % (X, len(S[0]), len(S[1])))
    for f in key:
        if f["r"] > 15.5: continue
        e = 1 if f["parity"] == "odd" else 0
        rows = S[e]
        near = sorted(rows, key=lambda d: abs(d["r"] - f["r"]))[:2]
        for d in near:
            dev = [d["mean_ap"][i] - f["ap"][str(p)] for i, p in enumerate(d["primes"][:4])]
            print("  r_j=%.4f %s grid r=%.2f t_low=% .5f  max|a_p - true| (p<=7) = %.3f  (%s)" % (
                f["r"], f["parity"], d["r"], d["t_low"], max(abs(v) for v in dev), " ".join("%+.3f" % v for v in dev)))
    # dominance on common grid points
    ev = {round(d["r"], 4): d["t_low"] for d in S[0]}
    viol = [(d["r"], d["t_low"], ev[round(d["r"], 4)]) for d in S[1] if round(d["r"], 4) in ev and d["t_low"] < ev[round(d["r"], 4)] - 1e-6]
    print("  points with t_odd < t_even - 1e-6: %d (of %d common)%s" % (len(viol), len(ev), (" e.g. %s" % (viol[:3],)) if viol else ""))
    # feasible set components at -0.002
    for e in (0, 1):
        rs = np.array([d["r"] for d in S[e]]); ts = np.array([d["t_low"] for d in S[e]])
        feas = ts >= -0.002
        comps = []; i = 0
        while i < len(rs):
            if feas[i]:
                j = i
                while j + 1 < len(rs) and feas[j + 1]: j += 1
                comps.append((rs[i], rs[j])); i = j + 1
            else: i += 1
        print("  eps=%d feasible components (t_low >= -0.002): %s" % (e, ", ".join("[%.2f,%.2f]" % c for c in comps)))
