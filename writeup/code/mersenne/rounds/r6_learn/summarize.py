#!/usr/bin/env python3
"""Empirical p-values of the real cross-validated scores against the simulated S-null label sets.
Statistic per model type: the best configuration of its grid (max gain in bits; max dAUC), the same
maximum taken on every simulated set. Also the best over all configurations. Output: cv_report.txt."""
import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
cv = json.load(open(os.path.join(HERE, "cv_results.json")))
sims = cv["sim"]; N = len(sims)
lines = []
P = lambda s: (print(s), lines.append(s))
P(f"simulated label sets: {N}; hits per set {np.mean([s['hits'] for s in sims]):.1f} +- {np.std([s['hits'] for s in sims]):.1f} (real 280)")
P(f"S-null within-family AUC (hit-weighted, survivors): {cv['real']['all']['_base_auc'][1]:.4f}")
for fs in ("all", "compact"):
    P(f"\n== feature set '{fs}' ==")
    real = cv["real"][fs]
    for c, (g, a) in real.items():
        if c.startswith("_"):
            continue
        null = np.array([s[fs][c][0] for s in sims])
        P(f"  {c:16s} gain {g:+8.2f} bits (null median {np.median(null):+8.2f}, 90% {np.quantile(null, 0.9):+8.2f})  dAUC {a:+.4f}")
    for typ in ("L2", "L1", "GB", ""):
        cs = [c for c in real if c.startswith(typ) and not c.startswith("_")]
        for idx, nm in ((0, "gain"), (1, "dAUC")):
            r = max(real[c][idx] for c in cs)
            nl = np.array([max(s[fs][c][idx] for c in cs) for s in sims])
            p = (1 + np.sum(nl >= r)) / (N + 1)
            P(f"  best {typ or 'any':3s} {nm:5s}: real {r:+9.4f}  null median {np.median(nl):+9.4f} "
              f"[5%, 95%] [{np.quantile(nl, 0.05):+.4f}, {np.quantile(nl, 0.95):+.4f}]  p {p:.3f}")
cs = [(fs, c) for fs in ("all", "compact") for c in cv["real"][fs] if not c.startswith("_")]
r = max(cv["real"][fs][c][0] for fs, c in cs)
nl = np.array([max(s[fs][c][0] for fs, c in cs) for s in sims])
P(f"\nbest gain over both feature sets and all {len(cs)} configurations: real {r:+.3f}, p {(1 + np.sum(nl >= r)) / (N + 1):.3f}")
open(os.path.join(HERE, "cv_report.txt"), "w").write("\n".join(lines) + "\n")
