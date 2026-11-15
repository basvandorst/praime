#!/usr/bin/env python3
"""Leave-families-out CV on the real labels and on NSIM label sets simulated from the S null
(y ~ Bernoulli(w), w = min(c_f lnQ/lnR, 0.95) with the real c_f). The identical pipeline
(c_f refitted, all configurations, same folds) runs on every label set.
Feature sets: 'all' (every column) and 'compact' (all groups except the residues n mod m).
Output: cv_results.json (real and simulated scores per configuration), oof_real.npz."""
import json, os, sys, time
import numpy as np, scipy.sparse as sp
import models as M

HERE = os.path.dirname(os.path.abspath(__file__))
NSIM = int(sys.argv[1]) if len(sys.argv) > 1 else 100
SEED = 20261116
d = np.load(os.path.join(HERE, "train.npz"))
X = d["X"].astype(np.float64); types = d["types"]; groups = d["groups"]
y, wraw, fam = d["y"], d["wraw"], d["fam"]
fold_of = M.folds(fam)
sets = {"all": np.ones(X.shape[1], bool), "compact": groups != "res"}
mats = {k: (sp.csr_matrix(X[:, m]), sp.csr_matrix(X[:, m & (types == "bin")])) for k, m in sets.items()}
out = {"folds": fold_of, "grids": dict(L2=M.L2_GRID, L1=M.L1_GRID, GB=M.GB_ROUNDS), "real": {}, "sim": []}
t0 = time.time()
oof = {}
for k, (Xs, Xb) in mats.items():
    r, pred, w = M.pipeline(Xs, Xb, y, wraw, fam, fold_of, return_preds=True)
    out["real"][k] = r
    for c, v in pred.items():
        oof[f"{k}|{c}"] = v
np.savez(os.path.join(HERE, "oof_real.npz"), **oof)
print("real", json.dumps(out["real"]), f"{time.time() - t0:.0f}s", flush=True)
w_real = M.s_null(y, wraw, fam)
rng = np.random.default_rng(SEED)
for s in range(NSIM):
    ys = (rng.random(len(y)) < w_real).astype(float)
    rec = {k: M.pipeline(Xs, Xb, ys, wraw, fam, fold_of) for k, (Xs, Xb) in mats.items()}
    rec["hits"] = int(ys.sum())
    out["sim"].append(rec)
    if s % 10 == 9 or s == NSIM - 1:
        with open(os.path.join(HERE, "cv_results.json"), "w") as fh:
            json.dump(out, fh)
        print(f"sim {s + 1}/{NSIM} {time.time() - t0:.0f}s", flush=True)
