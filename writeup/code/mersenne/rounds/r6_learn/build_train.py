#!/usr/bin/env python3
"""Training table: every sieve survivor n (log2 R >= 64) of the 79 round-5 families.
Rows with s(n) = 0 have S-null probability 0 and label 0: they add nothing to any likelihood,
so only survivors are kept. Saves train.npz (CSR features, labels, raw S weight lnQ/lnR, family)."""
import json, math, os, sys, time
import numpy as np, scipy.sparse as sp
import feats as F
import r5_features as FT

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "r5_transfer", "data.jsonl")

t0 = time.time()
rows = []
for fi, line in enumerate(open(DATA)):
    d = json.loads(line)
    b, kind, K = d["b"], d["kind"], d["K"]
    g = b - 1 if kind == "m" else b + 1
    recs = [r for r in d["recs"] if r[1] != "ndiv"]
    for r in recs:
        n = r[0]
        log2R = n * math.log2(b) - math.log2(g)
        if log2R < 64 or r[1] not in ("PRP", "fermat"):
            continue
        assert r[1] != "mr"
        rows.append((fi, kind, b, n, 1 if r[1] == "PRP" else 0, math.log(2 * K * n + 1) / (log2R * math.log(2))))
fam = np.array([r[0] for r in rows]); kinds = [r[1] for r in rows]
bs = np.array([r[2] for r in rows]); ns = np.array([r[3] for r in rows], np.int64)
y = np.array([r[4] for r in rows], float); wraw = np.array([r[5] for r in rows])
famname = {}
blocks = []
for fi in np.unique(fam):
    m = fam == fi
    kind, b = kinds[np.flatnonzero(m)[0]], int(bs[m][0])
    famname[int(fi)] = f"{kind}{b}"
    X, names, groups, types = F.features(kind, b, ns[m])
    # the vectorised H5/H6 must equal round 5's definitions
    assert np.array_equal(F.glide_vec(ns[m]), FT.collatz_glide(ns[m]).astype(np.int64))
    assert np.allclose(F.cf_vec(ns[m]), FT.cf_maxq_norm(ns[m]))
    blocks.append((np.flatnonzero(m), X))
X = np.zeros((len(rows), blocks[0][1].shape[1]))
for idx, B in blocks:
    X[idx] = B
X, names, groups, types, thr = F.add_thresholds(X, names, groups, types)
# standardisation of continuous columns (unsupervised; stored and reused for Mersenne)
mu = np.zeros(X.shape[1]); sd = np.ones(X.shape[1])
for j, t in enumerate(types):
    if t == "cont":
        mu[j] = X[:, j].mean(); sd[j] = X[:, j].std() or 1.0
Xs = (X - mu) / sd
const = np.flatnonzero(Xs.std(axis=0) == 0)
print("constant columns dropped:", [names[j] for j in const])
keep = np.setdiff1d(np.arange(X.shape[1]), const)
Xs = Xs[:, keep]; names = [names[j] for j in keep]; groups = [groups[j] for j in keep]
types = [types[j] for j in keep]; mu = mu[keep]; sd = sd[keep]
np.savez_compressed(os.path.join(HERE, "train.npz"), X=sp.csr_matrix(Xs).toarray().astype(np.float32),
                    y=y, wraw=wraw, fam=fam, n=ns, b=bs, kind=np.array(kinds), names=np.array(names),
                    groups=np.array(groups), types=np.array(types), mu=mu, sd=sd,
                    famname=json.dumps(famname), thr=json.dumps(thr))
from collections import Counter
print(f"rows {len(y)} hits {int(y.sum())} families {len(famname)} columns {Xs.shape[1]}",
      dict(Counter(groups)), dict(Counter(types)), f"nnz/row {np.count_nonzero(Xs) / len(y):.0f}",
      f"{time.time() - t0:.0f}s")
