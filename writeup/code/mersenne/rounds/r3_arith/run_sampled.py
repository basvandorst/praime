"""Sampled group (real quadratic class number and regulator: ~0.6 ms per prime): discovery in full,
validation on a uniform sample of 20,000 primes (expected counts scaled by 1/f) plus every validation
positive of every label set (real and synthetic)."""
import os, pickle, time
import numpy as np
import harness as H
import features as Fe

sets = H.label_sets()
rng = np.random.default_rng(4242)
vidx = np.flatnonzero(H.VMASK)
samp = rng.choice(vidx, 20000, replace=False)
f = 20000 / len(vidx)
smask = np.zeros(len(H.P), dtype=bool); smask[samp] = True
need = np.zeros(len(H.P), dtype=bool)
need[H.DMASK] = True; need[samp] = True
for s in sets:
    need[s[H.VMASK[s]]] = True
need[H.REAL] = True
idx = np.flatnonzero(need & (H.P >= 11))
t = time.time()
feats = Fe.g_realquad(H.P[idx])
print(f"computed {len(idx)} primes in {time.time()-t:.0f}s")
res = {}
for nm, vals in feats.items():
    v = np.full(len(H.P), np.nan); v[idx] = vals
    r = H.evaluate(v, sets, sample_mask=smask, sample_frac=f)
    r["real_vals_exps"] = v[H.REAL]
    res[nm] = r
    print("  %-22s %-14s D %d/%.2f  V %d/%.2f  p_v %.3g" % (nm, H.describe(r["classes"][int(r["cls"][0])]), r["Od"][0], r["Ed"][0], r["Ov"][0], r["Ev"][0], r["pv"][0]))
pickle.dump(res, open(os.path.join(H.HERE, "results", "realquad_sampled.pkl"), "wb"))
