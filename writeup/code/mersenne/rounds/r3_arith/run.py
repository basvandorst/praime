"""Compute one feature group on the primes < 1e7 and score every feature with the frozen harness."""
import sys, time, pickle, os
import numpy as np
import harness as H
import features as Fe

grp = sys.argv[1]
t = time.time()
feats = Fe.GROUPS[grp](H.P.copy())
tf = time.time() - t
sets = H.label_sets()
res = {}
for nm, v in feats.items():
    v = np.asarray(v, dtype=np.float64)
    assert len(v) == len(H.P), nm
    r = H.evaluate(v, sets)
    if r is None:
        print("degenerate", nm); continue
    r["real_vals_exps"] = v[H.REAL]
    res[nm] = r
os.makedirs(os.path.join(H.HERE, "results"), exist_ok=True)
pickle.dump(res, open(os.path.join(H.HERE, "results", grp + ".pkl"), "wb"))
best = sorted(res.items(), key=lambda kv: kv[1]["pv"][0])[:5]
print(f"{grp}: {len(res)} features, compute {tf:.0f}s, total {time.time()-t:.0f}s")
for nm, r in best:
    print("  %-28s %-14s D %d/%.2f  V %d/%.2f  p_v %.3g bits %.2f" % (nm, H.describe(r["classes"][int(r["cls"][0])]),
          r["Od"][0], r["Ed"][0], r["Ov"][0], r["Ev"][0], r["pv"][0], r["bits"][0]))
