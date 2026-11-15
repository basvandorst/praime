"""Final exam on 1e7 <= p < B for the 5 pre-registered features (frozen classes from prereg.json).
Expected counts from a uniform random sample of NSAMP final-range primes, scaled by 1/f; the 13 exponents
are scored exactly."""
import json, os, sys, time
import numpy as np
import harness as H
import features as Fe

NSAMP = 200000
pre = json.load(open(os.path.join(H.HERE, "prereg.json")))
ap = Fe.all_primes()
F = ap[(ap >= H.S2) & (ap < H.B)]
rng = np.random.default_rng(777)
samp = np.sort(rng.choice(F, NSAMP, replace=False))
f = NSAMP / len(F)
ex = np.array([e for e in H.EXPS if H.S2 <= e < H.B], dtype=np.int64)
k = np.union1d(samp, ex)
in_s = np.isin(k, samp)
pos = np.flatnonzero(np.isin(k, ex))
w = H.null_w(k)
print(f"final range: {len(F)} primes, sample {NSAMP} (f = {f:.5f}), {len(ex)} exponents; E_total(all F) = {H.null_w(F).sum():.3f}")


def gp_load_list(kk):
    g = Fe.pari()
    g("P=[" + ",".join(str(int(x)) for x in kk) + "]")
    return g


Fe.gp_load = gp_load_list
cache = {}
for t in pre["top5"]:
    g = t["group"]
    if g not in cache:
        t0 = time.time()
        cache[g] = Fe.GROUPS[g](k.copy())
        print(f"  computed group {g} on {len(k)} primes in {time.time()-t0:.0f}s")
    v = np.asarray(cache[g][t["feature"]], dtype=np.float64)
    cls = (t["class"][0], t["class"][1])
    st = H.class_stats(v[in_s], w[in_s], [cls], 1.0 / f)
    sf = H.survival(H.pmf_from_sums(st[:, 1], st[:, 2:]))
    O = int(H.in_class(v[pos], cls).sum())
    E = st[0, 0]
    print("%-28s %-14s D %d/%.2f  V %d/%.2f (p %.3g)  FINAL %d/%.2f  p_final %.3g" % (
        t["feature"], t["class_text"], t["D"][0], t["D"][1], t["V"][0], t["V"][1], t["V"][2], O, E, sf[0, min(O, H.K)]))
