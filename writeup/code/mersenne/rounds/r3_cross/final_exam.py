"""Final exam on [1e7, B) for the 5 pre-registered features (prereg.json, sha256 checked).
Thresholds and directions are read from prereg.json; nothing is re-chosen here."""
import os, sys, json, hashlib, math, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from scipy.stats import binom
from common3 import EXPS, B, V_HI, D_HI, wnull, primes_upto, pb_tail_from_powersums
import features as FT

HERE = os.path.dirname(os.path.abspath(__file__))
s = open(os.path.join(HERE, "prereg.json")).read()
h = open(os.path.join(HERE, "prereg.sha256")).read().split()[0]
assert hashlib.sha256(s.encode()).hexdigest() == h, "prereg.json changed after registration"
pre = json.loads(s)
T0 = time.time()
P = primes_upto(B + 60000, "B").astype(np.int64)
zl = open("/home/bas/praime/writeup/data/riemann/zeros_hp2419.txt").read().split()
ctx = FT.Ctx(P, B + 10, np.array([float(z) for z in zl])); ctx.zeros_str = zl
i0, i1 = int(np.searchsorted(P, V_HI)), int(np.searchsorted(P, B))
idx_f = np.arange(i0, i1); nf = len(idx_f)
w = wnull(P[idx_f]); r = w / (1 - w); l1 = np.log1p(-w)
isexp = np.isin(P[idx_f], np.array(EXPS))
Otot = int(isexp.sum()); Eall = w.sum()
print(f"final range [1e7, B): {nf} primes, {Otot} exponents, E = {Eall:.3f} (O/E {Otot/Eall:.2f})")
rng = np.random.default_rng(4242)
pis = np.minimum(1.0, (60000 / Eall) * w)
samp = np.flatnonzero((rng.random(nf) < pis) | isexp)
lines = []
group_of = {}
for g in FT.GROUPS:
    group_of[g[0]] = g
for f in pre:
    # find the group that yields this feature name (probe on 3 points)
    grp = None
    for gname, dom, mode, fn in FT.GROUPS:
        if dom != f["domain"] or mode != f["mode"]:
            continue
        names = [n for n, _, _ in fn(FT.Ctx(P[:200000], 3 * 10**5, ctx.zeros), np.arange(5000, 5003))] if False else None
        grp = (gname, mode, fn)
        try:
            ctx_small = FT.Ctx(P[:300000], 4 * 10**6, ctx.zeros); ctx_small.zeros_str = zl
            if any(n == f["name"] for n, _, _ in fn(ctx_small, np.arange(5000, 5003))):
                break
        except Exception:
            pass
        grp = None
    assert grp is not None, f["name"]
    gname, mode, fn = grp
    t = np.array(f["thresholds"]); j = f["chosen_j"]
    pts = idx_f if mode == "full" else idx_f[samp]
    vals = np.empty(len(pts))
    CH = 400000 if mode == "full" else 3000
    for s0 in range(0, len(pts), CH):
        for n, kind, v in fn(ctx, pts[s0:s0 + CH]):
            if n == f["name"]:
                vals[s0:s0 + CH] = v
                break
    b = np.searchsorted(t, vals, side="left"); defined = ~np.isnan(vals)
    inc = defined & ((b <= j) if f["direction"] == "<=" else (b > j))
    if mode == "full":
        u = np.ones(nf); pos = np.arange(nf)
    else:
        u = 1.0 / pis[samp]; pos = samp
    E = float((u * w[pos] * inc).sum())
    c0 = float((u * l1[pos] * inc).sum()); S = [float((u * r[pos] ** k * inc).sum()) for k in range(1, 6)]
    O = int((isexp[pos] & inc).sum())
    p = pb_tail_from_powersums(c0, S, 60)[O]
    pc = binom.sf(O - 1, Otot, E / Eall)
    exps_in = [int(x) for x in P[idx_f[pos]][isexp[pos] & inc]]
    ln = (f"{f['name']} [{f['domain']}] class {f['direction']} {f['threshold']:.6g} | disc {f['disc_O']}/{f['disc_E']:.2f} | "
          f"val {f['val_O']}/{f['val_E']:.2f} p {f['val_p']:.2e} | FINAL {O}/{E:.2f} p {p:.3e} pcond(given {Otot}) {pc:.3f} | in class: {exps_in}")
    print(ln, flush=True); lines.append(ln)
    ctx._cache.clear()
open(os.path.join(HERE, "final_exam.txt"), "w").write(
    f"final range [1e7,B): {nf} primes, {Otot} exponents, E = {Eall:.3f}\n" + "\n".join(lines) + f"\n({time.time()-T0:.0f}s)\n")
