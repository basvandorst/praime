"""Frozen harness: discovery (odd p < 1e5) chooses threshold+direction, validation [1e5, 1e7)
scores, NS synthetic null label sets calibrate the family-wise bar. Never touches p >= 1e7.

usage: harness.py NS [groups...]   (writes out/*.npz, results.txt, prereg.json)
"""
import os, sys, json, time, math, hashlib
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from common3 import EXPS, D_HI, V_HI, wnull, primes_upto, pb_tail_from_powersums
import features as FT

T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out"); os.makedirs(OUT, exist_ok=True)
NS = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
ONLY = sys.argv[2:]
SAMPLE_TARGET = 60000
KPOW = 5
LD = 45      # discovery pmf length
LV = 60      # validation pmf length


def log(*a):
    print(f"[{time.time()-T0:7.1f}s]", *a, flush=True)


# ---------------- data ----------------
P = primes_upto(3 * 10**7, "3e7").astype(np.int64)
zl = open("/home/bas/praime/writeup/data/riemann/zeros_hp2419.txt").read().split()
zeros = np.array([float(z) for z in zl])
ctx = FT.Ctx(P, V_HI, zeros); ctx.zeros_str = zl
i_d0, i_d1 = 1, int(np.searchsorted(P, D_HI))          # odd primes 3..<1e5
i_v1 = int(np.searchsorted(P, V_HI))
idx_d = np.arange(i_d0, i_d1); idx_v = np.arange(i_d1, i_v1)
nd, nv = len(idx_d), len(idx_v)
w_d = wnull(P[idx_d]); w_v = wnull(P[idx_v])
isexp_all = np.isin(P, np.array(EXPS))
lab_d_real = isexp_all[idx_d]; lab_v_real = np.flatnonzero(isexp_all[idx_v])     # positions in idx_v
log(f"discovery primes {nd}, exps {lab_d_real.sum()}, E {w_d.sum():.2f}; validation primes {nv}, exps {len(lab_v_real)}, E {w_v.sum():.3f}")

# ---------------- null label sets (seeded, drawn once) ----------------
rng = np.random.default_rng(20261115)
lab_d = np.zeros((nd, NS + 1), dtype=bool); lab_d[:, 0] = lab_d_real
lab_d[:, 1:] = rng.random((nd, NS)) < w_d[:, None]
wmax = w_v.max()
vsets = [lab_v_real]
for s in range(NS):
    g = np.cumsum(rng.geometric(wmax, size=int(nv * wmax * 1.6 + 200))) - 1
    assert g[-1] >= nv
    g = g[g < nv]
    keep = rng.random(len(g)) < w_v[g] / wmax
    vsets.append(g[keep])
U = np.unique(np.concatenate(vsets))                    # positions in idx_v that carry a label somewhere
v_set_id = np.concatenate([np.full(len(v), s) for s, v in enumerate(vsets)])
v_upos = np.searchsorted(U, np.concatenate(vsets))
d_pairs_i, d_pairs_s = np.nonzero(lab_d)
Otot_v = np.array([len(v) for v in vsets])
log(f"null sets {NS}: disc labels mean {lab_d[:,1:].sum(0).mean():.2f}, val labels mean {Otot_v[1:].mean():.2f}; |U| = {len(U)}")

# importance sample of validation for 'sample' mode: pi = min(1, c w)
c = SAMPLE_TARGET / w_v.sum()
pis = np.minimum(1.0, c * w_v)
srng = np.random.default_rng(777)
samp = np.flatnonzero(srng.random(nv) < pis)
log(f"importance sample: {len(samp)} validation primes (c = {c:.0f})")

# ---------------- feature pass ----------------
meta = []; DB = []; VB_U = []; PS = []          # per feature
r_v = w_v / (1 - w_v); l1w_v = np.log1p(-w_v)


def add_feature(name, domain, kind, mode, vd, vval, vpos, weights, vU):
    """vd: values at discovery; vval: values at the E-set of validation (positions vpos, HT weights)."""
    if kind == "bin":
        t = np.array([0.5])
    else:
        x = vd[~np.isnan(vd)]
        t = np.unique(np.quantile(x, np.arange(1, 10) / 10.0)) if len(x) else np.array([])
        # drop thresholds that leave a side empty on discovery
        t = t[(t >= x.min()) & (t < x.max())] if len(x) else t
    if len(t) == 0:
        meta.append(dict(name=name, domain=domain, kind=kind, mode=mode, skip=True)); return

    def binit(v):
        b = np.searchsorted(t, v, side="left").astype(np.int8)
        b[np.isnan(v)] = -1
        return b
    bd = binit(vd)
    if len(np.unique(bd[bd >= 0])) < 2:
        meta.append(dict(name=name, domain=domain, kind=kind, mode=mode, skip=True)); return
    bv = binit(vval)
    nb = len(t) + 1
    ok = bv >= 0
    ww = weights[ok]; pp = vpos[ok]; bb = bv[ok]
    ps = np.zeros((KPOW + 2, nb))
    ps[0] = np.bincount(bb, weights=ww * w_v[pp], minlength=nb)
    ps[1] = np.bincount(bb, weights=ww * l1w_v[pp], minlength=nb)
    for k in range(1, KPOW + 1):
        ps[k + 1] = np.bincount(bb, weights=ww * r_v[pp] ** k, minlength=nb)
    meta.append(dict(name=name, domain=domain, kind=kind, mode=mode, skip=False, thr=t.tolist()))
    DB.append(bd); VB_U.append(binit(vU)); PS.append(ps)


groups = [g for g in FT.GROUPS if not ONLY or g[0] in ONLY]
# planted features (tests of the harness, label-dependent by construction)
def planted(ctx, idx):
    k = ctx.P[idx]
    u = np.array([int.from_bytes(hashlib.sha256(b"u" + str(int(v)).encode()).digest()[:6], "big") / 2**48 for v in k])
    h = np.array([hashlib.sha256(b"h" + str(int(v)).encode()).digest()[0] & 1 for v in k])
    e = isexp_all[idx]
    yield ("PLANT strong: exponents in bottom 5%", "cont", np.where(e, 0.05 * u, u))
    yield ("PLANT medium: half the exponents in bottom 10%", "cont", np.where(e & (h == 1), 0.1 * u, u))
    yield ("PLANT weak: exponents in bottom 40%", "cont", np.where(e, 0.4 * u, u))
    yield ("PLANT discovery-only", "cont", np.where(e & (k < D_HI), 0.05 * u, u))
    yield ("PLANT none (pure noise)", "cont", u)
groups = groups + [("planted", "PLANTED", "full", planted)]

for gname, domain, mode, fn in groups:
    tg = time.time()
    if mode == "full":
        pts = np.concatenate([idx_d, idx_v])
        vpos = np.arange(nv); weights = np.ones(nv)
    else:
        vs = np.union1d(samp, U)
        pts = np.concatenate([idx_d, idx_v[vs]])
        insamp = np.isin(vs, samp)
        vpos_all = vs
    n_before = len(meta)
    for name, kind, vals in fn(ctx, pts):
        vals = np.asarray(vals, dtype=np.float64)
        vd = vals[:nd]; vv = vals[nd:]
        if mode == "full":
            add_feature(name, domain, kind, mode, vd, vv, vpos, weights, vv[U])
        else:
            vU = vv[np.searchsorted(vpos_all, U)]
            add_feature(name, domain, kind, mode, vd, vv[insamp], vpos_all[insamp], 1.0 / pis[vpos_all[insamp]], vU)
    log(f"group {gname}: {len(meta)-n_before} features in {time.time()-tg:.1f}s")
    ctx._cache.clear()

# ---------------- candidates ----------------
F = [i for i, m in enumerate(meta) if not m["skip"]]
cand_f, cand_j, cand_dir = [], [], []
for fi, mi in enumerate(F):
    nT = len(meta[mi]["thr"])
    for j in range(nT):
        for dr in (0, 1):
            cand_f.append(fi); cand_j.append(j); cand_dir.append(dr)
cand_f = np.array(cand_f); cand_j = np.array(cand_j); cand_dir = np.array(cand_dir)
NC = len(cand_f)
DBm = np.stack(DB)                                  # [nF, nd] int8
log(f"{len(F)} usable features ({len(meta)-len(F)} degenerate), {NC} candidate classes")

# discovery pmf per candidate: exact DP over the discovery primes
pmf = np.zeros((NC, LD + 1)); pmf[:, 0] = 1.0
for i in range(nd):
    b = DBm[cand_f, i]
    mem = (b >= 0) & np.where(cand_dir == 0, b <= cand_j, b > cand_j)
    wi = w_d[i] * mem
    sh = pmf[:, :-1] * wi[:, None]
    pmf *= (1 - wi)[:, None]
    pmf[:, 1:] += sh
tail_d = np.cumsum(pmf[:, ::-1], axis=1)[:, ::-1]
E_d = np.zeros(NC)
for i in range(nd):
    b = DBm[cand_f, i]
    E_d += w_d[i] * ((b >= 0) & np.where(cand_dir == 0, b <= cand_j, b > cand_j))
log("discovery DP done")

# validation: per-candidate power sums -> exact Poisson-binomial tails
tail_v = np.zeros((NC, LV + 1)); E_v = np.zeros(NC)
for c_ in range(NC):
    ps = PS[cand_f[c_]]
    nb = ps.shape[1]
    sel = np.arange(nb) <= cand_j[c_] if cand_dir[c_] == 0 else np.arange(nb) > cand_j[c_]
    s = ps[:, sel].sum(1)
    E_v[c_] = s[0]
    tail_v[c_] = pb_tail_from_powersums(s[1], s[2:], LV)
log("validation tails done")
Eall_v = w_v.sum()

# ---------------- scoring all label sets ----------------
nF = len(F)
NSET = NS + 1
res = {k: np.zeros((nF, NSET)) for k in ("cand", "Od", "pd", "Ov", "Ev", "pv", "bits", "pcond", "Ed")}
from scipy.stats import binom
for fi in range(nF):
    nb = len(meta[F[fi]]["thr"]) + 1
    cs = np.flatnonzero(cand_f == fi)
    # counts per (set, bin) in discovery and validation
    bdis = DBm[fi][d_pairs_i]
    cd = np.zeros((NSET, nb + 1)); np.add.at(cd, (d_pairs_s, bdis + 1), 1)
    bval = VB_U[fi][v_upos]
    cv = np.zeros((NSET, nb + 1)); np.add.at(cv, (v_set_id, bval + 1), 1)
    cd = cd[:, 1:]; cv = cv[:, 1:]                   # drop NaN bin
    cdl = np.cumsum(cd, 1); cvl = np.cumsum(cv, 1)
    Od = np.where(cand_dir[cs][None, :] == 0, cdl[:, cand_j[cs]], cdl[:, -1:] - cdl[:, cand_j[cs]]).astype(int)
    Ov = np.where(cand_dir[cs][None, :] == 0, cvl[:, cand_j[cs]], cvl[:, -1:] - cvl[:, cand_j[cs]]).astype(int)
    pdm = tail_d[cs[None, :], np.minimum(Od, LD)]
    best = np.argmin(pdm, axis=1)                    # first minimum = deterministic tie-break
    cb = cs[best]
    ov = Ov[np.arange(NSET), best]; od = Od[np.arange(NSET), best]
    res["cand"][fi] = cb; res["Od"][fi] = od; res["Ed"][fi] = E_d[cb]
    res["pd"][fi] = pdm[np.arange(NSET), best]
    res["Ov"][fi] = ov; res["Ev"][fi] = E_v[cb]
    res["pv"][fi] = tail_v[cb, np.minimum(ov, LV)]
    # bits: rates fitted on discovery (p > 7 only; +0.5 smoothing), scored on validation (Poisson form)
    det = (P[idx_d] <= 7)
    od_inf = od - 3 * 0  # p<=7 contribute 1 to O and E in either class; remove them below
    b37 = DBm[fi][:3]
    m37 = np.array([((b37 >= 0) & ((b37 <= cand_j[c_]) if cand_dir[c_] == 0 else (b37 > cand_j[c_]))).sum() for c_ in cb])
    lam_in = (od - m37 + 0.5) / (E_d[cb] - m37 + 0.5)
    Ed_tot = w_d.sum() - 3; Od_tot = lab_d[:, :].sum(0) - 3
    lam_out = (Od_tot - (od - m37) + 0.5) / (Ed_tot - (E_d[cb] - m37) + 0.5)
    ev = E_v[cb]; eo = Eall_v - ev; oo = Otot_v - ov
    res["bits"][fi] = (ov * np.log2(lam_in) - (lam_in - 1) * ev / math.log(2) + oo * np.log2(lam_out) - (lam_out - 1) * eo / math.log(2))
    res["pcond"][fi] = binom.sf(ov - 1, Otot_v, np.minimum(ev / Eall_v, 1.0))
log("scoring done")

names = [meta[m]["name"] for m in F]; doms = [meta[m]["domain"] for m in F]
is_plant = np.array([d == "PLANTED" for d in doms])
np.savez_compressed(os.path.join(OUT, "scores.npz"), **res, names=np.array(names), doms=np.array(doms),
                    cand_f=cand_f, cand_j=cand_j, cand_dir=cand_dir, E_v=E_v, E_d=E_d)
json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"))

# ---------------- report ----------------
L = []
def out(s=""):
    print(s); L.append(s)
real = ~is_plant
pv = res["pv"]; bits = res["bits"]
out(f"Harness: discovery odd p<1e5 ({nd} primes, {int(lab_d_real.sum())} exponents, E={w_d.sum():.2f}); validation [1e5,1e7) ({nv} primes, {len(lab_v_real)} exponents, E={Eall_v:.3f}); final [1e7,B) untouched.")
out(f"Features: {len(meta)} generated, {nF} usable ({int(real.sum())} real + {int(is_plant.sum())} planted); candidate classes {NC}; null sets {NS}")
from collections import Counter
out("Per domain (usable): " + ", ".join(f"{d}: {n}" for d, n in sorted(Counter(np.array(doms)[real]).items())))
out("Degenerate (skipped): " + "; ".join(m["name"] for m in meta if m["skip"]))
out("")
minp_null = pv[real][:, 1:].min(0); minp_real = pv[real][:, 0].min()
maxb_null = bits[real][:, 1:].max(0); maxb_real = bits[real][:, 0].max()
q = np.quantile(minp_null, [0.01, 0.05, 0.10, 0.5])
out(f"Family-wise (best validation p over {int(real.sum())} real features, {NS} null sets): 1% {q[0]:.2e}, 5% {q[1]:.2e}, 10% {q[2]:.2e}, median {q[3]:.2e}")
out(f"Real best validation p = {minp_real:.3e}; fraction of null sets with a smaller best p = {(minp_null <= minp_real).mean():.3f}")
qb = np.quantile(maxb_null, [0.5, 0.9, 0.95, 0.99])
out(f"Best bits gain: null median {qb[0]:.2f}, 90% {qb[1]:.2f}, 95% {qb[2]:.2f}, 99% {qb[3]:.2f}; real {maxb_real:.2f} (frac null >= real {(maxb_null >= maxb_real).mean():.3f})")
# per-feature calibration: share of null sets with pv <= 0.05
sh = (pv[real][:, 1:] <= 0.05).mean(1)
out(f"Per-feature null calibration: share of null sets with p_v <= 0.05: mean {sh.mean():.4f}, max {sh.max():.4f} (should be <= ~0.05)")
out("")
def row(fi):
    cb = int(res["cand"][fi, 0]); m = meta[F[fi]]
    thr = m["thr"][cand_j[cb]]
    cl = ("<=" if cand_dir[cb] == 0 else ">") + f" {thr:.6g}" if m["kind"] == "cont" else ("F=0" if cand_dir[cb] == 0 else "F=1")
    return (f"{names[fi][:62]:62s} {doms[fi]:10s} {cl:16s} disc {int(res['Od'][fi,0]):2d}/{res['Ed'][fi,0]:5.2f}  "
            f"val {int(res['Ov'][fi,0]):2d}/{res['Ev'][fi,0]:5.2f}  p_v {res['pv'][fi,0]:.2e}  bits {res['bits'][fi,0]:+6.2f}  pcond {res['pcond'][fi,0]:.3f}")
out("Planted tests (same pipeline; family-wise bar from real features):")
for fi in np.flatnonzero(is_plant):
    out("  " + row(fi) + f"  beats 5% bar: {res['pv'][fi,0] < q[1]}")
out("")
order = np.argsort(pv[:, 0])
order = [fi for fi in order if real[fi]]
out("All real features by validation p (class chosen on discovery; O/E in each range):")
for r_, fi in enumerate(order):
    out(f"{r_+1:3d} " + row(fi))
open(os.path.join(HERE, "results.txt"), "w").write("\n".join(L) + "\n")

# pre-registration of the top 5 (frozen before the final exam is computed)
pre = []
for fi in order[:5]:
    cb = int(res["cand"][fi, 0]); m = meta[F[fi]]
    pre.append(dict(name=names[fi], domain=doms[fi], kind=m["kind"], mode=m["mode"], thresholds=m["thr"],
                    chosen_j=int(cand_j[cb]), direction=("<=" if cand_dir[cb] == 0 else ">"),
                    threshold=m["thr"][cand_j[cb]], disc_O=int(res["Od"][fi, 0]), disc_E=float(res["Ed"][fi, 0]),
                    val_O=int(res["Ov"][fi, 0]), val_E=float(res["Ev"][fi, 0]), val_p=float(res["pv"][fi, 0]),
                    val_pcond=float(res["pcond"][fi, 0]), bits=float(res["bits"][fi, 0])))
s = json.dumps(pre, indent=1)
open(os.path.join(HERE, "prereg.json"), "w").write(s)
open(os.path.join(HERE, "prereg.sha256"), "w").write(hashlib.sha256(s.encode()).hexdigest() + f"  prereg.json  {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
log("done")
