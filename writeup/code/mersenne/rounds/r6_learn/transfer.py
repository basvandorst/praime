#!/usr/bin/env python3
"""Transfer to the Mersenne exponents.

Final models: for each model type (L2, L1, GB) on the 'all' feature set, the configuration with the
largest real cross-validated gain (cv_results.json), refitted on all 79 families.
Mersenne baseline: w(k) = s(k) * min(e^gamma log2(2Kk)/k, 0.95), K = 10^4, s = trial-factoring survival.
Model: logit p = logit w + score(k).  Population: the survivors among the primes 2 <= k < B (weighted sample).
Statistics: mean percentile of the 51 exponents' scores in the population (ties half), and the gain in bits
sum [y log2(p/w) + (W - y) log2((1-p)/(1-w))]. Null: NSIM label sets drawn from w on the population
(Bernoulli where W = 1, Poisson(W w) in the sampled octaves), scored by the same fixed model.
The 52nd exponent 136,279,841: percentile of its score. Output: transfer.txt, final_models.npz."""
import json, math, os, sys
import numpy as np, scipy.sparse as sp
import models as M

HERE = os.path.dirname(os.path.abspath(__file__))
EG = math.exp(0.5772156649015329)
K = 10000
NSIM = 4000
SEED = 777


def pct(pop_s, pop_W, s):
    o = np.argsort(pop_s)
    ss, ww = pop_s[o], pop_W[o]
    cw = np.concatenate(([0.0], np.cumsum(ww)))
    lo = np.searchsorted(ss, s, "left"); hi = np.searchsorted(ss, s, "right")
    return 100 * (cw[lo] + 0.5 * (cw[hi] - cw[lo])) / cw[-1]


def main():
    d = np.load(os.path.join(HERE, "train.npz"))
    names = list(d["names"]); types = d["types"]
    X = d["X"].astype(np.float64); y, wraw, fam = d["y"], d["wraw"], d["fam"]
    cv = json.load(open(os.path.join(HERE, "cv_results.json")))
    real = cv["real"]["all"]
    w = M.s_null(y, wraw, fam); act = w > 0
    off = np.full(len(y), -np.inf); off[act] = M.logit(w[act])
    Xs = sp.csr_matrix(X[act]); bincols = np.flatnonzero(types == "bin")
    Xb = sp.csr_matrix(X[act][:, bincols])
    me = np.load(os.path.join(HERE, "mersenne.npz"))
    XM = me["X"].astype(np.float64); k = me["k"]; W = me["W"]; st = me["stratum"]; surv = me["surv"]
    base = np.where(surv, np.minimum(EG * np.log2(2.0 * K * k) / k, 0.95), 0.0)
    popm = (st != -2) & surv
    isexp = st == -1
    newi = np.flatnonzero(st == -2)[0]
    lines = []
    P = lambda s: (print(s), lines.append(s))
    P(f"population: {int(round(W[st != -2].sum()))} primes (sample rows {int((st != -2).sum())}), "
      f"survivors (weighted) {W[popm].sum():.0f}; exponents {isexp.sum()} all survive {surv[isexp].all()}")
    P(f"baseline expected count sum W w = {np.sum(W[popm] * base[popm]):.2f} (observed 51)")
    scores, chosen, saved = {}, {}, {}
    for typ in ("L2", "L1", "GB"):
        cands = {c: v[0] for c, v in real.items() if c.startswith(typ)}
        c = max(cands, key=cands.get); chosen[typ] = c
        par = float(c.split("=")[1])
        if typ == "L2":
            th = M.fit_l2(Xs, y[act], off[act], par)
            s = th[0] + XM @ th[1:]; saved["L2"] = th
        elif typ == "L1":
            th = M.fit_l1(Xs, y[act], off[act], par)
            s = th[0] + XM @ th[1:]; saved["L1"] = th
        else:
            trees, _ = M.fit_gb(Xb, y[act], off[act], int(par))
            s = M.gb_apply(trees, sp.csr_matrix(XM[:, bincols]), int(par))
            import pickle
            pickle.dump(dict(trees=trees, rounds=int(par), bincols=bincols), open(os.path.join(HERE, "final_gb.pkl"), "wb"))
        scores[typ] = s
    np.savez(os.path.join(HERE, "final_models.npz"), **saved, names=np.array(names), bincols=bincols,
             chosen=json.dumps(chosen), mscore_L2=scores["L2"], mscore_L1=scores["L1"], mscore_GB=scores["GB"])
    for typ in ("L1", "L2"):
        th = saved[typ]
        o = np.argsort(-np.abs(th[1:]))[:12]
        P(f"{typ} ({chosen[typ]}) intercept {th[0]:+.3f}, nonzero {int(np.sum(np.abs(th[1:]) > 1e-8))}; largest: " +
          ", ".join(f"{names[j]} {th[1 + j]:+.3f}" for j in o))
    rng = np.random.default_rng(SEED)
    pop = np.flatnonzero(popm)
    Wp, wp = W[pop], base[pop]
    bern = (st[pop] == 15) | (st[pop] == -1)
    rows_, cols_, vals_ = [], [], []
    for i in range(NSIM):
        c = np.zeros(len(pop))
        c[bern] = rng.random(int(bern.sum())) < wp[bern]
        c[~bern] = rng.poisson(Wp[~bern] * wp[~bern])
        nz = np.flatnonzero(c)
        rows_.append(np.full(len(nz), i)); cols_.append(nz); vals_.append(c[nz])
    sims_y = sp.csr_matrix((np.concatenate(vals_), (np.concatenate(rows_), np.concatenate(cols_))),
                           shape=(NSIM, len(pop)))
    res = {}
    for typ in ("L2", "L1", "GB"):
        s = scores[typ]
        sp_ = s[pop]
        p = M.sigmoid(M.logit(wp) + sp_)
        pctl = pct(sp_, Wp, sp_)
        yreal = isexp[pop].astype(float)
        t1 = np.log2(p) - np.log2(wp); t0 = np.log2(1 - p) - np.log2(1 - wp)
        g_real = float(np.sum(yreal * t1 + (Wp - yreal) * t0))
        mp_real = float(pctl[yreal == 1].mean())
        g_sim = sims_y @ (t1 - t0) + np.sum(Wp * t0)
        cnt = np.asarray(sims_y.sum(1)).ravel()
        mp_sim = (sims_y @ pctl) / np.maximum(cnt, 1)
        p_g = (1 + np.sum(g_sim >= g_real)) / (NSIM + 1)
        p_m = (1 + np.sum(mp_sim >= mp_real)) / (NSIM + 1)
        new_pct = float(pct(sp_, Wp, np.array([s[newi]]))[0])
        res[typ] = dict(g_real=g_real, mp_real=mp_real, p_g=p_g, p_m=p_m, g_sim=g_sim, mp_sim=mp_sim)
        ex_p = np.sort(pctl[yreal == 1])
        P(f"\n{typ} [{chosen[typ]}]: score sd in population {np.sqrt(np.average((sp_ - np.average(sp_, weights=Wp)) ** 2, weights=Wp)):.4f}")
        P(f"  mean percentile of the 51 exponents {mp_real:.1f} (null {mp_sim.mean():.1f} +- {mp_sim.std():.1f}), p {p_m:.3f}")
        P(f"  gain {g_real:+.3f} bits (null mean {g_sim.mean():+.3f}, sd {g_sim.std():.3f}), p {p_g:.3f}")
        P(f"  exponents above the 90th percentile: {int(np.sum(ex_p > 90))} (5.1 expected), below the 10th: {int(np.sum(ex_p < 10))}")
        P(f"  52nd exponent 136,279,841: score {s[newi]:+.4f}, percentile {new_pct:.1f}")
        P("  exponent percentiles: " + " ".join(f"{int(kk)}:{v:.0f}" for kk, v in zip(k[pop][yreal == 1], pctl[yreal == 1])))
    # family-wise over the three models: min p, calibrated on the same simulations
    for stat in ("g", "m"):
        rp = min(res[t][f"p_{stat}"] for t in res)
        ps = []
        for t in res:
            arr = res[t][f"{'g' if stat == 'g' else 'mp'}_sim"]
            srt = np.sort(arr)
            ps.append((1 + NSIM - np.searchsorted(srt, arr, "left")) / (NSIM + 1))
        simp = np.min(ps, axis=0)
        fw = (1 + np.sum(np.array(simp) <= rp)) / (NSIM + 1)
        P(f"family-wise over L2/L1/GB ({'gain' if stat == 'g' else 'mean percentile'}): smallest p {rp:.3f}, family-wise {fw:.3f}")
    open(os.path.join(HERE, "transfer.txt"), "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
