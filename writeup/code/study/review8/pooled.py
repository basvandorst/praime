"""(a) calibration-free version of round 40: within groups of rows whose law share P_a(g) is equal (same unit pattern of
a+2..a+g mod 210), subtract the group mean; compare with v_g (n_grp - 1). (b) Monte Carlo of the sd of the pooled mean u^2
with rows that share gaps (multinomial per row), against the independent-row formula."""
import json, math
import numpy as np

D = '/home/bas/praime/writeup/data/study/'
r40 = json.load(open(D + 'r40.json')); r41 = json.load(open(D + 'r41.json'))
Q = 210; GAPS = list(range(2, 39, 2))
wins = {'train': ('train.npz', r40['train']['per_gap']), 'A': ('A.npz', r40['A']['per_gap']), 'w89': ('w89.npz', r41['per_gap'])}
res = {}
rng = np.random.default_rng(1)
for w, (f, per) in wins.items():
    z = np.load(f); N = z['N'].astype(float); Y = z['Y'].astype(float); U = list(z['U']); Yb = Y.mean()
    S = 0.0; E = 0.0; V = 0.0; repo_ok = []
    perg = {}
    for g in GAPS:
        rows = [i for i, a in enumerate(U) if math.gcd(a + g, Q) == 1]
        grp = {}
        for i in rows:
            a = U[i]; key = tuple(math.gcd(a + j, Q) == 1 for j in range(2, g, 2))
            grp.setdefault(key, []).append(i)
        v = per[str(g)]['pred_u2']; s = 0.0; dof = 0
        for rr in grp.values():
            if len(rr) < 2:
                continue
            X = N[rr, g // 2 - 1]; Pg = X.mean() / Yb
            uu = (X - Pg * Yb) / math.sqrt(Pg * Yb); s += float(np.sum((uu - uu.mean()) ** 2)); dof += len(rr) - 1
        perg[g] = (len(rows), len(grp), dof, s / dof if dof else None, v)
        S += s; E += v * dof; V += 2 * v * v * dof
    res[w] = {'within_group_sum': S, 'expected': E, 'ratio_meas_over_pred': S / E, 'z': (S - E) / math.sqrt(V),
              'dof_total': sum(p[2] for p in perg.values()),
              'per_gap(n,groups,dof,u2_within,pred)': {str(k): v for k, v in perg.items()}}
    # (b) Monte Carlo: calibrated multinomial rows, vY from the repo, P = empirical shares (so mean u^2 ~ 1-P+P vY)
    vY = r41['v2']['meas'] if w == 'w89' else r40[w]['summary']['v2']['meas']
    ph = N / Y[:, None]
    ns = []; vs = []
    for g in GAPS:
        rows = [i for i, a in enumerate(U) if math.gcd(a + g, Q) == 1]
        ns.append(len(rows)); vs.append(float(np.mean(1 - ph[rows, g // 2 - 1] + ph[rows, g // 2 - 1] * vY)))
    ns = np.array(ns); vs = np.array(vs)
    sims = []
    for _ in range(1500):
        Ys = np.maximum(1, np.round(rng.normal(Yb, math.sqrt(vY * Yb), 48))).astype(np.int64)
        Xs = np.array([rng.multinomial(Ys[i], ph[i]) for i in range(48)], float)
        Ybs = Ys.mean(); tot = 0.0
        for g in GAPS:
            rows = [i for i, a in enumerate(U) if math.gcd(a + g, Q) == 1]
            P = ph[rows, g // 2 - 1]; X = Xs[rows, g // 2 - 1]
            tot += float(np.sum((X - P * Ybs) ** 2 / (P * Ybs)))
        sims.append(tot / ns.sum())
    sims = np.array(sims)
    res[w]['mc_null_mean'] = float(sims.mean()); res[w]['mc_expected'] = float(np.sum(ns * vs) / ns.sum())
    res[w]['mc_sd'] = float(sims.std()); res[w]['indep_formula_sd'] = float(math.sqrt(2 * np.sum(ns * vs ** 2)) / ns.sum())
    print(w, {k: v for k, v in res[w].items() if not k.startswith('per_gap')}, flush=True)
json.dump(res, open('pooled.json', 'w'), indent=1)
