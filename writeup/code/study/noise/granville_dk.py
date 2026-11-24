"""Short-range comparator for round 83's band statistic: a sieve model (independent 'primes' on the residues coprime to
Q = 30030, each with probability (Q/phi(Q))/lam, lam = ln(1.05e9)); this is a first-order chain on the classes mod 30030
with the sieve's own gap law (no Hardy-Littlewood correction beyond 13) and no long-range rigidity.
For D_k = p_{n+k} - p_n (k = 1..4), bands [lo, lo+10], lo = 2, 14, ..., 230 (D <= 240), the Pearson chi^2 of each run
against the pooled band law, its mean over runs (x R/(R-1) for the pooled law), its sd, and the part along the mean
('scale') direction. usage: granville_dk.py runs span seed"""
import sys, math, json
import numpy as np

R = int(sys.argv[1]); SPAN = int(float(sys.argv[2])); SEED = int(sys.argv[3])
LAM = math.log(1.05e9); Q = 30030
adm = np.array([a for a in range(Q) if math.gcd(a, Q) == 1], np.int64)
p = (Q / len(adm)) / LAM
edges = np.arange(2, 241, 12)
rng = np.random.default_rng(SEED)
cnts = {k: [] for k in (1, 2, 3, 4)}
Ns = []; cls30 = []
U30 = [a for a in range(30) if math.gcd(a, 30) == 1]
for r in range(R):
    base = (np.arange(SPAN // Q, dtype=np.int64)[:, None] * Q + adm[None, :]).ravel()
    P = base[rng.random(len(base)) < p]
    Ns.append(len(P))
    cls30.append(np.bincount(P % 30, minlength=30)[U30])
    for k in (1, 2, 3, 4):
        d = P[k:] - P[:-k]; d = d[d <= 240]
        b = (d - 2) // 12
        cnts[k].append(np.bincount(b, minlength=20)[:20])
X30 = np.array(cls30, float); X30 = X30 - X30.mean(1, keepdims=True)   # class contrasts mod 30 (sum zero)
X30 = X30[:, :-1]
out = {'lam': LAM, 'p': p, 'runs': R, 'span': SPAN, 'N_mean': float(np.mean(Ns)), 'N_var_over_mean': float(np.var(Ns, ddof=1) / np.mean(Ns)), 'k': {}}
for k in (1, 2, 3, 4):
    C = np.array(cnts[k], float); tot = C.sum(1)
    Pb = C.sum(0) / C.sum()
    keep = Pb * tot.mean() >= 50                      # bands kept as in round 83 (expected >= 50); here all for k>=2
    E = tot[:, None] * Pb[None, :]
    z = (C - E) / np.sqrt(E)
    chi = (z[:, keep] ** 2).sum(1) * R / (R - 1)
    cen = edges + 5.0; m = (Pb * cen).sum()
    u = np.sqrt(Pb) * (cen - m); u /= np.linalg.norm(u)
    sc = (z @ u) ** 2 * R / (R - 1)
    # chi^2 explained by the class contrasts mod 30 (regression across runs, bias-corrected by permutation)
    def expl(Z, X):
        Zc = Z - Z.mean(0); Xc = X - X.mean(0)
        Cxz = Xc.T @ Zc / (R - 1); Vx = Xc.T @ Xc / (R - 1)
        return float(np.trace(Cxz.T @ np.linalg.solve(Vx, Cxz)))
    e30 = expl(z[:, keep], X30)
    perm = [expl(z[:, keep], X30[rng.permutation(R)]) for _ in range(20)]
    out_e30 = e30 - float(np.mean(perm))
    out['k'][k] = {'chi2_mean': float(chi.mean()), 'chi2_sd': float(chi.std(ddof=1)), 'scale_part_mean': float(sc.mean()),
                   'bands_kept': int(keep.sum()), 'dof': int(keep.sum()) - 1,
                   'explained_by_mod30_contrasts': out_e30}
    print(k, {a: round(b, 3) if isinstance(b, float) else b for a, b in out['k'][k].items()}, flush=True)
print('N mean %.0f, Var/mean of N %.3f (1 - p = %.3f)' % (out['N_mean'], out['N_var_over_mean'], 1 - p))
json.dump(out, open('granville_dk_%d.json' % SEED, 'w'), indent=1)
