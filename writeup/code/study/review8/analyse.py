"""Twin (g=2) and g=6 u^2 from own counts; P(2) from the Hardy-Littlewood product (own derivation);
null errors for mean u^2 (per gap and pooled); joint significance of the twin shortfall."""
import json, math
import numpy as np
from scipy import stats, integrate
import sympy

D = '/home/bas/praime/writeup/data/study/'
r40 = json.load(open(D + 'r40.json')); r41 = json.load(open(D + 'r41.json'))
Q = 210

# C' = prod_{q>=11} (1-2/q)/(1-1/q)^2 ; twin constant C2 / local factors at 3,5,7 (tail beyond 1e7 by 1 - 1/q^2 approx)
lp = 0.0
for q in sympy.primerange(11, 10 ** 7):
    lp += math.log1p(-2 / q) - 2 * math.log1p(-1 / q)
lp += -1 / (10 ** 7 * math.log(10 ** 7))   # tail sum_{q>1e7} -1/q^2
Cp = math.exp(lp)
print('C prime (q>=11) =', Cp)


def P2_law(L, W):
    f = lambda u: math.exp(u) / u
    num = integrate.quad(lambda u: f(u) / u, math.log(L), math.log(L + W))[0]
    den = integrate.quad(f, math.log(L), math.log(L + W))[0]
    return (Q / 48) * Cp * num / den


wins = {'train': ('train.npz', r40['train']['per_gap'], 1e7, 9e7),
        'A': ('A.npz', r40['A']['per_gap'], 1e10, 5e7),
        'w89': ('w89.npz', r41['per_gap'], 1e8, 9e8)}
out = {}
pvals = []
for w, (f, per, L, W) in wins.items():
    z = np.load(f); N = z['N'].astype(float); Y = z['Y'].astype(float); U = list(z['U'])
    Yb = Y.mean()
    rows2 = [i for i, a in enumerate(U) if math.gcd(a + 2, Q) == 1]
    X = N[rows2, 0]
    P = P2_law(L, W)
    ratio = X.sum() / (P * Y[rows2].sum())
    u = (X - P * Yb) / np.sqrt(P * Yb)
    u2 = float(np.mean(u * u))
    # data-only minimum over a common P (any common calibration offset removed)
    Pm = X.mean() / Yb
    umin = float(np.mean((X - Pm * Yb) ** 2 / (Pm * Yb)))
    # conditional on Y_a version (round-13 style): (X - P Y_a)^2/(P Y_a)
    ucond = float(np.mean((X - P * Y[rows2]) ** 2 / (P * Y[rows2])))
    pred = per['2']['pred_u2']; meas = per['2']['meas_u2']; n = len(rows2)
    # null: u_a ~ N(0, pred) independent -> n*mean/pred ~ chi2_n (and chi2_{n-1} if one common mode is lost)
    p_n = stats.chi2.cdf(n * meas / pred, n); p_n1 = stats.chi2.cdf(n * meas / pred, n - 1)
    sd_gauss = pred * math.sqrt(2 / n)
    pvals.append(p_n)
    out[w] = dict(P2_law=P, measured_over_law=ratio, u2_mine=u2, u2_repo=meas, u2_min_dataonly=umin,
                  u2_cond_on_Ya=ucond, pred=pred, rows=n, p_lower_chi2_n=p_n, p_lower_chi2_nm1=p_n1,
                  z_gauss=(meas - pred) / sd_gauss, sd_gauss=sd_gauss, eps2PY=(ratio - 1) ** 2 * P * Yb)
    # g = 6: group rows by which interior points a+2, a+4 are units (local pattern); data-only min per group
    rows6 = [i for i, a in enumerate(U) if math.gcd(a + 6, Q) == 1]
    grp = {}
    for i in rows6:
        a = U[i]; grp.setdefault((math.gcd(a + 2, Q) == 1, math.gcd(a + 4, Q) == 1), []).append(i)
    tot = 0.0
    for k, rr in grp.items():
        X6 = N[rr, 2]; Pg = X6.mean() / Yb
        tot += float(np.sum((X6 - Pg * Yb) ** 2 / (Pg * Yb)))
    out[w]['g6_u2_min_dataonly_by_interior_group'] = tot / len(rows6)
    out[w]['g6_groups'] = {str(k): len(v) for k, v in grp.items()}
    out[w]['g6_u2_repo'] = per['6']['meas_u2']; out[w]['g6_pred'] = per['6']['pred_u2']
    out[w]['primes'] = int(Y.sum()); out[w]['twins'] = int(N[:, 0].sum()); out[w]['gap6'] = int(N[:, 2].sum())

    # pooled null SE of mean u^2 over gaps 2..38 (rows independent, Gaussian): sqrt(2 sum n v^2)/sum n
    ns = np.array([per[g]['rows'] for g in per], float); vs = np.array([per[g]['pred_u2'] for g in per])
    out[w]['pooled_null_se_indep'] = float(math.sqrt(2 * np.sum(ns * vs ** 2)) / ns.sum())

# joint: Fisher over the three disjoint windows
chi = -2 * sum(math.log(p) for p in pvals)
out['twins_fisher'] = {'pvals': pvals, 'chi2_6': chi, 'p': float(stats.chi2.sf(chi, 6))}
# Stouffer on Gaussian z
zs = [out[w]['z_gauss'] for w in wins]
out['twins_stouffer_gauss'] = {'z': zs, 'Z': sum(zs) / math.sqrt(3), 'p_one_sided': float(stats.norm.cdf(sum(zs) / math.sqrt(3)))}
json.dump(out, open('analyse.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
