"""Measured spreads Var/mean (ddof 1) per modulus from the recount; within/between decomposition mod 30030."""
import math, json
import numpy as np

wins = ['train', '1e8', '1e9', 'A', 'B', 'C']
res = {}
for w in wins:
    c = np.load(f'c_{w}.npy')
    res[w] = {}
    for Q in (210, 2310, 30030):
        cq = np.bincount(np.arange(30030) % Q, weights=c, minlength=Q)
        adm = [a for a in range(Q) if math.gcd(a * (a + 2), Q) == 1]
        bad = sum(cq[a] for a in range(Q) if a not in set(adm))
        X = cq[adm]
        res[w][Q] = float(np.var(X, ddof=1) / X.mean())
        res[w][f'n{Q}'] = len(adm); res[w][f'outside{Q}'] = float(bad)
    # within/between mod 30030 split by the class mod 2310
    adm = np.array([a for a in range(30030) if math.gcd(a * (a + 2), 30030) == 1])
    X = c[adm].astype(float); mu = X.mean()
    grp = adm % 2310
    ug, inv = np.unique(grp, return_inverse=True)
    gm = np.bincount(inv, weights=X) / np.bincount(inv)
    within = float(np.sum((X - gm[inv]) ** 2)); between = float(np.sum(np.bincount(inv) * (gm - mu) ** 2))
    res[w]['within_per_dof'] = within / (len(adm) - len(ug)) / mu
    res[w]['between_dof'] = len(ug) - 1
    res[w]['check_between_vs_2310'] = between * 11 / (len(ug) - 1) / (11 * mu) / res[w][2310]
    print(w, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in res[w].items()})
