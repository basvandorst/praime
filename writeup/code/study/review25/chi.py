import json, math, numpy as np
from scipy import stats
r = np.load('rlaw.npz'); out = r['out']; ds = r['d']; lams = r['lams']
cnt = json.load(open('counts.json'))
u = np.linspace(lams[0], lams[-1], 401); wt = np.exp(u) / u ** 2
def avg(A):   # A: (3, nd)
    return np.array([np.sum(np.interp(u, lams, A[:, j]) * wt) / wt.sum() for j in range(A.shape[1])])
BANDS = [(lo, lo + 24) for lo in range(6, 600, 30)] + [(lo, lo + 54) for lo in range(606, 1200, 60)]
assert sorted(sum([list(range(lo, hi + 1, 6)) for lo, hi in BANDS], [])) == list(range(6, 1201, 6))
def chi2(c, law, minexp=50):
    law = law / law.sum(); N = sum(c.values()); x = 0; dof = 0; zs = []
    for lo, hi in BANDS:
        o = sum(c.get(str(d), 0) for d in range(lo, hi + 1, 6)); e = N * law[(ds >= lo) & (ds <= hi)].sum()
        if e >= minexp:
            x += (o - e) ** 2 / e; dof += 1; zs.append((o - e) / math.sqrt(e))
    return x, dof - 1, max(map(abs, zs)), stats.chi2.sf(x, dof - 1)
A = {k: avg(out[k - 1]) for k in (1, 2, 3)}
one = A[1] / A[1].sum()
for k in (1, 2, 3):
    c = cnt[f'{k}_in']
    print('k', k, 'pairs', sum(c.values()), 'mass', round(A[k].sum(), 4), 'law chi2 %.2f dof %d maxz %.2f p %.3f' % chi2(c, A[k]))
    print('   ext-edge version: %.2f %d %.2f %.3f' % chi2(cnt[f'{k}_ext'], A[k]))
    print('   single mid lam: %.2f %d %.2f %.3f' % chi2(c, out[k - 1, 1]), ' lam0: %.2f' % chi2(c, out[k - 1, 0])[0], ' lam2: %.2f' % chi2(c, out[k - 1, 2])[0])
    # independent comparator by explicit enumeration (own indexing): distribution of the sum of k iid gaps
    ind = np.zeros(len(ds))
    if k == 1: ind = one.copy()
    elif k == 2:
        for i, a in enumerate(one):
            for j, b in enumerate(one):
                s = ds[i] + ds[j]
                if s <= 1200: ind[s // 6 - 1] += a * b
    else:
        two = np.zeros(400)
        for i, a in enumerate(one):
            for j, b in enumerate(one): two[(ds[i] + ds[j]) // 6 - 1] += a * b
        for i, a in enumerate(two):
            for j, b in enumerate(one):
                s = 6 * (i + 1) + ds[j]
                if s <= 1200: ind[s // 6 - 1] += a * b
    print('   independent: %.1f %d %.2f' % chi2(c, ind)[:3])
    # data-driven independent comparator: shuffle actual consecutive gaps
