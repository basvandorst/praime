"""(a) Empirical noise: 90 disjoint windows of 1e8 in [1e9, 1e10); ratio measured/predicted spread per window,
mod 2310 and mod 30030; sd of ratios against sqrt(2/(K-1)).  (b) within/between split of the mod-30030 spread
for the six windows of round 43, with predictions (1484 v30030 - 134 v2310)/1350 (within) and v2310 (between)."""
import math, json
import numpy as np
import theory as T

segs = np.load('seg_1e9.npy')
out = {}
for Q in (2310, 30030):
    adm = [a for a in range(Q) if math.gcd(a * (a + 2), Q) == 1]; K = len(adm)
    rat = []
    for i in range(segs.shape[0]):
        cq = np.bincount(np.arange(30030) % Q, weights=segs[i], minlength=Q)
        X = cq[adm]
        meas = np.var(X, ddof=1) / X.mean()
        S = T.setup(Q, 10 ** 9 + i * 10 ** 8, 10 ** 8)
        r = S['r']; v = (float(np.sum(r * (1 - r))) + T.F(S, 0)[0]) / float(r.sum())
        rat.append(meas / v)
    rat = np.array(rat)
    out[Q] = dict(K=K, nwin=len(rat), mean_ratio=float(rat.mean()), se_mean=float(rat.std(ddof=1) / math.sqrt(len(rat))),
                  sd_ratio=float(rat.std(ddof=1)), sd_nominal=math.sqrt(2 / (K - 1)),
                  sd_over_nominal=float(rat.std(ddof=1) / math.sqrt(2 / (K - 1))),
                  z_mean_nominal=float((rat.mean() - 1) / (math.sqrt(2 / (K - 1)) / math.sqrt(len(rat)))))
    print(Q, out[Q], flush=True)

# (b)
R = json.load(open('/home/bas/praime/writeup/data/study/r43.json'))
wmap = {'train': 'train', '1e8-1e9': '1e8', '1e9-1e10': '1e9', 'A': 'A', 'B': 'B', 'C': 'C'}
adm = np.array([a for a in range(30030) if math.gcd(a * (a + 2), 30030) == 1])
grp = adm % 2310
ug, inv = np.unique(grp, return_inverse=True)
rows = {}
for w, f in wmap.items():
    c = np.load(f'c_{f}.npy'); X = c[adm].astype(float); mu = X.mean()
    gm = np.bincount(inv, weights=X) / np.bincount(inv)
    within = float(np.sum((X - gm[inv]) ** 2)) / 1350 / mu
    v30, v23 = R[w]['30030']['predicted'], R[w]['2310']['predicted']
    pw = (1484 * v30 - 134 * v23) / 1350
    rows[w] = dict(within_ratio=within / pw, se=math.sqrt(2 / 1349), between_ratio=R[w]['2310']['ratio'], se_b=math.sqrt(2 / 134))
    print(w, {k: round(v, 4) for k, v in rows[w].items()})
for name, ws in (('long3', ['train', '1e8-1e9', '1e9-1e10']), ('all6', list(wmap))):
    wr = np.mean([rows[w]['within_ratio'] for w in ws]); br = np.mean([rows[w]['between_ratio'] for w in ws])
    sw = math.sqrt(2 / 1349) / math.sqrt(len(ws)); sb = math.sqrt(2 / 134) / math.sqrt(len(ws))
    rows[name] = dict(within=wr, within_se=sw, within_z=(wr - 1) / sw, between=br, between_se=sb, between_z=(br - 1) / sb)
    print(name, {k: round(v, 4) for k, v in rows[name].items()})
for name, ws in (('2310_four_upto_1e10', ['train', '1e8-1e9', '1e9-1e10', 'A']),
                 ('2310_all6', list(wmap))):
    br = np.mean([R[w]['2310']['ratio'] for w in ws]); sb = math.sqrt(2 / 134) / math.sqrt(len(ws))
    rows[name] = dict(ratio=br, se=sb, z=(br - 1) / sb); print(name, {k: round(v, 4) for k, v in rows[name].items()})
out['split'] = rows
json.dump(out, open('subwin.json', 'w'), indent=1)
