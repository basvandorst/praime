"""Nested split of the spreads: level 1 = across the 15 classes mod 210 (14 dof), level 2 = across the 9 subclasses
mod 2310 inside each class mod 210 (120 dof), level 3 = across the 11 subclasses mod 30030 inside each class mod 2310
(1350 dof). Each level: measured SS per dof / mean, against the HL prediction from v210, v2310, v30030."""
import math, json
import numpy as np

R = json.load(open('/home/bas/praime/writeup/data/study/r43.json'))
v210_1e9 = json.load(open('v210_1e9.out'))['v'] if False else 0.44920251720633914
wmap = {'train': 'train', '1e8-1e9': '1e8', '1e9-1e10': '1e9', 'A': 'A', 'B': 'B', 'C': 'C'}


def ss(c, Q, coarse):
    adm = np.array([a for a in range(Q) if math.gcd(a * (a + 2), Q) == 1])
    X = np.bincount(np.arange(30030) % Q, weights=c, minlength=Q)[adm]
    if coarse is None:
        return float(np.sum((X - X.mean()) ** 2)), X.mean()
    _, inv = np.unique(adm % coarse, return_inverse=True)
    gm = np.bincount(inv, weights=X) / np.bincount(inv)
    return float(np.sum((X - gm[inv]) ** 2)), X.mean()


out = {}
acc = {1: [], 2: [], 3: []}
for w, f in wmap.items():
    c = np.load(f'c_{f}.npy').astype(float)
    v210 = R[w]['210']['predicted'] if '210' in R[w] else v210_1e9
    v2310, v30030 = R[w]['2310']['predicted'], R[w]['30030']['predicted']
    s1, m1 = ss(c, 210, None); s2, m2 = ss(c, 2310, 210); s3, m3 = ss(c, 30030, 2310)
    r1 = s1 / 14 / m1 / v210
    r2 = s2 / 120 / m2 / ((134 * v2310 - 14 * v210) / 120)
    r3 = s3 / 1350 / m3 / ((1484 * v30030 - 134 * v2310) / 1350)
    out[w] = dict(level1_210=r1, level2_2310_in_210=r2, level3_30030_in_2310=r3)
    for i, r in ((1, r1), (2, r2), (3, r3)):
        if w in ('train', '1e8-1e9', '1e9-1e10'): acc[i].append(r)
    print(w, {k: round(v, 3) for k, v in out[w].items()})
dof = {1: 14, 2: 120, 3: 1350}
for i in (1, 2, 3):
    m = float(np.mean(acc[i])); se = math.sqrt(2 / (dof[i] - 0)) / math.sqrt(3)
    out[f'long3_level{i}'] = dict(mean=m, se=se, z=(m - 1) / se)
    print('long3 level', i, round(m, 3), '+-', round(se, 3), 'z', round((m - 1) / se, 2))
json.dump(out, open('levels.json', 'w'), indent=1)
