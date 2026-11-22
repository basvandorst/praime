"""Review 26: own class-chain repeat probabilities vs r82c.json (1e12, 1e15) and vs round 90 data (1e7)."""
import json, math
import numpy as np
D7 = json.load(open('/home/bas/praime/writeup/data/study/r90.json'))
RC = json.load(open('/home/bas/praime/writeup/data/study/r82c.json'))
R = {'r82': json.load(open('/home/bas/praime/writeup/data/study/r82.json')), 'r84': json.load(open('/home/bas/praime/writeup/data/study/r84.json'))}
QA = (3, 5, 7, 11, 13, 17, 19, 23)


def rep(d, Ds, q):
    return d[Ds % q == 0].sum() / d.sum()


for name, f, G in (('r82', 'ch_1e12.npy', 320), ('r84', 'ch_1e15.npy', 400)):
    ch = np.load('/home/bas/scratch/study-review26/' + f); Ds = np.arange(2, G + 1, 2)
    mxd = 0; zc = {'lo': [], 'hi': []}; zl = {'lo': [], 'hi': []}
    for k in range(1, 7):
        for q in QA:
            c = rep(ch[k - 1], Ds, q); v = R[name]['data'][str(k)][str(q)]; rc = RC[name]['per_k'][str(k)][str(q)]
            mxd = max(mxd, abs(c - rc['chain']) / v['se'])
            if k >= 2:
                g = 'lo' if q <= 13 else 'hi'
                zc[g].append((c - v['p']) / v['se']); zl[g].append(v['z'])
    print(name, 'own chain vs r82c chain: max |diff|/se %.3f' % mxd)
    for g in ('lo', 'hi'):
        a = np.array(zc[g]); b = np.array(zl[g])
        print('  %s (%d values): chain max|z| %.1f sum z2 %.0f | law max|z| %.2f sum z2 %.1f' % (g, len(a), np.abs(a).max(), (a * a).sum(), np.abs(b).max(), (b * b).sum()))
    a = np.array(zc['lo']).reshape(5, 5); print('  chain z (k=2..6 x q<=13):\n', np.round(a, 1))
    print('  chain z q=17..23:\n', np.round(np.array(zc['hi']).reshape(5, 3), 1))
# 1e7 power check
ch = np.load('/home/bas/scratch/study-review26/ch_1e7.npy'); Ds = np.arange(2, 201, 2); zz = []
for k in range(1, 5):
    row = []
    for q in (3, 5, 7, 11, 13):
        v = D7['data'][str(k)][str(q)]; z = (rep(ch[k - 1], Ds, q) - v['p']) / v['se']; row.append(z)
        if k >= 2:
            zz.append(z)
    print('1e7 k=%d chain z' % k, np.round(row, 1))
zz = np.array(zz); print('1e7 chain k>=2: max|z| %.1f sum z2 %.0f (15 values)' % (np.abs(zz).max(), (zz * zz).sum()))
