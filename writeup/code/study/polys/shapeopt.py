"""Hill-climb root configurations of a monic sextic (sum of roots = -B5) to maximise the heuristic log-chance that a
window of L consecutive values |N(y)|/K is all prime (values < 67^2 count as certain; others e^gamma ln 67 / ln v).
usage: shapeopt.py K NSTART SEED OUT"""
import sys, math, json
import numpy as np
K = int(sys.argv[1]); NST = int(sys.argv[2]); rng = np.random.default_rng(int(sys.argv[3])); OUT = sys.argv[4]
B5 = 3 if K == 72 else 0
g = 0.5772156649; kap = math.exp(g) * math.log(67); P2 = 67 * 67; L = 56
yy = np.arange(-70, 71).astype(float)

def roots_of(x, typ):
    nr = 6 - 2 * typ; r = list(x[:nr]); j = nr
    for _ in range(typ): r += [complex(x[j], abs(x[j + 1])), complex(x[j], -abs(x[j + 1]))]; j += 2
    r = np.array(r, complex); r = r - (r.sum() + B5) / 6
    return r

def score(x, typ):
    r = roots_of(x, typ)
    v = np.ones_like(yy, dtype=complex)
    for z in r: v = v * (yy - z)
    v = np.abs(v.real) / K
    lp = np.where(v < P2, 0.0, np.log(np.minimum(1, kap / np.log(np.maximum(v, 3)))))
    c = np.concatenate([[0], np.cumsum(lp)])
    w = c[L:] - c[:-L]
    return w.max() / math.log(10)

res = []
for s in range(NST):
    typ = int(rng.integers(0, 3))
    nr = 6 - 2 * typ
    x = np.concatenate([rng.uniform(-25, 25, nr), np.ravel([[rng.uniform(-25, 25), rng.uniform(0.3, 6)] for _ in range(typ)])])
    best = score(x, typ); step = 3.0
    for it in range(400):
        xn = x + rng.normal(0, step, len(x))
        sn = score(xn, typ)
        if sn > best: x, best = xn, sn
        if it % 100 == 99: step *= 0.6
    r = roots_of(x, typ)
    c = np.real(np.poly(r))           # c[0] = 1, c[1] = b5, ...
    res.append({'score': best, 'type': typ, 'roots': [[z.real, z.imag] for z in r],
                'b': [int(round(t)) for t in c[2:6]]})
res.sort(key=lambda d: -d['score'])
json.dump(res, open(OUT, 'w'))
print('best scores:', [round(d['score'], 2) for d in res[:10]], 'types', [d['type'] for d in res[:10]])
print('top shape:', res[0]['b'], [tuple(round(v, 2) for v in z) for z in res[0]['roots']])
