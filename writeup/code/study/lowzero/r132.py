"""Round 132 (preregister/round132.txt): best run on ln p13 and S for the 600 filtered discriminants of round 126.
usage: r132.py PART NPARTS | combine"""
import os, sys, json, math, random
import numpy as np, cypari2, gmpy2
pari = cypari2.Pari(); pari.allocatemem(800 * 10 ** 6)
H = os.path.dirname(os.path.abspath(__file__))
def best_run(d):
    best = 0
    for A in range(1, 301):
        for B in range(0, 2 * A):
            if (B * B - d) % (4 * A): continue
            C = (B * B - d) // (4 * A); cur = 0
            for y in range(-80, 81):
                v = abs(A * y * y + B * y + C)
                if v > 1 and gmpy2.is_prime(v):
                    cur += 1
                    if cur > best: best = cur
                else: cur = 0
    return best
if sys.argv[1] != 'combine':
    part, nparts = int(sys.argv[1]), int(sys.argv[2])
    rng = random.Random(126); ds = set()
    while len(ds) < 600:
        x = rng.randint(10 ** 5, 10 ** 6)
        if x % 8 == 5 and all(int(pari.kronecker(x, p)) == -1 for p in (3, 5, 7, 11)) and int(pari.isfundamental(x)): ds.add(x)
    rows = []
    for d in sorted(ds)[part::nparts]:
        p13 = next(int(p) for p in pari.primes(200) if int(p) >= 13 and int(pari.kronecker(d, p)) != -1)
        L = pari.lfuncreate(d); S = 2 * float(pari.lfun(L, 1, 1)) / float(pari.lfun(L, 1)) + math.log(d / math.pi) + float(pari.psi(0.5))
        rows.append({'d': d, 'p13': p13, 'S': S, 'run': best_run(d)})
    json.dump(rows, open(os.path.join(H, 'r132_part%d.json' % part), 'w'))
else:
    rows = sum((json.load(open(os.path.join(H, 'r132_part%d.json' % i))) for i in range(2)), [])
    y = np.array([r['run'] for r in rows], float)
    z = lambda a: (a - a.mean()) / a.std()
    x1 = z(np.log([r['p13'] for r in rows])); x2 = z(np.array([r['S'] for r in rows]))
    def fit(X):
        X = np.column_stack([np.ones(len(y))] + X); c = np.linalg.lstsq(X, y, rcond=None)[0]; res = y - X @ c
        s2 = res @ res / (len(y) - X.shape[1]); se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X))); r2 = 1 - res @ res / np.sum((y - y.mean()) ** 2)
        return c, se, r2
    c, se, r2 = fit([x1, x2]); _, _, r2a = fit([x1]); _, _, r2b = fit([x2])
    print('n %d; ln p13: %+.3f +- %.3f (t %.1f); S: %+.3f +- %.3f (t %.1f); R2 p13 %.3f, S %.3f, both %.3f' % (len(y), c[1], se[1], c[1] / se[1], c[2], se[2], c[2] / se[2], r2a, r2b, r2))
    rd = [bool(c[1] / se[1] > 3), bool(c[2] / se[2] > 3)]
    print('readings', rd)
    json.dump({'n': len(y), 'coef': [float(v) for v in c], 'se': [float(v) for v in se], 'R2': {'p13': r2a, 'S': r2b, 'both': r2}, 'readings': rd}, open(os.path.join(H, 'r132.json'), 'w'), indent=1)
