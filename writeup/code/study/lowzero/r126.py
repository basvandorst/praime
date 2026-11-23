"""Round 126 (preregister/round126.txt): small-prime filter, then rank by the zero. usage: r125.py"""
import os, json, math, random
import numpy as np, cypari2, gmpy2
from scipy.stats import mannwhitneyu
pari = cypari2.Pari(); pari.allocatemem(900 * 10 ** 6)
H = os.path.dirname(os.path.abspath(__file__))
rng = random.Random(126)
ds = set()
while len(ds) < 600:
    x = rng.randint(10 ** 5, 10 ** 6)
    if x % 8 == 5 and all(int(pari.kronecker(x, p)) == -1 for p in (3, 5, 7, 11)) and int(pari.isfundamental(x)): ds.add(x)
ds = sorted(ds)
def g1(d):
    Lf = pari.lfuncreate(d); T = 2.0; z = []
    while not z: z = [float(t) for t in pari.lfunzeros(Lf, T) if float(t) > 0]; T *= 2
    return min(z) * math.log(d) / (2 * math.pi)
G = {d: g1(d) for d in ds}
order = sorted(ds, key=lambda d: G[d]); low = order[:15]; rnd = order[-15:]
def best_run(d):
    best = 0; Y = range(-80, 81)
    for A in range(1, 301):
        for B in range(0, 2 * A):
            if (B * B - d) % (4 * A): continue
            C = (B * B - d) // (4 * A); cur = 0
            for y in Y:
                v = abs(A * y * y + B * y + C)
                if v > 1 and gmpy2.is_prime(v):
                    cur += 1
                    if cur > best: best = cur
                else: cur = 0
    return best
out = {'low': [], 'random': []}  # 'random' holds the HIGH group here
for name, grp in (('low', low), ('random', rnd)):
    for d in grp:
        r = best_run(d); out[name].append({'d': d, 'g': G[d], 'best_run': r}); print(name, d, round(G[d], 3), r, flush=True)
a = [x['best_run'] for x in out['low']]; b = [x['best_run'] for x in out['random']]
p = mannwhitneyu(a, b, alternative='greater').pvalue
r1 = np.median(a) - np.median(b) >= 4; r2 = p < 0.01; r3 = np.median(a) > 6
print('median LOW %.1f, HIGH %.1f, Mann-Whitney p %.2g; readings %s %s %s' % (np.median(a), np.median(b), p, r1, r2, r3))
out.update({'median_low': float(np.median(a)), 'median_random': float(np.median(b)), 'p': float(p), 'readings': [bool(r1), bool(r2), bool(r3)], 'g_all_median': float(np.median(list(G.values())))})
json.dump(out, open(os.path.join(H, 'r126.json'), 'w'), indent=1)
