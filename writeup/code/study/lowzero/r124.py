"""Round 124 (preregister/round124.txt): scaled first zero of L(s, chi_d) for the fundamental discriminants of the
quadratic search's long runs (>= 36), as a percentile among 150 random fundamental discriminants of similar size.
usage: r124.py"""
import os, json, math, random
import numpy as np, cypari2
pari = cypari2.Pari(); pari.allocatemem(900 * 10 ** 6)
H = os.path.dirname(os.path.abspath(__file__))
polys = []
for line in open(os.path.join(H, '..', 'polys', 'dist300.txt')):
    A, B, C, s, L = map(int, line.split())
    if L >= 36: polys.append((A, B, C, L))
def g1(d):
    Lf = pari.lfuncreate(d); T = 2.0; z = []
    while not z: z = [float(t) for t in pari.lfunzeros(Lf, T) if float(t) > 0]; T *= 2
    return min(z) * math.log(abs(d)) / (2 * math.pi)
rng = random.Random(124); out = {}
for A, B, C, L in polys:
    d = int(pari.coredisc(B * B - 4 * A * C))
    if d in out: out[d]['runs'].append(L); continue
    gd = g1(d); ref = []
    lo, hi = int(abs(d) / 1.5), int(abs(d) * 1.5); sg = 1 if d > 0 else -1
    while len(ref) < 150:
        x = sg * rng.randint(lo, hi)
        if x != d and int(pari.isfundamental(x)): ref.append(g1(x))
    pct = float(np.mean(np.array(ref) < gd))
    out[d] = {'runs': [L], 'g': gd, 'percentile': pct, 'ref_median': float(np.median(ref)), 'L1': float(pari.lfun(pari.lfuncreate(d), 1))}
    print(d, 'run', L, 'g %.4f' % gd, 'percentile %.3f' % pct, flush=True)
pcts = [v['percentile'] for v in out.values()]
r1 = float(np.median(pcts)) <= 0.20; r2 = out.get(7213, {}).get('percentile', 1) < 0.05; r3 = np.mean(np.array(pcts) < 0.5) >= 0.75
print('n', len(out), 'median percentile %.3f, share below median %.2f; readings %s %s %s' % (np.median(pcts), np.mean(np.array(pcts) < 0.5), r1, r2, r3))
json.dump({'discs': {str(k): v for k, v in out.items()}, 'median_percentile': float(np.median(pcts)), 'readings': [bool(r1), bool(r2), bool(r3)]}, open(os.path.join(H, 'r124.json'), 'w'), indent=1)
