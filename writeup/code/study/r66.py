"""Round 66 (preregister/round66.txt): gap entropy at 1e18 from round 17's primes against the asymptote and the full
law at lam = 41.4 (r57.json). usage: r66.py [path]   (default: ../../data/study/r17_primes.txt, one prime per line)"""
import sys, json, math, os
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(OUT, 'r17_primes.txt')
P = [int(x) for x in open(src).read().split()] if src.endswith('.txt') else [int(x) for x in np.load(src, allow_pickle=True)]
g = np.array([b - a for a, b in zip(P[:-1], P[1:])], dtype=np.int64)
v, c = np.unique(g, return_counts=True); N = int(c.sum()); p = c / N
H = float(-np.sum(p * np.log2(p)) + (len(c) - 1) / (2 * N * math.log(2)))
se = float(np.sqrt(np.sum(p * np.log2(p) ** 2) - np.sum(p * np.log2(p)) ** 2) / math.sqrt(N))
lam = math.log(P[0] + (P[-1] - P[0]) / 2)
r57 = json.load(open(os.path.join(OUT, 'r57.json')))['41.4']
K = json.load(open(os.path.join(OUT, 'r52.json')))['K_nats']
asym = (math.log(math.e * lam) - K) / math.log(2)
out = {'primes': len(P), 'gaps': N, 'first': P[0], 'lambda': lam, 'H_data': H, 'se': se, 'H_law_lam41.4': r57['H_law'],
       'H_asym': asym, 'data_minus_law': H - r57['H_law'], 'data_minus_asym': H - asym, 'z_law': (H - r57['H_law']) / se}
print({k: (round(x, 4) if isinstance(x, float) else x) for k, x in out.items()})
json.dump(out, open(os.path.join(OUT, 'r66.json'), 'w'), indent=1)
