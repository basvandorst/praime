"""Round 110 (preregister/round110.txt): unfolded nearest-neighbour spacing variance of the zeros of zeta and of six
L(s, chi_D), in bins of L = ln(q t / 2 pi), with block-bootstrap errors; GUE (Gaudin) 0.180. usage: r110.py"""
import os, json, math
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(H, '..', '..', '..', 'data', 'riemann')
sets = {'zeta': (1, np.array([float(l) for l in open(os.path.join(DATA, 'zeros_10000_flint.txt')) if l.strip()]))}
for D in (-3, -4, 5, -163, 1009, -3635):
    sets['chi%d' % D] = (abs(D), np.array([float(l) for l in open(os.path.join(DATA, 'lzeros', 'D%d.txt' % D)) if l.strip()]))
rng = np.random.default_rng(110); out = {}
for name, (q, g) in sets.items():
    g = np.sort(g)
    N = g / (2 * math.pi) * np.log(q * g / (2 * math.pi * math.e))
    s = np.diff(N); mid = 0.5 * (g[1:] + g[:-1]); L = np.log(q * mid / (2 * math.pi))
    rows = {}
    for b in range(0, 14):
        m = (L >= b) & (L < b + 1)
        if m.sum() < 800: continue
        x = s[m]; nb = len(x) // 50
        blocks = x[:nb * 50].reshape(nb, 50)
        boot = [np.var(blocks[rng.integers(0, nb, nb)].ravel()) for _ in range(400)]
        rows[b] = {'n': int(m.sum()), 'mean_spacing': float(x.mean()), 'var': float(np.var(x)), 'se': float(np.std(boot)), 'L_mid': float(L[m].mean())}
    out[name] = rows
    print(name, {b: '%.4f±%.4f (n=%d, mean %.3f)' % (r['var'], r['se'], r['n'], r['mean_spacing']) for b, r in rows.items()})
# readings
r1 = all(r['var'] < 0.180 for rows in out.values() for r in rows.values())
comp = []
for b in range(14):
    if b not in out['zeta']: continue
    z = out['zeta'][b]
    for name in out:
        if name == 'zeta' or b not in out[name]: continue
        o = out[name][b]; zz = (o['var'] - z['var']) / math.hypot(o['se'], z['se'])
        comp.append({'bin': b, 'L': name, 'z': zz})
print('reading 1 (all below 0.180):', r1)
print('reading 2 comparisons:', [(c['bin'], c['L'], round(c['z'], 2)) for c in comp], 'all |z| < 3:', all(abs(c['z']) < 3 for c in comp))
json.dump({'bins': out, 'reading1': r1, 'comparisons': comp}, open(os.path.join(H, 'r110.json'), 'w'), indent=1)
