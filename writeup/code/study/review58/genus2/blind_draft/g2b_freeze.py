"""DRAFT. Seeded draw and frozen predictions for the blind genus-2 test. Stratum P: eps = +1, generic Sato-Tate
(|zM4| <= 4), c >= 1.5; draw n = 72 by numpy default_rng(SEED).permutation of the rows sorted by (N, f, h) (all if fewer).
Writes predictions_g2b.json with per-L: spec, Q, mu, eps, c, T2 = 150, the law (sharp (0,150]: -eps c T/2pi; Tukey(0.1):
0.9 of it; band (10,150]: -eps c 140/2pi), the alternative DEFICIT = 0.929 x law, the four lines; prints its sha256.
usage: .venv/bin/python -I g2b_freeze.py SEED"""
import sys, os, json, math, hashlib
import numpy as np, cypari2
HERE = os.path.dirname(os.path.abspath(__file__)); SEED = int(sys.argv[1])
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6)
rows = [json.loads(l) for l in open(os.path.join(HERE, 'candidates_g2b.jsonl'))]
P = sorted([r for r in rows if r['eps'] > 0 and abs(r['zM4']) <= 4 and r['c'] >= 1.5], key=lambda r: (r['N'], r['f'], r['h']))
idx = np.random.default_rng(SEED).permutation(len(P))[:72]; draw = [P[i] for i in sorted(idx)]
def lines(f, h, eps, c):
    L = pari('lfungenus2([%s,%s])' % (f, h)); out = {}
    for p in (2, 3):
        Pp = 1 / pari.lfuneuler(L, p); co = [float(x) for x in pari.Vecrev(Pp)] if str(pari.type(Pp)) == 't_POL' else [float(Pp)]
        co = [x / co[0] * p ** (-j / 2) for j, x in enumerate(co)]; J = 60; lam = [1.0] + [0.0] * J
        for j in range(1, J + 1): lam[j] = -sum(co[i] * lam[j - i] for i in range(1, min(j, len(co) - 1) + 1))
        F = sum(lam[j] * co[j] * p ** (-j) for j in range(len(co)))
        for v, nm in ((1, '%d' % p), (-1, '1/%d' % p)):
            S = sum(lam[k + v] * co[k] * p ** (-(2 * k + v) / 2) for k in range(len(co)) if 0 <= k + v <= J); out[nm] = -eps * c * S / F * 150 / (2 * math.pi)
    return out
L = {}
for i, r in enumerate(draw):
    law = -r['eps'] * r['c'] * 150 / (2 * math.pi)
    L['B_N%d_%d' % (r['N'], i)] = dict(spec=['g2', r['f'], r['h']], Q=float(r['N']), mu=[0.5, 0.5, 1.5, 1.5], eps=r['eps'], c=r['c'], T2=150,
        pred=law, pred_tukey=0.9 * law, pred_band=-r['eps'] * r['c'] * 140 / (2 * math.pi), alt_deficit=0.929 * law, lines=lines(r['f'], r['h'], r['eps'], r['c']))
out = json.dumps(dict(L=L, seed=SEED, n_stratum=len(P), rule=__doc__), indent=1, sort_keys=True)
open(os.path.join(HERE, 'predictions_g2b.json'), 'w').write(out)
print('stratum P: %d; drawn %d; sum c^2 %.1f; expected se(b) = 10/(23.87 sqrt(sum c^2)) = %.4f' % (len(P), len(draw), sum(r['c'] ** 2 for r in draw), 10 / (23.87 * math.sqrt(sum(r['c'] ** 2 for r in draw)))))
print('sha256', hashlib.sha256(out.encode()).hexdigest())
