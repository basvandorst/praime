"""T4b redone (genus2red's conductor omits the 2-part when PARI cannot determine it): for every genus-2 L-function found by
b0_sources.py (all_found list), build lfungenus2, read the conductor N_L from lfunparams, require lfuncheckfeq <= -20,
N_L < 1000, N_L not in {169, 249, 277}; one per (N_L, a_1..a_60); primitive if mean |lambda_p|^2 over good p <= 10^4 is in
[0.8, 1.2]. Then the seeded draw (SEED + 5) of 15 test curves and 2 proxies. Writes T4b into testset.json."""
import json, math, time, sys, os
import numpy as np, cypari2
pari = cypari2.Pari(); pari.allocatemem(800 * 10 ** 6)
SEED = 20261125
d = json.load(open('testset.json'))
t0 = time.time()
seen = {}; rows = []
for N0, f, h, m2 in d['T4b']['all_found']:
    try:
        L = pari('lfungenus2([%s,%s])' % (f, h))
        NL = int(pari.lfunparams(L)[0])
        if NL >= 1000 or NL in (169, 249, 277): continue
        feq = float(pari.lfuncheckfeq(L))
    except Exception as e:
        continue
    if feq > -20: continue
    an = pari.lfunan(L, 10000)
    key = (NL,) + tuple(int(an[i]) for i in range(60))
    if key in seen: continue
    ps = [int(p) for p in pari.primes([2, 10000]) if NL % int(p) != 0]
    m2b = float(np.mean([(float(an[p - 1]) / math.sqrt(p)) ** 2 for p in ps]))
    seen[key] = 1
    rows.append([NL, f, h, round(m2b, 4), round(feq)])
rows.sort()
prim = [r for r in rows if 0.8 <= r[3] <= 1.2]
rng = np.random.default_rng(SEED + 5)
perm = rng.permutation(len(prim))
sel = [prim[i] for i in perm]
d['T4b'] = dict(n_candidates=len(prim), test=sel[:15], proxy=sel[15:17], all_valid=rows,
                rule='lfungenus2 conductor from lfunparams, lfuncheckfeq <= -20, N < 1000, not 169/249/277, unique (N, a_1..a_60), mean |lambda_p|^2 in [0.8,1.2]')
json.dump(d, open('testset.json', 'w'), indent=1)
print('valid genus-2 L-functions with N < 1000:', len(rows), ' primitive:', len(prim), ' conductors:', sorted(set(r[0] for r in rows)))
print('test:', [(r[0], r[3], r[4]) for r in sel[:15]]); print('proxy:', [(r[0], r[3]) for r in sel[15:17]]); print('secs', round(time.time() - t0, 1))
