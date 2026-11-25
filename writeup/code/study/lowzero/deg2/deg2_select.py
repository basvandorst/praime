"""Development and fresh curve sets for the degree-2 round, drawn by a fixed seed from the semistable isogeny classes of
deg2_curves.py with NMIN <= N <= NMAX. Strata: g2 = |a_2| in {0, 1, 2} or 'b' (2 | N), g3 = |a_3| in {0, 1, 2, 3} or
'b' (3 | N): 20 cells. DEV: one class per cell. FRESH: up to K classes per cell, disjoint from DEV.
usage: deg2_select.py CANDIDATES.json NMIN NMAX K LO HI SEED"""
import sys, json
import numpy as np
C = json.load(open(sys.argv[1]))['classes']; NMIN, NMAX, K = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
LO, HI, SEED = float(sys.argv[5]), float(sys.argv[6]), int(sys.argv[7])
C = [c for c in C if NMIN <= c['N'] <= NMAX]
def g(c, p): return 'b' if p in c['bad'] else str(abs(c['ap'][str(p)]))
cells = {}
for c in C:
    cells.setdefault((g(c, 2), g(c, 3)), []).append(c)
rng = np.random.default_rng(SEED); dev, fresh = [], []
for key in sorted(cells):
    L = sorted(cells[key], key=lambda c: (c['N'], c['ainv'])); perm = rng.permutation(len(L))
    dev.append(L[perm[0]]); fresh += [L[i] for i in perm[1:1 + K]]
def ident(c): return ('N%d_%s' % (c['N'], '_'.join(str(a) for a in c['ainv']))).replace('-', 'm')
for S in (dev, fresh):
    for c in S: c['id'] = ident(c); c['g2'] = g(c, 2); c['g3'] = g(c, 3)
for name, S in (('dev', dev), ('fresh', fresh)):
    json.dump({'band': [LO, HI], 'seed': SEED, 'NMIN': NMIN, 'NMAX': NMAX, 'curves': S}, open('%s_set.json' % name, 'w'), indent=0)
    cnt = {}
    for c in S: cnt[(c['g2'], c['g3'])] = cnt.get((c['g2'], c['g3']), 0) + 1
    print(name, len(S), 'curves; cells', dict(sorted(cnt.items())))
print('cell sizes in the population', {k: len(v) for k, v in sorted(cells.items())})
