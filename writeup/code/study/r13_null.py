"""Round 13 check: the chi^2 machinery on multinomial samples drawn from the law itself (same row totals)."""
import os, json, math, pickle
import numpy as np
import r7, r13
import theory_r3 as R3
import common as C

cache = pickle.load(open(r7.CACHE7, 'rb'))
rng = np.random.default_rng(13)
out = {}
for win in ('train', 'C'):
    L, W = R3.WIN[win]
    M = C.primes_mask(L, W); P = np.nonzero(M)[0].astype(np.int64) + L
    for Q in (30, 210):
        U = r7.units(Q); idx = {a: i for i, a in enumerate(U)}
        rows = np.bincount([idx[int(a)] for a in (P[:-1] % Q)], minlength=len(U))
        T = r7.trans(cache, Q, L, W); Pa = T.sum(2).T; Pa = Pa / Pa.sum(1, keepdims=True)
        vals = []
        for _ in range(20):
            N = np.array([rng.multinomial(int(n), Pa[i]) for i, n in enumerate(rows)], float)
            c, d = r13.chi2(N, Pa); vals.append(c / d)
        out['%s_Q%d' % (win, Q)] = {'mean': float(np.mean(vals)), 'sd': float(np.std(vals, ddof=1))}
        print(win, Q, out['%s_Q%d' % (win, Q)], flush=True)
json.dump(out, open(os.path.join(R3.OUT, 'r13_null.json'), 'w'), indent=1)

# null 2: each of the 10 blocks drawn from its own block-level law, compared with the window-level law
out2 = {}
for win in ('train', 'A', 'C'):
    L, W = R3.WIN[win]
    M = C.primes_mask(L, W); P = np.nonzero(M)[0].astype(np.int64) + L
    Q = 210; U = r7.units(Q); idx = {a: i for i, a in enumerate(U)}
    Tw = r7.trans(cache, Q, L, W); Pw = Tw.sum(2).T; Pw = Pw / Pw.sum(1, keepdims=True)
    blocks = []
    for b in range(10):
        Lb, Wb = L + b * (W // 10), W // 10
        sel = (P[:-1] >= Lb) & (P[:-1] < Lb + Wb)
        rows = np.bincount([idx[int(a)] for a in (P[:-1][sel] % Q)], minlength=len(U))
        Tb = r7.trans(cache, Q, Lb, Wb); Pb = Tb.sum(2).T; Pb = Pb / Pb.sum(1, keepdims=True)
        blocks.append((rows, Pb))
    vals = []
    for _ in range(20):
        N = sum(np.array([rng.multinomial(int(n), Pb[i]) for i, n in enumerate(rows)], float) for rows, Pb in blocks)
        c, d = r13.chi2(N, Pw); vals.append(c / d)
    out2[win] = {'mean': float(np.mean(vals)), 'sd': float(np.std(vals, ddof=1))}
    print('drift null', win, out2[win], flush=True)
json.dump({'fixed_law': out, 'drifting_law': out2}, open(os.path.join(R3.OUT, 'r13_null.json'), 'w'), indent=1)
