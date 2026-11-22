"""Round 13 diagnostic: per gap class, the mean over the residue rows (mod 210) of the squared standardised residual
((N - E)/sqrt(E))^2, cells with E >= 5 only; about 1 for independent draws."""
import os, json, pickle
import numpy as np
import theory_r3 as R3
import common as C
import r7

cache = pickle.load(open(r7.CACHE7, 'rb'))
Q = 210; U = r7.units(Q); idx = {a: i for i, a in enumerate(U)}
out = {}
for win in ('train', 'A'):
    L, W = R3.WIN[win]
    M = C.primes_mask(L, W); P = np.nonzero(M)[0].astype(np.int64) + L
    cl = np.minimum(np.diff(P) // 2, 20) - 1
    st = np.array([idx[int(a)] for a in (P[:-1] % Q)])
    N = np.zeros((len(U), 20)); np.add.at(N, (st, cl), 1)
    T = r7.trans(cache, Q, L, W); Pa = T.sum(2).T; Pa = Pa / Pa.sum(1, keepdims=True)
    E = N.sum(1, keepdims=True) * Pa
    r2 = np.where(E >= 5, (N - E) ** 2 / np.where(E > 0, E, 1), np.nan)
    out[win] = {'class_%d' % (c + 1): (float(np.nanmean(r2[:, c])), int(np.sum(E[:, c] >= 5))) for c in range(20)}
    print(win, ' '.join('%d:%.2f(%d)' % (c + 1, *out[win]['class_%d' % (c + 1)]) for c in range(20)), flush=True)
json.dump(out, open(os.path.join(R3.OUT, 'r13e.json'), 'w'), indent=1)
