"""Round 13, addendum 3: is there a length constraint per residue row?"""
import os, json, math, pickle
import numpy as np
import theory_r3 as R3
import common as C
import r7

cache = pickle.load(open(r7.CACHE7, 'rb'))
Q = 210; U = r7.units(Q); idx = {a: i for i, a in enumerate(U)}
rng = np.random.default_rng(3)
out = {}


def analyse(N, Pa, hbar):
    E = N.sum(1, keepdims=True) * Pa
    proj2, chi_rest, dof_rest = [], 0.0, 0
    for i in range(N.shape[0]):
        m = E[i] >= 5
        e, o, hb = E[i][m], N[i][m], hbar[m]
        e = e * o.sum() / e.sum()
        r = (o - e) / np.sqrt(e)
        v = np.sqrt(e) * (hb - np.sum(e * hb) / e.sum()); v /= np.linalg.norm(v)
        pr = float(r @ v)
        proj2.append(pr * pr)
        chi_rest += float(r @ r) - pr * pr; dof_rest += len(o) - 2
    return float(np.mean(proj2)), chi_rest / dof_rest


for win, (L, W) in R3.WIN.items():
    M = C.primes_mask(L, W); P = np.nonzero(M)[0].astype(np.int64) + L
    g = np.diff(P); cl = np.minimum(g // 2, 20) - 1
    hbar = np.array([np.mean(g[cl == c]) for c in range(20)], float)        # mean gap per class (class 20 = >= 40)
    st = np.array([idx[int(a)] for a in (P[:-1] % Q)])
    N = np.zeros((len(U), 20)); np.add.at(N, (st, cl), 1)
    T = r7.trans(cache, Q, L, W); Pa = T.sum(2).T; Pa = Pa / Pa.sum(1, keepdims=True)
    p_data, rest_data = analyse(N, Pa, hbar)
    S = np.array([rng.multinomial(int(n), Pa[i]) for i, n in enumerate(N.sum(1))], float)
    p_null, rest_null = analyse(S, Pa, hbar)
    G = np.array([np.sum(g[st == i]) for i in range(len(U))], float)
    out[win] = {'mean_proj2_data': p_data, 'chi2_rest_per_dof_data': rest_data,
                'mean_proj2_null': p_null, 'chi2_rest_per_dof_null': rest_null,
                'G_a_rel_sd': float(G.std() / G.mean()), 'G_a_rel_sd_if_independent': float(math.sqrt(np.var(g) * len(g) / len(U)) / G.mean())}
    print(win, {k: round(v, 4) for k, v in out[win].items()}, flush=True)
json.dump(out, open(os.path.join(R3.OUT, 'r13d.json'), 'w'), indent=1)
