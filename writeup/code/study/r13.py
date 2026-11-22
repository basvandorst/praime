"""Round 13 (preregister/round13.txt): calibration of the residue-conditioned gap law, cell by cell."""
import os, json, math, pickle
import numpy as np
import common as C
import r7
import theory_r3 as R3

OUT = R3.OUT


def chi2(N, P):
    """N, P: (states, classes); expected = row total * P. Cells with exp < 5 pooled per row."""
    E = N.sum(1, keepdims=True) * P
    chi = 0.0; cells = 0
    for i in range(N.shape[0]):
        big = E[i] >= 5
        o = list(N[i][big]); e = list(E[i][big])
        if E[i][~big].sum() >= 5:
            o.append(N[i][~big].sum()); e.append(E[i][~big].sum())
        o = np.array(o); e = np.array(e)
        e = e * o.sum() / e.sum()
        chi += float(np.sum((o - e) ** 2 / e)); cells += len(o)
    dof = cells - N.shape[0]
    return chi, dof


def main():
    cache = pickle.load(open(r7.CACHE7, 'rb'))
    th3 = json.load(open(os.path.join(OUT, 'theory_r3.json')))
    res = {}
    for win, (L, W) in R3.WIN.items():
        M = C.primes_mask(L, W)
        P = np.nonzero(M)[0].astype(np.int64) + L
        g = np.diff(P); cl = np.minimum(g // 2, 20) - 1
        for Q in (30, 210):
            U = r7.units(Q); idx = {a: i for i, a in enumerate(U)}
            st = np.array([idx[int(a)] for a in (P[:-1] % Q)])
            N = np.zeros((len(U), 20)); np.add.at(N, (st, cl), 1)
            T = r7.trans(cache, Q, L, W)                          # (20, states, states)
            Pa = T.sum(2).T                                       # (states, classes)
            Pa = Pa / Pa.sum(1, keepdims=True)
            c1, d1 = chi2(N, Pa)
            J = np.clip(np.array(th3[win]['-1']['_joint']), 0, None); m = J.sum(1); m /= m.sum()
            c0, d0 = chi2(N, np.tile(m, (len(U), 1)))
            res['%s_Q%d' % (win, Q)] = {'primes': int(N.sum()), 'chi2': c1, 'dof': d1, 'chi2_per_dof': c1 / d1,
                                         'z': (c1 / d1 - 1) / math.sqrt(2 / d1),
                                         'unconditioned_chi2_per_dof': c0 / d0}
            print(win, Q, {k: round(v, 3) for k, v in res['%s_Q%d' % (win, Q)].items()}, flush=True)
    json.dump(res, open(os.path.join(OUT, 'r13.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
