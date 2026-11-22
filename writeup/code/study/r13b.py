"""Round 13 addendum: the stratification test (preregister/round13.txt, addendum)."""
import os, sys, json, math, pickle, time
from math import gcd
import numpy as np
import theory2 as T2
import theory_r3 as R3
import common as C
import r7, r13

Q2, HM2 = 2310, 38


def build():
    cache = pickle.load(open(r7.CACHE7, 'rb'))
    t0 = time.time()
    for a in r7.units(Q2):
        for h in range(2, HM2 + 1, 2):
            key = (Q2, a, h)
            if key in cache:
                continue
            if gcd(a + h, Q2) != 1:
                cache[key] = (0.0, None); continue
            pt = T2.GPattern([0, h], 4, T2.S_cond(a, Q2))
            cache[key] = (pt.SB, np.array([pt.log_avoid(l, 4) for l in R3.LAMG]) if pt.SB > 0 else None)
        if time.time() - t0 > 450:
            break
    pickle.dump(cache, open(r7.CACHE7, 'wb'))
    print('cached', sum((Q2, a, h) in cache for a in r7.units(Q2) for h in range(2, HM2 + 1, 2)), 'of', 480 * 19)


def law2310(cache, L, W):
    U = r7.units(Q2)
    u = np.linspace(math.log(L), math.log(L + W), 401); w = np.exp(u) / u
    P = np.zeros((len(U), 20))
    for i, a in enumerate(U):
        S0 = float(T2.S_cond(a, Q2)(np.array([[0]]))[0])
        for h in range(2, HM2 + 1, 2):
            SB, la = cache[(Q2, a, h)]
            if SB:
                f = SB / S0 / u * np.exp(np.interp(u, R3.LAMG, la))
                P[i, h // 2 - 1] = np.trapezoid(w * f, u) / np.trapezoid(w, u)
        P[i, 19] = 1 - P[i, :19].sum()
    return U, P


def test():
    cache = pickle.load(open(r7.CACHE7, 'rb'))
    rng = np.random.default_rng(2310)
    out = {}
    for win, (L, W) in R3.WIN.items():
        M = C.primes_mask(L, W); Pr = np.nonzero(M)[0].astype(np.int64) + L
        cl = np.minimum(np.diff(Pr) // 2, 20) - 1
        U2, P2 = law2310(cache, L, W); i2 = {a: i for i, a in enumerate(U2)}
        U1 = r7.units(210); i1 = {a: i for i, a in enumerate(U1)}
        T = r7.trans(cache, 210, L, W); P1 = T.sum(2).T; P1 = P1 / P1.sum(1, keepdims=True)
        s2 = np.array([i2[int(a)] for a in (Pr[:-1] % Q2)])
        N2 = np.zeros((len(U2), 20)); np.add.at(N2, (s2, cl), 1)
        c2, d2 = r13.chi2(N2, P2)
        map21 = np.array([i1[a % 210] for a in U2])
        rows2 = N2.sum(1).astype(int)
        vals = []
        for _ in range(20):
            D2 = np.array([rng.multinomial(n, P2[i]) for i, n in enumerate(rows2)], float)
            D1 = np.zeros((len(U1), 20)); np.add.at(D1, map21, D2)
            c, d = r13.chi2(D1, P1); vals.append(c / d)
        out[win] = {'data_vs_law2310_chi2_per_dof': c2 / d2, 'dof2310': d2,
                    'stratified_null_mod210_mean': float(np.mean(vals)), 'stratified_null_mod210_sd': float(np.std(vals, ddof=1))}
        print(win, {k: round(v, 3) for k, v in out[win].items()}, flush=True)
    json.dump(out, open(os.path.join(R3.OUT, 'r13b.json'), 'w'), indent=1)


if __name__ == '__main__':
    {'build': build, 'test': test}[sys.argv[1]]()
