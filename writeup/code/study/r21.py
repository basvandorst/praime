"""Round 21 (preregister/round21.txt): the exact Hardy-Littlewood distribution of the number of primes in [n, n+30)."""
import os, json, math, itertools
import numpy as np
from scipy.stats import t as student_t
import theory2 as T2
import theory_r3 as R3
import common as C
from r20 import counts

OUT = R3.OUT
H, K = 30, 7


def T_sums():
    T = np.zeros(H // 2 + 2)
    T[0] = 1.0
    T[1] = H                                                   # S({m}) = 1
    for par in (0, 1):
        pts = np.arange(par, H, 2)
        for j in range(2, len(pts) + 1):
            sets = np.array(list(itertools.combinations(pts, j)), dtype=np.int64)
            T[j] += float(T2.S_many(sets).sum())
    return T


def pk(T, lam):
    js = np.arange(len(T))
    return np.array([sum(math.comb(j, k) * (-1) ** (j - k) * T[j] / lam ** j for j in range(k, len(T))) for k in range(K + 1)])


def law(T, L, W):
    u = np.linspace(math.log(L), math.log(L + W), 7)
    Pu = np.array([pk(T, x) for x in u]); w = np.exp(u)
    return np.trapezoid(w[:, None] * Pu, u, axis=0) / np.trapezoid(w, u)


def main():
    T = T_sums()
    print('T_j', [float('%.4g' % x) for x in T], flush=True)
    out = {'T': T.tolist()}
    for win, (L, W) in R3.WIN.items():
        M = C.primes_mask(L, W + H)
        emp, th = [], []
        for b in range(10):
            lo, n = b * W // 10, W // 10
            cnt = counts(M[lo:lo + n + H - 1], H)
            emp.append(np.bincount(cnt, minlength=K + 1)[:K + 1] / len(cnt)); th.append(law(T, L + lo, n))
        emp = np.array(emp); th = np.array(th)
        r = emp - th
        t9 = r.mean(0) / (r.std(0, ddof=1) / math.sqrt(10))
        res = {'data': emp.mean(0).tolist(), 'law_exact': th.mean(0).tolist(), 't9': t9.tolist(),
               'sum_law': float(law(T, L, W).sum())}
        out[win] = res
        sel = th.mean(0) >= 1e-3
        print(win, 'misses', int(np.sum(np.abs(t9[sel]) > 3)), 'of', int(sel.sum()), ' '.join(
            '%d:%.4f/%.4f(t%.1f)' % (k, emp.mean(0)[k], th.mean(0)[k], t9[k]) for k in range(K + 1)), flush=True)
    json.dump(out, open(os.path.join(OUT, 'r21.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
