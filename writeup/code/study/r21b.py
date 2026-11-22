"""Round 21 addendum: exact Hardy-Littlewood count distribution for H = 60 (admissible subsets by depth-first search)."""
import os, json, math, time
import numpy as np
from scipy.stats import t as student_t
import theory2 as T2
import theory_r3 as R3
import common as C
from r20 import counts

OUT = R3.OUT
H, K = 60, 10
PS = [int(p) for p in C.small_primes(H) if p > 2]          # odd primes that can be fully covered by <= 30 points


def admissible_sets(pts):
    """all subsets of pts (sorted list) with no odd prime p <= H whose residues are all covered."""
    out = {}
    occ = {p: [0] * p for p in PS}
    cov = {p: 0 for p in PS}
    cur = []

    def rec(start):
        k = len(cur)
        out.setdefault(k, []).append(tuple(cur))
        for i in range(start, len(pts)):
            x = pts[i]
            bad = False
            for p in PS:
                r = x % p
                if occ[p][r] == 0 and cov[p] == p - 1:
                    bad = True; break
            if bad:
                continue
            for p in PS:
                r = x % p
                if occ[p][r] == 0:
                    cov[p] += 1
                occ[p][r] += 1
            cur.append(x)
            rec(i + 1)
            cur.pop()
            for p in PS:
                r = x % p
                occ[p][r] -= 1
                if occ[p][r] == 0:
                    cov[p] -= 1
    rec(0)
    return out


def T_sums():
    T = np.zeros(32)
    T[0] = 1.0
    nsets = 0
    for par in (0, 1):
        sets = admissible_sets(list(range(par, H, 2)))
        for k, lst in sets.items():
            if k == 0:
                continue
            arr = np.array(lst, dtype=np.int64)
            for i in range(0, len(arr), 200000):
                T[k] += float(T2.S_many(arr[i:i + 200000], P=H).sum())
            nsets += len(lst)
    return T, nsets


def pk(T, lam):
    return np.array([sum(math.comb(j, k) * (-1) ** (j - k) * T[j] / lam ** j for j in range(k, len(T))) for k in range(K + 1)])


def main():
    t0 = time.time()
    T, nsets = T_sums()
    print('admissible sets', nsets, 'T_j', [float('%.4g' % x) for x in T if x], round(time.time() - t0), 's', flush=True)
    out = {'T': T.tolist(), 'admissible_sets': nsets}
    for win, (L, W) in R3.WIN.items():
        M = C.primes_mask(L, W + H)
        emp, th = [], []
        for b in range(10):
            lo, n = b * W // 10, W // 10
            cnt = counts(M[lo:lo + n + H - 1], H)
            emp.append(np.bincount(cnt, minlength=K + 1)[:K + 1] / len(cnt))
            u = np.linspace(math.log(L + lo), math.log(L + lo + n), 7)
            Pu = np.array([pk(T, x) for x in u]); w = np.exp(u)
            th.append(np.trapezoid(w[:, None] * Pu, u, axis=0) / np.trapezoid(w, u))
        emp = np.array(emp); th = np.array(th)
        r = emp - th
        t9 = r.mean(0) / (r.std(0, ddof=1) / math.sqrt(10))
        sel = th.mean(0) >= 1e-3
        out[win] = {'data': emp.mean(0).tolist(), 'law_exact': th.mean(0).tolist(), 't9': t9.tolist()}
        print(win, 'misses', int(np.sum(np.abs(t9[sel]) > 3)), 'of', int(sel.sum()), flush=True)
    json.dump(out, open(os.path.join(OUT, 'r21b.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
