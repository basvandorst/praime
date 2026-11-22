"""Round 22 (preregister/round22.txt): the exact Hardy-Littlewood share of gap h (h <= HM) from admissible interior
subsets, and the train tail P(g >= HM + 2) = 1 - sum of the shares."""
import os, sys, json, math, time
import numpy as np
from scipy.stats import t as student_t
import theory2 as T2
import theory_r3 as R3

OUT = R3.OUT
HM = int(os.environ.get('R22_HM', 52))


def U_sums(h):
    """U[j] = sum over admissible j-subsets A of the interior of S({0, h} u A)."""
    PS = [int(p) for p in T2.small_primes(h + 1) if p > 2]
    B = [0, h]
    if T2.S1(B) == 0:
        return np.zeros(1)
    occ = {p: [0] * p for p in PS}; cov = {p: 0 for p in PS}
    for x in B:
        for p in PS:
            if occ[p][x % p] == 0:
                cov[p] += 1
            occ[p][x % p] += 1
    pts = [m for m in range(2, h, 2)]                       # odd interior points are dead (parity)
    sets = {}
    cur = []

    def rec(start):
        sets.setdefault(len(cur), []).append(tuple(cur))
        for i in range(start, len(pts)):
            x = pts[i]
            if any(occ[p][x % p] == 0 and cov[p] == p - 1 for p in PS):
                continue
            for p in PS:
                if occ[p][x % p] == 0:
                    cov[p] += 1
                occ[p][x % p] += 1
            cur.append(x); rec(i + 1); cur.pop()
            for p in PS:
                occ[p][x % p] -= 1
                if occ[p][x % p] == 0:
                    cov[p] -= 1
    rec(0)
    U = np.zeros(max(sets) + 1)
    for j, lst in sets.items():
        arr = np.array([list(B) + list(a) for a in lst], dtype=np.int64)
        for i in range(0, len(arr), 200000):
            U[j] += float(T2.S_many(arr[i:i + 200000], P=max(h, 3)).sum())
    return U, sum(len(v) for v in sets.values())


def share(U, lam):
    return sum((-1) ** j * U[j] / lam ** (j + 1) for j in range(len(U)))


def main():
    t0 = time.time()
    Us = {}; nsets = 0
    for h in range(2, HM + 1, 2):
        r = U_sums(h)
        if isinstance(r, tuple):
            Us[h], n = r; nsets += n
        print('h', h, 'sets so far', nsets, round(time.time() - t0), 's', flush=True)
    rows = [json.loads(l) for l in open(os.path.join(OUT, 'results_r3.jsonl'))]
    import common as C
    L, W = R3.WIN['train']
    M = C.primes_mask(L, W); P = np.nonzero(M)[0].astype(np.int64) + L; g = np.diff(P)
    res = {'HM': HM, 'admissible_sets': nsets, 'blocks': []}
    stats = {}
    for b in range(10):
        Lb, Wb = L + b * (W // 10), W // 10
        u = np.linspace(math.log(Lb), math.log(Lb + Wb), 9); w = np.exp(u) / u
        sh = {h: float(np.trapezoid(w * np.array([share(U, x) for x in u]), u) / np.trapezoid(w, u)) for h, U in Us.items()}
        sel = (P[:-1] >= Lb) & (P[:-1] < Lb + Wb); gb = g[sel]
        for h in Us:
            stats.setdefault('S.gap%d' % h, []).append((float(np.mean(gb == h)), sh[h]))
        stats.setdefault('T1.tail%d' % (HM + 2), []).append((float(np.mean(gb >= HM + 2)), 1 - sum(sh.values())))
    out = {}
    for k, v in stats.items():
        v = np.array(v); r = v[:, 0] - v[:, 1]
        t9 = float(r.mean() / (r.std(ddof=1) / math.sqrt(10)))
        out[k] = {'data': float(v[:, 0].mean()), 'law_exact': float(v[:, 1].mean()), 't9': t9}
    res['stats'] = out
    mis = [k for k, d in out.items() if abs(d['t9']) > 3]
    print('misses', mis, 'of', len(out), flush=True)
    for k in ('S.gap6', 'S.gap30', 'S.gap50', 'S.gap52', 'T1.tail%d' % (HM + 2)):
        if k in out:
            print(k, out[k], flush=True)
    json.dump(res, open(os.path.join(OUT, 'r22_HM%d.json' % HM), 'w'), indent=1)


if __name__ == '__main__':
    main()
