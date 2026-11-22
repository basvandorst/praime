"""Independent exact Hardy-Littlewood gap law: U_j(g) = sum over j-subsets A of the (even) interior of S({0,g} u A).
Singular series: S(H) = prod_{p<=g} (1-nu_p/p)/(1-1/p)^k * T_k,  T_k = prod_{p>g} (1-k/p)(1-1/p)^-k from the prime zeta."""
import math, sys, json
import mpmath as mp
from sympy import primerange

mp.mp.dps = 30


def tail_T(k, g):
    """prod_{p > g} (1-k/p)/(1-1/p)^k via log = -sum_{m>=2} (k^m - k)/m * P_{>g}(m)."""
    small = list(primerange(2, g + 1))
    s = mp.mpf(0)
    for m in range(2, 400):
        Pm = mp.primezeta(m) - sum(mp.mpf(p) ** -m for p in small)
        term = (mp.mpf(k) ** m - k) / m * Pm
        s -= term
        if abs(term) < mp.mpf(10) ** -28 and m > 5:
            break
    return mp.e ** s


def U_sums(g):
    ps = list(primerange(2, g + 1))
    base = [0, g]
    occ0 = []
    for p in ps:
        occ0.append({0 % p, g % p})
    if any(len(o) == p for o, p in zip(occ0, ps)):
        return None
    cand = list(range(2, g - 1, 2))
    # prefilter: a point inadmissible with {0,g} alone
    cand = [m for m in cand if all(len(o | {m % p}) < p for o, p in zip(occ0, ps))]
    res = [[m % p for p in ps] for m in cand]
    W = {}
    # iterative DFS over increasing index
    nu0 = [len(o) for o in occ0]
    w0 = 1.0
    for p, n in zip(ps, nu0):
        w0 *= (1 - n / p)
    stack = [(0, [set(o) for o in occ0], w0, 0)]
    count = 0

    def rec(start, occ, w, j):
        nonlocal count
        W[j] = W.get(j, 0.0) + w
        count += 1
        for idx in range(start, len(cand)):
            r = res[idx]
            neww = w
            ok = True
            changed = []
            for t, p in enumerate(ps):
                o = occ[t]
                if r[t] not in o:
                    n = len(o)
                    if n + 1 == p:
                        ok = False
                        break
                    neww *= (p - n - 1) / (p - n)
                    changed.append(t)
            if not ok:
                continue
            for t in changed:
                occ[t].add(r[t])
            rec(idx + 1, occ, neww, j + 1)
            for t in changed:
                occ[t].discard(r[t])
    sys.setrecursionlimit(10000)
    rec(0, [set(o) for o in occ0], w0, 0)
    Cinv = mp.mpf(1)
    for p in ps:
        Cinv *= mp.mpf(1) / (1 - mp.mpf(1) / p)
    U = []
    for j in sorted(W):
        k = j + 2
        U.append(mp.mpf(W[j]) * Cinv ** k * tail_T(k, g))
    return U, count


if __name__ == '__main__':
    out = {}
    for g in [int(a) for a in sys.argv[1:]]:
        r = U_sums(g)
        U, cnt = r
        out[g] = [str(u) for u in U]
        print(g, 'sets', cnt, 'U0', mp.nstr(U[0], 15), 'U1', mp.nstr(U[1], 12) if len(U) > 1 else '', flush=True)
    json.dump(out, open('/home/bas/scratch/study-review4/U_%s.json' % '_'.join(sys.argv[1:]), 'w'))
