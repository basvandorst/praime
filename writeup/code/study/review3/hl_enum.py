"""Independent check of rounds 21/22: admissible subsets, T_j (H=30, 60), U_j(h), P(k), gap shares.
Singular series S(A) = prod_p (1 - nu_p(A)/p) / (1 - 1/p)^k, with p = 2 included explicitly (no parity shortcut).
Admissibility test via big-int bit fields with guard bits (all residues of p covered -> carry into guard bit)."""
import math, json, sys, time
import numpy as np
import sympy

OUT = '/home/bas/scratch/study-review3/'
PR_ALL = np.array(list(sympy.primerange(2, 10_000_000)), dtype=np.float64)


def tail_logG(k, P):
    """log prod_{p > P} (1 - k/p)/(1 - 1/p)^k, primes to 1e7 plus tail -C(k,2) sum_{p>1e7} 1/p^2."""
    p = PR_ALL[PR_ALL > P]
    s = float(np.sum(np.log1p(-k / p) - k * np.log1p(-1.0 / p)))
    X = 1e7
    return s - k * (k - 1) / 2 / (X * math.log(X))


def layout(primes):
    off = {}; o = 0; LOW = 0; GUARD = 0
    for p in primes:
        off[p] = o; LOW |= 1 << o; GUARD |= 1 << (o + p); o += p + 1
    return off, LOW, GUARD


def enum(points, base, primes, maxk=40):
    """DFS over subsets A of points (sorted) such that base u A is admissible for every p in primes.
    returns (count of sets incl. empty, sums[j] = sum over |A|=j of prod_p (1 - nu_p(base u A)/p))."""
    off, LOW, GUARD = layout(primes)
    fmask = {p: ((1 << p) - 1) << off[p] for p in primes}
    M = {x: sum(1 << (off[p] + x % p) for p in primes) for x in points}
    st0 = 0
    for x in base:
        st0 |= sum(1 << (off[p] + x % p) for p in primes)
    if ((st0 + LOW) & GUARD):
        return 0, np.zeros(1)
    sums = np.zeros(maxk); cnt = [0]
    pl = [(p, off[p], (1 << p) - 1) for p in primes]
    n = len(points)
    pts = list(points)

    def val(st):
        v = 1.0
        for p, o, m in pl:
            v *= 1.0 - ((st >> o) & m).bit_count() / p
        return v

    def rec(start, st, j):
        cnt[0] += 1
        sums[j] += val(st)
        for i in range(start, n):
            s2 = st | M[pts[i]]
            if (s2 + LOW) & GUARD:
                continue
            rec(i + 1, s2, j + 1)
    rec(0, st0, 0)
    return cnt[0], sums


def T_sums(H):
    primes = list(sympy.primerange(2, H))            # p >= H > span: all residues distinct
    cnt, sums = enum(list(range(H)), [], primes)
    T = np.zeros(len(sums))
    for j in range(len(sums)):
        if sums[j] == 0:
            continue
        T[j] = sums[j] * math.exp(tail_logG(j, H - 1) - j * sum(math.log(1 - 1 / p) for p in primes))
    return cnt - 1, T                                  # nonempty admissible sets


def U_sums(h):
    primes = list(sympy.primerange(2, h + 1))        # p > h > span
    cnt, sums = enum(list(range(1, h)), [0, h], primes)
    U = np.zeros(len(sums))
    for j in range(len(sums)):
        if sums[j] == 0:
            continue
        k = j + 2
        U[j] = sums[j] * math.exp(tail_logG(k, h) - k * sum(math.log(1 - 1 / p) for p in primes))
    return cnt, U


def pk(T, lam, K=12):
    return [sum(math.comb(j, k) * (-1) ** (j - k) * T[j] / lam ** j for j in range(k, len(T))) for k in range(K + 1)]


def share(U, lam):
    return sum((-1) ** j * U[j] / lam ** (j + 1) for j in range(len(U)))


if __name__ == '__main__':
    res = {}
    t0 = time.time()
    for H in (30, 60):
        n, T = T_sums(H)
        nz = [j for j in range(len(T)) if T[j] != 0]
        res['H%d' % H] = {'nonempty_admissible': n, 'T': T[:max(nz) + 1].tolist(), 'maxj': max(nz),
                          'P_at_17.8': pk(T, 17.8)}
        print('H', H, 'admissible nonempty', n, 'max j', max(nz), 'T', ['%.6g' % x for x in T[:max(nz) + 1]],
              round(time.time() - t0), 's', flush=True)
        print('  P(k) at 17.8', ['%.5f' % x for x in pk(T, 17.8)], 'sum', sum(pk(T, 17.8)),
              'mean', sum(k * x for k, x in enumerate(pk(T, 17.8))), 'H/lam', H / 17.8, flush=True)
    Us = {}; tot = 0
    for h in range(2, 71, 2):
        c, U = U_sums(h)
        Us[h] = U; tot += c
    print('U sets total h<=70', tot, round(time.time() - t0), 's', flush=True)
    res['U'] = {h: U.tolist() for h, U in Us.items()}; res['U_sets'] = tot
    for lam in (17.8,):
        sh = {h: share(U, lam) for h, U in Us.items()}
        print('lam', lam, 'share6 %.6f share30 %.6f share52 %.6f' % (sh[6], sh[30], sh[52]),
              'tail>=54 %.6f tail>=72 %.6f' % (1 - sum(v for h, v in sh.items() if h <= 52), 1 - sum(sh.values())))
    # odd h: S({0,h}) = 0 check is implicit (p=2). Telescoping identity check: tail(>=H+2) vs direct avoidance sum
    json.dump(res, open(OUT + 'enum.json', 'w'), indent=1)
    # direct avoidance P(no prime in 1..30 | prime at 0) vs 1 - sum of shares h <= 30
    primes = list(sympy.primerange(2, 31))
    c, sums = enum(list(range(1, 31)), [0], primes)
    V = [sums[j] * math.exp(tail_logG(j + 1, 30) - (j + 1) * sum(math.log(1 - 1 / p) for p in primes)) if sums[j] else 0
         for j in range(len(sums))]
    for lam in (17.8, 30.0):
        direct = sum((-1) ** j * V[j] / lam ** j for j in range(len(V)))
        tel = 1 - sum(share(Us[h], lam) for h in range(2, 31, 2))
        print('identity check lam', lam, 'direct avoid', direct, 'telescoped', tel)
