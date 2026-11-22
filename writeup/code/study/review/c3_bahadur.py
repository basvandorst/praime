"""Check 3: exact Hardy-Littlewood inclusion-exclusion over all subsets of the live interior versus my own
implementation of the Bahadur form (orders 1..4) and the distinct-index factorial-cumulant exponential (v1).

Exact: P(B)/(S(B)/lam^|B|) = sum_{A subset live} (-1)^|A| S(B u A) / (S(B) lam^|A|).
Bahadur: ln = sum ln(1-p_m) + sum_{|A|=2} k/prod(1-p) - sum_{|A|=3} k/prod(1-p) + sum_{|A|=4} k/prod(1-p),
k = joint cumulant of the X's from the moments mu(A) = S(B u A)/(S(B) lam^|A|) (generic set-partition formula).
"""
import math, itertools, sys
import numpy as np
from myS import S, logG, PRI

LAMS = [17.8, 23.0, 27.6, 34.5]


def set_partitions(seq):
    if not seq:
        yield []
        return
    first, rest = seq[0], seq[1:]
    for part in set_partitions(rest):
        for i in range(len(part)):
            yield part[:i] + [[first] + part[i]] + part[i + 1:]
        yield [[first]] + part


def cumulant(A, mom):
    tot = 0.0
    for pi in set_partitions(list(A)):
        b = len(pi)
        tot += (-1) ** (b - 1) * math.factorial(b - 1) * np.prod([mom(tuple(sorted(x))) for x in pi])
    return tot


def live_interior(B):
    Bs = set(B); s = max(B)
    return [m for m in range(1, s) if m not in Bs and S(list(B) + [m], P=max(s, 2)) > 0]


def exact_all_lams(B, live, lams):
    """vectorised exact inclusion-exclusion over all 2^n subsets, for each lam."""
    s = max(B); n = len(live); kB = len(B)
    P = max(s, kB + n, 2)
    ps = [int(p) for p in PRI[PRI <= max(s, 2)]]
    N = 1 << n
    size = np.zeros(N, dtype=np.int64)
    for i in range(n):
        size[(np.arange(N) >> i) & 1 == 1] += 1
    logS = np.zeros(N)
    dead = np.zeros(N, dtype=bool)
    for p in ps:
        base = 0
        for b in B:
            base |= 1 << (b % p)
        masks = np.full(1, base, dtype=np.uint64)
        for i in range(n):
            masks = np.concatenate([masks, masks | np.uint64(1 << (live[i] % p))])
        nu = np.bitwise_count(masks).astype(np.float64)
        k = kB + size
        f = 1 - nu / p
        dead |= f <= 0
        logS += np.log(np.where(f > 0, f, 1.0)) - k * math.log(1 - 1 / p)
    # primes in (span, P] and beyond: all residues distinct when p > span
    k = kB + size
    for kk in np.unique(k):
        sel = k == kk
        # primes between s and P where kk could exceed p are handled exactly:
        extra = 0.0
        for p in [int(q) for q in PRI[(PRI > max(s, 2)) & (PRI <= max(kk, s, 2))]]:
            extra += math.log(1 - kk / p) - kk * math.log(1 - 1 / p) if kk < p else -np.inf
        logS[sel] += extra + logG(int(kk), max(kk, s, 2))
    SB = S(B)
    Sv = np.where(dead, 0.0, np.exp(logS)) / SB
    sign = np.where(size % 2 == 0, 1.0, -1.0)
    return {lam: float(np.sum(sign * Sv / lam ** size)) for lam in lams}


def bahadur(B, live, lams, orders=(1, 2, 3, 4)):
    s = max(B); SB = S(B)
    cache = {}

    def rho(A):
        if A not in cache:
            cache[A] = S(list(B) + list(A), P=max(s, 2)) / SB
        return cache[A]
    out = {}
    for lam in lams:
        mom = lambda A: rho(A) / lam ** len(A)
        p = {m: rho((m,)) / lam for m in live}
        terms = [sum(math.log1p(-p[m]) for m in live)]
        for r in (2, 3, 4):
            t = 0.0
            for A in itertools.combinations(live, r):
                if r == 4 and len(live) > 16:
                    pass
                k = cumulant(A, mom)
                t += k / np.prod([1 - p[m] for m in A])
            terms.append(t if r != 3 else -t)
        out[lam] = {o: math.exp(sum(terms[:o])) for o in orders}
    return out


if __name__ == '__main__':
    hs = [int(x) for x in sys.argv[1:]] or [6, 8, 10, 12, 14, 18, 20, 24, 30]
    print('B = {0,h}: ratio (Bahadur order k)/(exact HL inclusion-exclusion) - 1, at lam =', LAMS)
    for h in hs:
        B = [0, h]
        live = live_interior(B)
        ex = exact_all_lams(B, live, LAMS)
        bh = bahadur(B, live, LAMS)
        for lam in LAMS:
            e = ex[lam]
            print('h=%3d n_live=%2d lam=%4.1f exact=%.8f  o1 %+.2e  o2 %+.2e  o3 %+.2e  o4 %+.2e' % (
                h, len(live), lam, e, bh[lam][1] / e - 1, bh[lam][2] / e - 1, bh[lam][3] / e - 1, bh[lam][4] / e - 1),
                flush=True)
    # h = 6 algebra: exact 1 - x - y; v1 exponential e^{-x-y-xy}; Bahadur order 2
    B = [0, 6]; live = live_interior(B); lam = 17.8
    x, y = S([0, 2, 6]) / S(B) / lam, S([0, 4, 6]) / S(B) / lam
    print('h=6 live', live, 'S(0,2,4,6)=', S([0, 2, 4, 6]))
    print('  exact 1-x-y          = %.10f' % (1 - x - y))
    print('  v1 exp(-x-y-xy)      = %.10f' % math.exp(-x - y - x * y))
    print('  Bahadur order 2      = %.10f   ((1-x)(1-y) exp(-xy/((1-x)(1-y))))' % ((1 - x) * (1 - y) * math.exp(-x * y / ((1 - x) * (1 - y)))))
    print('  true factorial-cumulant series (incl. diagonal) = log(1-x-y) exactly: exp ->', math.exp(math.log(1 - x - y)))
