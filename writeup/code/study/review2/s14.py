"""Section 13 check: Sophie Germain primes. Own singular series from the forms themselves (roots counted by brute
force mod q, no token encoding), exact inclusion-exclusion for the shares of gaps 6..24, and the counts."""
import math, json
from itertools import combinations
import numpy as np
from sieve import small_primes, mask

PR = small_primes(10 ** 6)
LN2 = math.log(2)


def S_forms(forms):
    """forms: list of (a, b) meaning a x + b."""
    k = len(forms); logS = 0.0
    span = max(abs(b) for a, b in forms) + 3
    for q in PR:
        q = int(q)
        if q <= 2 * span + 2:
            roots = set()
            for a, b in forms:
                if a % q:
                    roots.add((-b * pow(a, -1, q)) % q)
                elif b % q == 0:
                    return 0.0
            nu = len(roots)
        else:
            nu = k
        if nu >= q:
            return 0.0
        logS += math.log(1 - nu / q) - k * math.log(1 - 1 / q)
    return math.exp(logS)


def sg(m):
    return [(1, m), (2, 2 * m + 1)]


def share_exact(h, L, W):
    S0 = S_forms(sg(0))
    B = sg(0) + sg(h)
    SB = S_forms(B)
    inter = [m for m in range(1, h) if S_forms(B + sg(m)) > 0]
    terms = []            # (sign, S, number of extra events)
    for k in range(len(inter) + 1):
        for A in combinations(inter, k):
            F = B + [f for m in A for f in sg(m)]
            terms.append(((-1) ** k, S_forms(F), k))
    u = np.linspace(math.log(L), math.log(L + W), 2001)
    l1, l2 = u, u + LN2
    val = sum(s * Sv / S0 / (l1 * l2) ** (k + 1) for s, Sv, k in terms)
    wgt = np.exp(u) / (l1 * l2)
    return float(np.trapezoid(wgt * val, u) / np.trapezoid(wgt, u)), len(inter)


if __name__ == '__main__':
    out = {}
    out['S_{x,2x+1}'] = S_forms(sg(0))
    out['S_p,p+6'] = S_forms(sg(0) + sg(6))
    out['S_p,p+3'] = S_forms(sg(0) + sg(3)); out['S_p,p+15'] = S_forms(sg(0) + sg(15))
    out['S_p,p+30'] = S_forms(sg(0) + sg(30))
    L, W = 10 ** 7, 9 * 10 ** 7
    for h in (6, 12, 18, 24, 30):
        out['law_exact_gap%d' % h] = share_exact(h, L, W)
    print(out, flush=True)
    # data: SG primes in [1e7, 1e8)
    sgs = []
    for s in range(L, L + W, 10 ** 7):
        M = mask(s, 10 ** 7); M2 = mask(2 * s + 1, 2 * 10 ** 7)        # M2[i] = 2 s + 1 + i prime
        p = np.nonzero(M)[0]
        sgs.append(p[M2[2 * p]] + s)
    T = np.concatenate(sgs)
    out['count'] = int(len(T))
    out['mod6'] = sorted(set((T % 6).tolist()))
    g = np.diff(T)
    for h in (6, 12, 18, 24, 30):
        out['data_gap%d' % h] = float(np.mean(g == h))
    sm = g <= 300
    out['data_le300'] = float(np.mean(sm))
    out['same5'] = float(np.mean(g[sm] % 5 == 0)); out['same7'] = float(np.mean(g[sm] % 7 == 0))
    out['mod5_classes'] = {int(c): int(np.sum(T % 5 == c)) for c in range(5)}
    print(out, flush=True)
    json.dump(out, open('/home/bas/scratch/study-review2/s14.json', 'w'), indent=1)
