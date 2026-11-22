"""Independent recomputation of the round-39 prediction for g = 2, 4 (exact within the law) and 6 (heuristic).
rho(k) for all k by sieving over ALL primes 11 <= p <= 210 n + g (no cofactor trick), local factors from the
definition S(set) = prod_p (1 - nu_p/p)/(1 - 1/p)^|set|, conditioned on the class (p | 210 factors = 1).
P_a(g) per row taken from r39.json (the law's shares; not recomputed)."""
import json, math, sys
import numpy as np

Q = 210
WIN = {'train': (10**7, 9 * 10**7), 'A': (10**10, 5 * 10**7)}
R39 = json.load(open('/home/bas/praime/writeup/data/study/r39.json'))


def primes_upto(N):
    s = np.ones(N // 2 + 1, dtype=bool); s[0] = False           # s[i] <-> 2i+1
    for i in range(1, (math.isqrt(N) - 1) // 2 + 1):
        if s[i]:
            p = 2 * i + 1
            s[p * p // 2::p] = False
    return 2 * np.nonzero(s)[0] + 1


def loc(p, nu_big, n_big, parts):
    """local ratio S(big)/prod S(part); parts = list of (nu, n)."""
    v = (1 - nu_big / p) * (1 - 1 / p) ** (-n_big)
    for nu, n in parts:
        v = v / ((1 - nu / p) * (1 - 1 / p) ** (-n))
    return v


def cases(kind, p):
    # returns dict case -> local factor; nu from the configuration (g not divisible by p >= 11 here)
    if kind == 'YY':        # {0,d} / ({0}{d})
        return {'gen': loc(p, 2, 2, [(1, 1), (1, 1)]), 'k': loc(p, 1, 2, [(1, 1), (1, 1)])}
    if kind == 'XX':        # {0,g,d,d+g} / ({0,g}{d,d+g})
        return {'gen': loc(p, 4, 4, [(2, 2), (2, 2)]), 'k': loc(p, 2, 4, [(2, 2), (2, 2)]),
                'lin': loc(p, 3, 4, [(2, 2), (2, 2)])}
    # XY: {0,g,d} / ({0,g}{d})
    return {'gen': loc(p, 3, 3, [(2, 2), (1, 1)]), 'k': loc(p, 2, 3, [(2, 2), (1, 1)]),
            'lin': loc(p, 2, 3, [(2, 2), (1, 1)])}


def modinv_vec(a, p):
    """a^(p-2) mod p, vectorised, p < 2^31."""
    res = np.ones_like(p); base = a % p; e = p - 2
    while (e > 0).any():
        m = (e & 1) == 1
        res = np.where(m, (res * base) % p, res)
        base = (base * base) % p
        e >>= 1
    return res


def logrho(n, g, kind, PR):
    P = PR[PR >= 11]
    pf = P.astype(np.float64)
    c = cases(kind, pf)
    L = np.full(n, float(np.log(c['gen']).sum()))
    # p | k (only p < n matter)
    small = P[P < n]
    ck = np.log(c['k'] / c['gen'])
    for i, p in enumerate(small):
        L[int(p)::int(p)] += ck[i]
    if kind == 'YY':
        return L
    # p | 210 k + s g : 210 k = -s g mod p
    sg = {'XX': (1, -1), 'XYp': (-1,), 'XYm': (1,)}[kind]     # XYp: {0,g,d}, d = g mod p -> 210k - g; XYm: d = -g
    cl = np.log(c['lin'] / c['gen'])
    inv = modinv_vec(np.full_like(P, Q), P)
    for s in sg:
        k0 = ((-s * g) % P * inv) % P
        sm = P < n
        for i in np.nonzero(sm)[0]:
            L[int(k0[i])::int(P[i])] += cl[i]
        big = (~sm) & (k0 < n) & (k0 >= 1)
        np.add.at(L, k0[big], cl[big])
    return L


def main(win):
    L0, W = WIN[win]
    n = W // Q
    PR = primes_upto(Q * n + 64).astype(np.int64)
    r = (Q / 48) / np.log(L0 + Q * np.arange(n) + Q / 2)
    N2 = 1 << int(np.ceil(np.log2(2 * n)))
    F = np.fft.rfft(r, N2)
    w = np.fft.irfft(F * np.conj(F), N2)[:n]                     # w[k] = sum_j r_j r_{j+k}
    m = r.sum(); s1 = np.sum(r * (1 - r)); s2 = np.sum(r * r)
    rYY = np.exp(logrho(n, 0, 'YY', PR))
    dY = w[1:] * (rYY[1:] - 1)
    varY = s1 + 2 * dY.sum()
    out = {'win': win, 'n': n, 'm': m, 'vY': varY / m, 'vY_r39': R39[win]['vY']}
    cs = np.cumsum(dY)
    out['YY_partial'] = {K: float(2 * cs[min(K, n - 1) - 1] / m) for K in (10, 100, 1000, 10000, 100000, n)}
    np.save('/home/bas/scratch/study-review7/rhoYY_%s.npy' % win, rYY[:3001])
    for g in (2, 4, 6):
        rXX = np.exp(logrho(n, g, 'XX', PR)); rXp = np.exp(logrho(n, g, 'XYp', PR)); rXm = np.exp(logrho(n, g, 'XYm', PR))
        if win == 'train':
            np.save('/home/bas/scratch/study-review7/rho_g%d_train.npy' % g, np.vstack([rXX[:3001], rXp[:3001], rXm[:3001]]))
        SXX = 2 * np.sum(w[1:] * (rXX[1:] - 1)); SXY = np.sum(w[1:] * (rXp[1:] - 1 + rXm[1:] - 1))
        vals = []
        for P in R39[win][str(g)]['P']:
            varX = P * m - P * P * s2 + P * P * SXX
            cov = P * s1 + P * SXY
            vals.append((varX - 2 * P * cov + P * P * varY) / (P * m))
        # partial sums of the combined long-range term for the first row
        P = R39[win][str(g)]['P'][0]
        comb = P * P * (2 * w[1:] * (rXX[1:] - 1) - 2 * w[1:] * (rXp[1:] - 1 + rXm[1:] - 1) + 2 * w[1:] * (rYY[1:] - 1)) / (P * m)
        cc = np.cumsum(comb)
        out['g%d' % g] = {'pred_mine': float(np.mean(vals)), 'pred_r39': R39[win][str(g)]['pred'],
                          'maxdiff_rows': float(np.max(np.abs(np.array(vals) - np.array(R39[win][str(g)]['per_row'])))),
                          'null': float(np.mean([1 - p for p in R39[win][str(g)]['P']])),
                          'longrange_partial': {K: float(cc[min(K, n - 1) - 1]) for K in (10, 100, 1000, 10000, 100000, n)},
                          'deficit_over_P': float((1 - P - vals[0]) / P)}
    print(json.dumps(out, indent=1, default=float))
    json.dump(out, open('/home/bas/scratch/study-review7/mine_%s.json' % win, 'w'), indent=1, default=float)


if __name__ == '__main__':
    main(sys.argv[1])
