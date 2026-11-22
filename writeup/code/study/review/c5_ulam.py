"""Check 5: Ulam diagonals 4n^2 + bn + c.
(a) omega rules by brute force: omega_d(p) = |R u (R - d)|; rule: p | d -> w; p | 16 d^2 - D -> 3, with D = b^2 - 16c
    (the discriminant of 4x^2 + bx + c) versus the write-up's D = b^2 - 64c.
(b) C_1/C^2 for 4n^2 + 1 by brute force over p < 3000.
(c) predicted chi^2/curve on the train window [1e7, 1e8) for all 195 curves, own implementation:
    Var/E = 1 - C sum w^2 / sum w + 2 C sum_d (R_d - 1) sum_n w_n w_{n+d} / sum w,  w_n = 1/ln f(n),
    R_d = C_d/C^2 = prod_p (1 - omega_d(p)/p)/(1 - w(p)/p)^2: brute force for p <= 3600 (all d), generic
    (1 - 2w/p)/(1 - w/p)^2 for 3600 < p < 1e6 (p | 16 d^2 - D for p > 3600 ignored: factor < 1 + 1/p).
"""
import math, sys
import numpy as np
from myS import PRI


def roots(f, p):
    a, b, c = f
    return [x for x in range(p) if (a * x * x + b * x + c) % p == 0]


def curves():
    out = []
    for b in (-4, -2, 0, 2, 4):
        for c in range(1, 81, 2):
            D = b * b - 16 * c
            if D >= 0 and math.isqrt(D) ** 2 == D:
                continue
            out.append((4, b, c))
    return out


def part_a():
    bad16 = bad64 = checked = 0
    for f in [(4, 0, 1), (4, 2, 3), (4, -4, 5), (4, 4, 41), (4, 2, 79), (4, 0, 7)]:
        a, b, c = f
        D16 = b * b - 16 * c
        D64 = b * b - 64 * c
        for p in [int(q) for q in PRI[(PRI > 2) & (PRI < 400)]]:
            R = set(roots(f, p)); w = len(R)
            if w != 2:
                continue
            for d in range(1, 3 * p):
                om = len(R | {(r - d) % p for r in R})
                if d % p == 0:
                    continue
                checked += 1
                pred16 = 3 if (16 * d * d - D16) % p == 0 else 4
                pred64 = 3 if (16 * d * d - D64) % p == 0 else 4
                bad16 += om != pred16
                bad64 += om != pred64
    # p | d rule
    badd = 0
    for f in [(4, 0, 1), (4, 2, 3)]:
        for p in [int(q) for q in PRI[(PRI > 2) & (PRI < 400)]]:
            R = set(roots(f, p))
            om = len(R | {(r - p) % p for r in R})
            badd += om != len(R)
    print('(a) w=2 cases checked %d: rule with D=b^2-16c wrong in %d; with D=b^2-64c wrong in %d; p|d rule wrong in %d'
          % (checked, bad16, bad64, badd))


def ratio_table(f, dmax, pbrute=3600, pgen=1_000_000):
    a, b, c = f
    D = b * b - 16 * c
    d = np.arange(dmax + 1)
    logR = np.zeros(dmax + 1)
    for p in [int(q) for q in PRI[PRI <= pbrute]]:
        R = roots(f, p); w = len(R)
        tab = np.empty(p)
        for t in range(p):
            om = len(set(R) | {(r - t) % p for r in R})
            tab[t] = math.log((1 - om / p) / (1 - w / p) ** 2) if om < p else -np.inf
        logR += tab[d % p]
    ps = PRI[(PRI > pbrute) & (PRI < pgen)]
    leg = np.array([pow(int(D % p), (int(p) - 1) // 2, int(p)) for p in ps])
    w = np.where(leg == 1, 2, np.where(leg == 0, 1, 0)).astype(float)
    logR += float(np.sum(np.log1p(-2 * w / ps) - 2 * np.log1p(-w / ps)))
    return np.exp(logR)


def bh_C(f, pbrute=3600, pgen=1_000_000):
    a, b, c = f
    D = b * b - 16 * c
    s = 0.0
    for p in [int(q) for q in PRI[PRI <= pbrute]]:
        w = len(roots(f, p))
        s += math.log((1 - w / p) / (1 - 1 / p))
    ps = PRI[(PRI > pbrute) & (PRI < pgen)]
    leg = np.array([pow(int(D % p), (int(p) - 1) // 2, int(p)) for p in ps])
    w = np.where(leg == 1, 2, np.where(leg == 0, 1, 0)).astype(float)
    s += float(np.sum(np.log1p(-w / ps) - np.log1p(-1.0 / ps)))
    return math.exp(s)


def part_b():
    f = (4, 0, 1)
    r = 1.0
    for p in [int(q) for q in PRI[PRI < 3000]]:
        R = set(roots(f, p)); w = len(R)
        om = len(R | {(x - 1) % p for x in R})
        r *= (1 - om / p) / (1 - w / p) ** 2
    print('(b) C_1/C^2 for 4n^2+1, product over p < 3000: %.5f' % r)


def part_c(L=10**7, W=9 * 10**7, ncurves=None):
    vals = []
    cs = curves()
    if ncurves:
        cs = cs[:ncurves]
    for f in cs:
        a, b, c = f
        n = np.arange(1, int(math.isqrt((L + W) // 4)) + 3, dtype=np.int64)
        v = a * n * n + b * n + c
        n = n[(v >= L) & (v < L + W)]
        v = (a * n * n + b * n + c).astype(np.float64)
        wt = 1 / np.log(v)
        N = len(n)
        C = bh_C(f)
        R = ratio_table(f, N)
        cor = np.correlate(wt, wt, mode='full')[N:]           # lags 1..N-1
        var_over_E = 1 - C * np.sum(wt ** 2) / np.sum(wt) + 2 * C * np.sum((R[1:N] - 1) * cor) / np.sum(wt)
        vals.append(var_over_E)
    print('(c) curves %d, n per curve %d, predicted chi^2/curve at train = %.4f' % (len(vals), N, np.mean(vals)))
    return vals




def big_factor_corr(f, dmax, pbrute=3600):
    """log-correction for primes p > pbrute dividing 16 d^2 - D (omega = 3 instead of 4), all such primes."""
    a, b, c = f
    D = b * b - 16 * c
    d = np.arange(dmax + 1, dtype=np.int64)
    m = 16 * d * d - D                      # positive, < 2e8
    corr = np.zeros(dmax + 1)
    rem = m.copy()
    for p in [int(q) for q in PRI[(PRI > pbrute) & (PRI < 15000)]]:
        hit = rem % p == 0
        if hit.any():
            corr[hit] += math.log((1 - 3 / p) / (1 - 4 / p))
            while True:
                h2 = rem % p == 0
                if not h2.any():
                    break
                rem[h2] //= p
    # strip small primes (<= pbrute) from the remainder too
    for p in [int(q) for q in PRI[PRI <= pbrute]]:
        while True:
            h2 = rem % p == 0
            if not h2.any():
                break
            rem[h2] //= p
    big = rem > 1                            # leftover cofactor is a prime > 15000 (m < 15000^2)
    corr[big] += np.log((1 - 3 / rem[big]) / (1 - 4 / rem[big]))
    return corr


def part_c2(L=10**7, W=9 * 10**7):
    vals = []
    for f in curves():
        a, b, c = f
        n = np.arange(1, int(math.isqrt((L + W) // 4)) + 3, dtype=np.int64)
        v = a * n * n + b * n + c
        n = n[(v >= L) & (v < L + W)]
        v = (a * n * n + b * n + c).astype(np.float64)
        wt = 1 / np.log(v)
        N = len(n)
        C = bh_C(f)
        R = ratio_table(f, N) * np.exp(big_factor_corr(f, N))
        cor = np.correlate(wt, wt, mode='full')[N:]
        vals.append(1 - C * np.sum(wt ** 2) / np.sum(wt) + 2 * C * np.sum((R[1:N] - 1) * cor) / np.sum(wt))
    vals = np.array(vals)
    print('(c2) all primes | 16d^2-D included: predicted chi^2/curve at train = %.4f (198 curves); '
          'without the 39 b=+4 duplicates of b=-4: %.4f' % (vals.mean(), vals[[i for i, f in enumerate(curves()) if f[1] != 4]].mean()))


if __name__ == '__main__':
    print('number of curves:', len(curves()))
    part_a()
    part_b()
    part_c()
    part_c2()
