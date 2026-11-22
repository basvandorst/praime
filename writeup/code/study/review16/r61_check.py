"""Independent check of round 61.
(1) brute force: for small q, enumerate residues mod q of the current prime (0) and k earlier and the next prime,
    weight = local Hardy-Littlewood factor (q - nu) (uniform base measure, which is the lam -> inf limit of the
    e^(-g/lam)-weighted residues at fixed q), and compute I(next residue; past residues) directly.
(2) the coupon formula, own code: exact DP for small q, Taylor (mean/variance of the occupancy) for large q.
(3) k-table, ln ln k growth, the k -> inf sum."""
import itertools, math, sys
import numpy as np
from sympy import primerange


def H(p):
    p = np.asarray(p, float); p = p[p > 0]; return float(-(p * np.log(p)).sum())


def brute(q, k):
    # positions: x0 = 0, x1..xk the earlier primes, xn the next. weight (q - nu)
    J = {}
    for past in itertools.product(range(q), repeat=k):
        for nx in range(q):
            nu = len(set((0,) + past + (nx,)))
            w = q - nu
            if w > 0:
                J[(past, nx)] = J.get((past, nx), 0) + w
    tot = sum(J.values())
    P = {key: v / tot for key, v in J.items()}
    Pp, Pn = {}, {}
    for (pa, nx), v in P.items():
        Pp[pa] = Pp.get(pa, 0) + v; Pn[nx] = Pn.get(nx, 0) + v
    return H(list(Pp.values())) + H(list(Pn.values())) - H(list(P.values()))


def Hpred(q, m):
    m = np.asarray(m, float)
    a = 1.0 / (q - 1)
    c = (q - m - 1) / ((q - 1) * (q - m))
    with np.errstate(divide='ignore', invalid='ignore'):
        t = np.where(c > 0, (q - m) * c * np.log(np.where(c > 0, c, 1)), 0.0)
    return -(m * a * np.log(a) + t)


def coupon_exact(q, k):
    # m = distinct residues among {0} U k draws uniform on q-1 allowed residues (0 is one of them)
    n = q - 1
    dist = np.zeros(n + 1); dist[1] = 1
    for _ in range(k):
        m = np.arange(n + 1)
        stay = dist * m / n; up = dist * (n - m) / n
        dist = stay.copy(); dist[1:] += up[:-1]
    m = np.arange(n + 1)
    return float(Hpred(q, 1) - np.sum(dist[1:] * Hpred(q, m[1:])))


def coupon_taylor(qs, k):
    qs = np.asarray(qs, float); n = qs - 1
    # U = unseen allowed coupons among the n - 1 not pre-seen
    r1 = (1 - 1 / n) ** k; r2 = (1 - 2 / n) ** k
    EU = (n - 1) * r1
    VU = (n - 1) * r1 + (n - 1) * (n - 2) * r2 - EU ** 2
    mu = n - EU
    h = 1e-3 * np.maximum(1, np.sqrt(np.maximum(VU, 0)))
    H0 = Hpred_vec(qs, mu); Hp = Hpred_vec(qs, mu + h); Hm = Hpred_vec(qs, mu - h)
    d2 = (Hp - 2 * H0 + Hm) / h ** 2
    return Hpred_vec(qs, np.ones_like(qs)) - (H0 + 0.5 * d2 * VU)


def Hpred_vec(q, m):
    a = 1.0 / (q - 1); c = (q - m - 1) / ((q - 1) * (q - m))
    c = np.maximum(c, 1e-300)
    return -(m * a * np.log(a) + (q - m) * c * np.log(c))


if __name__ == '__main__':
    print('(1) brute force vs coupon, nats')
    for q, kmax in ((3, 6), (5, 5), (7, 4), (11, 3), (13, 3)):
        for k in range(1, kmax + 1):
            b = brute(q, k); c = coupon_exact(q, k)
            print(' q=%2d k=%d brute %.10f coupon %.10f diff %.1e' % (q, k, b, c, b - c))
    # round 56's local I_p formula
    for p in (3, 5, 7):
        fA = (1 - 2 / p) / (1 - 1 / p) ** 3; fB = (1 - 3 / p) / (1 - 1 / p) ** 3
        Ip = (p - 1) / p ** 2 * fA * math.log((1 - 1 / p) / (1 - 2 / p)) + ((p - 1) * (p - 2) / p ** 2 * fB * math.log((1 - 3 / p) * (1 - 1 / p) / (1 - 2 / p) ** 2) if p > 3 else 0)
        print(' round56 I_p p=%d %.10f coupon k=1 %.10f' % (p, Ip, coupon_exact(p, 1)))
    sys.stdout.flush()
    print('(2) Taylor vs exact DP')
    for q in (101, 1009, 3001):
        for k in (10, 1000, 10000):
            print(' q=%d k=%d exact %.4e taylor %.4e' % (q, k, coupon_exact(q, k), coupon_taylor([q], k)[0]))
    sys.stdout.flush()
    QX = 3000
    small = [int(p) for p in primerange(3, QX)]
    import gmpy2
    big = np.array(list(primerange(QX, 2 * 10 ** 7)), float)
    print('(3) I_k table')
    res = {}
    for k in (1, 2, 3, 5, 10, 30, 100, 1000, 10000):
        s = sum(coupon_exact(q, k) for q in small)
        tb = coupon_taylor(big, k)
        sb = float(tb.sum())
        # tail beyond 2e7: term ~ c/q^2 -> sum_{q>Q} c/q^2 ~ c/(Q ln Q)
        cq = tb[-1000:].mean() * big[-1000:].mean() ** 2
        Q = big[-1]; tail = cq / (Q * math.log(Q))
        # contribution of q in (2e5, 2e7] (r61 stops at 2e5)
        beyond = float(tb[big > 2e5].sum())
        res[k] = s + sb + tail
        print(' k=%5d I_k %.6f nats %.5f bits  [q>2e5 part %.2e nats, tail>2e7 %.1e]' % (k, res[k], res[k] / math.log(2), beyond, tail))
        sys.stdout.flush()
    print(' I2-I1 = %.8f nats' % (res[2] - res[1]))
    print(' growth k=1e3->1e4: %.4f nats; lnln diff %.4f' % (res[10000] - res[1000], math.log(math.log(1e4)) - math.log(math.log(1e3))))
    print(' growth k=100->1e3: %.4f nats; lnln diff %.4f' % (res[1000] - res[100], math.log(math.log(1e3)) - math.log(math.log(1e2))))
    allp = np.array(list(primerange(3, 10 ** 5 + 1)), float)
    t = ((allp - 2) / (allp - 1)) * np.log((allp - 1) / (allp - 2))
    for Q in (100, 1000, 10 ** 4, 10 ** 5):
        print(' k->inf sum up to %d: %.4f bits' % (Q, t[allp <= Q].sum() / math.log(2)))
