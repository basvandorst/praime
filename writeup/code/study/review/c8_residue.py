"""Check 8: (a) S_a(B) = S(B) prod_{p|Q} 1[a+b != 0 mod p for all b]/(1 - nu_p(B)/p): average over a mod Q = S(B),
and S_a({0,h}) equals a direct product with conditional local factors 1[...]/(1-1/p)^k.
(b) measured H1 = H(c_{n+1}|c_n), H2 = H(c_{n+1}|c_n,c_{n-1}) at train [1e7, 1e8), c = min(g/2, 20), Miller-Madow,
from my own odd-only sieve."""
import math
import numpy as np
from myS import S, PRI


def nu(B, p):
    return len({b % p for b in B})


def S_a(B, a, Q):
    out = S(B)
    for p in [int(q) for q in PRI[PRI <= Q] if Q % q == 0]:
        ok = all((a + b) % p for b in B)
        out *= (1.0 if ok else 0.0) / (1 - nu(B, p) / p)
    return out


def S_a_direct(B, a, Q):
    """direct: product over p of local factors, conditional at p | Q."""
    k = len(B)
    span = max(B)
    out = S(B)
    for p in [int(q) for q in PRI[PRI <= Q] if Q % q == 0]:
        uncond = (1 - nu(B, p) / p) / (1 - 1 / p) ** k
        cond = (1.0 if all((a + b) % p for b in B) else 0.0) / (1 - 1 / p) ** k
        out = out / uncond * cond
    return out


def odd_sieve_primes(lo, hi):
    n = hi
    s = np.ones(n // 2, dtype=bool)          # s[i] <-> 2i+1
    s[0] = False
    for i in range(1, (int(n ** 0.5) - 1) // 2 + 1):
        if s[i]:
            p = 2 * i + 1
            s[p * p // 2::p] = False
    P = 2 * np.nonzero(s)[0] + 1
    return P[(P >= lo) & (P < hi)]


def H(cnt):
    c = cnt[cnt > 0].astype(float)
    p = c / c.sum()
    return float(-np.sum(p * np.log2(p)))


def Hmm(cnt):
    N = cnt.sum()
    return H(cnt) + (np.count_nonzero(cnt) - 1) / (2 * N * math.log(2))


if __name__ == '__main__':
    worst = 0.0
    for B in ([0, 2], [0, 6], [0, 2, 6], [0, 4, 10, 12], [0, 30], [0, 6, 12, 18]):
        for Q in (6, 30, 210):
            avg = np.mean([S_a(B, a, Q) for a in range(Q)])
            worst = max(worst, abs(avg / S(B) - 1))
            for a in range(Q):
                d = abs(S_a(B, a, Q) - S_a_direct(B, a, Q))
                worst = max(worst, d / S(B))
    print('(a) max relative deviation of avg_a S_a(B) from S(B), and of S_a from the direct product: %.1e' % worst)
    P = odd_sieve_primes(10**7, 10**8)
    print('(b) primes in [1e7,1e8):', len(P))
    g = np.diff(P)
    c = np.minimum(g // 2, 20) - 1
    del P, g
    pair = np.bincount(c[:-1] * 20 + c[1:], minlength=400)
    trip = np.bincount((c[:-2] * 20 + c[1:-1]) * 20 + c[2:], minlength=8000)
    H1 = Hmm(pair) - Hmm(pair.reshape(20, 20).sum(1))
    H2 = Hmm(trip) - Hmm(pair)
    # self-consistent variant: both terms from the triple counts
    t3 = trip.reshape(400, 20)
    H2b = Hmm(trip) - Hmm(t3.sum(1))
    H1p = H(pair) - H(pair.reshape(20, 20).sum(1)); H2p = H(trip) - H(t3.sum(1))
    print('    Miller-Madow  H1 = %.4f  H2 = %.4f  H1-H2 = %.4f   (H2 with triple marginal: %.4f, H1-H2 = %.4f)' % (H1, H2, H1 - H2, H2b, H1 - H2b))
    print('    plug-in       H1 = %.4f  H2 = %.4f  H1-H2 = %.4f' % (H1p, H2p, H1p - H2p))
