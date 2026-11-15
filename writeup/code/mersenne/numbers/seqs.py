"""Real sequence (A000043, 52 terms) and two synthetic nulls with the same length and growth.

null A ('poisson'): t_n = ln a_n is a Poisson process with the empirical rate
    lam = 51 / (ln a_52 - ln a_1); a_1 = 2, a_{k+1} = max(round(e^t), a_k + 1).
    (log log M_p = ln p + ln ln 2, so the points log log M_p are the ln a_n up to a shift.)
null B ('wagstaff'): each prime p is taken with probability min(1, e^gamma*log2(2p)/p).
    Exact Bernoulli over all primes below B = 10^7; above B the hits form a Poisson process of
    intensity e^gamma*log2(2x)/(x ln x) dx (prime density 1/ln x times the probability), sampled by
    inverting its cumulative intensity, each point mapped to the next prime (gmpy2.next_prime).
    Stop at 52 terms.
Output: seqs.npz with real (52,), nullA (NA,52), nullB (NB,52) int64.
"""
import numpy as np, gmpy2, math, sys

REAL = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423,
        9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433,
        1257787, 1398269, 2976221, 3021377, 6972593, 13466917, 20996011, 24036583, 25964951, 30402457,
        32582657, 37156667, 42643801, 43112609, 57885161, 74207281, 77232917, 82589933, 136279841]
assert len(REAL) == 52 and REAL == sorted(REAL)
L = 52
EG = math.exp(0.5772156649015329)


def null_poisson(n, rng):
    lam = (L - 1) / (math.log(REAL[-1]) - math.log(REAL[0]))
    out = np.zeros((n, L), np.int64)
    for i in range(n):
        t = math.log(2.0); a = [2]
        while len(a) < L:
            t += rng.exponential(1 / lam)
            a.append(max(int(round(math.exp(t))), a[-1] + 1))
        out[i] = a
    return out, lam


def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]: s[p * p::p] = False
    return np.nonzero(s)[0]


def Lam(U, UB):  # cumulative intensity from B to x, U = ln x
    return EG / math.log(2) * ((U - UB) + math.log(2) * math.log(U / UB))


def null_wagstaff(n, rng, B=10 ** 7):
    P = primes_upto(B)
    prob = np.minimum(1.0, EG * np.log2(2.0 * P) / P)
    seqs = [[] for _ in range(n)]
    CH = 50000
    for s in range(0, len(P), CH):
        pc, qc = P[s:s + CH], prob[s:s + CH]
        hit = rng.random((n, len(pc))) < qc
        for i in range(n):
            if len(seqs[i]) < L:
                seqs[i].extend(int(v) for v in pc[hit[i]])
    UB = math.log(B)
    out = np.zeros((n, L), np.int64)
    for i in range(n):
        a = seqs[i][:L]
        U = UB; c = 0.0
        while len(a) < L:
            c += rng.exponential(1.0)  # next point of unit-rate process in Lambda
            # solve Lam(U) = c by Newton
            for _ in range(60):
                f = Lam(U, UB) - c
                df = EG / math.log(2) * (1 + math.log(2) / U)
                U -= f / df
            q = int(gmpy2.next_prime(int(math.exp(U))))
            if q > a[-1]:
                a.append(q)
        out[i] = a[:L]
    return out


if __name__ == '__main__':
    NA = int(sys.argv[1]) if len(sys.argv) > 1 else 250
    NB = int(sys.argv[2]) if len(sys.argv) > 2 else 250
    rng = np.random.default_rng(20261115)
    A, lam = null_poisson(NA, rng)
    Bq = null_wagstaff(NB, rng)
    np.savez('seqs.npz', real=np.array(REAL, np.int64), nullA=A, nullB=Bq)
    print('lambda (per unit ln a) =', lam)
    for nm, M in (('nullA', A), ('nullB', Bq)):
        print(nm, M.shape, 'median a_52 = %.3g' % np.median(M[:, -1]),
              'quartiles a_52 = %.3g %.3g' % tuple(np.quantile(M[:, -1], [.25, .75])))
        print('  example:', list(M[0, :14]))
        print('  share starting 2,3,5,7: %.2f' % np.mean((M[:, :4] == [2, 3, 5, 7]).all(1)))
    print('real a_52 =', REAL[-1])
