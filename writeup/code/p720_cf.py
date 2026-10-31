# Puzzle 720, odd k: Z^2 - m^2 = 10^k s^2, m = s - eps. A = Z - m = 2^a1 5^b1 U^2, B = Z + m = 2^a2 5^b2 V^2,
# V/U ~ sqrt(R), R = 2^x 5^y (x, y odd), |V/U - sqrt R| < 1/U^2  =>  V/U is a convergent or an adjacent
# intermediate fraction of sqrt(R) (Fatou). Each candidate gives s exactly; test it.
from fractions import Fraction as Fr
from math import isqrt
import sys

def cf_sqrt_candidates(R, qmax):
    """convergents p/q of sqrt(R) (R = Fraction, non-square) with q <= qmax, plus the two neighbouring mediants"""
    # continued fraction of sqrt(num/den) = sqrt(num*den)/den
    n, d = R.numerator, R.denominator
    N = n * d; r = isqrt(N)
    # expansion of (sqrt(N) + P)/Q with Q = d, P = 0
    P, Q = 0, d
    if (N - P * P) % Q:            # make Q | N - P^2
        N *= Q * Q; P *= Q; Q *= Q; r = isqrt(N)
    p0, q0, p1, q1 = 1, 0, (P + r) // Q, 1
    a = p1
    out = [(p1, q1)]
    for _ in range(100000):
        P = a * Q - P; Q = (N - P * P) // Q
        a = (r + P) // Q
        p0, q0, p1, q1 = p1, q1, a * p1 + p0, a * q1 + q0
        if q1 > qmax: break
        out.append((p1, q1))
        out.append((p1 + p0, q1 + q0)); out.append((p1 - p0, q1 - q0))
    return out

def solve_k(k):
    sols = set()
    T2 = 10**k
    smax = isqrt(10**k) + 2          # m^2 < 10^k  =>  s <= m + 1 <= sqrt(10^k) + 1
    from math import log
    xm = k + 2 * int(log(smax, 2)) + 2     # |a1 - a2| <= a1 + a2 = k + 2 v2(s)
    ym = k + 2 * int(log(smax, 5)) + 2     # |b1 - b2| <= k + 2 v5(s)
    for x in range(-xm, xm + 1):
        if (x + k) % 2: continue
        for y in range(-ym, ym + 1):
            if (y + k) % 2: continue
            R = Fr(2)**x * Fr(5)**y
            q = R * T2                                   # 10^k R must be a square, hence parity of x, y
            sq = Fr(isqrt(q.numerator), isqrt(q.denominator))
            if sq * sq != q: continue
            c = Fr(T2) / R; rc = Fr(isqrt(c.numerator), isqrt(c.denominator))   # sqrt(10^k / R)
            assert rc * rc == c
            for V, U in cf_sqrt_candidates(R, smax):
                if U <= 0 or V <= 0: continue
                den = 2 - rc * (V * V - R * U * U) / (U * V)
                if den == 0: continue
                for eps in (1, -1):
                    s = Fr(2 * eps) / den
                    if s.denominator != 1 or s <= 0: continue
                    s = int(s); m = s - eps
                    if m <= 0: continue
                    if len(str(m * m)) != k: continue
                    N = s * s * 10**k + m * m; z = isqrt(N)
                    if z * z == N: sols.add((s, m, N))
    return sols

K = int(sys.argv[1])
for k in range(1, K + 1, 2):
    print(k, sorted(solve_k(k)), flush=True)
