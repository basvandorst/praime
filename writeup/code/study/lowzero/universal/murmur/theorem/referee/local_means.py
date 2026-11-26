# Referee: exact local means of F_l(D_p), D_p = N(4p - r^2 N), p uniform over l-adic units.
# F_l(D) = (1 - chi(l)/l)^(-1) * l^(-v) [sigma(l^v) - chi(l) sigma(l^(v-1))],
# where D = D0 m^2, -D0 fundamental, v = v_l(m), chi = kron(-D0, .).
# Product over l of F_l(D) = pi H(D)/sqrt(D).
from fractions import Fraction as Fr
from math import gcd
import sympy

def sig(l, v):
    if v < 0:
        return 0
    return sum(l**i for i in range(v + 1))

def kron(a, l):
    return int(sympy.jacobi_symbol(a % l, l)) if l > 2 else (0 if a % 2 == 0 else (1 if a % 8 in (1, 7) else -1))

def F_local(D, l):
    # D > 0, D = 0,3 mod 4 ; returns exact Fraction
    if l == 2:
        v = 0
        while D % 4 == 0 and (-(D // 4)) % 4 in (0, 1):
            D //= 4
            v += 1
        d0 = -D
        chi = 0 if d0 % 2 == 0 else (1 if d0 % 8 == 1 else -1)
    else:
        a = 0
        while D % l == 0:
            D //= l
            a += 1
        v = a // 2
        chi = 0 if a % 2 else kron(-D, l)
    return Fr(l, l - chi) * Fr(sig(l, v) - chi * sig(l, v - 1), l**v)

def mean_F(l, N, r, J):
    M = l**J
    tot = Fr(0)
    cnt = 0
    for a in range(1, M):
        if a % l == 0:
            continue
        p = a + M * (10 * r * r * N + 7)  # positive lift, 4p > r^2 N
        D = N * (4 * p - r * r * N)
        tot += F_local(D, l)
        cnt += 1
    return tot / cnt

def mean_F_r0(l, N, J):
    M = l**J
    tot = Fr(0); cnt = 0
    for a in range(1, M):
        if a % l == 0:
            continue
        p = a + M * 5
        tot += F_local(4 * p * N, l)
        cnt += 1
    return tot / cnt

if __name__ == "__main__":
    for N in (10007, 10009):
        print("N", N, "N mod 8", N % 8)
        for l, J in ((2, 14), (3, 9), (5, 6)):
            print(" l", l, "r=0 mean:", float(mean_F_r0(l, N, J)))
            for r in range(1, 7):
                m = mean_F(l, N, r, J)
                if l == 2:
                    pred = "?"
                elif r % l == 0:
                    pred = Fr(l * l, l * l - 1)
                else:
                    pred = 1 - Fr(1, (l - 1) * (l * l - 1))
                print("   r", r, "mean", float(m), "pred", float(pred) if pred != "?" else pred, "approx frac", m.limit_denominator(2000))
