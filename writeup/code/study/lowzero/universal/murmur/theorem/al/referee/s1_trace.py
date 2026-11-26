"""Referee step 1: Martin Prop 5.5 (trace of T_p W_Q on S_k^new(N)) against PARI matrices, own code."""
import math, itertools, sys
from fractions import Fraction as Fr
import cypari2
from sympy import primerange, factorint, divisors
pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)

def H(n):
    return Fr(str(pari.qfbhclassno(n)))

def xi(Delta, l):
    if Delta % (l * l) == 0:
        return l * H(-Delta // (l * l)) / H(-Delta) - 1
    return int(pari.kronecker(Delta, l)) - 1

def mu(n):
    f = factorint(n)
    return 0 if any(e > 1 for e in f.values()) else (-1) ** len(f)

def martin(N, Q, k, p):
    M = N // Q
    tot = Fr(0)
    s = -int(math.isqrt(4 * p // Q)) - 1
    while s * s * Q <= 4 * p + 4 * Q * 10:
        if s * s * Q < 4 * p:
            D = 4 * p * Q - s * s * Q * Q
            x = Fr(1)
            for l in factorint(M):
                x *= xi(-D, l)
            poly = 1 if k == 2 else (s * s * Q - p)   # p^(k/2-1) U_{k-2}(s sqrt(Q)/(2 sqrt p))
            tot += poly * x * H(D)
        s += 1
    sig = sum(divisors(p))
    return -tot / 2 + (mu(M) * sig if k == 2 else 0)

bad = 0; cnt = 0; cnt1 = 0; cnt_eig = 0
for N in (6, 10, 14, 15, 21, 30, 42, 66):
    for k in (2, 4):
        pari('mf=mfinit([%d,%d],0)' % (N, k))
        d = int(pari('mfdim(mf)'))
        for Q in divisors(N):
            if Q == 1:
                continue
            if d == 0:
                # empty space: formula must give 0
                for p in [1] + [p for p in primerange(2, 101) if N % p]:
                    v = martin(N, Q, k, p); cnt += 1
                    if v != 0:
                        bad += 1; print('EMPTY-FAIL', N, k, Q, p, v)
                continue
            pari('A=mfatkininit(mf,%d); W=A[2]/A[3]' % Q)
            assert pari('W^2==matid(%d)' % d)
            v = martin(N, Q, k, 1); w = Fr(str(pari('trace(W)'))); cnt1 += 1
            if v != w:
                bad += 1; print('TRW-FAIL', N, k, Q, v, w)
            for p in primerange(2, 101):
                if N % p == 0:
                    continue
                v = martin(N, Q, k, p)
                w = Fr(str(pari('trace(W*mfheckemat(mf,%d))' % p)))
                cnt += 1
                if v != w:
                    bad += 1; print('FAIL', N, k, Q, p, v, w)
        # eigenform route: a_f(l) = -l^(k/2-1) w_l(f) at l | N; eigenspace traces
        if d:
            pari('L=mfeigenbasis(mf)')
            n = int(pari('#L'))
            for i in range(n):
                for l in factorint(N):
                    # compare W_l eigenvalue sign with a_f(l)
                    al = pari('mfcoef(L[%d],%d)' % (i + 1, l))
                    cnt_eig += 1
print('step1: %d (N,k,Q,p) cases p<=100 incl. empty spaces, %d traces of W_Q (p=1), bad = %d' % (cnt, cnt1, bad))
