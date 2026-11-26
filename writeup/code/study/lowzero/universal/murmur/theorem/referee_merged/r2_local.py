# Referee: main-term double series e_l(i,j) by exact enumeration; constants of the note.
from fractions import Fraction as Fr
import math, sys
sys.path.insert(0, '.')
from sympy import primerange, factorint
import mpmath as mp

def kron2(a):
    if a % 2 == 0: return 0
    return 1 if a % 8 in (1, 7) else -1

def leg(a, l):
    a %= l
    if a == 0: return 0
    return 1 if pow(a, (l - 1) // 2, l) == 1 else -1

def T_loc(l, i, j, D):
    # T_l(i, j; D) of Step 2(iv)
    if D % (l ** (2 * i)): return 0
    a = -D // l ** (2 * i)
    if l == 2:
        if a % 4 not in (0, 1): return 0
        s = kron2(a)
    else:
        s = leg(a, l)
    return s ** j if j else 1

def e_loc(l, i, j, N, r):
    K = 2 * i + j + (3 if l == 2 else 1)
    mod = l ** K
    tot = Fr(0); cnt = 0
    for p in range(1, mod):
        if p % l == 0: continue
        D = N * (4 * p - r * r * N)
        # D may be negative as an integer; T depends on D mod l^K only, shift by multiple
        Dm = D % (l ** (K + 2)) + l ** (K + 2) * 10
        # make it a true residue representative: T_loc only uses divisibility and residues
        tot += T_loc(l, i, j, D); cnt += 1
    return tot / cnt

def e_closed(l, N, r):
    if N % l == 0: return Fr(1)
    if r % l == 0: return Fr(l * l, l * l - 1)
    return 1 - Fr(1, (l - 1) ** 2 * (l + 1))

print("double series sum_{i<=I, j<=J} e_l(i,j) l^(-i-j) against E_l^(N)(r)")
worst = 0
for N in [2 * 10007, 3 * 10007, 6 * 3343, 30 * 1009]:
    for r in range(4):
        row = []
        for l in [2, 3, 5, 7]:
            I = 4 if l == 2 else (3 if l == 3 else 2)
            J = 9 if l == 2 else (6 if l == 3 else 4)
            s = Fr(0)
            for i in range(I + 1):
                for j in range(J + 1):
                    if (2 * i + j) > (12 if l == 2 else (9 if l == 3 else 6)): continue
                    s += e_loc(l, i, j, N, r) * Fr(1, l ** (i + j))
            ex = e_closed(l, N, r)
            # tail bound: sum over (i > I or j > J or cut) of |e| l^{-i-j}: crude
            row.append((l, float(s), float(ex)))
            worst = max(worst, abs(float(s - ex)))
        print(N, "r=%d" % r, " ".join("l=%d: %.5f/%.5f" % t for t in row))
print("max |truncated series - closed form| =", worst)

# structure at l | N: e_l(i, j) = [i = j = 0]
bad = 0
for N in [2 * 10007, 3 * 10007, 6 * 3343, 30 * 1009]:
    for l in [2, 3, 5]:
        if N % l: continue
        for r in range(5):
            for i in range(3):
                for j in range(5):
                    if 2 * i + j > 7: continue
                    v = e_loc(l, i, j, N, r)
                    if v != (1 if i == j == 0 else 0): bad += 1
print("l | N: e_l(i,j) = [i=j=0] violations:", bad)

# closed forms for l not dividing 2rN
bad = 0
for N in [2 * 10007, 3 * 10007]:
    for r in [1, 2, 4]:
        for l in [5, 7]:
            if N % l == 0 or r % l == 0: continue
            for j in range(5):
                v = e_loc(l, 0, j, N, r)
                w = Fr(1) if j == 0 else (Fr(-1, l - 1) if j % 2 else Fr(l - 2, l - 1))
                if v != w: bad += 1
            for i in range(1, 3):
                if e_loc(l, i, 0, N, r) != Fr(l ** (1 - 2 * i) if False else 1, l ** (2 * i - 1) * (l - 1)): bad += 1
print("closed forms e_l(0,j), e_l(i,0) for l not | 2rN: violations", bad)

# constants
mp.mp.dps = 30
P = list(primerange(2, 10 ** 6))
C0 = mp.mpf(1); a = mp.mpf(1); b = mp.mpf(1); g = mp.mpf(1)
for l in P:
    l = mp.mpf(l)
    C0 *= 1 - 1 / ((l - 1) ** 2 * (l + 1))
    b *= (l ** 3 + l ** 2 - 1) / (l * (l ** 2 + l - 1))
    g *= l * (l + 1) / (l ** 2 + l - 1)
    a *= (l ** 4 - 2 * l ** 2 - l + 1) / (l ** 4 - 2 * l ** 2 + l)
print("C0 = %.9f  alpha = %.6f  beta = %.6f  gamma = %.6f" % (C0, 2 * mp.pi * a, 2 * mp.pi * b, 12 * g))
# cross-check alpha factor from (1 - l^-2) * l (l e + 1)/(l^2 + l - 1), e = 1 - u
a2 = mp.mpf(1)
for l in P[:5000]:
    l = mp.mpf(l); e = 1 - 1 / ((l - 1) ** 2 * (l + 1))
    f1 = (1 - l ** -2) * l * (l * e + 1) / (l ** 2 + l - 1)
    f2 = (l ** 4 - 2 * l ** 2 - l + 1) / (l ** 4 - 2 * l ** 2 + l)
    a2 = max(a2, abs(f1 - f2))
print("max |alpha factor from mixture - closed factor| over first 5000 primes:", mp.nstr(a2 - 1 if a2 >= 1 else a2, 3))
c = mp.zeta(2) * mp.zeta(3) / mp.zeta(6)
print("zeta(2)zeta(3)/zeta(6) = %.6f, c0 = 2x = %.6f" % (c, 2 * c))
# w_d bound: w_d = sum_{Q<q<=2Q, d|q} q^-1/2/phi(q) <= c0 Q^-1/2/phi(d), all q (not only squarefree)
import numpy as np
X = 2 ** 15
phi = np.arange(X + 1, dtype=np.int64)
for l in primerange(2, X + 1):
    phi[l::l] -= phi[l::l] // l
worst = 0; arg = None
for e in range(1, 15):
    Q = 2 ** e
    for d in range(1, 2 * Q + 1):
        qs = np.arange(((Q // d) + 1) * d, 2 * Q + 1, d)
        if len(qs) == 0: continue
        w = np.sum(qs.astype(float) ** -0.5 / phi[qs])
        val = phi[d] * Q ** 0.5 * w
        if val > worst: worst, arg = val, (Q, d)
print("max phi(d) Q^1/2 w_d over Q = 2^1..2^14, all d <= 2Q: %.4f at (Q,d) = %s (bound c0 = %.4f)" % (worst, arg, 2 * c))
# sum_{q <= X} q/phi(q) <= c X for all X <= 2^15
s = np.cumsum(np.r_[0, np.arange(1, X + 1) / phi[1:]])
ratio = s[1:] / np.arange(1, X + 1)
print("max_{X<=2^15} (sum_{q<=X} q/phi(q))/X = %.5f (c = %.5f)" % (ratio.max(), c))
# epsilon for N = 5*7*...*127
eps = sum(1.0 / l for l in primerange(5, 128))
print("eps for 5*7*...*127 = %.4f" % eps)
