# Puzzle 1142: k = sum (p_m - 1)/p_(m-1)#; p_n# k = S_n + T_n with 1 < T_n < 2, so
# floor(p_n# k - p_n floor(p_(n-1)# k - 1)) = p_(n+1). Checked with exact rational tails.
from fractions import Fraction as Fr
from sympy import primerange
P = list(primerange(2, 12000))                 # tails truncated far beyond the n checked
N = 300
def tail(n):                                   # T_n with p_n = P[n-1] (1-based n), exact up to the truncation
    t = Fr(0); den = 1
    for j in range(1, len(P) - n - 1):
        den *= P[n + j - 1]                    # p_(n+1) ... p_(n+j)
        t += Fr(P[n + j] - 1, den)
        if den > 10**200: break
    return t
k = Fr(0); pr = 1
for p in P[:400]: k += Fr(p - 1, pr); pr *= p
for n in range(1, N):
    T = tail(n)
    assert 1 < T < 2, n
prim = 1; prims = [1]
for p in P[:N + 1]: prim *= p; prims.append(prim)   # prims[n] = p_n#
for n in range(2, N):
    pn = P[n - 1]
    v = (prims[n] * k - pn * ((prims[n - 1] * k - 1).__floor__())).__floor__()
    assert v == P[n], n
print(f"1 < T_n < 2 for n < {N}; formula gives p_(n+1) for 2 <= n < {N}; k = {float(k):.15f}")
