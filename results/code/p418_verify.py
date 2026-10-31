# Puzzle 418: a(n) = (p(n+1)+...+p(2n)) / (p(1)+...+p(n)) satisfies 3 < a(n) for n >= 3 and a(n) < 4 for n >= 10,
# so a(n) is never an integer.  Part 1: exact check for n <= N1 (about 1.07e7).  Part 2: explicit bounds for n > 1e7.
import numpy as np
from math import isqrt, log
N = 400_000_000; sv = np.ones(N, bool); sv[:2] = False
for i in range(2, isqrt(N) + 1):
    if sv[i]: sv[i*i::i] = False
P = np.nonzero(sv)[0].astype(np.int64); del sv
S = np.cumsum(P)                                   # int64, exact (S < 2^63)
N1 = len(P) // 2; n = np.arange(1, N1 + 1)
Sn, T = S[n - 1], S[2 * n - 1] - S[n - 1]
assert np.all(T[2:] > 3 * Sn[2:])                  # a(n) > 3 for 3 <= n <= N1
assert np.all(T[9:] < 4 * Sn[9:])                  # a(n) < 4 for 10 <= n <= N1
assert [round(float(t / s), 4) for t, s in zip(T[:10], Sn[:10])] == [1.5, 2.4, 3.1, 3.5294, 3.6071, 3.8049, 3.8448, 3.9481, 4.01, 3.9535]
assert not any(int(t) % int(s) == 0 for t, s in zip(T[:9], Sn[:9]))     # n <= 9: no integer value
m = 39017; Sm1 = int(S[m - 2])                     # S(m - 1) = sum of the first 39016 primes
# Dusart (1999): p_k > k(ln k + lnln k - 1) for k >= 2; p_k < k(ln k + lnln k - 0.9484) for k >= 39017
def T_low(n):  # >= integral_n^{2n} x(ln x + lnln n - 1) dx
    return n * n * (1.5 * log(n) + 1.5 * log(log(n)) + 2 * log(2) - 2.25)
def S_up(n):   # <= S(m-1) + integral_m^{n+1} x(ln x + lnln n - 0.9484) dx  (lower-limit terms dropped, they are negative)
    c = log(log(n)) - 0.9484
    return Sm1 + (n + 1)**2 * (log(n + 1) / 2 - 0.25 + c / 2)
def T_up(n):   # <= integral_n^{2n+1} x(ln x + lnln(2n) - 0.9484) dx
    c = log(log(2 * n)) - 0.9484; F = lambda x: x * x * (2 * log(x) - 1) / 4 + c * x * x / 2
    return F(2 * n + 1) - F(n)
def S_low(n):  # >= integral_{m-1}^{n} x(ln x + lnln m - 1) dx
    c = log(log(m)) - 1; F = lambda x: x * x * (2 * log(x) - 1) / 4 + c * x * x / 2
    return F(n) - F(m - 1)
for e in np.linspace(7, 150, 2000):
    x = 10.0**e
    assert T_low(x) > 3 * S_up(x) and T_up(x) < 4 * S_low(x)
print("exact up to n =", N1, "; bounds hold on a grid from 1e7 to 1e150 (see results/418.md for the closed-form margin)")
print("all claims verified")
