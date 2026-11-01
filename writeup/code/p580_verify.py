# Puzzle 580: P(n) - P(n-k) - (n-k)P(k) = 0 iff P(n), P(k) are twin primes (writeup/580.md).
# With m = n - k: E(m, k) = P(m+k) - P(m) - m P(k). m = 1 is the twin case; for m >= 2 we show E(m, k) < 0.
# Run: python3 p580_verify.py   (about 10 seconds; needs numpy and sympy)
from math import log
import numpy as np
from sympy import primerange, primepi

X0 = 396738                          # Dusart (2010): a prime in [x, x + x/(25 ln^2 x)] for x >= 396738
M0 = int(primepi(X0 - 1))            # P(m) < X0  <=>  m <= M0
P = np.array([0] + list(primerange(2, 14_500_000)), dtype=np.int64)   # P[i] = i-th prime, i < 930000
assert P[1] == 2 and P[M0] < X0 <= P[M0 + 1]
print('M0 = pi(396737) =', M0)

# 1. m = 1: E(1, k) = P(k+1) - P(k) - 2, zero exactly for twin primes; and P(n) - P(k) = 2 forces n = k + 1.
for k in range(1, 100000):
    assert (P[k + 1] - P[k] - 2 == 0) == (P[k + 1] - P[k] == 2)
    assert P[k + 2] - P[k] != 2

# 2. Finite part: E(m, k) < 0 for all 2 <= m <= M0 and 1 <= k <= M0 (this covers case A with P(m) < X0 and
#    case B with P(k) < X0).
ks = np.arange(1, M0 + 1)
worst = -10**18
for m in range(2, M0 + 1):
    E = P[m + ks] - P[m] - m * P[ks]
    mx = int(E.max())
    assert mx < 0, (m, int(ks[E.argmax()]))
    worst = max(worst, mx)
print('finite part: E(m, k) < 0 for 2 <= m <= %d, 1 <= k <= %d (largest value %d)' % (M0, M0, worst))

# 3. The constants of the analytic part.  G(p) = 1 + (p+1)/(25 ln^2 (p+1)) bounds the gap after a prime p >= X0;
#    U(n) = n(ln n + ln ln n) >= p_n (n >= 6, Rosser-Schoenfeld), L(n) = n(ln n + ln ln n - 1) <= p_n (n >= 3500).
U = lambda n: n * (log(n) + log(log(n)))
L = lambda n: n * (log(n) + log(log(n)) - 1)
G = lambda p: 1 + (p + 1) / (25 * log(p + 1) ** 2)
for n in [6, 100, 3500, 10**4, 10**5, 10**6]:
    assert P[n] <= U(n) if n < len(P) else True
for n in range(3500, 900000, 997):
    assert L(n) <= P[n] <= U(n)
# Case A (k <= m, m > M0): the k gaps after P(m) are each <= G(P(2m)) <= 1 + (U(2m)+1)/(25 ln^2(2m)) < m,
# so P(m+k) - P(m) < k*m <= m*P(k).  The bound A(m) = (1 + (U(2m)+1)/(25 ln^2 (2m)))/m is decreasing in m:
grid = np.unique(np.geomspace(M0 + 1, 10**15, 4000).astype(np.int64))
A = [(1 + (U(2 * int(m)) + 1) / (25 * log(2 * int(m)) ** 2)) / int(m) for m in grid]
assert all(a2 <= a1 for a1, a2 in zip(A, A[1:])) and A[0] < 0.01
print('case A: G(P(2m))/m <= %.4f for m > M0' % A[0])
# Case B (k > m >= 2, k > M0): P(m+k) - P(m) < P(k) + m*G(P(2k)), and m*G <= (m-1)*P(k) holds when G <= P(k)/2.
#    P(2k) <= U(2k) <= 2.3*L(k) <= 2.3*P(k) for k > M0 (the ratio U(2k)/L(k) is decreasing), and then
#    G(P(2k)) <= 1 + 2.3*P(k)/(25 ln^2 X0) < P(k)/2 for P(k) >= X0.
R = [U(2 * int(k)) / L(int(k)) for k in grid]
assert all(r2 <= r1 for r1, r2 in zip(R, R[1:])) and R[0] <= 2.3
assert 1 + 2.3 * X0 / (25 * log(X0) ** 2) < X0 / 2
print('case B: U(2k)/L(k) <= %.3f for k > M0' % R[0])
print('all checks passed')
