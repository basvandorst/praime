"""Puzzle 156 (one million doors): door i is toggled once per distinct prime factor of i, so it ends open
iff omega(i) is odd.  m(n) = +1 (open) / -1 (closed), M(n) = sum of m(i), i <= n.

Checks: the count of open doors among 10^6 (page: 500954) and M(10^6) = open - closed = +1908 (the page's sign is for closed - open); the ratio open/n at powers
of ten up to 10^8; the zeros of M(n) (the page conjectures there are infinitely many); the identity
(-1)^omega = mu * g with g supported on squareful numbers, g(p^k) = -(k-1) for k >= 2, on all n <= 10^5;
and the values of H(s) = G(s) zeta(2s) near s = 1/2 and at the first zeta zero (G(1/2+) = 0, G(rho_1) != 0).
"""
import numpy as np
from sympy import factorint, primerange, mobius, divisors
from math import log

N = 10**8
omega = np.zeros(N + 1, dtype=np.uint8)
for p in primerange(2, N + 1):
    omega[p::p] += 1
m = np.where(omega % 2 == 1, 1, -1).astype(np.int8)       # m(1) = -1: door 1 is never touched
m[0] = 0
M = np.cumsum(m, dtype=np.int32)                          # M(n) = open - closed among the first n doors
open_1e6 = int(np.count_nonzero(omega[1:10**6 + 1] % 2 == 1))
assert open_1e6 == 500954 and M[10**6] == 1908            # the page writes M(10^6) = -1908 but computes closed - open = -1908
print("open doors among 10^6:", open_1e6, " M(10^6) =", M[10**6])
for e in range(1, 9):
    n = 10**e
    op = (n + M[n]) // 2
    print(f"n = 10^{e}: open = {op}, open/n = {op / n:.6f}, M(n) = {M[n]}, M(n)/sqrt(n) = {M[n] / n**0.5:+.3f}")
zeros = np.flatnonzero(M[1:] == 0) + 1
print("number of n <= 10^8 with M(n) = 0:", len(zeros), " first:", zeros[:12].tolist(), " last:", zeros[-3:].tolist())
pos = np.flatnonzero(M[1:] > 0) + 1
print("M(n) > 0 for", len(pos), "values of n <= 10^8; last such n:", int(pos[-1]) if len(pos) else None,
      "; max M =", int(M.max()), "at", int(M.argmax()), "; min M =", int(M.min()), "at", int(M.argmin()))

# (-1)^omega = mu * g, g multiplicative, g(p^k) = -(k-1) for k >= 2, g(p) = 0; note m(n) = -(-1)^omega(n)
def g(n):
    r = 1
    for p, k in factorint(n).items():
        if k == 1:
            return 0
        r *= -(k - 1)
    return r
for n in range(1, 10**5 + 1):
    s = sum(mobius(n // d) * g(d) for d in divisors(n))
    assert s == (-1 if omega[n] % 2 == 1 else 1), n
print("(-1)^omega = mu * g checked for n <= 10^5")
# sum |g(d)|/d converges: it equals prod_p (1 + 1/(p-1)^2) (Euler product of |g|, which is multiplicative)
S = sum(abs(g(d)) / d for d in range(1, 10**5))
Pr = 1.0
for p in primerange(2, 10**5):
    Pr *= 1 + 1 / (p - 1)**2
assert S < Pr < 3
print(f"sum |g(d)|/d over d < 10^5 = {S:.4f}; prod_(p < 10^5) (1 + 1/(p-1)^2) = {Pr:.4f}")

# G(s) = prod_p (1 - 2p^-s)/(1 - p^-s)^2 = H(s)/zeta(2s) with H(s) = prod_p (1 - 2p^-s)/((1 - p^-s)^2 (1 - p^-2s)),
# whose factors are 1 - 2p^-3s + O(p^-4s): H converges absolutely for Re s > 1/3.  So G(sigma) -> 0 as sigma -> 1/2+
# (pole of zeta(2 sigma)), and G(rho_1) = H(rho_1)/zeta(2 rho_1) != 0 at the first zeta zero rho_1 = 1/2 + 14.1347i
# (a factor of H vanishes only where p^-s = 1/2, impossible on Re s = 1/2; zeta(1 + 28.27i) != 0).
import mpmath as mp
mp.mp.dps = 20
def H(s, P=2*10**5):
    return mp.fprod((1 - 2*mp.power(p, -s)) / ((1 - mp.power(p, -s))**2 * (1 - mp.power(p, -2*s))) for p in primerange(2, P))
for sig in (mp.mpf('0.5'), mp.mpf('0.51'), mp.mpf('0.55'), mp.mpf('0.6')):
    h = H(sig)
    gs = h / mp.zeta(2*sig) if sig > 0.5 else mp.mpf(0)
    print(f"sigma = {mp.nstr(sig, 3)}: H = {mp.nstr(h, 6)}, G = H/zeta(2 sigma) = {mp.nstr(gs, 6)}")
rho1 = mp.mpc(0.5, mp.zetazero(1).imag)
h1 = H(rho1); z1 = mp.zeta(2*rho1)
print("rho_1 =", mp.nstr(rho1, 12), " H(rho_1) =", mp.nstr(h1, 6), " zeta(2 rho_1) =", mp.nstr(z1, 6),
      " |G(rho_1)| =", mp.nstr(abs(h1 / z1), 6), " zeta(1/2) =", mp.nstr(mp.zeta(0.5), 8))
assert abs(h1) > 0.1 and abs(z1) > 0.1
