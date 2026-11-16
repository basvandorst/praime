# Item C: primes with prime order of 2; closed form, finite form, direct count.
import numpy as np, math
from fractions import Fraction
from sympy import primerange, n_order, isprime
import mpmath as mp

primes = np.array(list(primerange(3, 10**7)), dtype=np.float64)
r = primes
C2 = math.exp(np.sum(np.log1p(-1 / (r - 1) ** 2)) - 1 / (1e7 * math.log(1e7)))
print('C2 = %.10f' % C2)

# closed form: sum_k m g / k^2 = (1/2 + sum_{a>=2} 4^-a) * prod_{r odd} (1 + (r-1)/(r-2) * 1/(r^2-1))
# = 7/12 * prod_r (1 + 1/((r-2)(r+1)))
logprod = np.sum(np.log1p(1 / ((r - 2) * (r + 1)))) + 1 / (1e7 * math.log(1e7))
Ssum = 7 / 12 * math.exp(logprod)
A = C2 * Ssum
print('A (Euler product) = %.6f' % A)
# share of 3|k: factor at r=3 is 1 + 2/(1*8) = 5/4, of which the 3|k part is 2/8 = 1/4
f3 = 1 + Fraction(2, 1) * Fraction(1, 3**2 - 1)
print('share 3|k exactly =', (f3 - 1) / f3)

# direct partial sum check of A
K = 2_000_000
k = np.arange(1, K + 1)
g = np.ones(K + 1)
for p in primerange(3, K + 1):
    g[p::p] *= (p - 1) / (p - 2)
g = g[1:]
m = np.where(k % 2 == 1, 0.5, np.where(k % 4 == 0, 1.0, 0.0))
terms = m * g / k.astype(float) ** 2
print('A partial sum k<=%d: %.6f ; 3|k share: %.6f' % (K, C2 * terms.sum(), terms[k % 3 == 0].sum() / terms.sum()))

# finite form N(x)
xg, wg = np.polynomial.legendre.leggauss(30)
def N_of(x):
    kmax = int((x - 1) / 6)
    kk = k[:kmax].astype(float)
    T = (x - 1) / (2 * kk)
    ok = T > 3
    kk, T = kk[ok], T[ok]
    a0 = math.log(3); b = np.log(T); l2k = np.log(2 * kk)
    I = np.zeros_like(kk)
    # split [ln3, lnT] into pieces of length <= 0.5
    npieces = int(np.ceil((b.max() - a0) / 0.5))
    for i in range(npieces):
        lo = a0 + i * 0.5
        hi = np.minimum(lo + 0.5, b)
        act = hi > lo
        if not act.any():
            break
        h = np.where(act, hi - lo, 0.0)
        for xi, wi in zip(xg, wg):
            u = lo + (xi + 1) / 2 * h
            I += wi * h / 2 * np.exp(u) / (u * (u + l2k))
    w = 2 * C2 * m[:kmax][ok] * g[:kmax][ok] / kk * I
    three = (kk % 3 == 0)
    return w.sum(), w[three].sum() / w.sum(), (kk, I)

for x in (10**6, 10**7):
    Nx, sh, (kk, I) = N_of(x)
    # spot-check the integral with mpmath
    for kt in (1, 7, 1000):
        Tt = (x - 1) / (2 * kt)
        ref = mp.quad(lambda t: 1 / (mp.log(t) * mp.log(2 * kt * t)), [3, 10, 100, 1000, 10**4, 10**5, Tt] if Tt > 10**5 else [3, Tt])
        assert abs(ref - I[kt - 1]) < 1e-8 * abs(ref), (kt, ref, I[kt - 1])
    print('x = %d: N = %.1f, 3|k share = %.4f' % (x, Nx, sh))

# direct count of primes q <= 1e7 with ord_q(2) an odd prime
X = 10**7
spf = np.zeros(X + 1, np.int32)
for p in primerange(2, int(X**0.5) + 1):
    idx = np.arange(p, X + 1, p)
    sel = spf[idx] == 0
    spf[idx[sel]] = p
spf[spf == 0] = np.arange(X + 1)[spf == 0]
isp = (spf == np.arange(X + 1)) ; isp[:2] = False
print('pi(1e7) from spf:', int(isp.sum()))
spf_l = spf.tolist()
cnt = 0; cnt1 = 0; cnt6 = 0
qs = np.nonzero(isp)[0].tolist()
for q in qs:
    if q == 2:
        continue
    n = q - 1
    while n % 2 == 0:
        n //= 2
    found = False
    while n > 1:
        p = spf_l[n]
        while n % p == 0:
            n //= p
        if pow(2, p, q) == 1:
            found = True
            break
    if found:
        cnt += 1
        if q % 3 == 1:
            cnt1 += 1
        if q <= 10**6:
            cnt6 += 1
print('direct: primes q <= 1e7 with ord_q(2) odd prime: %d, q = 1 mod 3: %d; (q <= 1e6: %d)' % (cnt, cnt1, cnt6))
# second method up to 1e6: sympy n_order
c2 = sum(1 for q in primerange(3, 10**6) if isprime(n_order(2, q)) and n_order(2, q) > 2)
print('check via sympy n_order, q <= 1e6:', c2)
