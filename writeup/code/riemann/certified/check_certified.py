"""Independent check (mpmath, sympy) of the certified chain.

1. The certified primes equal sympy's primes below X (none missing, none extra).
2. The identity and the error bound at random sample points: zero side from the 28-digit zero strings,
   evaluated in mpmath with mpmath's own Bessel functions; prime side exact; |difference| <= delta.
3. Two certificates rebuilt from scratch in mpmath (own window logic): every feasible 0/1 pattern agrees
   on the candidate, and the true pattern is feasible.
"""
import json, math, random, sys
import os
import mpmath as mp
from sympy import primerange, isprime, factorint, perfect_power
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cert import sample   # only for delta (the rigorous budget under test)

D = os.environ.get('WORK', os.getcwd()) + '/'
mp.mp.dps = 30
Zs = [mp.mpf(l.strip()) for l in open(D + 'hp2419.txt')]
Gf = np.load(D + 'hp2419.npy')

def Kmp(v, a, beta, m):
    t = 1 - (v / a) ** 2
    return t ** (mp.mpf(m) / 2) * mp.besseli(m, beta * mp.sqrt(t)) if t > 0 else mp.mpf(0)

def Khmp(w, a, beta, m):
    nu = m + mp.mpf(1) / 2; s = abs(a * w)
    if s < beta:
        q = mp.sqrt(beta ** 2 - s ** 2); z = mp.besseli(nu, q) / q ** nu if q > 0 else 1 / (2 ** nu * mp.gamma(nu + 1))
    else:
        y = mp.sqrt(s ** 2 - beta ** 2); z = mp.besselj(nu, y) / y ** nu if y > 0 else 1 / (2 ** nu * mp.gamma(nu + 1))
    return a * mp.sqrt(2 * mp.pi) * beta ** m * z

def Pmp(x, a, beta, m):
    nu = m + mp.mpf(1) / 2; r = a / (2 * x); q = mp.sqrt(beta ** 2 + r ** 2)
    return a * mp.sqrt(2 * mp.pi) * beta ** m * mp.besseli(nu, q) / q ** nu

def lam(n):
    f = factorint(n); return mp.log(next(iter(f))) if len(f) == 1 else mp.mpf(0)

def zero_side(x, a, beta, m, T):
    x = mp.mpf(x); L = mp.log(x); s = mp.mpf(0)
    for g in Zs:
        if g > T: break
        s += mp.cos(g * L) * Khmp(g / x, a, beta, m)
    triv = mp.quad(lambda v: Kmp(v, a, beta, m) * mp.sqrt(x * mp.e ** (v / x)) / (x * ((x * mp.e ** (v / x)) ** 3 - x * mp.e ** (v / x))), [-a, 0, a])
    return Pmp(x, a, beta, m) - 2 / mp.sqrt(x) * s - mp.sqrt(x) * triv

def prime_side(x, a, beta, m):
    x = mp.mpf(x); lo = int(mp.floor(x * mp.e ** (-a / x))); hi = int(mp.ceil(x * mp.e ** (a / x)))
    return sum(lam(n) * mp.sqrt(x / n) * Kmp(x * mp.log(n / x), a, beta, m) for n in range(max(2, lo), hi + 1))

out = []
def say(s):
    print(s); out.append(s)

# 1
for f in ('chain10000_low.json', 'chain10000_hp.json'):
    if not os.path.exists(D + f): continue
    J = json.load(open(D + f))
    got = J['primes']; true = list(primerange(3, J['X']))
    say(f'{f}: {len(got)} certified odd primes below {J["X"]}; equal to sympy list: {got == true}; stop = {J["result"]["stop"]}; zeros up to {J["result"]["maxT"]:.2f}')

# 2
random.seed(1)
pts = [(random.uniform(30, 190), 1.9, 6, 2), (random.uniform(300, 990), 12, 12, 2), (random.uniform(1200, 2900), 60, 14, 2), (random.uniform(3200, 9900), 130, 16, 2), (random.uniform(3200, 9900), 130, 16, 2)]
for x, a, beta, m in pts:
    T = 1.2 * beta * x / a
    zs = zero_side(x, a, beta, m, T); ps = prime_side(x, a, beta, m)
    _, d, _ = sample(float(x), Gf, T, a, beta, m, 1e-12)
    say(f'x = {x:.3f} (a = {a}, beta = {beta}, T = {T:.1f}): zero side {mp.nstr(zs, 15)}, prime side {mp.nstr(ps, 15)}, |diff| = {mp.nstr(abs(zs - ps), 3)}, delta = {d:.3e}, inside: {abs(zs - ps) <= d}, diff/delta = {float(abs(zs-ps))/d:.1e}')

# 3
def certificate(c, a, beta, m, k):
    def unknown(n): return n % 2 == 1 and not perfect_power(n) and n >= c
    U = []; n = c
    while len(U) < k + 1:
        if unknown(n): U.append(n)
        n += 2
    U, w = U[:-1], U[-1]
    xs = [x for x in range(int(c - a - 1), w) if x * math.exp(-a / x) > 1 and abs(x * math.log(c / x)) < a and x * math.log(w / x) >= a]
    A, y, d = [], [], []
    for x in xs:
        T = 1.2 * beta * x / a
        zs = zero_side(x, a, beta, m, T)
        xm = mp.mpf(x); lo = int(mp.floor(xm * mp.e ** (-a / xm))); hi = int(mp.ceil(xm * mp.e ** (a / xm)))
        known = sum(lam(n) * mp.sqrt(xm / n) * Kmp(xm * mp.log(n / xm), a, beta, m) for n in range(max(2, lo), hi + 1) if n not in U)
        A.append([float(mp.log(u) * mp.sqrt(xm / u) * Kmp(xm * mp.log(u / xm), a, beta, m)) for u in U])
        y.append(float(zs - known)); d.append(sample(float(x), Gf, T, a, beta, m, 1e-12)[1])
    A = np.array(A); y = np.array(y); d = np.array(d)
    sols = []
    def dfs(j, b, acc):
        if np.any(acc > y + d) or np.any(acc + A[:, j:].sum(axis=1) < y - d): return
        if j == len(U): sols.append(tuple(b)); return
        for v in (0, 1): dfs(j + 1, b + [v], acc + v * A[:, j])
    dfs(0, [], np.zeros(len(xs)))
    truth = tuple(1 if isprime(u) else 0 for u in U)
    return len(xs), len(sols), {s[0] for s in sols}, truth in sols
for (c, a, beta, m, k) in [(1019, 12, 12, 2, 6), (9973, 130, 16, 2, 65)]:
    ns, nsol, firsts, ok = certificate(c, a, beta, m, k)
    say(f'certificate for {c} (a = {a}, k = {k}): {ns} samples, {nsol} feasible patterns, values of the first bit: {sorted(firsts)} (truth: {1 if isprime(c) else 0}); true pattern feasible: {ok}')
open(D + 'check_certified.txt', 'w').write('\n'.join(out) + '\n')
