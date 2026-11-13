"""Part A: run the classical prime formulas and measure what each one needs and where it breaks."""
import math, time
import os
from fractions import Fraction
import mpmath as mp
from sympy import isprime, prime, primerange, divisors, mobius, primorial, factorial, nextprime

out = []
def say(s):
    print(s); out.append(s)

# 1. Wilson: (n-1)! = -1 mod n iff n prime; cost n-2 multiplications mod n
for n in (101, 1009, 10007):
    t = time.time(); f = 1
    for j in range(2, n): f = f * j % n
    say(f'Wilson n={n}: (n-1)! mod n = {f} (= n-1: {f == n-1}), {n-2} multiplications, {1e3*(time.time()-t):.1f} ms')

# 2. Willans (1964): p_n = 1 + sum_{i=1}^{2^n} floor((n / sum_{j<=i} F(j))^(1/n)),  F(j) = floor(cos^2(pi ((j-1)!+1)/j))
def F(j):
    return 1 if (math.factorial(j - 1) + 1) % j == 0 else 0   # exact value of floor(cos^2(...))
for n in range(1, 6):
    s, acc = 0, 0
    for i in range(1, 2 ** n + 1):
        acc += F(i)
        s += 1 if acc <= n else 0   # floor((n/acc)^(1/n)) = 1 iff acc <= n, else 0 (1 <= acc, n^(1/n) < 2)
    say(f'Willans n={n}: p_n = 1 + {s} = {1 + s} (true {prime(n)}); {2**n} outer terms, largest factorial ({2**n - 1})! has {len(str(math.factorial(2**n - 1)))} digits')

# 3. Mills (1947): floor(A^(3^n)) prime; A from the primes themselves (Caldwell-Cheng, assuming RH)
mp.mp.dps = 60
A = mp.mpf('1.3063778838630806904686144926026057129167845851567136443680537599664340537668265988215014037011973957072')
for n in range(1, 5):
    v = mp.floor(A ** (3 ** n))
    say(f'Mills n={n}: floor(A^(3^{n})) = {mp.nstr(v, 30)} ({len(mp.nstr(v, 100).split(".")[0])} digits), prime: {isprime(int(v))}')
for digits in (10, 20, 30):
    mp.mp.dps = 60
    Ad = mp.mpf(mp.nstr(A, digits))
    vals = [int(mp.floor(Ad ** (3 ** n))) for n in range(1, 5)]
    ok = [isprime(v) and v == int(mp.floor(A ** (3 ** n))) for n, v in zip(range(1, 5), vals)]
    say(f'Mills with A to {digits} digits: correct terms n=1..4: {ok}')

# 4. Gandhi (1971): p_n = floor(1 - log2(-1/2 + sum_{d | P} mu(d)/(2^d - 1))), P = p_1 ... p_{n-1}
for n in range(2, 8):
    P = int(primorial(n - 1))
    ds = divisors(P)
    s = sum(Fraction(int(mobius(d)), 2 ** d - 1) for d in ds)
    val = 1 - math.log2(float(-Fraction(1, 2) + s)) if float(-Fraction(1, 2) + s) > 0 else float('nan')
    bits = max(d for d in ds)
    say(f'Gandhi n={n}: p_n = {math.floor(val)} (true {prime(n)}); {len(ds)} divisors of {P}; exact arithmetic needs 2^{bits}-sized denominators')

# 5. Rowland (2008): a(n) = a(n-1) + gcd(n, a(n-1)), a(7) = 7; differences are 1 or prime
a, n, firsts, steps = 7, 1, {}, 0
while len(firsts) < 12 and n < 10 ** 6:
    n += 1; g = math.gcd(n, a); a += g
    if g > 1 and g not in firsts: firsts[g] = n
say('Rowland: first appearance of each new prime difference (prime: step): ' + ', '.join(f'{p}: {s}' for p, s in sorted(firsts.items(), key=lambda t: t[1])))

# 6. Euler n^2 + n + 41; any polynomial P has P(n + k P(n)) = 0 mod P(n)
bad = next(n for n in range(100) if not isprime(n * n + n + 41))
say(f'Euler n^2+n+41: first composite at n = {bad}: {bad*bad+bad+41} = 41^2')

# 7. Prime-representing constant (Fridman et al. 2019): f = sum (p_n - 1)/(p_1...p_{n-1}); f_{k+1} = floor(f_k)(f_k - floor(f_k) + 1)
mp.mp.dps = 400
ps = list(primerange(2, 600))
f = mp.mpf(0); P = mp.mpf(1)
for p in ps:
    f += (p - 1) / P; P *= p
say(f'Fridman constant = {mp.nstr(f, 25)}: built from the primes it outputs')
for digits in (20, 50, 100, 200):
    mp.mp.dps = digits + 10
    x = mp.mpf(mp.nstr(f, digits)); got = []
    for _ in range(80):
        q = int(mp.floor(x)); got.append(q); x = q * (x - q + 1)
    good = 0
    for g, p in zip(got, ps):
        if g == p: good += 1
        else: break
    say(f'  with {digits} digits: first {good} primes correct (up to {ps[good-1]}); then {got[good:good+3]}')

# 8. Golomb (1970), see Part I of the explainer: needs every prime <= p and ~p log2 p bits
say('Golomb: p = 1499 needs s = 1515 and 16,587-bit integers, and all primes <= 1499 as input (explainer, Part I)')

# 9. Truncated explicit formula: psi(x) from the first K zeros; error at x = 100.5 and 1000.5
import numpy as np
G = np.fromfile(os.environ.get('WORK', os.getcwd()) + '/z2m.f64')
def psi_true(x):
    s = 0.0
    for p in primerange(2, int(x) + 1):
        q = p
        while q <= x: s += math.log(p); q *= p
    return s
for x in (100.5, 1000.5):
    u = math.log(x); line = []
    for K in (10, 100, 1000, 10000, 100000):
        g = G[:K]
        acc = float(np.sum((0.5 * np.cos(g * u) + g * np.sin(g * u)) / (0.25 + g * g)))
        v = x - math.log(2 * math.pi) - 0.5 * math.log(1 - 1 / x ** 2) - 2 * math.sqrt(x) * acc
        line.append(f'K={K}: {v - psi_true(x):+.4f}')
    say(f'truncated explicit formula, psi({x}) error: ' + ', '.join(line))

open(os.environ.get('WORK', os.getcwd()) + '/existing_formulas.txt', 'w').write('\n'.join(out) + '\n')
