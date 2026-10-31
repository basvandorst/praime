# Conjecture 99: Why strong semiprimes outnumber primes by a factor ln 3

> **In short:** Alain Rochelli counted semiprimes p·q with p ≤ q < p³ ("strong semiprimes", A362910) and conjectured that near any large number there are ln 3 = 1.0986 of them for every prime. It is a consequence of the prime number theorem: the condition q < p³ means p lies between n^(1/4) and n^(1/2), and ∫ du/(u(1−u)) over that range is exactly ln 3. The same calculation done with the actual primes predicts his measured ratios: 1.0826 at 10⁹ (he found 1.0833) and 1.0929 at 10¹² (he found 1.0937). The values sit slightly below ln 3 because of Mertens' remainder at p ≈ n^(1/4), and the gap closes as n grows.

**Conjecture:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_099.htm)

**Status:** proved (Q2, for growing intervals, from the prime number theorem) · fixed-length windows explained, not provable · Q1 counts reproduced exactly

## The question

A strong semiprime is n = p·q with primes p ≤ q < p³. Let SP(A, B) count them in [A, B] and PI(A, B) count the primes. Rochelli conjectures that SP/PI tends to ln 3 = 1.098612…, and measured it on windows of length 2·10⁶:

```
A        SP       PI       R = SP/PI
10⁹      104445   96417    1.08326
10¹²     79196    72413    1.09367
10¹⁸     52850    48427    1.09133
10³⁰     31729    29054    1.09207
```

Q1 asks for checks, Q2 for a proof.

## What the page had vs ours

| | On the page | Ours |
|---|---|---|
| Q1 | four measurements, no contributions | rows 10⁹ and 10¹² reproduced exactly |
| Q2 | open | **proved** for intervals [A, cA] and for counts up to x; the exact prime sum predicts the measured values |

## How we found it

Fix n and ask which splittings n = p·q count. The condition p ≤ q means p ≤ √n. The condition q < p³ means n = pq < p⁴, so p > n^(1/4). So p runs over the primes between n^(1/4) and n^(1/2), and for each such p the cofactor q = n/p must be prime. That has chance about 1/ln(n/p), while n itself is prime with chance 1/ln n. The ratio is therefore a sum over p, and it turns into a simple integral.

## The math

The density of strong semiprimes near n, relative to the density of primes:

```
R(n) = ln n · Σ_{n^(1/4) < p ≤ n^(1/2)} 1 / (p · ln(n/p))
```

With t = n^u, the prime number theorem turns the sum into an integral:

```
Σ 1/(p ln(n/p)) ≈ ∫ dt / (t · ln t · ln(n/t)) = (1/ln n) ∫_{1/4}^{1/2} du / (u(1 − u)) = (1/ln n) · [ln(u/(1−u))]_{1/4}^{1/2} = ln 3 / ln n
```

so R(n) → ln 3. (For ordinary semiprimes the same integral runs from 0 to 1/2 and diverges like ln ln n. That is why all semiprimes outnumber primes by a growing factor, while strong semiprimes keep a fixed ratio.)

**A proof for counts up to x.** Let S(x) count strong semiprimes ≤ x. Primes p ≤ x^(1/4) contribute at most Σ π(p³) = O(x/ln² x), and the cases q = p add only π(√x). For x^(1/4) < p ≤ x^(1/2), the count is π(x/p) − π(p) + 1, and Σ π(p) = O(x/ln² x). The prime number theorem and partial summation give

```
Σ_{x^(1/4) < p ≤ x^(1/2)} π(x/p) = (ln 3 + o(1)) · x / ln x
```

So S(x)/π(x) → ln 3. The same argument works for intervals [A, cA] with fixed c > 1. For windows of fixed length like 2·10⁶, nothing can be proved, since that would need primes in short intervals. But the local formula R(n) is the natural prediction there too.

**Why the data sit below ln 3.** The sum over actual primes differs from the integral by Mertens' remainder at the lower end p ≈ n^(1/4), where the primes are smallest and the discreteness is largest. At n = 10⁹ the cutoff is p ≈ 178, and the exact sum gives:

```
A        measured R   exact prime sum R(A)   ln 3      noise (1 s.d.)
10⁹      1.08326      1.08263                1.09861   ±0.005
10¹²     1.09367      1.09290                1.09861   ±0.006
10¹⁸     1.09133      1.09810                1.09861   ±0.007
10³⁰     1.09207      1.09861                1.09861   ±0.009
```

At 10⁹ the measured deficit of 0.015 is 3 noise widths below ln 3, but it matches the exact sum to 0.0006. At 10¹⁸ and 10³⁰ the cutoff is so large that the sum equals ln 3. There the measured values are about one noise width below, which is ordinary scatter for counts of 30,000–50,000.

## Verify it yourself

This script recounts the window at 10⁹ exactly, checks the integral, and evaluates the exact prime sum at 10⁹ and 10¹². It needs sympy and gmpy2, and runs in a few seconds.

```python
import math, gmpy2
from sympy import primerange

def window(A, L=2*10**6):
    B = A + L
    PI = sum(1 for n in range(A | 1, B + 1, 2) if gmpy2.is_prime(n))
    SP = sum(1 for p in primerange(2, math.isqrt(B) + 1)
               for q in range(max(p, -(-A // p)), B // p + 1)
               if q < p**3 and gmpy2.is_prime(q))
    return SP, PI

def R(n):        # exact prime sum: ln n * sum over n^(1/4) < p <= n^(1/2) of 1/(p ln(n/p))
    return math.log(n) * sum(1/(p*math.log(n/p)) for p in primerange(int(n**0.25) + 1, math.isqrt(int(n)) + 1))

# Q1: Rochelli's first row, recounted.
assert window(10**9) == (104445, 96417)

# The integral: int_{1/4}^{1/2} du/(u(1-u)) = ln 3.
f = lambda u: math.log(u/(1 - u))
assert abs((f(0.5) - f(0.25)) - math.log(3)) < 1e-12

# The exact prime sum predicts the measured ratios (1.08326 and 1.09367).
assert abs(R(10**9 + 10**6) - 1.08263) < 1e-4 and abs(R(10**12 + 10**6) - 1.09290) < 1e-4
print("window 1e9 reproduced; predicted R = %.5f (1e9), %.5f (1e12); ln 3 = %.5f"
      % (R(10**9 + 10**6), R(10**12 + 10**6), math.log(3)))
```
