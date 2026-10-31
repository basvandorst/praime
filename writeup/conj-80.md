# Conjecture 80: Why the numerator always ends in 25 (proof)

> **In short:** It was claimed that the reduced numerator of ∏(p² + 1)/(p − 1) over the first n primes ends in 25 for n > 2, but shown only for large n. The proof: (1) the numerator is ≡ 1 (mod 4), so 25 must divide it; (2) primes ≡ ±2 (mod 5) add factors 5, primes ≡ 1 remove them; (3) explicit bounds settle p_n ≥ 8·10⁹, a sieve the rest. This proves it for every n ≥ 3.

**Conjecture:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_080.htm) (from [Puzzle 900](https://www.primepuzzles.net/puzzles/puzz_900.htm))

**Status:** proved (for every n ≥ 3; explicit bounds from the literature plus an exhaustive computation to 8·10⁹)

## The question

```
x(n) = ∏_{i=1}^{n} (p_i² + 1) / (p_i − 1)
n = 3:  325/2      n = 4:  8125/6      n = 5:  99125/6      n = 6:  8425625/36     n = 7:  1221715625/288
```

Why does the numerator of x(n), in lowest terms, always end in 25 for n > 2? Vantieghem checked n < 500,000 and called it a conjecture. Q: give a complete proof or a counterexample.

## What the page had vs ours

| | On the page | Ours |
|---|---|---|
| Numerator odd, ≡ 1 (mod 4) | shown (Vantieghem) | used, with the details below |
| Divisible by 25 | "certainly true for n big enough", checked to n = 500,000 | **proved for all n ≥ 3** |

## How we found it

Write A = ∏(p² + 1) and B = ∏(p − 1). The reduced numerator is A/gcd(A, B), and its power of 5 is W(n) = v₅(A) − v₅(B), when that is positive.

- A prime ending in 3 or 7 (p ≡ ±2 mod 5) has 5 | p² + 1, so it adds to W.
- A prime ending in 1 (p ≡ 1 mod 5) has 5 | p − 1, so it subtracts.
- Primes ending in 9 change nothing.

Primes are spread evenly over the four classes mod 5. Counting higher powers of 5 too, the average step is 5/8 − 5/16 = +5/16 per prime. So W climbs steadily: W(3) = 2, W(4) = 4, W(5) = 3, W(6) = 4, W(7) = 5, …, W(10⁸) ≈ 1.8 million. Proving it for every n needs two things: an explicit bound that makes "evenly spread" precise from some point on, and a computation up to that point.

## The math

**Step 1: the numerator is odd and ≡ 1 (mod 4).** For odd p, p² + 1 ≡ 2 (mod 8), so v₂(A) = n − 1 (the factor for p = 2 is 5). Every odd prime contributes at least one factor 2 to B, and p₃ = 5 contributes two, so v₂(B) ≥ n for n ≥ 3. Hence the reduced numerator N is odd. Every odd prime q dividing p² + 1 has −1 as a square mod q, so q ≡ 1 (mod 4). As a product of such primes, N ≡ 1 (mod 4).

**Step 2: reduction to the prime race.** A number ≡ 1 (mod 4) ends in 25 exactly when it is divisible by 25. So the claim is equivalent to

```
W(n) = Σ_{i≤n, p_i ≡ ±2 (mod 5)} v₅(p_i² + 1)  −  Σ_{i≤n, p_i ≡ 1 (mod 5)} v₅(p_i − 1)  ≥  2      for all n ≥ 3.
```

**Step 3: small n, by computation.** W(3) = 2 (the primes 2 and 3 each give one factor 5). From p = 7 on, a segmented sieve adds up the steps of all 367,783,651 primes up to 8·10⁹. The running total never drops below its starting value, so W(n) ≥ 2 whenever p_n ≤ 8·10⁹. The sieve's prime count matches π(8·10⁹) = 367,783,654. At 8·10⁹, W = 114,937,685.

**Step 4: large n, by explicit bounds.** Three published inequalities are used:

```
(BMOR)  |θ(x; q, a) − x/φ(q)| < x/(160 ln x)     for 3 ≤ q ≤ 10⁵, x ≥ 8·10⁹
        (Bennett, Martin, O'Bryant, Rechnitzer, Illinois J. Math. 62 (2018))
(RS)    π(y) < 1.25506 y/ln y                    for y > 1   (Rosser–Schoenfeld 1962)
(MV)    π(x; q, a) < 2x/(φ(q) ln(x/q))           for q < x   (Montgomery–Vaughan 1973)
```

Let x = p_n ≥ 8·10⁹. Each prime ≡ ±2 (mod 5) adds at least 1. Each prime ≡ 1 (mod 5) subtracts v₅(p − 1) = Σ_k [5^k | p − 1]. So

```
W ≥ [π(x;5,2) + π(x;5,3)] − π(x;5,1) − Σ_{k≥2} π(x; 5^k, 1).
```

- **Plus side:** π(x;5,2) + π(x;5,3) ≥ (θ(x;5,2) + θ(x;5,3))/ln x > (x/2 − x/(80 ln x))/ln x.
- **Primes ≡ 1 (mod 5):** with y = x^0.8, π(x;5,1) ≤ π(y) + θ(x;5,1)/ln y < 1.25506·x^0.8/(0.8 ln x) + (x/4 + x/(160 ln x))/(0.8 ln x).
- **Higher powers of 5:** use (MV) for 2 ≤ k ≤ 6. For k ≥ 7 with 5^k ≤ √x, use (MV) with ln(x/5^k) ≥ ½ ln x. For 5^k > √x, use the trivial bound x/5^k < √x, with at most log₅ x such k.

Divided by x/ln x, the lower bound is

```
c(x) = ½ − 1/(80 L) − (¼ + 1/(160 L))/0.8 − 1.25506·x^(−0.2)/0.8 − Σ_{k=2}^{6} L/(2·5^(k−1)·(L − k ln 5)) − 1/(4·5⁵) − L² /(ln 5·√x),
```

with L = ln x. Every negative term shrinks as x grows, so c(x) is increasing, and c(8·10⁹) = 0.0176. Hence W ≥ 0.0176·x/ln x > 6·10⁶ for all x ≥ 8·10⁹.

Steps 3 and 4 together give W(n) ≥ 2 for every n ≥ 3, so the numerator always ends in 25. ∎

**Why it was hard to see.** The step +5/16 per prime is an average. Locally, primes ending in 1 can run ahead of those ending in 3 or 7 for a while. Before the explicit bound one could only say "for n large enough". With the bound, "large enough" becomes 8·10⁹, which a sieve covers in 12 seconds.

## Validation

This script checks Steps 1–2 directly on the reduced fractions for n ≤ 400 (numerator ≡ 1 mod 4, ends in 25, v₅ equals W). It computes W for all primes below 10⁷ (minimum 2, at n = 3), and evaluates the Step 4 bound and its monotonicity. It runs in about 15 seconds. The computation to 8·10⁹ is [code/conj-80.c](code/conj-80.c) (`./conj80 7 8000000001`, 12 seconds; it prints `minprefix=0`).

```python
import math
from fractions import Fraction
import numpy as np
from sympy import prime, multiplicity

def v5(m):
    return multiplicity(5, m)

# Steps 1-2 on the actual fractions: numerator = 1 (mod 4), ends in 25, and v5(numerator) = W(n).
x, W = Fraction(1), 0
for n in range(1, 401):
    p = prime(n)
    x *= Fraction(p*p + 1, p - 1)
    W += v5(p*p + 1) - v5(p - 1)
    if n >= 3:
        N = x.numerator
        assert N % 4 == 1 and N % 100 == 25 and v5(N) == W, n

# Step 3 (small part): W over all primes below 10^7.
M = 10**7
s = np.ones(M, dtype=bool); s[:2] = False
for i in range(2, int(M**0.5) + 1):
    if s[i]: s[i*i::i] = False
W, minW = 0, None
for n, p in enumerate(np.nonzero(s)[0].tolist(), start=1):
    r = p % 5
    if r == 1:
        W -= v5(p - 1)
    elif r in (2, 3):
        W += v5(p*p + 1)
    if n >= 3:
        minW = W if minW is None else min(minW, W)
assert minW == 2

# Step 4: the explicit lower bound c(x) for W / (x / ln x), positive at 8e9 and increasing.
def c(x):
    L = math.log(x)
    return (0.5 - 1/(80*L) - (0.25 + 1/(160*L))/0.8 - 1.25506*x**-0.2/0.8
            - sum(L/(2*5**(k-1)*(L - k*math.log(5))) for k in range(2, 7)) - 1/(4*5**5) - L*L/(math.log(5)*math.sqrt(x)))
assert c(8e9) > 0.0175 and all(c(8e9*1.01**(i+1)) > c(8e9*1.01**i) for i in range(5000))
print("reduction checked to n = 400; min W below 1e7 is 2; c(8e9) = %.4f" % c(8e9))
```
