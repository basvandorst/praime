# Conjecture 114: Rotondo's semiprime sequence, extended to 200 terms (attempt)

> **In short:** For f(n+1) = f(n)·R(n), with R(n) the least composite c making f(n) + c prime, it was conjectured from 48 terms that every R(n) is a semiprime, no prime repeats, and every prime appears. Proved here: a prime shared by c and f(n) divides f(n) + c, so no prime repeats. The other two claims stay open; 200 terms confirm them.

**Conjecture:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_114.htm)

**Status:** open · claim 2 proved · Q1: extended to 200 terms · claims 1 and 3 explained heuristically

## The question

```
f(1) = 1,   R(n) = the least composite c with f(n) + c prime,   f(n+1) = f(n) · R(n)
R = 4, 9, 25, 77, 323, 713, 169, 1189, 1961, 1849, …
```

Rotondo checked 48 terms and conjectured: (1) every R(n) is a semiprime; (2) the prime factors of different R(n) never repeat, so a prime that appears twice does so as a square p²; (3) every prime eventually appears. Q1 asks for more terms, Q2 for a proof.

## What the page had vs ours

| | On the page | Ours |
|---|---|---|
| Terms computed | 48 (all primes below 491 used) | **200** (all primes below 2309 used) |
| Claim 2 (no reuse) | observed | **proved** |
| Claims 1 and 3 | observed | true in all 200 terms; heuristic explanation below |

## How we found it

Claim 2 is a divisibility fact. If a prime q divides both c and f(n), it divides f(n) + c, and f(n) + c > q, so f(n) + c is not prime. So R(n) is always coprime to f(n), the product of all earlier R's. Its prime factors are therefore new, and a prime can only occur twice inside the same R(n), as a square.

That turns the search into a clean picture. Let q₁ < q₂ < … be the primes not yet used. The candidates are the composites built from these primes, in increasing order: q₁², q₁q₂, q₁q₃, …, q₂², …. Every such composite below q₁³ is a semiprime.

## The math

**Claim 1, heuristically.** The number f(n) + c has no prime factor dividing f(n), which covers all primes below q₁ apart from a few stragglers. So by Mertens it is prime with probability about e^γ·ln q₁ / ln f(n). At n = 200 that is roughly 1.78 · 7.7 / 2620 ≈ 1/190 per candidate. There are about q₁²/ln q₁, a few hundred thousand, semiprime candidates below q₁³. The chance that all of them fail is less than e^(−2000). Summed over all n this stays negligible, so a non-semiprime R(n) should never occur. A proof would need a prime in a very thin, explicit set of numbers, which is far beyond current methods.

**Claim 3, heuristically.** The first success is among the first ~190 candidates. Those are products (q₁ + a)(q₁ + b) with a + b small. A fraction of order 1/√190 of them contains the smallest unused prime q₁. So q₁ is used at each step with probability roughly 5–10%, falling only slowly as f grows. The sum of these probabilities diverges, so q₁ is eventually used, and with it every prime. The data agree: the smallest unused prime moves up steadily, to 491 after 48 terms and 2309 after 200.

**Squares.** R(n) = q₁² exactly when the very first candidate works. 26 of the 200 terms are squares: 2², 3², 5², 13², 43², 167², 241², 251², ….

## The numbers

```
terms 1–10:     4 9 25 77 323 713 169 1189 1961 1849            (as on the page)
terms 49–52:    198691 195911 277931 177241
term 200:       5947937 = 2411 · 2467                           f(200) has 1139 digits
all 200 terms semiprime, no prime factor repeated
smallest primes not yet used after 200 terms: 2309, 2333, 2477, 2503, …
```

## Validation

This script recomputes the first 80 terms with PRP tests. It checks that they start as on the page, are all semiprimes, and never reuse a prime. It also checks claim 2 directly: whenever c shares a factor with f, f + c is composite. It runs in a few seconds. The 200-term run uses the same loop.

```python
import gmpy2
from sympy import factorint

f, used, R = gmpy2.mpz(1), set(), []
for n in range(1, 81):
    c = 4
    while True:
        if not gmpy2.is_prime(c):
            if gmpy2.gcd(c, f) > 1:
                assert not gmpy2.is_prime(f + c)      # claim 2: a shared factor divides f + c
            elif gmpy2.is_prime(f + c, 30):
                break
        c += 1
    fac = factorint(int(c))
    assert sum(fac.values()) == 2                     # claim 1 holds for these terms
    assert not (set(fac) & used)                      # no reuse
    used |= set(fac); R.append(int(c)); f *= c
assert R[:10] == [4, 9, 25, 77, 323, 713, 169, 1189, 1961, 1849]
print("80 terms: all semiprimes, no prime reused; f(81) has", len(str(f)), "digits")
```
