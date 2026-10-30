# Conjecture 55: Sand's "p² + q + d" strengthening is false (first counterexample 305)

> **In short:** On the page of Conjecture 55, Werner Sand proposed a stronger form: every n > 5 is p² + q + d with p, q prime and d ∈ {0, 1, 2}, "checked until n = 10⁶". It fails already at **n = 305**, and 72,048 numbers below 10⁶ fail. The reason is a remainder mod 6. If n ≡ 5 (mod 6) and p ≥ 5, then n − p² − d is ≡ 4, 3 or 2 (mod 6), so it can only be prime if it equals 2 or 3. Such n can only be written using p = 2 or p = 3, which leaves just three chances to hit a prime. So almost every n ≡ 5 (mod 6) fails. Sand's original conjecture 55 (with any square m²) is not affected.

**Conjecture:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_055.htm)

**Status:** variant disproved (Sand's p² + q + d version, Oct 29 contribution) · infinitely many counterexamples (proof) · smallest counterexample 305 (exhaustive) · the main conjecture 55 stays open

## The question

Conjecture 55 says every n > 2 is m² + p + d with m ≥ 1, p prime, d ∈ {0, 1}. Later on the same page, Sand refined it, restricting the square to the square of a prime:

> Every natural number n > 5 is the sum of a square prime, a prime, and one of the values 0, 1, 2: n = p² + q + d, d out of (0, 1, 2). I checked this until n = 10⁶.

He also counted 47 such representations of n = 10⁶.

## What the page had vs ours

| | On the page | Ours |
|---|---|---|
| p² + q + d claim | "checked until n = 10⁶" | **false**: 305, 325, 329, 479, 521, …; 72,048 failures below 10⁶ |
| Why | none | n ≡ 5 (mod 6) leaves only three chances; almost all such n fail |
| Representations of 10⁶ | 47 | 47 (same reading of the rule) |

## How we found it

Squares of primes p ≥ 5 are all ≡ 1 (mod 6), because p ≡ ±1 (mod 6). So for n ≡ 5 (mod 6) every candidate n − p² − d lands in a fixed residue class. Working these out:

```
n ≡ 5 (mod 6), p ≥ 5:   n − p² ≡ 4,  so  n − p² − d ≡ 4, 3, 2  (mod 6)  for d = 0, 1, 2
```

A prime ≡ 4, 3 or 2 (mod 6) can only be 2 or 3. So p ≥ 5 is useless unless q ≤ 3, and everything rests on p = 2 and p = 3.

## The math

For n ≡ 5 (mod 6), n > 11, the candidates for q are:

```
p = 2:  n − 4, n − 5, n − 6    ≡ 1, 0, 5 (mod 6)   → only n − 4 and n − 6 can be prime
p = 3:  n − 9, n − 10, n − 11  ≡ 2, 1, 0 (mod 6)   → only n − 10 can be prime
p ≥ 5:  q = n − p² − d ∈ {2, 3}  → only n = p² + 4
```

So **n ≡ 5 (mod 6) is representable iff n − 4, n − 6 or n − 10 is prime, or n − 4 is the square of a prime.** By the prime number theorem the chance that none of n − 4, n − 6, n − 10 is prime tends to 1. So almost every n ≡ 5 (mod 6) is a counterexample, and there are infinitely many. At n = 305, the numbers 301 = 7·43, 299 = 13·23 and 295 = 5·59 are composite, and 301 is not a square, so 305 fails.

The other residue classes are fine for large n. For example, for n ≡ 1 (mod 6), d = 1 and any p ≥ 5 give a candidate ≡ 5 (mod 6). They only produce a few small failures (325, 615, …).

Sand's count of 47 representations of 10⁶ is exactly what we get with this reading of his rule, so the failures are not due to a different interpretation. His "checked until 10⁶" must have contained an error.

## Verify it yourself

This script finds all failures below 10⁶, confirms that 305 is the first, recounts the 47 representations of 10⁶, and checks the mod-6 criterion on every n ≡ 5 (mod 6) up to 10⁶. It runs in a few seconds.

```python
import numpy as np

LIM = 10**6
s = np.ones(LIM + 1, dtype=bool); s[:2] = False
for i in range(2, int(LIM**0.5) + 1):
    if s[i]: s[i*i::i] = False
psq = [p*p for p in range(2, int(LIM**0.5) + 1) if s[p]]

ok = np.zeros(LIM + 1, dtype=bool)              # ok[n]: n = p^2 + q + d for some primes p, q and d in {0,1,2}
for P2 in psq:
    for d in (0, 1, 2):
        ok[P2 + d:] |= s[:LIM + 1 - P2 - d]
fails = [n for n in range(6, LIM + 1) if not ok[n]]
assert fails[:5] == [305, 325, 329, 479, 521] and len(fails) == 72048

# Sand's own count for n = 10^6
assert sum(1 for P2 in psq for d in (0, 1, 2) if LIM - P2 - d >= 2 and s[LIM - P2 - d]) == 47

# n = 5 (mod 6): representable iff n-4, n-6 or n-10 is prime, or n-4 is a prime square
sq = set(psq)
assert all(bool(ok[n]) == bool(s[n-4] or s[n-6] or s[n-10] or (n - 4) in sq) for n in range(17, LIM + 1, 6))
print("305 is the first of", len(fails), "counterexamples below 10^6")
```
