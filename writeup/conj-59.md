# Conjecture 59: n! + p_n is never a square (attempt; odd powers stay open)

> **In short:** Alexander Povolotsky conjectured that n! + p_n is never a perfect power m^k (k > 1), where p_n is the n-th prime. Two AI "proofs" on the page were shown to be wrong. We prove the case of even k: n! + p_n is never a perfect square. A square would force every number up to n to be a quadratic residue mod p_n, but the least non-residue of a prime p is always below √p + 1, and p_n is far smaller than n². Odd exponents (cubes, fifth powers, …) remain open. We also checked that no n ≤ 3000 gives any perfect power.

**Conjecture:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_059.htm)

**Status:** open · even exponents proved · all k checked for n ≤ 3000

## The question

Prove that n! + p_n = m^k has no solution with k > 1, or find one. The page (and [conjecture 98](https://www.primepuzzles.net/conjectures/conj_098.htm)) discusses two ChatGPT "proofs". J. R. Howell and others showed they only rule out even m.

## What the page had vs ours

| | On the page | Ours |
|---|---|---|
| Even k (squares) | open; AI proofs invalid | **proved impossible** |
| Odd k | open | open; every prime factor of m exceeds n, and k < n |
| Search | n < 200,000 (Casini, conj. 98) | n ≤ 3000, every k (re-check) |

## How we found it

If n! + p_n = s², then reducing mod a small prime q ≤ n kills n!, so p_n must be a square mod q. That is a condition at *every* prime up to n at once. Quadratic reciprocity turns it around: every prime up to n is then a square mod p_n. But primes always have a small non-residue, and that clash gives the proof.

## The math

**Theorem.** n! + p_n is not a perfect square for any n ≥ 1. Hence it is not m^k for any even k.

*Proof.* For n ≤ 3: 3, 5 and 11 are not squares. Let n ≥ 4 and suppose n! + p_n = s².

1. 8 divides n!, and p_n is odd, so s is odd and s² ≡ 1 (mod 8). Hence p_n ≡ 1 (mod 8).
2. For every odd prime q ≤ n: q | n!, so p_n ≡ s² (mod q). Also q ∤ p_n, because p_n > n ≥ q. So the Legendre symbol (p_n/q) = 1.
3. Since p_n ≡ 1 (mod 4), reciprocity gives (q/p_n) = (p_n/q) = 1 for every odd prime q ≤ n. Since p_n ≡ 1 (mod 8), (2/p_n) = 1 too. So every prime ≤ n is a quadratic residue mod p_n, and therefore so is every integer 1 … n.
4. So the least quadratic non-residue n₀ of p_n is at least n + 1. Classical lemma: n₀(n₀ − 1) < p for every odd prime p. (Let m = ⌊p/n₀⌋ + 1, so 0 < m·n₀ − p < n₀. Then m·n₀ − p is a residue, hence m is a non-residue, hence m ≥ n₀, which gives p > n₀(n₀ − 1).) Therefore p_n > n(n + 1).
5. But p_n < n(n + 1) for every n ≥ 2. For n ≥ 6 this follows from p_n < n(ln n + ln ln n) (Rosser–Schoenfeld); 3, 5, 7 and 11 check it for n = 2…5. Contradiction. ∎

**What is left: odd k.** If n! + p_n = m^k, every prime q ≤ n has n! + p_n ≡ p_n ≢ 0 (mod q), so every prime factor of m exceeds n. Then m^k ≥ (n + 1)^k, while n! + p_n < n^n, so k < n. For odd k the same idea gives "p_n is a k-th power residue mod every prime q ≤ n with q ≡ 1 (mod k)". Primes q ≢ 1 (mod k) give no information. There is no reciprocity law strong enough to turn this into a contradiction, so cubes and higher odd powers remain open.

## Verify it yourself

This script checks every n ≤ 3000 for perfect powers of any exponent. It also checks the two numerical facts in the proof: p_n < n(n + 1) for n ≤ 10⁵, and n₀(n₀ − 1) < p for all primes below 2·10⁵. It runs in a few seconds.

```python
import gmpy2
from sympy import primerange, legendre_symbol

P = list(primerange(2, 1_400_000))
f = gmpy2.mpz(1)
for n in range(1, 3001):
    f *= n
    assert not gmpy2.is_power(f + P[n - 1]), n          # no m^k, k >= 2

assert all(P[n - 1] < n * (n + 1) for n in range(2, 100001))
for p in primerange(3, 200000):
    n0 = next(a for a in range(2, p) if legendre_symbol(a, p) == -1)
    assert n0 * (n0 - 1) < p
print("no perfect powers for n <= 3000; both lemmas hold on the tested ranges")
```
