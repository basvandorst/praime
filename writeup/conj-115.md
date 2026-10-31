# Conjecture 115: The prime 3917 fails Rotondo's triangle test

> **In short:** Davide Rotondo conjectured that an odd number passes his "triangle test" exactly when it is a prime or 9. That is false: the prime **3917** fails. For a prime p the test only depends on a few numbers modulo p², which makes it cheap to scan, and 3917 is the smallest odd prime that fails. A probability argument says infinitely many primes fail, but they are rare: about eight more below 10¹⁰⁰. We also prove all three of Rotondo's valuation formulas, and with them that no odd prime power other than 9 passes.

**Conjecture:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_115.htm)

**Status:** disproved (counterexample, proven prime) · Q2: valuation formulas proved, prime powers settled · smallest counterexample (exhaustive check of all odd primes below it)

## The question

Let G(m) = 1^m + 2^m + … + m^m. For n > 1 and 1 ≤ k ≤ n/2, reduce the fraction

```
T(n, k) = G(n) / (G(k) · G(n − k))
```

to lowest terms. The number n **passes** if n divides the numerator for every k. Rotondo's main conjecture: an odd n > 1 passes if and only if n is a prime or n = 9. He tested every n up to 3000 and every odd prime power up to 100,000. Q1 asks to check his ranges; Q2 asks for ideas: prove his p-adic valuation formulas, prove the conjecture, or find a faster way to run the test.

## What the page had vs ours

| | On the page | Ours |
|---|---|---|
| Main conjecture | tested to n = 3000, no contributions | **false: the prime 3917 fails** (at k = 500); every odd prime below 3917 passes, and none between 3917 and 30,000 fails |
| Valuation formulas | "partly conjectural", tested to 100,000 | **proved** for every odd prime power |
| Odd prime powers | only 9 passes, tested to 100,000 | **proved**: 9 is the only odd prime power that passes |
| Speed | linear in n, with numbers of thousands of digits | for a prime p: only arithmetic mod p², about p²/4 multiplications |

## How we found it

For a prime p we can compute exactly how often p divides G(p): **always exactly twice** (v_p(G(p)) = 2, proof below). So p divides the numerator of T(p, k) unless p² divides G(k) · G(p − k). The test for a prime therefore only needs G(m) modulo p², never the full numbers with thousands of digits.

That made it easy to look at the remainders. Two things stood out:

- For every prime p ≡ 1 (mod 4), p divides G((p − 1)/2). That is a law (Legendre symbols, below), not chance.
- All other cases where p divides some G(m) look random: about one such m per prime, as a random model predicts.

In a random model a prime fails when two of these chance events line up. Below 3000 that should happen a few times, yet Rotondo saw none. So we scanned further, with two independent C programs (one working mod p² throughout, one mod p with a lift to p² only at zeros), and the first failure turned up at p = 3917: p² divides G(3417), and 3417 = 3917 − 500.

## The counterexample

```
p = 3917  (prime)       k = 500,   p − k = 3417

v_p(G(3917)) = 2    (G(3917) has 14,074 digits)
v_p(G(500))  = 0    (G(500)  has  1,350 digits)
v_p(G(3417)) = 2    (G(3417) has 12,075 digits)

v_p(T(3917, 500)) = 2 − 0 − 2 = 0,   numerator ≡ 2362 (mod 3917)
```

So 3917 does not divide the numerator of T(3917, 500), and the prime 3917 fails the test. This was rechecked with exact integer arithmetic (Python fractions on the full numbers), see the script at the end.

Scan results (all odd primes up to 30,000):

```
primes p < 30000 that fail:   3917 only
cases p | G(m), 1 ≤ m < p:    4956  (1611 of them the law for p ≡ 1 mod 4; the rest ≈ 1.03 per prime)
cases p² | G(m):              1     (p = 3917, m = 3417)
```

## The math

**1. A prime p passes iff p² never divides G(k) · G(p − k).** Pair j with p − j in G(p):

```
j^p + (p − j)^p ≡ p² · j^(p−1) ≡ p²   (mod p³)        (binomial theorem, p odd)
G(p) = Σ_{j<p} j^p + p^p ≡ (p − 1)/2 · p² ≡ −p²/2   (mod p³)
```

So v_p(G(p)) = 2 for every odd prime, and p passes iff v_p(G(k)) + v_p(G(p − k)) ≤ 1 for all k ≤ (p − 1)/2.

**2. The upper half reduces to harmonic sums.** For 2 ≤ k ≤ p − 1, using j^(p−k) ≡ j^(1−k) (mod p) and Σ_{j=1}^{p−1} j^(1−k) ≡ 0:

```
G(p − k) ≡ (−1)^k · (1 + 1/2^(k−1) + 1/3^(k−1) + … + 1/(k−1)^(k−1))   (mod p)
```

So "p | G(p − k)" means p divides the numerator of a fixed fraction, and "p | G(k)" means p divides the fixed integer G(k). These are unrelated numbers, so their residues behave like random ones.

**3. The one law: p ≡ 1 (mod 4).** With k = (p − 1)/2, j^k ≡ (j/p) (mod p) (Euler's criterion), so G(k) ≡ Σ_{j ≤ (p−1)/2} (j/p). When p ≡ 1 (mod 4), (j/p) = ((p − j)/p), so this half sum is half of the full sum Σ_{j<p}(j/p) = 0. Hence p | G((p − 1)/2) for every p ≡ 1 (mod 4). (For p ≡ 3 (mod 4) the half sum is (2 − (2/p))·h(−p) for p > 3, a number between 0 and p, so it never vanishes mod p.)

**4. How many primes should fail.** In a random model, each m < p has a chance 1/p that p | G(m), and given that, a chance 1/p that p² | G(m). A prime fails when two zeros pair up (k and p − k) or one zero is a double zero:

```
P(p fails) ≈ 1.5/p   for p ≡ 3 (mod 4)
P(p fails) ≈ 3.5/p   for p ≡ 1 (mod 4)    (the law of step 3 adds two chances)
expected number of failing primes ≤ x  ≈  2.5 · ln ln x + constant
```

The sum diverges, so infinitely many primes should fail, but only very slowly (like Wieferich primes): the model expects about eight more between 30,000 and 10¹⁰⁰. Given the zeros actually found, the model expects 3.5 failures below 30,000, and we see one. A result this low has a chance of about 14% under the model, so the data fit it.

**5. Rotondo's valuation formulas are true.** Write N = p^f · m with p ∤ m and f ≥ 1.

```
Lemma A.  If (p − 1) | m:      v_p(G(N)) = f − 1.
Lemma B.  v_p(G(p^e)) = 3e − 1   for e ≥ 1.
Lemma C.  If (p − 1) ∤ m:      v_p(G(N)) ≥ f + 1.
```

*A:* for p ∤ j, j^N ≡ 1 (mod p^(f+1)) because φ(p^(f+1)) divides N; terms with p | j vanish. So G(N) ≡ #{j ≤ N : p ∤ j} = p^(f−1)·m·(p − 1) (mod p^(f+1)).
*B:* pairing j with N − j as in step 1 gives G(N) ≡ (N²/2)·S (mod p^(3e)) with S = Σ_{j<N} j^(N−1). The terms of S coprime to p are j^(p^(e−1)−1) (mod p^e); these powers run (p − 1) times over the subgroup 1 + pℤ, whose elements sum to p^(e−1) (mod p^e). So v_p(S) = e − 1 and v_p(G(N)) = 2e + e − 1.
*C:* for p ∤ j, j^N ≡ ω(j)^m (mod p^(f+1)) with ω(j) the Teichmüller root of unity, and Σ ω(a)^m over a = 1..p−1 is 0 when (p − 1) ∤ m.

Rotondo's three formulas are B for e and for e − 1, and A with f = e − 1, m = p − 1 (giving e − 2).

**6. Odd prime powers: only 9 passes.**

```
e ≥ 3:  k = p^(e−1):  v_p(T) = (3e − 1) − (3e − 4) − (e − 2) = 5 − e < e      (B, B, A)  → fails
e = 2, p ≥ 5:  k = 2p:  v_p(T) ≤ 5 − 2 − 2 = 1 < 2                              (B, C, C)  → fails
p = 3, e = 2:  2p = 6 > 9/2 is outside the triangle; k = 1..4 all pass        → 9 passes
```

This is exactly Rotondo's "geometric accident", now proved for all prime powers.

**7. On the other ideas in Q2.** For a prime p the test costs about p²/4 multiplications mod p², which is exponential in the number of digits of p. Nothing here suggests a way to make it polynomial, and the counterexample shows it cannot be a primality test anyway. The composite half of the conjecture ("only 9 among odd composites") is still open beyond prime powers. For composites with two or more distinct prime factors, Lemma C plus the same random model gives a convergent sum, so we expect no other composite to pass.

## Validation

This script checks every claim: the definition on small n, the exact counterexample with full integers, that every odd prime below 3917 passes (mod p² arithmetic with numpy), the lemmas on sample cases, and the prime powers. It runs in about 15 seconds. The C scanner for larger ranges is in [code/conj-115.c](code/conj-115.c).

```python
import sys
import numpy as np
from fractions import Fraction
from sympy import isprime, primerange, multiplicity as v
sys.set_int_max_str_digits(0)

def G(m):                       # G(m) = 1^m + 2^m + ... + m^m, exact
    return sum(j**m for j in range(1, m + 1))

def passes(n):                  # Rotondo's test, exactly as defined on the page
    Gn = G(n)
    return all(Fraction(Gn, G(k) * G(n - k)).numerator % n == 0 for k in range(1, n // 2 + 1))

# 1. The definition: up to 300, the odd numbers that pass are the odd primes and 9.
assert [n for n in range(3, 301, 2) if passes(n)] == [n for n in range(3, 301, 2) if isprime(n) or n == 9]

# 2. The counterexample: 3917 is prime, and T(3917, 500) has a numerator that 3917 does not divide.
p, k = 3917, 500
assert isprime(p)
T = Fraction(G(p), G(k) * G(p - k))
assert T.numerator % p == 2362
assert (v(p, G(p)), v(p, G(k)), v(p, G(p - k))) == (2, 0, 2)

# 3. 3917 is the smallest odd prime that fails. For a prime p, v_p(G(p)) = 2, so p passes iff
#    v_p(G(k)) + v_p(G(p-k)) <= 1 for every k <= (p-1)/2; that only needs G(m) mod p^2.
def fails(p):
    p2, J = p * p, np.arange(1, p, dtype=np.int64)
    pw, val = np.ones(p - 1, dtype=np.int64), [0] * p
    for m in range(1, p):
        pw = pw * J % p2                      # pw[j-1] = j^m mod p^2
        s = int(pw[:m].sum()) % p2            # G(m) mod p^2
        val[m] = 2 if s == 0 else (1 if s % p == 0 else 0)
    return [k for k in range(1, (p - 1) // 2 + 1) if val[k] + val[p - k] >= 2]
for q in primerange(3, 3918):
    assert (fails(q) != []) == (q == 3917), q
assert fails(3917) == [500]

# 4. Lemmas behind the proofs, on sample cases.
assert all(v(q, G(q**e)) == 3*e - 1 for q, e in [(3,1),(3,2),(3,3),(5,1),(5,2),(7,1),(7,2),(11,2)])   # Lemma B
assert all(v(q, G(q**f * m)) == f - 1 for q, f, m in [(3,1,2),(3,2,2),(5,1,4),(5,2,4),(7,1,6),(11,1,10)])  # Lemma A
assert all(v(q, G(q**f * m)) >= f + 1 for q, f, m in [(5,1,2),(5,1,3),(7,1,2),(7,1,5),(5,2,3),(11,1,7)])  # Lemma C
for q in primerange(5, 120):                                # step 2: harmonic congruence for G(q-k)
    for kk in range(2, q):
        H = sum(Fraction(1, i**(kk - 1)) for i in range(1, kk))
        assert G(q - kk) % q == (-1)**kk * H.numerator * pow(H.denominator, -1, q) % q
assert all((G((q - 1)//2) % q == 0) == (q % 4 == 1) for q in primerange(5, 400))   # step 3

# 5. Prime powers: 9 passes, the next ones fail (step 6 proves it for all of them).
assert passes(9) and not any(passes(n) for n in [25, 27, 49, 81, 121, 125])
print("all claims verified: 3917 is the smallest odd prime that fails the triangle test")
```
