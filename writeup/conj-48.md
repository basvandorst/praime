# Conjecture 48: P = A ± B

> **In short:** Daniele Degiorgi asked on the page whether one can prove that no exponents give |2ⁿ − 3^m| = 103. Yes, by congruences: 2ⁿ − 3^m = 103 forces n and m even (mod 3 and mod 4), and then 103 would factor as a difference of squares (2^(n/2) = 52 is impossible). And 3^m − 2ⁿ = 103 fails mod 8. So 103 is the least prime not representable with only the primes 2 and 3, as Degiorgi suspected. For Werner Sand's corrected conjecture (A·B made of exactly the primes up to Q, nothing else), 397 and 419 do have representations, with Q = 23 and Q = 29. But **433 and 26 more primes below 1000 have none for any admissible Q ≤ 31 with A ≤ 10¹⁶**, so the corrected conjecture is most likely false. A disproof would also need to rule out every larger Q.

**Puzzle:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_048.htm)

**Status:** answered (Degiorgi's question, proof) · corrected conjecture: 433 is the least prime without a representation for Q ≤ 31, A ≤ 10¹⁶ (exhaustive in that range); larger Q not excluded

## The question

Werner Sand: let A, B be coprime positive integers with A·B containing every prime from 2 to Q. Then P = A + B or |A − B| is prime when 1 < P < R², where R is the prime after Q. He conjectured that every prime P > 3 arises in this way.

On the page:

- Luke Pebody proved the conjecture as stated: take B = 2·3⋯Q (Q the prime before P) and A = P + B.
- Sand then explained that he meant "and no factor else": A·B must have exactly the primes up to Q.
- Degiorgi took R = 5 (only 2 and 3) and wrote: "The first candidate is 103 but I am not able to prove that no m,n exists such that abs(2^n − 3^m) = 103 … Does such a proof exist?"
- Sand noted 103 = 3·5·7 − 2 with Q = 7, and that the restricted question is related to Catalan's equation.

## Record vs ours

| | Status |
|---|---|
| Page | stated conjecture proved (Pebody); corrected version open; |2ⁿ − 3^m| = 103 open (Degiorgi) |
| **Ours** | **|2ⁿ − 3^m| ≠ 103 (proof); corrected version: 433 has no representation for Q ≤ 31, A ≤ 10¹⁶** |

## The math

**Degiorgi's question.** All four forms Degiorgi lists fail for 103:
- **2ⁿ + 3^m = 103:** 103 − 3^m for m ≤ 4 is 102, 100, 94, 76, 22, none a power of 2.
- **2ⁿ3^m ± 1 = 103:** 102 = 2·3·17 and 104 = 2³·13, so neither works.
- **2ⁿ − 3^m = 103:** mod 3, 2ⁿ ≡ 1, so n is even. Mod 4 (n ≥ 2), −3^m ≡ 3, so 3^m ≡ 1 and m is even. Then (2^(n/2) − 3^(m/2))(2^(n/2) + 3^(m/2)) = 103 is prime, so the first factor is 1 and 2^(n/2) = (1 + 103)/2 = 52. Not a power of 2.
- **3^m − 2ⁿ = 103:** mod 3, 2ⁿ ≡ 2, so n is odd. For n ≥ 3, mod 8 gives 3^m ≡ 103 ≡ 7, but 3^m mod 8 is 1 or 3. For n = 1, 3^m = 105, not a power of 3.

**Sand's corrected conjecture.** For a prime P, an admissible Q is any prime with nextprime(Q)² > P, that is Q ≥ the largest prime ≤ √P (Sand's own example: 103 = 3·5·7 − 2 with Q = 7). A representation is P = A + B or |A − B| with gcd(A, B) = 1 and rad(AB) = 2·3·5⋯Q. Any such value is automatically prime when P < R², since no prime ≤ Q divides it.

| # | Claim | How | Range |
|---|---|---|---|
| 1 | every prime below 397 has a representation (with its smallest admissible Q) | search over Q-smooth A | A ≤ 10⁹ |
| 2 | 397 = 8008462 − 8008065 (Q = 23) and 419 = 12128480 − 12128061 (Q = 29) | found | — |
| 3 | 433 has no representation for Q = 19, 23, 29, 31 | all Q-smooth A ≤ 10¹⁶, checking A + 433, A − 433 and 433 − A | exhaustive in that range |
| 4 | the same holds for 26 more primes below 1000: 541, 547, 563, 599, 631, 647, 659, 733, 743, 761, 769, 797, 823, 827, 839, 853, 859, 863, 877, 881, 941, 967, 977, 983, 991, 997 | same search | Q ≤ 31, A ≤ 10¹⁶ |

For fixed P and Q, the equation is a two-term S-unit equation, which has only finitely many solutions (Evertse). A representation of 433 beyond our range would need two Q-smooth numbers above 10¹⁶ that differ by exactly 433. Below 10¹⁶ there are only about 10⁷ 19-smooth numbers, so such coincidences are extremely rare. Larger Q make A·B ≥ Q# grow quickly while the number of usable smooth pairs does not keep up. So we expect 433 to be a genuine counterexample, but the claim is limited to the ranges above.

## Validation

[`writeup/code/conj48_verify.py`](code/conj48_verify.py) (about 20 seconds) checks every step of the proof for 103 and confirms it by brute force for small exponents. It verifies the representations of 397 and 419, and reruns the search for 433 with Q ∈ {19, 23, 29, 31}: up to A ≤ 10¹³ by default, and up to 10¹⁶ with `--full`.
