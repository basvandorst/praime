# Conjecture 19: A bound to the largest prime factor of certain Carmichael numbers

> **In short (attempt, not a record):** Problem 33, the follow-up page to this conjecture, already has far larger extreme triplets: Michael Bell's 3661-digit C3 (p = 5927222·1450# + 1, 2000) and Phil Carmody's C3s of up to 9534 digits (p = k·1996995²⁵⁰ + 1, 2004). The triplet below (1056 digits, all three primes proven) does **not** beat them. The rest of this note is kept for the method.
>
> For a Carmichael number pqr with p < q < r, Richard Pinch and Chris Nash proved on the page that r ≤ (p³ + p² − p + 1)/2. The bound is reached exactly when q = p² + p − 1 and r = (p³ + p² − p + 1)/2 are both prime, an "extreme" triplet. Rivera's record extreme triplet gives a 764-digit Carmichael number. We found **p = 25971283·2⁵⁶⁰ + 1**, for which p, q and r are all **proven** prime, giving an extreme Carmichael number of **1056 digits**. The key is that q − 1 = (p − 1)(p + 2) and r − 1 = (p − 1)(p + 1)²/2 both contain p − 1 = k·2ⁿ, so Proth, Pocklington and Brillhart–Lehmer–Selfridge prove all three.

**Puzzle:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_019.htm)

**Status:** open · not a record: problem 33 has extreme C3s up to 9534 digits (Carmody 2004); ours has 1056 digits, all three primes proven

## The question

A Carmichael number with three prime factors is C3 = pqr, p < q < r. Pinch and Nash showed (May 2000) that for given p, r ≤ (p³ + p² − p + 1)/2. Equality needs Nash's case I₂ = 2, that is q = p² + p − 1. Rivera then asked:

- a) are there infinitely many extreme triplets p, p² + p − 1, (p³ + p² − p + 1)/2, all prime?
- b) find larger and larger such triplets.

He started with p = 7352265·2⁴⁰⁰ + 1, which gives an extreme C3 of 764 digits. The search continued on [problem 33](https://www.primepuzzles.net/problems/prob_033.htm), where the records are far larger.

## Record vs ours

| | p | digits of C3 |
|---|---|---|
| Rivera (2000, this page) | 7352265·2⁴⁰⁰ + 1 | 764 |
| Bell (2000, problem 33) | 5927222·1450# + 1 | 3661 |
| Carmody (2004, problem 33) | 65506387867086·1996995²⁵⁰ + 1 | 9534 |
| Ours | 25971283·2⁵⁶⁰ + 1 | 1056 (all three proven) |

## How we found it

The chain of claims:

| # | Claim | Why | Validated |
|---|---|---|---|
| 1 | r is an integer for odd p | p³ + p² − p + 1 is even | algebra |
| 2 | q − 1 = (p − 1)(p + 2) and 2(r − 1) = (p − 1)(p + 1)² | expand | asserted in the script |
| 3 | with p = k·2ⁿ + 1, the parts k·2ⁿ of q − 1 and k·2ⁿ⁺¹ of r − 1 are fully factored | claim 2 | — |
| 4 | that factored part is ≈ √q and ≥ r^(1/3), so Pocklington (q) and BLS Theorem 5 (r) give proofs | p ≈ q^(1/2) ≈ (2r)^(1/3) | proofs run |
| 5 | k is bad mod a small prime ℓ when p ≡ 0, p is a root of x² + x − 1, or p is a root of x³ + x² − x + 1 (mod ℓ) | ℓ would divide p, q or 2r | sieve |
| 6 | pqr is a Carmichael number | Korselt: (p − 1), (q − 1), (r − 1) all divide pqr − 1 | asserted |

We fixed n = 560, sieved k < 3·10⁷ by every prime ℓ < 10⁵ (35552 survivors), and tested the survivors in order. The first triplet came at k = 25971283, about 10 seconds in. A Hardy–Littlewood estimate predicted roughly 3.7·10⁻⁵ hits per survivor, so about 1.3 hits in that range.

**The proofs.**
- p = 25971283·2⁵⁶⁰ + 1 with k < 2ⁿ: Proth's theorem, with a base a for which (a/p) = −1.
- q: q − 1 = k·2⁵⁶⁰·(p + 2), and p + 2 has the small factors 11·41·67·246527·747557. The factored part then exceeds √q, so Pocklington's theorem proves q.
- r: r − 1 = k·2⁵⁶¹·((p + 1)/2)². With k = 13·73·27367 the factored part F = k·2⁵⁶¹ satisfies r^(1/3) ≤ F < √r. Pocklington's conditions plus Brillhart–Lehmer–Selfridge (1975, Theorem 5) prove r: write r = c₂F² + c₁F + 1, then r is prime iff c₁² − 4c₂ is not a perfect square.

PARI/GP's isprime (APR-CL) proves all three independently.

## The math

Question a) is open. By the Bateman–Horn conjecture the three polynomials p, p² + p − 1 and (p³ + p² − p + 1)/2 take prime values together infinitely often, since they have no fixed prime divisor. The number of extreme triplets with p ≤ x should be about C·x/(ln x)³. The search above confirms the density: about one triplet per 27,000 sieved candidates at 176 digits.

## Validation

[`writeup/code/conj19_verify.py`](code/conj19_verify.py) (under a second) proves p, q and r with the theorems above and checks Korselt's criterion for pqr. With `--pari` it also proves all three with PARI's APR-CL, which takes longer.
