# Conjecture 28: Buss's conjecture on nextprime(n!^m + 1) − n!^m

> **In short:** Q(n, m) = nextprime(n!^m + 1) − n!^m is prime for all n when m = 1 … 5 (verified to n = 3734, 2274, 1733, 1439, 1240), while composite values exist for m ≥ 6; the page asks for a counterexample. The analysis: (1) every prime up to n divides n!^m, so Q is the first k > n without a prime factor up to n for which n!^m + k is prime, and such k is prime unless it passes the square of the next prime after n; (2) each candidate succeeds with probability F_n / ln(n!^m), F_n = ∏_{p ≤ n} (1 − 1/p)⁻¹, so the chance that Q is composite is about exp(−e^γ n / (2m ln n)) times the share of composites among the candidates; (3) summed over n this gives 0.00, 0.03, 0.16, 0.38, 0.66 expected counterexamples for m = 1 … 5 and 1.0, 1.3, 1.7, 2.2, 2.6 for m = 6 … 10, against 0, 0, 0, 0, 0 and 1, 2, 3, 4, 3 found for n ≤ 40. So no counterexample is expected for m ≤ 5 (probability about 10⁻⁴⁰ beyond the verified range), and the conjecture holds by a 29 % chance at small n, not by structure.

**Conjecture:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_028.htm)

**Status:** explained (heuristic, calibrated on m = 6 … 10): no counterexample is expected for m ≤ 5; for every m ≥ 6 the model predicts the composite values that are found

## The question

Let P = nextprime(n!^m + 1) and Q(n, m) = P − n!^m (OEIS A037153, A067362–A067365 for m = 1 … 5). Buss's conjecture: Q(n, m) is prime for all n when m = 1 … 5. Known composite values: Q(28, 6), Q(6, 7), Q(5, 8).

**Q.** Can you find a counterexample?

## What the page had vs ours

| | On the page | Ours |
|---|---|---|
| m = 1 … 5 | all prime to n = 3734, 2274, 1733, 1439, 1240 (Jacobsen, 2014) | the expected number of composite values is 0.00, 0.03, 0.16, 0.38, 0.66, almost all of it at n < 60; beyond the verified ranges less than 10⁻³⁰ |
| m ≥ 6 | three composite values known | all composite values for n ≤ 40, m ≤ 10 (13 of them), and the model's expected counts 1.0 … 2.6 that match them |

## How we found it

Q(n, m) is a Fortune-type number: the gap from n!^m to the next prime. The key is that the gap can only land on numbers coprime to n!, and those are primes until the square of the next prime p′ after n. So Q is composite only if the search runs past p′² ≈ n², while the expected length of the search is about m·n. The probability of such an overrun is exponentially small in n/(m ln n), which is why counterexamples appear only at small n and only for larger m.

## The math

**Step 1 (proved). The candidates.** If p ≤ n then p | n!^m, so p ∤ n!^m + k iff p ∤ k. Hence Q(n, m) is the least k ≥ 2 with gcd(k, n!) = 1 and n!^m + k prime. The admissible k are the integers above n with no prime factor ≤ n: in increasing order these are the primes in (n, p′²) and then also composites, the first being p′² (p′ the next prime after n). In particular Q(n, m) < p′² implies Q prime.

**Step 2 (heuristic). One candidate.** An admissible k already avoids the primes up to n, so n!^m + k is prime with probability about q = F_n / ln(n!^m) with F_n = ∏_{p ≤ n} (1 − 1/p)⁻¹ ≈ e^γ ln n (Mertens): q ≈ e^γ ln n / (m n ln n) = e^γ/(m n). The candidates are nearly independent.

**Step 3 (heuristic). Q composite.** With candidates c₁ < c₂ < … and r(c) the number of candidates below c, P(Q = c) = q(1 − q)^{r(c)}, so E(n, m) = Σ_{c composite} q(1 − q)^{r(c)}. The first composite candidate is p′², preceded by π(p′²) − π(n) ≈ n²/(2 ln n) primes, so E(n, m) ≈ exp(−e^γ n/(2m ln n)) × (share of composites among the next candidates). For n = 28, m = 6: exp(−1.78·28/(12·3.33)) = 0.29 before the composite share; the exact sum gives E(28, 6) = 0.016.

**Step 4 (computed). The totals.** Summing E(n, m) over 4 ≤ n ≤ 200 (the terms decrease like exp(−n/m) beyond; for n ≤ 3 the model's q exceeds 1 and those values are checked exactly):

| m | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| expected composite Q | 0.000 | 0.031 | 0.161 | 0.377 | 0.656 | 0.98 | 1.35 | 1.74 | 2.16 | 2.60 |
| found, n ≤ 40 | 0 | 0 | 0 | 0 | 0 | 1 | 2 | 3 | 4 | 3 |

The composite values found: Q(28, 6) = 29·41; Q(6, 7) = 7·11, Q(10, 7) = 11·19; Q(5, 8) = 7·11, Q(10, 8) = 11², Q(17, 8) = 23·59; Q(2, 9) = 3², Q(3, 9) = 5², Q(8, 9) = 11·13, Q(24, 9) = 31²; Q(8, 10) = 13·17, Q(9, 10) = 11·19, Q(28, 10) = 29·37. They are products of two primes just above n, as step 1 says they must be.

**Step 5 (conclusion).** For m ≤ 5 the expected number of counterexamples is 1.22 in total, so observing none has probability e^(−1.22) = 29 %: the conjecture is true by chance at small n. Beyond the verified ranges (n > 1240 for m = 5) the remaining expectation is below 10⁻³⁰ (E(200, 5) = 7·10⁻⁷ and falling exponentially), so no counterexample will be found; but nothing structural prevents one, as the m ≥ 6 cases show.

## Validation

[code/conj28_model.py](code/conj28_model.py) computes Q(n, m) exactly for n ≤ 40, m ≤ 10 and lists the composite values, then evaluates the model E(n, m) by an exact sum over sieved candidates and prints the expected totals and the 29 %. Run after `pip install sympy gmpy2 numpy`:

```
python3 conj28_model.py      # a few seconds
```
