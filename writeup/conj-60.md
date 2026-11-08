# Conjecture 60: A prime between n^k and (n+1)^k

> **In short:** The conjecture states that for n ≥ 1 and k ≥ K = log 127/log 16 = 1.74717… there is always a prime strictly between n^k and (n+1)^k; it generalises Legendre's conjecture (k = 2). (1) At k = K exactly the statement fails: 16^K = 127 and 15^K = 113.45, and no prime lies strictly between them. So the right form is k > K, as the author also noted on the page. (2) A failure at (n, k) needs consecutive primes p < p′ with p ≤ n^k and (n+1)^k ≤ p′. Then p′ − p ≥ k·n^(k−1) ≥ k·p^(1−1/k) ≥ K·p^(1−1/K), so a gap of at least 1.747·p^0.4276. (3) Below 4·10¹⁸ no gap exceeds 1476, which leaves only p < 6.98·10⁶. Among those primes only the gap 113 → 127 is that large, and its failing exponents end exactly at K. So the conjecture holds for every k > K and every n with n^k ≤ 4·10¹⁸. Beyond, it needs gaps below 1.747·p^0.4276, which is stronger than Legendre's conjecture and open.

**Puzzle:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_060.htm)

**Status:** proved for n^k ≤ 4·10¹⁸ (k > K; false at k = K itself); open in general (needs gaps below K·p^(1−1/K))

## The question

T. D. Noe: for n > 0 and k ≥ K = log(127)/log(16) = 1.74717117…, there is always a prime between n^k and (n+1)^k. On the page:

- the proposal K = √3, which Noe refutes (the statement is false for k between log 113/log 15 and log 127/log 16);
- the remark that Firoozbakht's conjecture implies it;
- Noe's reply that with k > K the gap 113 → 127 "appears" to be the hardest case.

There is no proof on the page.

## Record vs ours

| | Status |
|---|---|
| Page | checked numerically; implied by Firoozbakht's conjecture |
| **Ours** | **proved for all k > K and n^k ≤ 4·10¹⁸**; false at k = K (n = 15); open beyond |

## The math

**Step 1 (proved): what a failure is.** For n ≥ 2 and k > 1 let p be the largest prime ≤ n^k and p′ the next prime. The open interval (n^k, (n+1)^k) has no prime iff (n+1)^k ≤ p′, that is iff

```
ln p / ln n ≤ k ≤ ln p′ / ln(n + 1).
```

For n = 1 the interval (1, 2^k) contains 2.

**Step 2 (proved): a failure needs a large gap.** By the mean value theorem, (n+1)^k − n^k ≥ k·n^(k−1) = k·(n^k)^(1−1/k) ≥ k·p^(1−1/k). The function k ↦ k·p^(1−1/k) is increasing, because its logarithmic derivative is 1/k + ln p/k² > 0. So for k ≥ K a failure needs p′ − p ≥ K·p^(1−1/K) = 1.74717·p^0.42765.

**Step 3 (computed): the only candidate gap.** Every gap below 4·10¹⁸ is at most 1476 (complete computation of Oliveira e Silva, Herzog and Pardi; FACTS 3023). So Step 2 needs 1476 ≥ K·p^(1−1/K), that is p < (1476/K)^(K/(K−1)) = 6.98·10⁶. Among all pairs of consecutive primes below that, only 113 → 127 has a gap (14) of at least K·p^(1−1/K) (13.2). Its failing exponents are [ln 113/ln 15, ln 127/ln 16] = [1.745680, 1.747171] for n = 15, and two shorter intervals for other n end lower. So the largest failing exponent is exactly K: the statement fails at k = K and holds for every k > K with n^k ≤ 4·10¹⁸.

**Beyond 4·10¹⁸ (open).** The argument needs gaps below 1.747·p^0.428. Legendre's conjecture (k = 2) needs gaps below about 2√p, so the statement with K < k < 2 is stronger than Legendre's, and Firoozbakht's conjecture would imply it. Under Cramér's model a gap of length K·p^0.428 at p ≥ 4·10¹⁸ has probability about exp(−K·p^0.428/ln p) < exp(−10⁶), so no failure is expected.

## Validation

[`code/conj60_verify.py`](code/conj60_verify.py) (about 1 s) sieves the primes below 7·10⁶ and lists every pair with gap ≥ K·p^(1−1/K). For each such pair it enumerates every n with a nonempty interval of failing k. It checks with 50-digit arithmetic that the highest one ends exactly at K (n = 15, 113 → 127) and that all others end below K.

```
K = 1.747171171693; failures need a gap >= K p^(1-1/K); with gaps <= 1476 this needs p < 6.9801e+06
1 prime pairs below 6.980e+06 have a gap that large; the largest p among them is 113
3 failing intervals with k >= 1 from these pairs; the highest ends at k = 1.747171171693 (n = 15, gap 113 -> 127)
no failure with k > K and n^k <= 4e18; at k = K exactly the statement fails at n = 15
all claims verified
```
