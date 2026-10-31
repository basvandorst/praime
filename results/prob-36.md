# Problem 36 (Q4 and Liskovets' bias question): Liskovets–Gallot numbers

> **In short:**
> - **Q4: yes.** There are k with 3 | k such that k·2ⁿ + 1 and k·2ⁿ − 1 are both composite for every even n, and also k for every odd n. Write n = 2m, so the numbers are k·4^m ± 1, or 2k·4^m ± 1 for odd n. That makes it a Brier problem in base 4. Cover the + side with {5, 7, 17, 13, 241} (m mod 12) and the − side with 28 other primes (m mod 720). The Chinese remainder theorem then gives a 93-digit k in each case.
> - **The bias.** Primes with orders 3, 4, 8, 12 and 24 of 2, namely 7, 5, 17, 13 and 241, act together on n mod 24. For k such as 87 or 93 they leave only about 10 % of the odd n free of small factors, against about 80 % of the even n, a near-covering of one parity class. Computing this survival jointly mod 55440 predicts the odd share of primes for all eight k on the page, close to what n ≤ 3000 shows.

**Puzzle:** [primepuzzles.net](https://www.primepuzzles.net/problems/prob_036.htm)

**Status:** answered · Q4 constructed and proved (covering argument); bias explained quantitatively

## The questions

Liskovets noticed that for some k divisible by 3, the primes k·2ⁿ + 1 occur almost only for odd n (51, 231, 261) or almost only for even n (87, 93, 177, 213, 267). Gallot found k for which one parity class is fully covered. Two questions are open:

- **Q4.** Are there k with both k·2ⁿ ± 1 composite for all even n (or all odd n)? If not, why not?
- **Liskovets.** Can anybody explain these distribution biases quantitatively?

## Q4

**Even n.** Put n = 2m. Each covering prime p divides k·4^m + 1 for m in one class modulo the order of 4 mod p, and similarly for k·4^m − 1.

- **+ side:** 5, 7, 17, 13 and 241 have orders 2, 3, 4, 6 and 12 of 4. With the classes m ≡ 0, 0, 1, 1, 11 they cover every m mod 12.
- **− side:** a greedy search covers ℤ/720 with 28 other primes.
- **k:** the conditions k ≡ 0 (mod 3), k odd, k ≡ −4^(−a) (mod p) and k ≡ 4^(−b) (mod q) are congruences modulo distinct primes. The Chinese remainder theorem gives

```
k = 188217134092166051946880247390376872865453249058827201192162143997299056928084310516854597049
```

**Odd n.** The same with 2k in place of k (n = 2m + 1) gives

```
k = 181872797385868280825936047503437703923497956594762799373133260571546640129450351204820852217
```

Both are divisible by 3 and larger than every prime used, so every number in question is a proper multiple of a covering prime. We do not claim these are the smallest such k.

## The bias

Treat the primes as independent and the prediction comes out right in direction but too weak: for k = 87 it gives an odd share of 0.25. The primes with orders dividing 24 work together on n mod 24, and for these k they nearly cover one parity class. So the survival has to be computed jointly modulo L₀ for all primes whose order divides L₀, with the remaining primes treated as independent.

| k | odd n surviving mod 55440 | even n surviving | predicted odd share | observed, n ≤ 3000 |
|---|---|---|---|---|
| 51 | 0.780 | 0.129 | 0.83 | 24/26 = 0.92 |
| 231 | 0.823 | 0.254 | 0.76 | 35/41 = 0.85 |
| 261 | 0.876 | 0.307 | 0.73 | 35/42 = 0.83 |
| 87 | 0.077 | 0.417 | 0.14 | 1/17 = 0.06 |
| 93 | 0.102 | 0.787 | 0.12 | 0/27 = 0.00 |
| 177 | 0.191 | 0.693 | 0.21 | 6/38 = 0.16 |
| 213 | 0.057 | 0.455 | 0.12 | 1/12 = 0.08 |
| 267 | 0.112 | 0.713 | 0.13 | 4/26 = 0.15 |

The prediction gets the direction right in every case and the size roughly right. The counts are small, so the differences are within Poisson noise; for k = 93 the expected number of odd n is 3.3, and observing none has probability about 4 %. In six of the eight cases the observed bias is a little stronger than predicted, so the model may still miss some correlation.

The Liskovets–Gallot numbers are the limit of this effect: they are the k for which these small primes cover a parity class completely.

## Verify it yourself

| script | what it does |
|---|---|
| `results/code/prob36_q4.py` | builds both coverings and both k, and checks every class of m mod 12 and mod 720 |
| `results/code/prob36_bias.py 55440` | computes the joint survival and the predictions |
