# Problem 69: A Brier number with disjoint covering sets

> **In short:** Yes, such a Brier number exists. Cover k·2ⁿ + 1 with the classic set A = {3, 5, 7, 13, 17, 241} (n mod 24). Cover k·2ⁿ − 1 with 60 other primes B (n mod 5040), using only orders of 2 that are at least 5. All conditions on k are then congruences modulo distinct primes, and the Chinese remainder theorem gives a 230-digit k. Every n is covered on both sides, so k·2ⁿ ± 1 is composite for every n ≥ 0, and A ∩ B = ∅.

**Puzzle:** [primepuzzles.net](https://www.primepuzzles.net/problems/prob_069.htm)

**Status:** answered · constructed and proved (a covering argument)

## The question

A Brier number k makes k·2ⁿ + 1 and k·2ⁿ − 1 composite for every n. Its covering sets A (for +1) and B (for −1) known so far always share a prime: Cohen and Selfridge's example shares 3. The problem asks for a Brier number with A ∩ B = ∅, or a proof that none exists.

On the page, Vantieghem found Sierpiński and Riesel multiples of 3 with one shared covering set, but no disjoint pair.

## The construction

**+1 side.** The primes 3, 5, 7, 13, 17, 241 have orders 2, 4, 3, 12, 8, 24 modulo 2. With residues n ≡ 0, 1, 0, 7, 3, 23 they cover every n mod 24. Prime p divides k·2ⁿ + 1 for n ≡ a (mod ord_p(2)) exactly when k ≡ −2^(−a) (mod p).

**−1 side.** Without those six primes, the moduli 2, 3, 4, 6, 8, 12 and 24 are unavailable: 2^m − 1 has no other primitive prime factor for those m. A greedy search over the divisors m ≥ 5 of 5040 finds a covering of ℤ/5040 by 60 primes with ord_q(2) = m, from 31 (m = 5) to 8369281 (m = 720). Prime q divides k·2ⁿ − 1 for n ≡ b exactly when k ≡ 2^(−b) (mod q).

**k.** All 66 primes are distinct, so the 66 congruences plus k odd have a common solution by the Chinese remainder theorem:

```
k = 42471920436627110676037647161862126761669187385276115106852869973057000545121168659751880739506377592186583970737835284711016161019613201498432465176564066940871081454518063549613802622632214564326803107663668779501442676850345767
```

For every n mod 24 a prime of A divides k·2ⁿ + 1, and for every n mod 5040 a prime of B divides k·2ⁿ − 1. Since k is larger than every prime in A ∪ B, each k·2ⁿ ± 1 is a proper multiple of that prime, so it is composite.

| set | covers | primes |
|---|---|---|
| A | k·2ⁿ + 1, n mod 24 | 3, 5, 7, 13, 17, 241 |
| B | k·2ⁿ − 1, n mod 5040 | 11, 19, 29, 31, 37, 41, 43, 61, 71, 73, 97, 109, 113, 127, 151, 181, 211, 257, 281, 331, 337, 421, 433, 577, 631, 673, 1009, 1321, 1429, 2017, 2521, 3361, 4481, 5153, 5419, 14449, 21169, 23311, 29191, 38737, 54001, 61681, 86171, 92737, 106681, 122921, 152041, 557761, 649657, 664441, 736961, 870031, 983431, 1130641, 1564921, 1765891, 2627857, 7416361, 8369281, 29247661 |

The k found is not claimed to be the smallest such number.

## Validation

[`writeup/code/prob69_brier.py`](code/prob69_brier.py) [`writeup/data/prob69_cover.json`](data/prob69_cover.json) rebuilds k from the covering. It checks that A and B are disjoint sets of primes, that the orders are right, and that every residue n mod 24 and every residue n mod 5040 is covered. [`writeup/code/prob69_cover.py`](code/prob69_cover.py) finds the two coverings.
