# Conjecture 72: Ghaffor's twin-count formula is the sieve of Eratosthenes on 6m ± 1 (proved for all n)

> **In short:** A formula G_n = 6n² + 10n + 4 − |U_n|, with U_n a union of 4n arithmetic progressions, was claimed to count the twin pairs up to 36n² + 60n + 21; it was verified to n = 10000 without explanation. (1) The progressions for k are exactly the m with 6m ± 1 a proper multiple of 6k ± 1. (2) Since (6n + 5)² > 36n² + 60n + 21, the moduli with k ≤ n catch every composite 6m ± 1 in range. So G_n is the twin count for every n; the limit question is the twin prime conjecture itself.

**Conjecture:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_072.htm)

**Status:** proved (Q1, for all n) · explained (Q2: the Eratosthenes mesh on 6m ± 1, identical to A002822) · Q3 is equivalent to the twin prime conjecture, open

## The question

Islem Ghaffor (2014) gave an algorithm: the number of twin pairs (p, p + 2) with p + 2 ≤ 36n² + 60n + 21 equals 6n² + 10n + 4 minus the number of distinct terms in 4n arithmetic progressions. His paper shows the formula and three worked examples (n = 1, 2, 3) and withholds the proof. Rivera verified it to n = 44, Kamenetsky to 2356, Schneider to 10000 (twin primes to 3.6·10⁹), Heleen to 50, Vantieghem to 100. Vantieghem noted the similarity to the theorem "6n − 1 and 6n + 1 both prime iff n is not 6xy ± x ± y" (A002822, Balestrieri 2011) but found that "Islem's numbers seem to have no direct relation to these forms".

- Q1: verify beyond n = 10000.
- Q2: what are the foundations; does the algorithm hide an Eratosthenes mesh?
- Q3: prove that G_n → ∞.

## What the page had vs ours

| | On the page | Ours |
|---|---|---|
| Q1 | verified to n = 10000, no proof | **proved for every n** (so no further verification is needed) |
| Q2 | "no idea how it works"; Perry's forms "seem to have no direct relation" | **yes**: the 4n progressions are the proper multiples of 6k ± 1, k ≤ n, inside the two progressions 6m ± 1; in m they are exactly 6xy ± x ± y |
| Q3 | five ideas, no result | equivalent to the twin prime conjecture; nothing in the formula helps |

## The algorithm, reconstructed

The PDF's notation is garbled in extraction, so the rule below is read off the three worked examples (n = 1, 2, 3) and matched to them exactly (union sizes 10, 29, 63; counts 10, 19, 25). For k = 1 … n and b = 0, 1, 2, … the four progressions are

```
5k − 3 + b(6k − 1)      7k − 3 + b(6k − 1)      5k − 1 + b(6k + 1)      7k − 1 + b(6k + 1)
```

with every term at most 6n² + 10n + 1 (his length l(v; r) = ⌈(6n² + 10n + 2 − v)/r⌉). U_n is the union of their terms and G_n = 6n² + 10n + 4 − |U_n|. For n = 1: {2, 7, 12, 17} ∪ {4, 9, 14} ∪ {4, 11} ∪ {6, 13} has 10 elements and G_1 = 20 − 10 = 10, the number of twin pairs up to 117.

## The math

**Step 1, the index.** Put m = v + 2. Then v = 1 … 6n² + 10n + 1 runs over m = 3 … 6n² + 10n + 3, and 6m + 1 ≤ 36n² + 60n + 19 < N = 36n² + 60n + 21 = (6n + 5)² − 4 exactly for these m (the next pair, m = 6n² + 10n + 4, has 6m + 1 = N + 4). Every twin pair above (3, 5) is (6m − 1, 6m + 1). The pairs with m ≤ 2 are (5, 7) and (11, 13); with (3, 5) that is three pairs, all twin, and 6n² + 10n + 4 = (6n² + 10n + 1) + 3. So G_n = 3 + #{3 ≤ m ≤ 6n² + 10n + 3 : v = m − 2 ∉ U_n}.

**Step 2, the progressions.** Let q = 6k − 1. Then 6k ≡ 1 (mod q), so 6m − 1 ≡ 0 ⟺ m ≡ k, and 6m + 1 ≡ 0 ⟺ m ≡ −k ≡ 5k − 1 (mod q). Let q = 6k + 1. Then 6k ≡ −1, so 6m + 1 ≡ 0 ⟺ m ≡ k, and 6m − 1 ≡ 0 ⟺ m ≡ −k ≡ 5k + 1 (mod q). In the index v = m − 2 the four classes are v ≡ k − 2, 5k − 3 (mod 6k − 1) and v ≡ k − 2, 5k − 1 (mod 6k + 1). Ghaffor's progressions start at 5k − 3 and 5k − 1 (b = 0: 6m + 1 = 5(6k − 1), resp. 6m − 1 = 5(6k + 1), proper multiples) and at 7k − 3 = (k − 2) + (6k − 1) and 7k − 1 = (k − 2) + (6k + 1), i.e. with b ≥ 1 in the class k − 2: that skips exactly m = k, where 6m ∓ 1 = q itself. So

```
U_n = { v : 6(v + 2) − 1 or 6(v + 2) + 1 is a proper multiple of some 6k ± 1, k ≤ n }.
```

**Step 3, completeness.** Every element of U_n is the index of a composite 6m ± 1, hence of a non-twin pair. Conversely, if 6m − 1 or 6m + 1 (≤ N) is composite, it has a prime factor p ≤ √N < 6n + 5; p ≠ 2, 3, and the numbers 6k ± 1 below 6n + 5 are those with k ≤ n; the quotient exceeds 1; so v ∈ U_n. The progressions of a composite modulus 6k ± 1 add nothing (they lie inside the progressions of its prime factors). Therefore

```
G_n = 3 + #{3 ≤ m ≤ 6n² + 10n + 3 : 6m − 1 and 6m + 1 prime} = π₂(36n² + 60n + 21)   for every n ≥ 1.
```

**Step 4, the forms.** In m, the four progressions read k + b(6k − 1) (b ≥ 1), 5k − 1 + b(6k − 1), 5k + 1 + b(6k + 1), k + b(6k + 1) (b ≥ 1). With x = k and y = b or b + 1 they are 6xy + x − y, 6xy − x − y, 6xy − x + y, 6xy + x + y: Perry's forms (A002822; Balestrieri 2011). Ghaffor's only addition is the bound x ≤ n, which Step 3 shows is enough. So Vantieghem's doubt is settled: the numbers are these forms, written per modulus.

**Q3.** G_n → ∞ says there are infinitely many twin primes. The formula is the sieve written out; it carries no information beyond the sieve, and the parity problem (FACTS 2017) keeps the sieve from deciding Q3.

## Validation

[code/conj72_check.py](code/conj72_check.py) rebuilds U_n from the rule above, reproduces the paper's three examples, compares G_n with the twin count from a sieve for every n ≤ 2000 (twin primes to 1.44·10⁸), checks U_n against direct factorisation of 6(v + 2) ± 1 for n ≤ 60, and checks that the progressions and Perry's forms give the same set below 9047. Run time 40 s.

```
Ghaffor's worked examples: |U_1|, |U_2|, |U_3| = [10, 29, 63] (paper: 10, 29, 63)
G_1, G_2, G_3 = [10, 19, 25] (paper: 10, 19, 25)
G_n equals the twin count up to 36n^2+60n+21 for all n <= 2000: True
U_n (index v) = {v : 6(v+2)-1 or 6(v+2)+1 composite} for all n <= 60: True
Ghaffor's terms (as m = v + 2) and Perry's 6xy ± x ± y agree below 9047: True
```
