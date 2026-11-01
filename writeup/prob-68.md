# Problem 68: Brier numbers whose two covering sets share 4 or more primes

> **In short:** The page asked for the smallest product P of the primes in the covering sets A (for k·2ⁿ + 1) and B (for k·2ⁿ − 1) of a Brier number when A and B share s primes; values were known for s = 1, 2, 3 and none for s = 4. The proof: (1) a shared prime has even order m and covers n ≡ a on one side and n ≡ a + m/2 on the other; with 3 shared, the others need m ≡ 2 (mod 4), which for periods dividing 720 leaves only 11, 19, 331 and 18837001; (2) splitting n by parity turns both sides into covering problems linked only through these primes; (3) an exact integer program over all residue classes finds the cheapest coverings. This gives the first Brier numbers with s = 4 (P of 37 digits) up to s = 8, and proves the page's P for s = 1 and 2 minimal for periods dividing 720 with 3 shared.

**Puzzle:** [primepuzzles.net](https://www.primepuzzles.net/problems/prob_068.htm)

**Status:** answered (s = 4 to 8: first examples; s = 1 and 2 proved minimal for periods dividing 720) · proof (covering checked for every n mod the period)

## The question

A Brier number k makes k·2ⁿ + 1 and k·2ⁿ − 1 composite for every n. Let A be a covering set for k·2ⁿ + 1 and B one for k·2ⁿ − 1, and P the product of the primes in A ∪ B.

**Q.** Find the smallest P for which the number of primes in both A and B is 1, 2, 3, 4, etc.

## What the page had vs ours

| shared primes s | on the page | ours |
|---|---|---|
| 1 | P = 33241542480794255062005345795 (29 digits) | **the smallest possible** when 3 is shared and the period divides 720 |
| 2 | P = 14393587894183912441848314729235 (32 digits) | **the smallest possible**, same conditions |
| 3 | P = 2708540859502146859799487840589422795 (37 digits) | the best the integer program finds is this same P (not proved minimal) |
| 4 | ? | **P = 8515242199115959558869800655081023415 (37 digits)**, A ∩ B = {3, 11, 19, 331} |
| 5 | — | **P = 71515612535007008102463292424534086157745 (41 digits)**, A ∩ B = {3, 11, 19, 331, 18837001}; s = 5 is the most possible when 3 is shared and the period divides 720 |
| 6 | — | **P of 69 digits**, A ∩ B = {3, 11, 19, 43, 281, 331} |
| 7 | — | **P of 84 digits**, A ∩ B = {3, 11, 19, 43, 281, 331, 5419} |
| 8 | — | **P of 60 digits**, A ∩ B = {3, 11, 43, 211, 281, 331, 5419, 86171} |

The s = 4 example:

```
k = 12015917365937901464995049971524255311
A = {3, 7, 11, 13, 19, 37, 61, 73, 109, 151, 331, 1321}       covers k·2ⁿ + 1, period 360
B = {3, 5, 11, 17, 19, 31, 41, 181, 241, 331, 433, 61681}     covers k·2ⁿ − 1, period 360
```

## How we found it

Every prime of A ∪ B fixes k modulo itself, and the two covering conditions only see k mod p. So the search is over residue patterns, not over k. Two facts make it small. Only primes of order m ≡ 2 (mod 4) can be shared next to 3. And after splitting n by parity, both sides become ordinary covering problems on the same residues t, linked only through the shared primes. An integer program over all residue classes mod 720 then finds the cheapest pattern exactly; the Chinese remainder theorem gives k.

## The math

Notation: m = ord_p(2), the order of 2 mod p.

**Step 1 (proved). Which n a prime covers.** If k ≡ −2⁻ᵃ (mod p), then p | k·2ⁿ + 1 exactly when n ≡ a (mod m). If k ≡ 2⁻ᵇ (mod p), then p | k·2ⁿ − 1 exactly when n ≡ b (mod m). So p can cover the +1 side only when k lies in −⟨2⟩ mod p, and the −1 side only when k lies in ⟨2⟩.

**Step 2 (proved). Shared primes.** If p ∈ A ∩ B, then −2⁻ᵃ ≡ 2⁻ᵇ, so 2^(b−a) ≡ −1 (mod p). This needs m even, and then b ≡ a + m/2 (mod m). Conversely, every prime of even order with p ∤ k covers one class on each side. A covering set is minimal exactly when each member covers some n that no other member covers (a private n). The A and B on the page are minimal in this sense; [code/prob68_verify.py](code/prob68_verify.py) checks it.

**Step 3 (proved). With 3 shared, the other shared primes have m ≡ 2 (mod 4).** For p = 3 (m = 2), step 2 gives that 3 covers one parity of n on the +1 side and the other parity on the −1 side. A shared prime p needs a private n on the +1 side, so its class a has the parity that 3 leaves open there; likewise b has the parity 3 leaves open on the −1 side. These parities differ, so m/2 ≡ b − a is odd. The primes with m ≡ 2 (mod 4) and m | 720 are 3, 11, 19, 331 and 18837001 (m = 2, 10, 18, 30, 90; 2⁶ − 1 has no new prime). So with 3 shared and period dividing 720, s ≤ 5, and s = 4 forces {3, 11, 19, 331} unless 18837001 is used.

**Step 4 (proved). The reduced problem.** Take k ≡ 2 (mod 3), so 3 covers even n on the +1 side and odd n on the −1 side; the other choice is the shift n → n + 1. Then A must cover n = 2t + 1 and B must cover n = 2t. By step 1:

- m odd: one class t ≡ c (mod m), on one side of our choice;
- m ≡ 0 (mod 4): one class mod m/2 on one side (its class on the other side has the parity 3 already covers);
- m ≡ 2 (mod 4), h = m/2 odd: the class t ≡ c (mod h) on the +1 side **and** t ≡ c + ½ (mod h) on the −1 side, where ½ = (h + 1)/2. Indeed n = 2c + 1 on the +1 side moves to n + h = 2(c + (h+1)/2) on the −1 side.

So 11, 19 and 331 cover t ≡ c (mod 5, 9, 15) on one side and the same classes shifted by 3, 5 and 8 on the other.

**Step 5 (exact search). The integer program.** For t mod T, one 0/1 variable per prime, class and side, with at most one class per prime; every t on both sides must be covered; a shared prime must have a private t on both sides; minimise Σ log p. T = 360 covers every Brier number with 3 shared and period dividing 720, with the complete list of the 39 primes of order m | 720 (two of them, of 18 and 51 digits, proved by APRCL). Solved with HiGHS:

| s at least | T | result |
|---|---|---|
| 1 | 360 | optimum P = 33241542480794255062005345795, the page's set: minimal |
| 2 | 360 | optimum P = 14393587894183912441848314729235, the page's set: minimal |
| 3 | 360 | best found in 25 minutes: the page's P (not proved minimal) |
| 4 | 360 | best found in 25 minutes: P = 8515242199115959558869800655081023415 (not proved minimal; the solver's lower bound reached 33.8 digits) |
| 5 | 360 | best found in 25 minutes: P = 71515612535007008102463292424534086157745 (not proved minimal; lower bound 32.6 digits) |

With T = 1260 the factor 7 brings in more primes of order ≡ 2 (mod 4): 43 (m = 14), 5419 (42), 281 and 86171 (70), 211, 664441 and 1564921 (210), 1765891 (630). The program at T = 1260 is too large to solve to the end (primes below 10⁷ only), but its best solutions give the examples for s = 6, 7 and 8. They are far from optimal: the s = 8 example is cheaper than the s = 6 one.

**Step 6 (checked). Composite, not just divisible.** Each k is odd and larger than every prime it uses, so each k·2ⁿ ± 1 is a proper multiple of a prime of A or B, hence composite.

**What this does not settle.** Whether 37 digits is the least P for s = 4 (the solver's lower bound reached 33.8 digits for periods dividing 720), the least P for s ≥ 5, and Brier numbers in which 3 is not shared.

## Validation

[code/prob68_verify.py](code/prob68_verify.py) checks from scratch the page's three records and our examples. For every n mod the period, a prime of A divides k·2ⁿ + 1 and a prime of B divides k·2ⁿ − 1. It also checks that A and B are minimal, that the shared count and P are as claimed, steps 2 and 3, and the completeness of the prime list [data/prob68_primes.json](data/prob68_primes.json) for all orders m | 720. [code/prob68_ilp.py](code/prob68_ilp.py) is the integer program of step 5. Run after `pip install sympy numpy scipy` (cypari2 for the APRCL check):

```
python3 prob68_verify.py                 # about 1 second
python3 prob68_ilp.py 360 1              # s >= 1: optimum, a few minutes
python3 prob68_ilp.py 360 4 1500         # s >= 4: the 37-digit example within 25 minutes
```
