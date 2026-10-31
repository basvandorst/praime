# Problem 92: Multiple covering sets for Sierpinski and Riesel integers

> **In short:** Emilia Gurisatti found a Sierpinski number and a Riesel number with 4 distinct covering sets each, and asked for more (Q1). We give **a Sierpinski number and a Riesel number with 8 distinct minimal covering sets each** (253 digits), and a much smaller **27-digit Sierpinski number with 5**: k = 660664662814032650416080443. The idea is to keep the core {3, 5, 7, 13} that all of Gurisatti's sets share, and add tails built from new primes that each cover what the core leaves open. Every prime imposes one congruence on k, and disjoint tails never conflict, so the Chinese remainder theorem glues any number of them into one k. The 27-digit number is also an explicit upper bound for Rivera's Q4 on [puzzle 1279](https://www.primepuzzles.net/puzzles/puzz_1279.htm) (the least S or R integer with five covering sets).

**Puzzle:** [primepuzzles.net](https://www.primepuzzles.net/problems/prob_092.htm)

**Status:** answered (Q1) · 8 covering sets for a Sierpinski and a Riesel number, 5 for a 27-digit Sierpinski number; all verified from k alone

## The question

A Sierpinski number k makes k·2ⁿ + 1 composite for every n ≥ 1, a Riesel number makes k·2ⁿ − 1 composite. A covering set is a finite set of primes such that every n has one of them dividing k·2ⁿ ± 1. Emilia Gurisatti found numbers with 3 distinct covering sets (OEIS A394518, A394519), and then:

- Sierpinski 218164505156603 with four: {3, 5, 7, 13, 19, 37, 73}, {3, 5, 7, 13, 17, 241}, {3, 5, 7, 13, 19, 37, 109}, {3, 5, 7, 13, 97, 241, 257};
- Riesel 16107442173648029 with the same four sets.

Q1: "Can you find Sierpinski or Riesel integers with more than 4 covering sets?" Arkadiusz Wesołowski added a construction with two disjoint covering sets (153 digits), and Gurisatti found Brier numbers with two (Q2). Puzzle 1279 (August 2026) asks for the least S or R integer with five.

## Record vs ours

| | Sierpinski | Riesel |
|---|---|---|
| Gurisatti (2026) | 4 sets (15 digits) | 4 sets (17 digits) |
| **Ours** | **8 sets** (253 digits); **5 sets with 27 digits** | **8 sets** (253 digits) |

## The math

| # | Claim | Why | Validated |
|---|---|---|---|
| 1 | for fixed k, a prime p divides k·2ⁿ + s exactly for n ≡ r (mod ord_p 2), where 2^r ≡ −s/k (mod p), or for no n at all | 2ⁿ runs through ⟨2⟩ mod p | recomputed from k for every prime used |
| 2 | Gurisatti's core {3, 5, 7, 13} leaves exactly 1/12 of the classes n mod 720720 uncovered | classes 0 (mod 2), 3 (mod 4), 2 (mod 3), 9 (mod 12) for her Sierpinski number | computed |
| 3 | any set of primes, each with a chosen class, that covers those leftover classes is a "tail"; core + tail is a covering set | — | each set checked to cover all n mod 720720 |
| 4 | tails built from pairwise disjoint primes are compatible: each prime fixes k mod p, so CRT gives one k for all of them | Chinese remainder theorem | the k's below |
| 5 | each core + tail is minimal: dropping any prime uncovers some class | — | checked |
| 6 | among primes below 2·10⁷ whose order of 2 divides 720720, 8 disjoint tails fit (greedy) | the remaining "capacity" Σ 1/ord ≈ 2.8 against 1/12 per tail | constructed |

Our tool reproduces Gurisatti's 4 sets for her number. Among primes below 2·10⁷ with order dividing 720720 it finds no fifth, so her 4 is exact for that prime range.

**The eight tails** (Sierpinski; the Riesel tails are the same up to small changes in the third and sixth):
1. {17, 241};
2. {19, 37, 73};
3. {109, 257, 433, 673, 38737};
4. {11, 31, 41, 151, 331};
5. {29, 43, 113, 127, 337, 1429, 5419};
6. a 23-prime tail from 61 up to 18837001;
7. {23, 67, 89, 353, 397, 683, 2113, 7393, 20857, 312709, 599479, 4327489};
8. {53, 79, 157, 313, 937, 1249, 1613, 2731, 3121, 4993, 6553, 8191, 21841, 86113, 94849, 121369, 858001}.

The numbers k are in `results/data/prob92_k_S.txt` and `prob92_k_R.txt`.

**Five sets with 27 digits.** Keep Gurisatti's Sierpinski number modulo her 12 primes, add tail 4, and take the least odd k in the resulting class modulo the product of the 17 primes: k = 660664662814032650416080443. Its five minimal covering sets are her four plus {3, 5, 7, 13, 11, 31, 41, 151, 331}. This is not claimed to be the least such k (puzzle 1279, Q4), but it is far below the 253-digit construction.

**Q3 (significance).** Multiple covering sets cost almost nothing once the core leaves a small remainder. Any number of tails can be glued, so the number of covering sets of a Sierpinski number is unbounded. The interesting questions are the least k with m sets, and how many sets the known small Sierpinski and Riesel numbers have. The verification tool here answers the second question for any given k.

## Verify it yourself

[`results/code/prob92_verify.py`](code/prob92_verify.py) (under a second) recomputes, from each k alone, the class of every prime used. It checks that each of the 8 (Sierpinski), 8 (Riesel) and 5 (27-digit) sets covers all n mod 720720 and is minimal, and that k·2ⁿ ± 1 has a covering factor below itself for n < 3000. [`results/code/prob92_build.py`](code/prob92_build.py) rebuilds the 253-digit numbers.
