# Conjecture 52: ⌊(p² + q²)/(2q)⌋ = p for consecutive primes

> **In short:** The conjecture says ⌊(p² + q²)/(2q)⌋ = p for consecutive primes p < q. (1) (p² + q²)/(2q) = p + g²/(2q) with g = q − p, so the claim is g² < 2q. (2) This holds below 2·10⁸ (direct check). Above that every gap below 4·10¹⁸ is at most 1476, and 1476² < 4·10⁸. So the conjecture holds for every p ≤ 4·10¹⁸; the general case is open. The proof is that of formula 2 in puzzle 1157.

**Puzzle:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_052.htm)

**Status:** proved for p ≤ 4·10¹⁸ (direct check below 2·10⁸ plus the complete gap table); at most O(X^(0.07+ε)) exceptions up to X; open in general (needs gaps below √(2q))

## The question

S. M. Ruiz: for consecutive primes p < q, ⌊(p² + q²)/(2q)⌋ = p. "Can you obtain the SMR relationship from another conjecture or theorem about consecutive primes?" On the page:

- a derivation from Cramér's conjecture (L. Rodríguez);
- a derivation from Eismann's conjecture that there is a prime in (p, p + √p) except for p = 3, 7, 13, 23, 31, 113 (F. Firoozbakht);
- two proofs from Chebyshev's theorem and from limits, both shown invalid by Firoozbakht, since ⌊x⌋ = ⌊y⌋ does not give ⌊xz⌋ = ⌊yz⌋.

The invalid step really matters. For the non-consecutive primes 7 < 13 one has ⌊7/13 + 13/7⌋ = 2, but ⌊(49 + 169)/26⌋ = 8 ≠ 7.

## Record vs ours

| | Status |
|---|---|
| Page | conditional derivations (Cramér, Eismann); two invalid proofs |
| **Ours** | **proved for every p ≤ 4·10¹⁸**, from the complete gap computation; O(X^(0.07+ε)) possible exceptions up to X |

## The math

**Step 1 (proved): the reduction.** (p² + q²)/(2q) = ((q − g)² + q²)/(2q) = q − g + g²/(2q) = p + g²/(2q). Since g²/(2q) > 0, the floor is p iff g² < 2q.

**Step 2 (computed): p < 2·10⁸.** g² < 2q holds for every pair of consecutive primes below 2·10⁸ ([puzzle 1157](1157.md), formula 2; its script checks the floor itself exactly).

**Step 3 (proved from the literature): 2·10⁸ ≤ p ≤ 4·10¹⁸.** Every gap below 4·10¹⁸ is at most 1476 (complete computation of Oliveira e Silva, Herzog and Pardi; FACTS 3023), and 1476² = 2178576 < 4·10⁸ < 2q.

**Step 4 (proved from the literature): beyond.** A failure needs a gap ≥ √(2q) > √q. Gaps of length ≥ √x in [x, 2x] have total length O(x^(0.57+ε)) (Järviniemi 2022, arXiv:2212.10965). So there are at most O(x^(0.07+ε)) of them per dyadic block, and O(X^(0.07+ε)) possible exceptions up to X.

## Validation

[`code/p1157_verify.py`](code/p1157_verify.py) (puzzle 1157, a few seconds) checks the floor and g² < 2q on every pair of consecutive primes below 2·10⁸, and asserts 1476² < 4·10⁸. The counterexample to the invalid step is a one-line check: ⌊7/13 + 13/7⌋ = 2 and ⌊218/26⌋ = 8.
