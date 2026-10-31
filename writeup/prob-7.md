# Problem 7: Parity of x and its inverse mod p

> **In short:** Two counts for x and its inverse x′ mod p, Np(all) = p − 2 and Np(oo) = Np(ee) − 1, were found by computation, not proved. The proof: (1) x ↦ p − x sends x′ to p − x′ and flips both parities; (2) so even-even and odd-odd pairs are equally many; (3) the counts omit x = 1, an odd-odd pair. This proves both for every odd modulus, so no counterexamples exist.

**Puzzle:** [primepuzzles.net](https://www.primepuzzles.net/problems/prob_007.htm)

**Status:** answered · Q2 proved (both relationships); Q1 explained; no counterexamples exist

## The question

For 0 < x < p let x′ be the inverse of x mod p, with 0 < x′ < p. D. H. Lehmer asked how often x and x′ have opposite parity. Example: p = 13 gives Np = 6. On the page, Anton Vrba found by computation that Np(all) = p − 2, Np(op) ≈ p/2, Np(ee) ≈ p/4 and Np(oo) = Np(ee) − 1. The questions:
1. Explain the theory behind this.
2. Prove Np(all) = p − 2 and Np(oo) = Np(ee) − 1.
3. Find counterexamples.

## The answer

**The involution.** For odd p put σ(x) = p − x.

1. (p − x)(p − x′) ≡ x·x′ ≡ 1 (mod p), so the inverse of σ(x) is σ(x′).
2. Since p is odd, p − x has the opposite parity to x.

So σ maps the pair (x, x′) to (p − x, p − x′) and flips both parities. Even-even goes to odd-odd and back, and opposite parity stays opposite. σ is a bijection on 1, …, p − 1, so over all x:

#ee = #oo.

**What Vrba counted.** His totals are p − 2, for example 11 for p = 13. They come from leaving out x = 1. The pair (1, 1) is odd-odd, so over x = 2, …, p − 1:

- Np(all) = (p − 1) − 1 = p − 2;
- Np(oo) = #oo − 1 = #ee − 1 = Np(ee) − 1.

The pair (p − 1, p − 1) stays in the count, and it is σ's image of (1, 1). That is the one even-even pair left without an odd-odd partner. Nothing beyond oddness of the modulus is used, so the same identity holds for every odd n, with x running over the units mod n other than 1. There are no counterexamples, and p = 2 has no x ≥ 2 at all.

**More from the same symmetry.** Np(op) is even, because x ↔ x′ pairs the opposite-parity cases and no x = x′ has opposite parity. Then Np(ee) = (p − 1 − Np(op))/2 and Np(oo) = Np(ee) − 1. So "Np(ee) ≈ p/4" is equivalent to "Np(op) ≈ p/2".

**Q1: why Np(op) ≈ p/2.** Lehmer's question is whether the parity of x and the parity of its inverse behave like independent coins. Zhang Wenpeng proved that they do: Np(op) = (p − 1)/2 + O(p^(1/2)·log²p) ("On a problem of D. H. Lehmer and its generalization", Compositio Math., 1993 and 1994; cited from memory, with no copy in `references/` to check against). The proof writes the parity conditions as sums of additive characters. The error term then becomes a sum of Kloosterman sums Σ e((ax + bx′)/p), and Weil's bound |K(a, b; p)| ≤ 2√p controls those. For 1000 < p < 5000 the measured deviation |Np(op) − (p − 1)/2| stays below 2.9·√p.

Example p = 13: ee = 3, oo = 2, op = 6, and 3 + 2 + 6 = 11 = p − 2.

## Validation

[`writeup/code/prob7_verify.py`](code/prob7_verify.py) does the following:

- checks Np(all) = p − 2, Np(oo) = Np(ee) − 1 and that Np(op) is even, for every odd prime p < 5000;
- checks the involution facts for each x;
- checks Np(oo) = Np(ee) − 1 for every odd composite n < 2000;
- prints the largest |Np(op) − (p − 1)/2|/√p for 1000 < p < 5000.
