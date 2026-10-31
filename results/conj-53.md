# Conjecture 53: Floor function and consecutive primes

> **In short:** For consecutive primes r < q < p, Capelle conjectured ⌊pq/r⌋ = p + q − r for every r except 15 small ones, and later listed 18 more such identities. All 19 follow from one inequality: **the two gaps after r sum to less than √r**, that is (p − r)² < r. That inequality fails only for r ≤ 1327 below 2·10⁸, and above 2·10⁸ it holds up to 4·10¹⁸ because no gap there exceeds 1476. So all 19 identities are **proved for 1327 < r with p ≤ 4·10¹⁸**, and beyond that at most O(X^(0.13+ε)) values of r ≤ X can fail. Capelle also asked for a simple equation for ⌊pq/r − qr/p⌋. The answer: **⌊pq/r − qr/p⌋ = 2(p − r) − [p − q > q − r]**, which subtracts 1 exactly when the second gap is the larger one. It holds for every r > 113.

**Puzzle:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_053.htm)

**Status:** answered (Capelle's question: the identity for ⌊pq/r − qr/p⌋) · all 19 identities proved for r > 1327 up to 4·10¹⁸, and for all but O(X^(0.13+ε)) r beyond; the full conjecture is open (it needs gaps < √r)

## The question

Patrick Capelle: if r < q < p are consecutive primes and r is not one of 2, 3, 5, 7, 11, 13, 19, 23, 31, 47, 83, 89, 113, 199, 1327, then ⌊pq/r⌋ = p + q − r. Prove it or find a counterexample.

On the page:

- Giovanni Resta, Phil Carmody and Kermite Rose showed it is equivalent to (p − r)(q − r) < r, which is "dangerously similar to Legendre's conjecture".
- Farideh Firoozbakht derived it from her conjecture 30 (gaps below 3 ln² p) plus a check of r < 201770. That proof is conditional.
- Capelle (5/5/07) added 18 more identities, conjecturally true for r > 1327 (numbered 2–19 below), and asked: "Can you propose a simple equation with floor(p*q/r − q*r/p)?"
- Werner Sand noted the condition is not about primes as such.

## Record vs ours

| | Status |
|---|---|
| Page | identity 1 equivalent to (p − r)(q − r) < r; proved assuming Firoozbakht's conjecture; identities 2–19 conjectured; ⌊pq/r − qr/p⌋ asked |
| **Ours** | **all 19 identities reduce to (p − r)² < r; proved unconditionally for 1327 < r, p ≤ 4·10¹⁸; failures O(X^(0.13+ε)) beyond; ⌊pq/r − qr/p⌋ = 2(p − r) − [p − q > q − r]** |

## The math

Write d = q − r and e = p − q for the two gaps. Three exact identities of rational numbers carry everything:

```
pq/r = (p + q − r) + A,   A = d(d + e)/r
qr/p = (q + r − p) + B,   B = e(d + e)/p
rp/q = (r + p − q) − C,   C = de/q
```

For example, pq − r(p + q − r) = (p − r)(q − r) = (d + e)·d. A, B and C are positive, and none of them is an integer, because r, p and q are primes larger than the gaps.

The chain of claims, each validated on every prime triple before it is used:

| # | Claim | Why | Validated |
|---|---|---|---|
| 1 | identity 1 ⇔ A < 1; identity 2 ⇔ B < 1; identity 3 ⇔ C < 1 | the three expansions above | failures of each, r < 10⁶, from the definitions |
| 2 | identities 10–19 follow from 1–3 | the three floors are then p + q − r, q + r − p and r + p − q − 1; e.g. 17 is (p + q − r)(q + r − p) = q² − (p − r)² | same |
| 3 | identities 4–9 ⇔ A + B < 1, 0 ≤ A − C < 1, 0 ≤ B − C < 1, A + C < 1, B + C < 1, 0 ≤ A + B − C < 1 | add the expansions | same |
| 4 | A > C and B > C always | A > C ⇔ (d + e)q > er; B > C ⇔ (d + e)q > dp ⇔ er > 0 | algebra |
| 5 | (p − r)² < r implies all of 1–19 | A + B, A + C and B + C are each ≤ (d + e)²/r, since B ≤ e(d + e)/r and C ≤ de/r | algebra; no identity fails above r = 1327 |
| 6 | (p − r)² < r for every r > 1327 with p < 2·10⁸ | direct check | all prime triples below 2·10⁸ |
| 7 | (p − r)² < r for 2·10⁸ < r, p ≤ 4·10¹⁸ | every gap below 4·10¹⁸ is ≤ 1476, so (p − r)² ≤ 2952² = 8714304 | theorem from a complete computation (literature) |
| 8 | #{r ≤ X : (p − r)² ≥ r} ≪ X^(0.13+ε) | a failure needs a gap ≥ √r/2; Järviniemi: gaps ≥ x^0.45 in [x, 2x] have total length ≪ x^(0.63+ε) | theorem (literature + one line) |
| 9 | ⌊pq/r − qr/p⌋ = 2(p − r) − [e > d] whenever identities 1 and 2 hold | pq/r − qr/p = 2(p − r) + A − B, and A − B has the sign of (d − e)r + d(d + e) | failures r < 10⁶: only 2, 3, 5, 7, 13, 19, 23, 109, 113 |

**Claim 9 in detail.** A − B = (d + e)(dp − er)/(rp), and dp − er = (d − e)r + d(d + e).
- If d ≥ e, this is positive, so 0 < A − B < 1 and the floor is 2(p − r).
- If d < e, then (d − e)r + d(d + e) ≤ −r + d(d + e) < 0, using A < 1. So −1 < A − B < 0, and the floor is 2(p − r) − 1.

In words: ⌊pq/r − qr/p⌋ is twice the distance from r to p, minus 1 when the second gap is the larger one. It fails exactly where identity 1 or 2 fails in a way that tips the floor: r ∈ {2, 3, 5, 7, 13, 19, 23, 109, 113}.

**The exceptions.** Below 2·10⁸ the failures of each identity, computed from the definitions, are:
- identity 1: Capelle's 15 values;
- identity 2: 2, 3, 5, 11, 17, 19, 29, 43, 47, 83, 109, 199, 283;
- identity 3: 5, 19.

The largest failure of any identity is r = 1327, where the gap 34 after 1327 is a famous maximal gap. The condition (p − r)² ≥ r holds only for 34 primes r, from 2 to 1327 (the largest are 523, 1321 and 1327).

**What stays open.** A proof for all r needs every pair of consecutive gaps below √r, which is Legendre/Oppermann territory. It does not follow even from the Riemann hypothesis (gaps O(√r ln r)). Firoozbakht's conjecture would be enough, as she noted. Our range 4·10¹⁸ is where her conjecture itself has been verified (Kourbatov 2015), so both routes stop at the same place.

## Verify it yourself

[`results/code/conj53_verify.py`](code/conj53_verify.py) (a few seconds) checks claim 6 on all prime triples below 2·10⁸. It evaluates all 19 identities and the new one from their definitions with exact rationals for every r < 10⁶, and asserts the failure lists above and that nothing fails above 1327. Claims 7 and 8 rest on the published maximal gap 1476 below 4·10¹⁸ and on Järviniemi's theorem (arXiv:2212.10965).
