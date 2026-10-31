# Conjecture 44: A property of prime twins, only?

> **In short:** It was conjectured that for consecutive primes p < q, √((p² + q²)/2 − 1) is an integer only for twins; this reduces to a gap bound of about 2√(2p). The proof: (1) a direct check below 2·10⁸; (2) every gap below 4·10¹⁸ is at most 1476, far below the bound. This proves it for every p ≤ 4·10¹⁸; the general case stays open.

**Puzzle:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_044.htm)

**Status:** proved for p ≤ 4·10¹⁸ (direct check below 2·10⁸ plus the complete gap table) · at most O(X^(0.07+ε)) exceptions up to X · open in general

## The question

Sebastián Martín Ruiz: for consecutive primes p < q, z = √((p² + q²)/2 − 1) is an integer if and only if p and q are twin primes. Prove it or show it false.

On the page:

- Farideh Firoozbakht derived it from her conjecture 30 (gaps below 2 ln² p) plus a check below 1229.
- Danielle Degiorgi showed that a non-twin pair giving an integer needs q − p ≥ 2√(2p + 3) + 2.
- Mike Oakes wrote the same bound as (t − 1)² < 2p + 3 with q = p + 2t, and checked all p ≤ 10⁹.
- Robin García parametrised the possible exceptions.

All agree a full proof needs a gap bound of size √p.

## Record vs ours

| | Status |
|---|---|
| Page | reduced to gaps < 2(√(2p+3) + 1); checked p ≤ 10⁹ (Oakes); conditional proofs |
| **Ours** | **proved for every p ≤ 4·10¹⁸; exceptions up to X at most O(X^(0.07+ε))** |

## The math

| # | Claim | Why | Validated |
|---|---|---|---|
| 1 | with q = p + 2t: z² = (p + t)² + (t² − 1) | expand (p² + (p + 2t)²)/2 − 1 | algebra |
| 2 | t = 1 ⇒ z = p + 1; 1 < t < √(2p + 3) + 1 ⇒ (p + t)² < z² < (p + t + 1)², so z is not an integer | 0 < t² − 1 < 2(p + t) + 1 ⇔ (t − 1)² < 2p + 3 (Oakes) | exact square test on every pair below 2·10⁸ |
| 3 | the bound holds for every pair below 2·10⁸ | direct check | all consecutive odd primes below 2·10⁸ |
| 4 | the bound holds for 2·10⁸ < p ≤ 4·10¹⁸ | every gap there is ≤ 1476 (complete computation of Oliveira e Silva, Herzog and Pardi), so (t − 1)² ≤ 737² < 2p + 3 | theorem (literature) |
| 5 | #{p ≤ X : the bound fails} ≪ X^(0.07+ε) | a failure needs a gap ≥ 2√(2p) > √p; gaps ≥ √x in [x, 2x] have total length ≪ x^(0.57+ε) (Järviniemi, arXiv:2212.10965) | theorem (literature + one line) |

Claims 1–4 prove the conjecture for all consecutive primes with p ≤ 4·10¹⁸. Claim 5 bounds the possible exceptions beyond. This is the same pattern as [puzzle 997](997.md) and [conjecture 53](conj-53.md): a statement about consecutive primes that fails only if some gap reaches the size √p.

**Beyond 4·10¹⁸ (conditional).** The published list of 85 maximal gaps (OEIS A005250; the 85th is 1854 after 101412319996363309069, found in 2026) comes from distributed searches with no refereed statement of completeness (FACTS 3023). If that list is complete, every gap below 1.01·10²⁰ is at most 1854, and the same argument extends the result to p < 1.01·10²⁰, since 1854² ≈ 3.4·10⁶ and (2·1854)² ≈ 1.4·10⁷ are both far below 2·10⁸.

## Validation

[`writeup/code/conj44_verify.py`](code/conj44_verify.py) (about 2 seconds) runs the exact square test on every pair of consecutive odd primes below 2·10⁸. It checks that z is an integer exactly for the twin pairs, and that Oakes' bound holds for every non-twin pair. Claim 4 uses the published maximal gap 1476 below 4·10¹⁸.
