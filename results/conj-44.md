# Conjecture 44: A property of prime twins, only?

> **In short:** For consecutive primes p < q, z = √((p² + q²)/2 − 1) is an integer exactly when q = p + 2. Mike Oakes showed on the page that this is equivalent to q − p < 2(√(2p + 3) + 1) for every non-twin pair, and checked it to 10⁹. We prove it **for every p ≤ 4·10¹⁸**: a direct check below 2·10⁸, then the maximal gap 1476 is far below the bound. Beyond that, at most O(X^(0.07+ε)) primes p ≤ X can fail (Järviniemi's bound on large gaps). A full proof needs gaps below about 2√(2p), which nobody can prove, not even under the Riemann hypothesis.

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

## Verify it yourself

[`results/code/conj44_verify.py`](code/conj44_verify.py) (about 2 seconds) runs the exact square test on every pair of consecutive odd primes below 2·10⁸. It checks that z is an integer exactly for the twin pairs, and that Oakes' bound holds for every non-twin pair. Claim 4 uses the published maximal gap 1476 below 4·10¹⁸.
