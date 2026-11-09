# Conjecture 30: Firoozbakht's conjecture under the calibrated gap model (no violation expected below 10⁴⁰)

> **In short:** Firoozbakht's conjecture, p_(n+1) < p_n^(1+1/n), amounts to the gap bound g_n < ln² p_n − ln p_n − b and is verified below 4·10¹⁸. (1) A sieve-and-coin gap model, which predicted the record gaps to 4·10¹⁸, expects only 0.010 violating gaps between 4·10¹⁸ and 10⁴⁰. (2) Typical record gaps stay below the bound, so only a fluctuation can cross it. So the conjecture is safe for every computable range (99 % below 10⁴⁰), though a refinement of Cramér's model says it must fail eventually.

**Conjecture:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_030.htm)

**Status:** explained (heuristic): 0.010 violations expected below 10⁴⁰, rate falling · the median record gap is ln² x (1 − 2.2 ln ln x/ln x) for 10¹⁶ ≤ x ≤ 10⁴⁰ (model, 6 points within 0.001) · verified below 4·10¹⁸ (Kourbatov, literature) · the conjecture itself stays open

## The question

Firoozbakht (2002, on the page): (p_n)^(1/n) decreases. Q1: prove or disconfirm. Q2: consequences if true. The page shows (Bukh, Rodríguez, El Aidi) that it implies Cramér's g_n = O(ln² p_n), with Rodríguez's sharper form g_n < ln p (ln p − 1). Kourbatov (2015, [arXiv 1506.03042](https://arxiv.org/abs/1506.03042), on the page) proved: the conjecture implies g_n < ln² p_n − ln p_n − 1 for n > 9, and g_n < ln² p_n − ln p_n − 1.17 for all n > 9 implies the conjecture; and the conjecture holds for all p_n < 4·10¹⁸ (the complete gap table). Three claimed proofs on the page were withdrawn or are unrefereed.

## What the page had vs ours

| | On the page | Ours |
|---|---|---|
| Q1 | verified to 4·10¹⁸; no proof; withdrawn proofs | **model forecast**: 0.0012 violations expected below 10²⁰, 0.010 below 10⁴⁰; the rate per unit ln p falls (3.1·10⁻⁴ at 4·10¹⁸, 7.7·10⁻⁵ at 10⁴⁰) |
| Where records sit | Rodríguez: the 1132 gap "passes" the bound | the median record gap below x is ln² x (1 − 2.2 ln ln x/ln x); Firoozbakht's line is ln² x (1 − 1/ln x); the two never meet, only fluctuations cross |
| Q2 | Cramér, Rodríguez's bound | nothing new; Kourbatov's bounds are the sharp form |

## How we found it

The model is the one of FINDINGS ("one explicit law for the next prime") and [facts_maxgap_is.py](code/facts_maxgap_is.py): integers with no prime factor ≤ y = x^0.405 survive, each survivor is prime with probability ρ = 1/(ln x ∏_(q≤y)(1 − 1/q) e^γ ω(1/0.405)) = 0.722, independently. P(gap ≥ G after a prime near x) = E[(1 − ρ)^S] with S the number of survivors in the window, computed by importance sampling over the residues of the window start modulo the primes below G and an exact one-kill-per-prime process for the primes between G and 2·10⁵. Calibration (FINDINGS): gaps ≥ 100, 200, 240 in [10⁸, 4·10⁹] within 1-8 %; the record 1476 below 4·10¹⁸ at the model's 45th percentile where Cramér's model expects 56 such gaps.

Two changes for this page ([conj30_model.py](code/conj30_model.py)): the Mertens product is continued past the prime table (3·10⁷) with ∏(1 − 1/q) ∝ 1/ln y, and the death process runs in log space, because (1 − ρ)^G underflows for G > 2000. That takes the model to x = 10⁴⁰ at 2-15 s per probability.

Before trusting it near the Firoozbakht line (gap ratio g/ln² p ≈ 0.98) we checked the only data there is, the Cramér-Shanks-Granville ratios of the maximal gaps below 4·10¹⁸ ([conj-7](conj-7.md)): the model expects 12.3 gaps with ratio ≥ 0.80 (the table has 11 maximal ones, plus any non-maximal), 1.4 with ratio ≥ 0.85 (real: one, the gap 1132 at 0.9206), and puts that record ratio at its 93rd percentile. So the model is right in count where the data reach, and the record is a 1-in-15 event, not a sign of a fatter tail.

## The math

**Violations.** A gap g after a prime p violates the conjecture only if g ≥ ln² p − ln p − 1 (Kourbatov's necessary condition), so the expected number of violations among primes in (4·10¹⁸, X] is at most

```
N_F(X) = ∫ P(gap ≥ ln² t − ln t − 1 after a prime near t) · dt / ln t ,   t from 4·10¹⁸ to X.
```

Model values of the integrand per unit ln t: 3.1·10⁻⁴ (4·10¹⁸), 4.2·10⁻⁴ (10²⁰), 2.6·10⁻⁴ (10²²), 2.2·10⁻⁴ (10²⁵), 2.0·10⁻⁴ (10³⁰), 1.6·10⁻⁴ (10³⁵), 7.7·10⁻⁵ (10⁴⁰); sampling error 9-18 % per point. Integrated: N_F(10²⁰) = 0.0012, N_F(10²²) = 0.0027, N_F(10²⁵) = 0.0044, N_F(10³⁰) = 0.0068, N_F(10³⁵) = 0.0089, N_F(10⁴⁰) = 0.0102. With a falling integrand the total beyond 10⁴⁰ grows slower than ln ln X; no range anyone can compute reaches an expected count of 1.

**Why the typical record stays below the line.** The median maximal gap below x in the model (6 points, [conj-7](conj-7.md)) is

```
G_med(x) = ln² x · (1 − 2.2 ln ln x / ln x)        (10¹⁶ ≤ x ≤ 10⁴⁰, each point within 0.001 of the ratio)
```

against Cramér's own model, where the record is ln² x (1 − ln ln x / ln x). Firoozbakht's line is ln² x − ln x = ln² x (1 − 1/ln x). The record is below the line by (2.2 ln ln x − 1) ln x, which grows; a violation therefore needs a gap far out in the tail of the record distribution, and the model's tail is steeper than Cramér's (the slope d ln N / dG is 1.24 times Cramér's at 10¹⁶, 1.05 times at 10⁴⁰). That is why the expected count stays at 10⁻⁴ per unit ln p.

**What this does not say.** Granville's refinement (FACTS 3007) and the Banks-Ford-Tao random sieve model ([arXiv 1908.08613](https://arxiv.org/abs/1908.08613), Inventiones 2023) argue that limsup g_n / ln² p_n ≥ 2e^(−γ) = 1.1229, in which case Firoozbakht's conjecture is false infinitely often. Our model contains the same mechanism (the exact sieve, with Maier-type windows at their natural frequency), and its ratio does rise, but only like 1 − 2.2 ln ln x/ln x: it reaches 0.95 near 10¹⁰⁰. The two statements are compatible: the conjecture is expected to fail, and nowhere anyone can look.

## Validation

[code/conj30_model.py](code/conj30_model.py): `python3 writeup/code/conj30_model.py firoozbakht` (1 minute), `trend` (8 minutes), `csg` (1 minute), `T <x> <G>` for one probability. Importance sampling with 200 samples per probability; reruns differ by the quoted 9-18 % per point and by about 0.003 in the ratios.

```
expected number of gaps g >= ln^2 p - ln p - 1 among primes above 4e18:
  t=4e+18 G=1790 (ratio 0.9757) P(gap>=G)=3.33e-21 (rel 0.09) dN/dln t=3.11e-04
  t=1e+20 G=2073 (ratio 0.9775) P(gap>=G)=1.91e-22 (rel 0.11) dN/dln t=4.15e-04
  t=1e+22 G=2514 (ratio 0.9797) P(gap>=G)=1.33e-24 (rel 0.13) dN/dln t=2.62e-04
  t=1e+25 G=3255 (ratio 0.9823) P(gap>=G)=1.28e-27 (rel 0.18) dN/dln t=2.23e-04
  t=1e+30 G=4701 (ratio 0.9852) P(gap>=G)=1.37e-32 (rel 0.16) dN/dln t=1.98e-04
  t=1e+35 G=6413 (ratio 0.9874) P(gap>=G)=1.27e-37 (rel 0.14) dN/dln t=1.57e-04
  t=1e+40 G=8389 (ratio 0.9889) P(gap>=G)=7.05e-43 (rel 0.13) dN/dln t=7.65e-05
  cumulative N_F: 10^20 1.17e-3, 10^22 2.73e-3, 10^25 4.40e-3, 10^30 6.82e-3, 10^35 8.87e-3, 10^40 1.02e-2
```

The trend and the CSG calibration are in [conj-7.md](conj-7.md).
