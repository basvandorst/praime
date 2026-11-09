# Conjecture 66: the maximal gap between twin pairs is not 0.45 ln³ p but (1/(2C₂)) ln³ p times a slowly rising factor (Kourbatov 2013)

> **In short:** The conjecture says the record gap between consecutive twin pairs grows like 0.45 (ln p₁)³ and asks for a justification or a better formula. (1) The Hardy–Littlewood twin count gives the record gap near p as a(ln(p/a) − b), with a = ln² p/(2C₂) the mean twin gap (literature). (2) So gap/ln³ p tends to 1/(2C₂) = 0.757, very slowly; 0.45 is only its value in the computed range. (3) The 22 OEIS records above 1000 follow this formula with mean ratio 1.036.

**Conjecture:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_066.htm)

**Status:** explained (literature: Kourbatov 2013, in references/2013_maximal-gaps-between-prime-k-tuples) · the constant is not 0.45 · the records to 2.6·10⁸ checked against the formula

## The question

Luis Rodríguez: the maximal gap between consecutive twin pairs near p₁ is about k (ln p₁)³ with k ≈ 0.45. Q1: justify it or suggest a better one. Carlos Rivera added the champion gaps below 2³² and their ratios gap/(ln p₁)³ (6.0, 0.53, 0.35, 0.39, 0.19, 0.36, 0.55, …, 0.55).

## What the page had vs ours

| | On the page | Ours |
|---|---|---|
| Q1 | an empirical 0.45 with ratios scattered between 0.19 and 0.57 | Kourbatov's a(ln(p/a) − b), a = ln² p/(2C₂), b ≈ 1: the ratio is (1/(2C₂))(1 − (2 ln ln p + c)/ln p + …) and tends to 0.757 |
| Data | champions below 2³² | the 22 OEIS records (A113274, A113275) above 1000: gap/formula has mean 1.036, sd 0.15 |

## The math

Twin pairs near p have mean spacing a = ln² p/(2C₂) (Hardy-Littlewood). If the spacings behaved like exponential variables with mean a, the largest of the N ≈ p/a spacings below p would be about a ln N = a ln(p/a), the Cramér argument transferred to twins. Kourbatov (2013) writes the record near p as

```
G(p) ≈ a (ln(p/a) − b),    a = ln² p / (2C₂),    b ≈ 1 (b = 1.26 makes the formula median-unbiased on the records to 10¹⁵)
```

Expanding ln(p/a) = ln p − 2 ln ln p + ln(2C₂):

```
G(p)/ln³ p ≈ (1/(2C₂)) · (1 − (2 ln ln p − ln(2C₂) + b)/ln p).
```

With 1/(2C₂) = 0.7574 this ratio is 0.44 at p ≈ 10⁷, 0.51 at 10⁹, 0.56 at 10¹², 0.59 at 10¹⁵, 0.62 at 10²⁰ and 0.66 at 10³⁰. Rodríguez's 0.45 is the value of a slowly rising function where the data happened to be. Our sieve-and-coin model for prime gaps ([conj-7](conj-7.md)) shows the same mechanism for single primes: the record ratio to ln² p rises as 1 − 2.2 ln ln p/ln p.

## Validation

[code/conj66_check.py](code/conj66_check.py) holds the 30 OEIS records (A113274 gaps, A113275 starting primes) and prints, for the 22 above 1000, the ratio gap/ln³ p₁, Kourbatov's formula with b = 1 and their quotient.

```
Rodriguez's ratio gap/ln^3 p1: min 0.321, max 0.570, mean 0.440 over 22 records
Kourbatov's formula with b = 1: gap/formula mean 1.036, sd 0.150
ratio implied by Kourbatov's formula at p = 1e9, 1e12, 1e15, 1e20, 1e30: 0.509, 0.556, 0.586, 0.620, 0.657
```
