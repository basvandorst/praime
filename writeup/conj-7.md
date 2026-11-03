# Conjecture 7: Cramér's conjecture under the calibrated gap model (the record ratio is 1 − 2.2 ln ln x / ln x)

> **In short:** The page states Cramér's conjecture as limsup g_n / ln² p_n = 1 and holds one withdrawn proof. The sieve-and-coin gap model, calibrated below 4·10⁹ and checked out of sample on the record gap 1476 below 4·10¹⁸, gives: (1) among the maximal gaps below 4·10¹⁸, eleven have Cramér-Shanks-Granville ratio g/ln² p ≥ 0.80 and one has ≥ 0.85 (1132 after 1693182318746371, ratio 0.9206); the model expects 12.3 and 1.4, puts the record ratio at its 93rd percentile, and gives a ratio ≥ 1 below 4·10¹⁸ a chance of 0.2 %; (2) the median maximal gap below x is ln² x · (1 − 2.2 ln ln x / ln x) at every tested x from 10¹⁶ to 10⁴⁰ (ratios 0.785, 0.807, 0.817, 0.845, 0.865, 0.892, each within 0.001), against 1 − ln ln x / ln x in Cramér's own model; (3) the model's tail is steeper than Cramér's (1.24 times at 10¹⁶, 1.05 times at 10⁴⁰). So the ratio climbs toward 1 only like ln ln x / ln x (0.95 near 10¹⁰⁰), fewer than 0.01 gaps above ln² p are expected below 10⁴⁰, and the limsup itself (1 by Cramér, at least 2e^(−γ) = 1.1229 by Granville) is not decided by any computable range.

**Conjecture:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_007.htm)

**Status:** explained (heuristic): the record ratio below x follows 1 − 2.2 ln ln x/ln x for 10¹⁶ ≤ x ≤ 10⁴⁰ (model) · the real records 1132 (below 10¹⁶) and 1476 (below 4·10¹⁸) sit at the model's 90th and 48th percentiles · the record CSG ratio 0.9206 at the 93rd · the conjecture stays open

## The question

Cramér (1936): g_n = O(ln² p_n), and his probabilistic model suggests limsup g_n / ln² p_n = 1. Granville (1995, FACTS 3007) corrected the model for divisibility by small primes and found limsup ≥ 2e^(−γ) = 1.1229. Banks, Ford and Tao ([arXiv 1908.08613](https://arxiv.org/abs/1908.08613), Inventiones 2023) analysed a random-sieve model rigorously and expressed the size of the largest gap through the extremal bounds of the interval sieve problem. The page itself has no contribution besides a withdrawn claim.

## What the page had vs ours

| | On the page | Ours |
|---|---|---|
| Statement | limsup = 1 (Cramér) | not decided; the model's record ratio rises 0.785 → 0.892 from 10¹⁶ to 10⁴⁰ as 1 − 2.2 ln ln x/ln x |
| Data | none | eleven maximal gaps below 4·10¹⁸ with ratio ≥ 0.80 (model 12.3), one ≥ 0.85 (model 1.4), record 0.9206 at the model's 93rd percentile |
| Ratio ≥ 1 | the conjectured limit | 0.2 % chance below 4·10¹⁸; fewer than 0.01 expected below 10⁴⁰ |

## How we found it

The model and its extension to 10⁴⁰ are described in [conj-30.md](conj-30.md) and implemented in [conj30_model.py](code/conj30_model.py): survivors of the sieve to y = x^0.405, a fair coin ρ = 0.722 per survivor, the gap tail P(gap ≥ G) = E[(1 − ρ)^S] by importance sampling over the window's residues. Nothing is fitted; the only choice, the depth 0.405, is the crossing e^γ ω(u) = 1 of Buchstab's function (FINDINGS).

For the median record below x we solve N(G, x) = ln 2 for G, where N(G, x) = ∫ P(gap ≥ G near t) dt/ln t is the expected number of gaps ≥ G among primes in [x/10⁴, x] (below that the contribution is negligible). The 10 % and 90 % points are N = 2.30 and N = 0.105.

## The math

**Records in the data.** The ratios g/ln² p of the maximal gaps (Wikipedia's table, archived in references/web) above 0.80 are: 766 (0.818), 906 (0.831), 1132 (0.921), 1184 (0.806), 1198 (0.806), 1220 (0.805), 1328 (0.814), 1356 (0.825), 1370 (0.832), 1442 (0.848), 1476 (0.845); beyond 4·10¹⁸ the distributed searches add 1676 (0.847), 1724 (0.827), 1854 (0.874). Model counts of gaps with ratio ≥ r among primes in [10¹², 4·10¹⁸] (below 10¹² the ratios are smaller):

```
r        0.80    0.85    0.90    0.9206   0.95    1.00
model   12.3    1.43    0.162   0.067    0.019   0.002
real    ≥ 11    ≥ 1     1       1        0       0        (maximal gaps; non-maximal gaps can add to the first two)
```

P(no gap with ratio ≥ 0.9206) = e^(−0.067) = 0.935: the famous 1132 is a 1-in-15 event under the model, not a sign of a fatter tail.

**The median record.** With ln² x and the model's median G_med:

```
x        ln² x   G_med   G_med/ln² x   10 %–90 %      1 − 2.2 ln ln x/ln x   tail slope / Cramér's
10¹⁶     1357    1066    0.785         0.759–0.827    0.785                  1.24
4·10¹⁸   1835    1480    0.807         0.783–0.843    0.807                  1.20
10²⁰     2121    1733    0.817         0.795–0.852    0.817                  1.17
10²⁵     3314    2801    0.845         0.827–0.875    0.845                  1.12
10³⁰     4772    4127    0.865         0.849–0.890    0.865                  1.07
10⁴⁰     8483    7567    0.892         0.880–0.911    0.892                  1.05
```

The real records: 1132 below 10¹⁶ has ratio 0.834 to ln² 10¹⁶ and sits at the model's 90th percentile; 1476 below 4·10¹⁸ (ratio 0.805) at the 48th. The tail slope is d ln N/dG; Cramér's model gives −1/ln x, the sieve model is steeper by the factor in the last column, which fades as x grows.

**Reading.** In Cramér's model the record below x is ln² x (1 − ln ln x/ln x + …); the sieve pushes the coefficient of ln ln x/ln x from 1 to 2.2 over this range. Extrapolating the fitted form, the ratio reaches 0.95 near 10¹⁰⁰ and never exceeds 1; but the form is an interpolation over 10¹⁶–10⁴⁰, and Granville's mechanism (windows after multiples of primorials, present in the model at their natural frequency) is what should eventually carry the ratio past 1. The model shows how slow that is, not whether it happens.

## Validation

[code/conj30_model.py](code/conj30_model.py): `python3 writeup/code/conj30_model.py trend` (8 minutes) and `csg` (1 minute); importance sampling with 200 samples per probability, so reruns move the ratios by about 0.003 and the counts by 10-20 %.

```
x=1e+16 ln^2x=1357: median max gap 1066 (ratio 0.785), 10-90% 0.759-0.827, dlnN/dG -0.0336 (Cramer -0.0271)  real 1132 (ratio 0.834) at model quantile 0.90
x=4e+18 ln^2x=1835: median max gap 1480 (ratio 0.807), 10-90% 0.783-0.843, dlnN/dG -0.0279 (Cramer -0.0233)  real 1476 (ratio 0.805) at model quantile 0.48
x=1e+20 ln^2x=2121: median max gap 1733 (ratio 0.817), 10-90% 0.795-0.852, dlnN/dG -0.0253 (Cramer -0.0217)
x=1e+25 ln^2x=3314: median max gap 2801 (ratio 0.845), 10-90% 0.827-0.875, dlnN/dG -0.0195 (Cramer -0.0174)
x=1e+30 ln^2x=4772: median max gap 4127 (ratio 0.865), 10-90% 0.849-0.890, dlnN/dG -0.0155 (Cramer -0.0145)
x=1e+40 ln^2x=8483: median max gap 7567 (ratio 0.892), 10-90% 0.880-0.911, dlnN/dG -0.0114 (Cramer -0.0109)
ratio >= 0.80: 12.284 expected below 4e18 (from 1e12); 0.85: 1.433; 0.90: 0.162; 0.9206: 0.067; 0.95: 0.019; 1.00: 0.002
```
