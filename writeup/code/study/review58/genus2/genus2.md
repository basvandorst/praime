# The genus-2 deficit of the one-point law: a search for the cause

Object: fact 2055, `Σ_(0<γ≤T) e^(2iθ_L(γ)) = −ε·c_L·T/2π + o(T)`, `c_L = Π_p G_p(1/p)` (fact 2057), on the 100 genus-2
L-functions of rounds 151 (15), 152 (40) and the review-58 replication (45). Registered estimator throughout:
`b = Σ Re(G·conj pred)/Σ|pred|²`, jackknife error over L-functions. All numbers below come from the scripts named with them
(own code; the repo's zeros and constants were loaded as data only). Labels: **proved**, **computed**, **heuristic**,
**post hoc** (a reading chosen after seeing the data).

## Verdict

**No cause found among H1–H5; each is rejected by a number.** What the data do show (computed):

* The deficit is robust inside the sample: 0.929 ± 0.020 (bootstrap sd 0.0195, P(b ≥ 1) = 0.0005 in 4,000 resamples);
  without round 151 (the sample that first showed it) 0.931 ± 0.022; without the one non-USp(4) curve 0.932 ± 0.022.
* It is **uniform in height and in window shape**: thirds of (0, 150] read 0.903, 0.968, 0.915 (± 0.033–0.044); Hann
  0.935 ± 0.031, Tukey(0.1) 0.946 ± 0.024; the r = 1 line keeps its sinc shape (0.943, 0.932 at u = ±π/T).
* It **sits in the r = 1 line only**: the lines r = 2, ½, 3, ⅓ read 0.982 ± 0.027 pooled (root number +1: 1.013 ± 0.041),
  1.6 errors above r = 1 (2.1 for root number +1).
* It **sits in the arithmetic part of the constant**: `Re G = b₁(−εT/2π) + b₂(−ε(c − 1)T/2π)` gives b₁ = 0.983 ± 0.048 (the
  all-integers term "1" of fact 2053) and **b₂ = 0.895 ± 0.028** (the part c − 1). A proportional excess fits better than an
  additive one on (10, 150] (δ = 0.067 ± 0.021, AIC 482.2 against 486.4).
* The per-curve shape of c is right: alternative constants c·L(1, Λ²₀)^(±½, ±1) raise the residual sum of squares from 4,010
  to 6,935–18,918.

Two premises of the brief need correcting (computed):

1. "Nothing else reads low" is not established. Sym³ and the Rankin–Selberg products pooled read **0.989 ± 0.043**, which is
   **1.3 errors** from genus 2 (difference −0.061 ± 0.048); all 195 degree-4 L-functions read 0.938 ± 0.018. Degree 3
   (round 151) is 1.003 ± 0.040. Only degrees 1 and 2 exclude a 7 % deficit (elliptic curves at r = 1: 0.989 ± 0.011 sharp,
   1.002 ± 0.014 Hann; round 148 `s2_out.txt`). So genus 2 shows a degree-4 reading most sharply (it carries the weight),
   and whether the effect belongs to genus 2 alone or to degree 4 is **undecided**.
2. In the analytic normalisation the three classes have three different gamma types: genus 2 μ = (½, ½, 3/2, 3/2),
   Sym³ μ = (½, 3/2, 3/2, 5/2), products μ = (0, 1, 1, 2) (from the frozen prediction files).

Surviving: **S1** chance (p ≈ 5·10⁻⁴ for the sharp reading, 2.1–2.3 errors for smooth windows); **S2** a real deficit of
about 7 % of the r = 1 reading in degree 4 at these heights, in the part c − 1, with no derivation (the ratios recipe predicts
none, below); **S2′** (post hoc) an effective truncation of the singular series at R ≈ 100. A registration draft that
separates S1 from S2 on new curves is at the end.

## H1. The constant c is biased: rejected (computed)

Recomputed for 51 curves: all 37 with root number +1 and 14 with −1 spread over c from −0.94 to 1.59; **93.8 % of the slope
weight Σc²** ([c_own.py](c_own.py), [h1_analyze.py](h1_analyze.py) → [h1_out.txt](h1_out.txt), table at the end).
Own code: exact local factors for p ≤ 1000 (hyperellcharpoly at odd good p, lfuneuler at p = 2 and the bad primes), own point
counts on (2y + h)² = 4f + h² for 1000 < p ≤ 10⁵ (equal to lfunan at every p ≤ 2·10⁴ for all 51 curves: 0 mismatches), and the
first-order factor 1 − (e₁² − 1)/p beyond p = 1000.

| check | result |
|---|---|
| registered rule (log-Cesàro over [5·10³, 5·10⁴], ν = 2) on own data vs the registered c | max 2.1·10⁻⁴ relative, mean −1·10⁻⁵ |
| new c (same rule at 10⁵) / old c − 1 | mean +0.00005, sd 0.0012, range −0.0030 … +0.0024 |
| out of sample: rule at 5·10⁴ vs raw product to 10⁵ | rms 0.22 %, max 0.61 % |
| rule at 5·10⁴ vs rule at 10⁵ | rms 0.12 %, max 0.30 % |
| ν = 0, 1, 3 instead of 2 (at 10⁵) | ≤ 0.12 % |
| tail size K = (ln c(10⁵) − ln c(P))·√P·ln P on [10³, 3·10⁴] | rms 3.3, max 5.9 ⇒ error of ln c at 10⁵: 9·10⁻⁴ rms, 1.6·10⁻³ max |
| second-order terms neglected beyond p = 1000 | ≤ 5·10⁻⁵ |
| slope, all 100 / root number +1 / rounds 151+152 / replication | 0.9286 → 0.9288, 0.9174 → 0.9176, 0.9137 → 0.9140, 0.9480 → 0.9481 |
| correlation of the change with c / ln N / the size of the tail correction | −0.18 / −0.17 / −0.87 |

A 7 % bias would need a tail about 23 times the largest change seen between P = 5·10⁴ and 10⁵. The correlation −0.87 says the
small changes are the tail correction itself, not a trend with c.

## H2. The form of the constant is wrong for GSp(4): rejected (proved algebra, computed checks)

[h2_local.py](h2_local.py) → [h2_local_out.txt](h2_local_out.txt), 20,000 random Satake sets each (USp(4) Haar by the Weyl
density, Sym³ of SU(2), SU(2) × SU(2)):

* F_p(1) = 0 to 2·10⁻¹⁴; G_p = F_p/(1 − X) has degree 3 with g₁ = 1 − |λ|² and g₃ = −conj h₄; for USp(4)
  **g₂ = 1 − e₁² + e₁²e₂ − e₂²** exactly (7·10⁻¹⁵). This is fact 2057's polynomial; the algebra is a theorem (Lemma 6 of the
  degree-d note) and the check agrees.
* |λ|² − 1 = χ₁₀ + χ₅ (Sym² ⊕ Λ²₀) for USp(4), S₂ + S₄ + S₆ for Sym³, S₂ + S₂′ + S₂S₂′ for products (to 2·10⁻¹⁴).
* `B_p = F_p/Π_(i,j)(1 − α_i ᾱ_j X) = 1 + O(X²)` for all three (coefficient of X below 7·10⁻¹⁵), so Π G_p(1/p) = B(1)/L(1, Ad)
  with **Ad the 15-dimensional L × L̃/ζ**, which for USp(4) is Sym²(std) ⊕ Λ²₀ (factorisation checked to 2·10⁻¹⁵). The
  derivation needs only this Euler product; neither L(1, Sym²) nor L(1, Λ²₀) is ever computed, for any class.
* Haar: χ₁₀ and χ₅ are irreducible and orthogonal (norms 1.009, 1.016; Frobenius–Schur indicators +0.995, +0.998), so the
  prime-square bias is ν = 2; Sym³ and products give ν = 3.02 and 3.01. Genus 2 and Sym³ are both symplectic
  (E tr Λ² = 1.000 and 0.997, E tr g² = −0.99 and −1.02); products are orthogonal (0.010 and +0.99). "Λ² contains the
  trivial representation" holds for Sym³ exactly as for genus 2, and it enters only the sign of the prime-square terms.

What is the same and what differs (read from `c6_constants.py`): one function `direct()` computes c for every class as
Π F_p(1/p)/(1 − 1/p) from lfuneuler; the classes differ only in P0 (1000 vs 3000), the end of the first-order range (5·10⁴
vs 10⁵) and ν (2 vs 3). By H1 these change c by at most 0.3 %.

Alternative constants (post hoc, [h6_models.py](h6_models.py) → [h6_out.txt](h6_out.txt)): with L5 = Π_(p≤1000) L_p(1, Λ²₀)
(range 0.28–1.85 over the 51 curves), the constant c·L5^k gives residual sums of squares 18,918 / 7,953 / **4,010** / 6,935 /
16,494 for k = −1, −½, 0, ½, 1. So a missing or extra degree-5 factor is rejected.

Bad primes ([h2_badfactor.py](h2_badfactor.py) → [h2_bad_out.txt](h2_bad_out.txt)): the factors are PARI's lfuneuler for
objects that pass lfuncheckfeq at −56 to −65 (all 85 box curves in `genus2.json`, and the 15 of round 151 in `blind/testset.json`). The check is sensitive to the factor at 2: for
G2_N2862_17 (2 ‖ N, factor 1/(1 + X² + 2X³)), removing it or flipping its top coefficient takes lfuncheckfeq from −64 to −5.
Curves with root number +1 and 2 | N read 0.896 ± 0.039 (n = 7) against 0.924 ± 0.023 (n = 30). The product of the
bad-prime factors of c has median 1.061 (range 0.50–2.59).

## H3. A finite-height term: rejected as the cause (heuristic theory, computed data)

*Theory (heuristic: the ratios recipe).* On Re s = ½ + δ, X(s)⁻¹ = ε·L̃(1 − s)/L(s) exactly, so the right-line integral is
ε∫h·L̃(1 − s)L′(s)/L(s)² ds/(2πi); one integration by parts gives ε∫[h′L̃(1 − s) − hL̃′(1 − s)]/L(s). The h′ term is exactly
∫h′X⁻¹ (bounded, the edges). For ⟨L̃(1 − s + α)/L(s)⟩ the recipe keeps the diagonal A(α) = Σ conj λ(n)μ_L(n)n^(−1−α) =
B(1 + α)/L(1 + α, L × L̃) and has no swap term (a single X-factor oscillates). A(0) = 0 and A′(0) = c, so the main term is
exactly −ε·c·∫h/2π, with no term in ln T, no 1/ln(QT⁴) and no dependence on the gamma factor (the same answer as Route B of
the degree-d note, step 5, where F_L(1) = 0 removes the 2θ′(t) term). The same recipe written with two numerators and two
denominators has a swap term that carries the whole main term, again with no t-dependence.

*The kernel expansion.* The degree-d note keeps the closed-form leading kernel (error ≪ T^((d−1)/2)) or K ≥ ⌈d/2⌉ terms; its
first correction is the factor 1 + iκ/t, κ = (d² − 1)/(24d) + ½Σ(μ_j − μ_j²): genus 2 **−0.344**, Sym³ **−2.344**,
products **−0.844**. It belongs to the closed-form prime-side kernel only; the zero-side main term does not use it. Even
there it is a phase at first order, and for genus 2 seven times smaller than for Sym³ (|κ|/t = 0.005 at t = 75).

*Data* ([h4_windows.py](h4_windows.py) → [h4_windows_out.txt](h4_windows_out.txt)): (0, 50] 0.903 ± 0.033, (50, 100]
0.968 ± 0.042, (100, 150] 0.915 ± 0.044; (0, 75] 0.934 ± 0.026, (75, 150] 0.923 ± 0.037; (10, 150] 0.933 ± 0.021. A relative
term ∝ T^(−½) would make the (0, 50] deficit 1.73 times the full one (observed 1.37; the two readings share zeros, so the
ratio is imprecise): not seen, not excluded either.
The lowest zeros carry a bounded offset, +1.58 ± 0.13 per curve on (0, 10] (root number +1: +1.09, −1: +1.86); it moves the
slope by −0.014 only.

## H4. Missing or spurious zeros, windowing: rejected (computed)

* Zeros recomputed at realprecision 57 ([h4_zerocheck.py](h4_zerocheck.py) → [h4_zerocheck_out.txt](h4_zerocheck_out.txt))
  for the heaviest curve on (120, 150] (96 zeros) and three heavy low-reading curves on (0.01, 150] (G2_N2862_17: 398,
  G2_N2575_6: 396, T4b_N909_10: 370): identical as doubles; the sign changes of the Hardy function on a 0.02 grid equal the
  zero counts in all four.
* Window shape ([h4_smooth.py](h4_smooth.py) → [h4_smooth_out.txt](h4_smooth_out.txt)): sharp 0.929 ± 0.020, Hann
  0.935 ± 0.031, Tukey 0.1 / 0.25 / 0.5: 0.946 ± 0.024 / 0.939 ± 0.029 / 0.935 ± 0.031; five bands of 30: 0.940, 0.895,
  0.895, 1.003, 0.910 (± 0.048–0.078).
* Line shape at r = 1 ([h8_profile.py](h8_profile.py) → [h8_out.txt](h8_out.txt)): slope at u = ±π/T 0.943 ± 0.035 and
  0.932 ± 0.043, at ±3π/T 1.056 ± 0.119 and 0.889 ± 0.125: the amplitude drops, the shape stays.
* Outliers and errors ([h5_null.py](h5_null.py) → [h5_out.txt](h5_out.txt)): bootstrap sd 0.0195 against jackknife
  0.0198, 95 % interval [0.891, 0.968]; weighted median ratio 0.910; without the 5 heaviest curves 0.929 ± 0.025; residuals
  clipped at ±22: 0.928. No phase rotation (imaginary slope +0.014, [h7_lines.py](h7_lines.py) → [h7_out.txt](h7_out.txt)).
  Far twists r = e^(−u), u ∈ [6, 9], give a residual variance of 191 per curve against 120 at r = 1, so they are not a valid
  null for r = 1 (the zeros are rigid at the Gram frequency); the jackknife stands.

## H5. Atypical curves: rejected as the cause, one atypical curve found (computed)

Moments of a_p/√p ([h5_null.py](h5_null.py), own counts to 10⁵ for 51 curves, lfunan to 3·10⁴ for 49; USp(4): M2 = 1,
M4 = 3, M6 = 14): 99 curves have |z(M4)| ≤ 2.15 and |z(M2)| ≤ 2.06. **One curve is not USp(4):** R_N2187_28,
y² + (1 + x³)y = −x³ − x⁶, N = 3⁷, c = 5.43, **9.1 % of the slope weight**: M4 = 5.80 (z = 31.6), M6 = 52.4, M2 = 0.983.
End_Q is Z (a simple pole of L × L̃, so the L-function is primitive), but the Sato–Tate group is smaller than USp(4); I did
not identify it. The registered filter (mean |λ_p|² in [0.75, 1.25]) cannot see this. Its ratio G/pred is 0.900, the sample's
mean; without it the slope is 0.932 ± 0.022. Among root number +1, the weighted regression of G/pred on z(M4), z(M2), z(M6)
gives −0.0005 ± 0.013, +0.011 ± 0.033, −0.0001 ± 0.0027; the 8 lowest ratios are all generic. Mean e₂ over 30 < p ≤ 1000 is
0.933 ± 0.007 (USp(4): 1), as expected from the prime-square bias of L(s, Λ²₀) at x = 1000.

## Other readings (post hoc)

* Covariates, root number +1 ([h5_out.txt](h5_out.txt)): c > 3: 0.936 ± 0.022 vs 0.880 ± 0.039; N > 1500: 0.901 ± 0.028 vs
  0.933 ± 0.028; deg f = 5: 0.929 ± 0.027 vs 0.894 ± 0.022; two or more bad primes 0.923 ± 0.025 vs 0.902 ± 0.026. No
  covariate separates. All 37 curves with root number +1 have no central zero.
* Residual scatter falls with c: sd 12.0 (c ≤ 1), 11.2 (1–2.5), 8.2 (> 2.5) on (10, 150].
* Comparison class split (h9, [h9_out.txt](h9_out.txt)): Sym³ + products give b₁ = 0.784 ± 0.110, b₂ = 1.146 ± 0.080 (weak,
  opposite to genus 2); elliptic curves (round 148) 0.988 ± 0.021 and 0.991 ± 0.023.
* Truncated singular series c_R = Σ_(q≤R) μ(q)Π_(p|q)(1 − c_p) ([h10_trunc.py](h10_trunc.py) → [h10_out.txt](h10_out.txt),
  51 curves): (c_R − 1)/(c − 1) has median 0.848, 0.931, 0.960, 0.983 at R = 30, 100, 300, 1000; slopes on c_R: 1.034 ± 0.030,
  0.979 ± 0.023, 0.951 ± 0.021, 0.934 ± 0.020; per-curve residual sum of squares 5,454, 4,118, 3,968, 3,998 against 4,010
  for c. So R = 30 is rejected by the shape, R ≈ 100 fits the level and is not rejected by the shape (Δχ² ≈ 1). This has
  one free parameter and no derivation; the recipe above says there is no truncation.

## Surviving hypotheses

* **S1, chance.** 3.6 errors for the sharp reading (P(b ≥ 1) = 0.0005, bootstrap), 3.1 without the discovery sample,
  2.1–2.3 for smooth windows. Unlikely, not excluded.
* **S2, a real deficit of the r = 1 reading in degree 4**, about 7 %, constant over 50 ≤ T ≤ 150, in the part c − 1,
  absent from the lines r = 2, 3 (1.6 errors). Genus-2-only or all of degree 4: undecided (1.3 errors). No mechanism is
  derived; it contradicts the recipe's prediction of no lower-order term, so if real it is a property of the recipe's error
  at analytic conductor QT⁴/(2π)⁴ ≈ 10⁸–10⁹, or a failure of the conjecture in degree 4.
* **S2′ (post hoc), an effective truncation of the singular series at R ≈ 100**: makes the same level prediction as S2 on
  curves like these, and predicts that the deficit shrinks with T (a truncation that grows with T) and follows c − c_R
  curve by curve.

## Registration draft (not run)

Scripts in [blind_draft/](blind_draft/). Only the selection ran, for 115 s as a dry run of its code (4 accepted L-functions
in the first 22.7 % of the new box at N ≤ 3000; that is why the conductor range below is widened). The dry run's
candidate file was deleted. [g2b_test.py](blind_draft/g2b_test.py) `--dry` on the 100 old curves gives b = 0.9286 ± 0.0198,
Tukey 0.9464, band 0.9325: the readings it will make.

```
Round 15x (genus 2, blind): is the one-point law's r = 1 reading in genus 2 equal to 1 (the law) or 0.929 (the deficit)?
POPULATION. y^2 + h y = f, h one of the 16 0/1-polynomials of degree <= 3; f in {-2..2}^7 with deg f in {5, 6} and some
  |coefficient| = 2 (box B'), or f in {-1,0,1}^7 (the old box B) above its old conductor bound. Squarefree 4f + h^2;
  genus2red N_red < 10000; one per fingerprint; lfungenus2 conductor 1000 <= N_L <= 8000, lfuncheckfeq <= -20, one per
  (N_L, a_1..a_60); every L-function of rounds 1, 151, 152 and review 58 (tests, proxies, all valid box rows) excluded by
  that key; primitive (mean |lambda_p|^2, good p <= 3000, in [0.75, 1.25]); generic Sato-Tate (|z(M4)| <= 4, good
  100 < p <= 3e4); eps = lfunrootres.
STRATUM P (the test): eps = +1, generic, c >= 1.5 (pre-zero quantities; maximises power). Draw n = 72 by
  numpy default_rng(SEED).permutation of the rows sorted by (N, f, h); all rows if fewer.
CONSTANT. c by the H1 rule (own code): exact local factors p <= 1000, own point counts to 1e5, first order beyond 1000,
  log-Cesaro over [1e4, 1e5], nu = 2 (accuracy <= 0.3 %).
ZEROS. bz4_lib.zeros of round 152 unchanged (completeness rule unchanged), T2 = 150, root number set from lfunrootres
  before lfuninit.
READINGS. b = sum Re(G conj law)/sum law^2, law = -eps c T/2pi on (0, 150], jackknife se over COMPLETE L.
  D1 (the only decision): LAW iff |b - 1| <= 2.5 se and b - 0.929 >= 2.5 se; DEFICIT iff |b - 0.929| <= 2.5 se and
  1 - b >= 2.5 se; otherwise UNDECIDED.
  Reported, no decision: Tukey(0.1); band (10, 150]; the four lines r = 2, 1/2, 3, 1/3 pooled and their difference from
  r = 1; additive against proportional excess on (10, 150]; for S2', the slope on c_100 and the cumulative slopes at
  T = 50, 100, 150.
POWER. Residual sd at T = 150 for c > 1.5: 8-11 (this sample); with sd 10, se = 10/(23.87 sqrt(sum c^2)): 0.0175 needs
  sum c^2 >= 573 (72 curves of mean c^2 = 8). The two hypotheses are then 4.06 se apart. The expected se is computed from
  the drawn constants at the freeze and written into the registration; with fewer curves the rule stays and UNDECIDED is
  the likely outcome.
COST. About 110-180 s of zeros per curve (113 s mean in round 152 at N <= 3000), about 3-4 CPU-hours for 72.
FALSIFIER of the deficit claim: b >= 0.99. Falsifier of the law at this height: b <= 0.95 with se <= 0.02.
```

Commands (from `blind_draft/`; each job under the memory cap, repeated until done):

```
(ulimit -v 1500000; timeout 600 nice /home/bas/praime/.venv/bin/python -I g2b_select.py 540)   # until state i = total
/home/bas/praime/.venv/bin/python -I g2b_freeze.py SEED      # prints the stratum size, the expected se and the sha256
#   commit predictions_g2b.json + candidates_g2b.jsonl + this text; put the sha256 into g2b_test.py (REG) and commit
(ulimit -v 1500000; timeout 600 nice /home/bas/praime/.venv/bin/python -I g2b_zeros.py 540)    # until ALL DONE
(ulimit -v 1500000; timeout 600 nice /home/bas/praime/.venv/bin/python -I g2b_test.py)         # refuses another hash
```

## Not checked

* No exact cross-check of c through L-values: L(1, Sym²) (degree 10) and L(1, Λ²₀) (degree 5) of a genus-2 curve were not
  built (their bad factors and conductors are not at hand). The Euler product is the definition of c and converges to
  0.1–0.3 % (H1), so this cannot hide 7 %.
* Zeros of 96 of the 100 curves were not recomputed (only 4, the heaviest and three low-reading ones).
* The Sato–Tate group of R_N2187_28 is not identified.
* No theory of the size of the recipe's error term at these conductors; no zeros beyond T = 150 (the T-dependence that would
  separate S2 from S2′ needs new zeros).
* No web source was used.

Compute used: about 32 minutes (c to 10⁵ for 51 curves 7 min; zero checks 5 min; the rest seconds to a minute each).

## Files

All in `this folder: `:
[g2data.py](g2data.py) (own loader, θ, sums, slope), [pick.py](pick.py) → [sample.txt](sample.txt) (the 51 curves), [c_own.py](c_own.py) + [c_own_run1.txt](c_own_run1.txt),
[c_own_run2.txt](c_own_run2.txt), `cown/*.npz` (local data, 51 curves), [h1_analyze.py](h1_analyze.py) →
[h1_out.txt](h1_out.txt), [c_new.json](c_new.json), [c_table.md](c_table.md); [h2_local.py](h2_local.py),
[h2_badfactor.py](h2_badfactor.py); [h4_windows.py](h4_windows.py), [h4_smooth.py](h4_smooth.py),
[h4_zerocheck.py](h4_zerocheck.py); [h5_null.py](h5_null.py) (+ [h5_table.json](h5_table.json)); [h6_models.py](h6_models.py);
[h7_lines.py](h7_lines.py); [h8_profile.py](h8_profile.py); [h9_split.py](h9_split.py); [h10_trunc.py](h10_trunc.py);
[h11_extra.py](h11_extra.py); every `*_out.txt`; [blind_draft/](blind_draft/) (g2b_select.py, g2b_freeze.py, g2b_zeros.py,
g2b_test.py).

## Appendix: old and new c (51 curves; 93.8 % of the weight)

"c new" is the registered rule (log-Cesàro, ν = 2) applied at P = 10⁵ to own local data; "raw" is the plain partial product
to 10⁵. Slopes: 0.9286 ± 0.0198 with the old c, 0.9288 ± 0.0198 with the new.

| key | sample | N | eps | c old | c new (rule at 1e5) | raw product to 1e5 | new/old - 1 |
|---|---|---|---|---|---|---|---|
| R_N2187_28 | rep | 2187 | +1 | 5.4305 | 5.4337 | 5.4335 | +0.0006 |
| R_N2348_31 | rep | 2348 | +1 | 4.2246 | 4.2151 | 4.2041 | -0.0022 |
| G2_N2862_17 | 152 | 2862 | +1 | 4.1826 | 4.1799 | 4.1770 | -0.0007 |
| G2_N2613_21 | 152 | 2613 | +1 | 3.9644 | 3.9526 | 3.9403 | -0.0030 |
| G2_N1239_8 | 152 | 1239 | +1 | 3.9013 | 3.9027 | 3.8977 | +0.0003 |
| R_N1055_0 | rep | 1055 | +1 | 3.6896 | 3.6900 | 3.6892 | +0.0001 |
| G2_N2575_6 | 152 | 2575 | +1 | 3.5565 | 3.5577 | 3.5557 | +0.0003 |
| R_N1137_6 | rep | 1137 | +1 | 3.4525 | 3.4583 | 3.4553 | +0.0017 |
| G2_N1532_38 | 152 | 1532 | +1 | 3.4468 | 3.4493 | 3.4492 | +0.0007 |
| G2_N1473_39 | 152 | 1473 | +1 | 3.3359 | 3.3322 | 3.3299 | -0.0011 |
| R_N1164_7 | rep | 1164 | +1 | 3.2729 | 3.2748 | 3.2768 | +0.0006 |
| G2_N1811_37 | 152 | 1811 | +1 | 3.1701 | 3.1684 | 3.1635 | -0.0006 |
| T4b_N797_12 | 151 | 797 | +1 | 3.1037 | 3.1014 | 3.1048 | -0.0008 |
| T4b_N603_3 | 151 | 603 | +1 | 3.0402 | 3.0380 | 3.0298 | -0.0007 |
| R_N1109_3 | rep | 1109 | +1 | 2.9253 | 2.9214 | 2.9183 | -0.0013 |
| G2_N1051_14 | 152 | 1051 | +1 | 2.7697 | 2.7765 | 2.7811 | +0.0024 |
| R_N1695_16 | rep | 1695 | +1 | 2.7581 | 2.7581 | 2.7599 | +0.0000 |
| R_N1989_23 | rep | 1989 | +1 | 2.6076 | 2.6064 | 2.6067 | -0.0004 |
| T4b_N834_6 | 151 | 834 | +1 | 2.5801 | 2.5802 | 2.5789 | +0.0000 |
| R_N1269_10 | rep | 1269 | +1 | 2.5640 | 2.5611 | 2.5559 | -0.0012 |
| T4b_N909_10 | 151 | 909 | +1 | 2.5349 | 2.5347 | 2.5309 | -0.0001 |
| T4b_N997_1 | 151 | 997 | +1 | 2.5345 | 2.5314 | 2.5305 | -0.0012 |
| G2_N2261_33 | 152 | 2261 | +1 | 2.2930 | 2.2904 | 2.2891 | -0.0011 |
| G2_N1047_0 | 152 | 1047 | +1 | 2.2495 | 2.2489 | 2.2466 | -0.0003 |
| T4b_N713_0 | 151 | 713 | +1 | 2.1242 | 2.1234 | 2.1232 | -0.0004 |
| R_N1345_11 | rep | 1345 | +1 | 1.9971 | 2.0011 | 2.0015 | +0.0020 |
| T4b_N970_11 | 151 | 970 | +1 | 1.9883 | 1.9931 | 1.9968 | +0.0024 |
| G2_N2722_18 | 152 | 2722 | +1 | 1.9753 | 1.9719 | 1.9686 | -0.0017 |
| R_N2611_35 | rep | 2611 | +1 | 1.8514 | 1.8520 | 1.8532 | +0.0003 |
| T4b_N691_2 | 151 | 691 | +1 | 1.8022 | 1.8056 | 1.8070 | +0.0019 |
| G2_N1109_4 | 152 | 1109 | +1 | 1.7127 | 1.7137 | 1.7138 | +0.0006 |
| R_N2949_42 | rep | 2949 | -1 | 1.5924 | 1.5896 | 1.5892 | -0.0018 |
| T4b_N523_9 | 151 | 523 | +1 | 1.4581 | 1.4599 | 1.4592 | +0.0013 |
| T4b_N353_4 | 151 | 353 | +1 | 1.3409 | 1.3399 | 1.3404 | -0.0007 |
| R_N1123_4 | rep | 1123 | +1 | 1.2225 | 1.2251 | 1.2263 | +0.0022 |
| R_N2511_33 | rep | 2511 | -1 | 1.1041 | 1.1038 | 1.1034 | -0.0002 |
| G2_N2139_5 | 152 | 2139 | -1 | 0.7697 | 0.7714 | 0.7726 | +0.0021 |
| G2_N2697_3 | 152 | 2697 | +1 | 0.7558 | 0.7571 | 0.7560 | +0.0017 |
| T4b_N349_8 | 151 | 349 | +1 | 0.5737 | 0.5730 | 0.5717 | -0.0011 |
| R_N2101_27 | rep | 2101 | -1 | 0.5191 | 0.5200 | 0.5192 | +0.0018 |
| G2_N2031_28 | 152 | 2031 | -1 | 0.4850 | 0.4847 | 0.4849 | -0.0006 |
| G2_N1929_32 | 152 | 1929 | -1 | 0.4473 | 0.4475 | 0.4482 | +0.0005 |
| G2_N1069_26 | 152 | 1069 | +1 | 0.4452 | 0.4450 | 0.4447 | -0.0005 |
| R_N1961_22 | rep | 1961 | -1 | 0.3922 | 0.3922 | 0.3917 | -0.0001 |
| R_N1383_12 | rep | 1383 | -1 | 0.2163 | 0.2164 | 0.2164 | +0.0004 |
| R_N1198_8 | rep | 1198 | -1 | 0.1724 | 0.1724 | 0.1724 | -0.0001 |
| G2_N1051_9 | 152 | 1051 | -1 | 0.1273 | 0.1273 | 0.1270 | +0.0004 |
| R_N1403_13 | rep | 1403 | -1 | 0.1062 | 0.1063 | 0.1064 | +0.0009 |
| R_N1127_5 | rep | 1127 | -1 | -0.5188 | -0.5190 | -0.5190 | +0.0005 |
| G2_N1253_23 | 152 | 1253 | -1 | -0.5352 | -0.5356 | -0.5364 | +0.0006 |
| R_N2813_41 | rep | 2813 | -1 | -0.9357 | -0.9338 | -0.9314 | -0.0020 |
