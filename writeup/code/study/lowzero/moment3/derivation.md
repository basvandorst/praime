# Mean square of a product of three different Dirichlet L-functions: the recipe, its checks, and a registered test

Scratch directory: `writeup/code/study/lowzero/moment3/`
(all files named below are in it). Status labels: **proved** (a theorem in print, cited), **derived** (worked out
here, elementary), **checked** (agreement verified numerically by a separate script), **conjecture**, **computed**.

## 1. The question

For real primitive characters χ_j = (D_j/·) (D_j fundamental discriminants, q_j = |D_j|, distinct), put

  M_k(T1, T2) = ∫_{T1}^{T2} |L(1/2+it, χ_1) ··· L(1/2+it, χ_k)|² dt,   k = 1, 2, 3.

* k = 1: the second moment; theorem with all lower-order terms (classical; the formula is (2.9) below).
* k = 2: Topacogullari, "The fourth moment of individual Dirichlet L-functions on the critical line",
  Math. Z. 298 (2021), arXiv 1909.11517, **Theorem 1.3** (two different primitive characters):
  M_2(1, T) = ∫_1^T P_{χ1,χ2}(log t) dt + O((q1* q1 + q2* q2)^{1/2} (q1q2)^{3/4−3θ/2} T^{1/2+θ+ε} + (…)^{1/3}(q1q2)^{1/3} T^{2/3+ε}),
  P quadratic, explicit residue (his eq. (5576) for q1 ≠ q2; (5583) with an extra constant for q1 = q2, which vanishes
  for opposite parity, the only way two real primitive characters share a conductor: D and −D). Theorem 1.5 is the
  case ζ·L(χ_D) = ζ_K (improving Motohashi 1970, Müller 1989). **Proved, but only as T → ∞ with the q fixed:** at our
  sizes (q ≈ 10²–10³, T ≤ 1000) the error bound is far larger than the main term, so at the test's heights the pair
  formula is a theorem's main term with an unproved error size.
* k = 3: no asymptotic is known (degree 3, like the sixth moment of ζ). Known: lower bound ≫ T(log T)³ (unconditional by the
  Rudnick–Soundararajan short-polynomial argument, see §4); upper bound ≪ T(log T)³ under GRH
  (Hagen, "Sharp conditional moment bounds for products of L-functions", arXiv 2409.19780, GRH + GRC; earlier
  T(log T)^{3+ε}, Milinovich–Turnage-Butterbaugh 2014). So the ORDER is theorem-grade under GRH; the constant and the
  lower-order polynomial are **conjecture**: Heap, arXiv 1303.6119 (only version v1, 2013; Math. Proc. Cambridge
  Philos. Soc. 170 (2021) 191–219), **Conjecture 4** (numbering of v1): for L = Π L_j^{e_j} with distinct primitive L_j,
  (1/T)∫_0^T |L(1/2+it)|^{2k} dt ~ a_L(k) Π_j G²(e_j k+1)/G(2 e_j k+1) (log(Q_j T^{d_j}))^{(e_j k)²}, a_L(k) = Π_p (1−1/p)^{n_L k²} Σ_n |α_{L,k}(p^n)|²/p^n.
  For three distinct Dirichlet L (e_j = d_j = k = 1, G(2) = G(3) = 1): M_3 ~ a_L(1) T Π_j log(q_j T/π), and a_L(1)
  contains Π_{i<j} L(1, χ_iχ_j)² (section 2.3). Heap derives the leading term from a modified CFKRS recipe
  (§6 and §7 of his paper: the "first step" of the recipe applied to each L_j separately) and from the hybrid model.
  I could not see the published 2021 text; arXiv has only v1.

## 2. The recipe for k distinct real characters (derived; checked in §3)

### 2.1 Shifted moment and the swap rule
I(α, β) = ∫ Π_j L(1/2+α_j+it, χ_j) L(1/2+β_j−it, χ_j) dt (real χ: L(1/2−it) = conj L(1/2+it)).
Approximate functional equation per factor: L(s, χ) = Σ χ(m) m^{−s} + κ_χ(s) Σ χ(n) n^{s−1}, κ_χ(1/2+it) = ε_χ (q t/2π)^{−it} e^{it+iπ/4}(1+O(1/t)),
ε_χ = 1 for real χ. Heap's modification (his §6–7): keep only the terms with the same number of κ from L_j(s) and from
L_j(1−s) **for each j separately** (cross terms carry (q_i/q_j)^{±it}; their diagonal is empty, because
n m' q_j = m n' q_i with χ_i(n m') χ_j(m n') ≠ 0 is impossible for q_i ≠ q_j: a prime of q_i/(q_i,q_j) or of q_j/(q_i,q_j)
would have to divide a character argument. derived). So with k = 1 per L_j there are 2^k terms, indexed by the set S of
swapped L-functions: for j ∈ S, (α_j, β_j) → (−β_j, −α_j) and a factor X_j^{−α_j−β_j}, X_j = q_j t/2π. Then keep the
diagonal of the whole product, m_1···m_k = n_1···n_k, extended to all integers:

  P(t; α, β) = Σ_{S ⊆ {1..k}} Π_{j∈S} X_j^{−α_j−β_j} R(α^S, β^S),
  R(α, β) = Σ_{m_1···m_k = n_1···n_k} Π_j χ_j(m_j n_j) m_j^{−1/2−α_j} n_j^{−1/2−β_j}.        (2.1)

Conjecture (CFKRS 2005 §2; Heap): M_k = ∫ P(t; 0, 0) dt + O(T^{1/2+ε}).

### 2.2 The Euler product
R = Π_p R_p, R_p = Σ_{n≥0} h_n(x) h_n(y), x_j = χ_j(p) p^{−1/2−α_j}, y_j = χ_j(p) p^{−1/2−β_j}, h_n the complete
homogeneous symmetric polynomial (Σ_n h_n(x) u^n = Π_j (1 − x_j u)^{−1}; the condition Σe = Σf is the constant term in u).
Closed form (residues of the u-integral, distinct y): R_p = Σ_k y_k^{k−1} / (Π_i (1 − x_i y_k) Π_{j≠k} (y_k − y_j)).
(derived; checked against the series to 1e-15, check_recipe.py.) Code uses the recurrence h_n = e1 h_{n−1} − e2 h_{n−2} + e3 h_{n−3}.

### 2.3 Polar part and the arithmetic factor
R(α, β) = Π_{i,j} Λ_ij(1+α_i+β_j) · A(α, β), Λ_ij(s) = Σ_n χ_i(n)χ_j(n) n^{−s},
A = Π_p A_p, A_p = R_p · Π_{i,j} (1 − x_i y_j).
* i = j: Λ_jj(s) = ζ(s) c_j(s−1), c_j(x) = Π_{p|q_j}(1 − p^{−1−x}).
* i ≠ j: Λ_ij(s) = L(s, χ_{D*}) Π_{p | D_iD_j}(1 − χ_{D*}(p) p^{−s}), D* = coredisc(D_i D_j) (the primitive character
  inducing χ_iχ_j); imprimitive whenever (q_i, q_j) > 1 or the 2-parts differ. Local factors are exact (no prime is dropped).
* log A_p = −e2(x) e2(y) + O(p^{−3}) (derived: the degree-4 part of log R_p + Σ log(1 − x_i y_j) is
  h2(x)h2(y) − p1(x)²p1(y)²/2 − p2(x)p2(y)/2 = −e2(x)e2(y)), so A converges absolutely for Re shifts > −1/4.
* k = 2: A_p = 1 − x1x2y1y2 exactly (Ramanujan's identity), i.e. A = Π_{p∤q1q2}(1 − p^{−2−α1−α2−β1−β2}). (checked to 6e-14.)
* At zero shift (real characters, α = β = 0):
  h(0) = a_L(1) = Π_j φ(q_j)/q_j · Π_{i<j} Λ_ij(1)² · A(0)
  — Heap's a_L(1) (the residue of Σ |α_L(n)|² n^{−s}), with the factor Π_{i<j} |L(1, χ_iχ_j)|² of the question.

**The three-body factor.** For k = 3, A_3(0) ≠ Π_{i<j} A_ij(0). The ratio is a local effect at small primes. Closed form at
p = 2 when χ_1(2) = χ_2(2) = χ_3(2) = ±1 (x_j = ±2^{−1/2}): R_2 = Σ C(n+2,2)² 2^{−n} = (1+4y+y²)/(1−y)^5 at y = 1/2 = 104,
A_{3,2} = 104/2^9 = 0.2031, Π A_{ij,2} = (3/4)³ = 0.4219, ratio 0.4815 (derived; sizes.py). Over the 56 development
triples A_3(0)/Π A_ij(0) ranges 0.284 … 1.701 (median 0.936). This factor is what separates the recipe from the
"pairwise" hypothesis (§5).

### 2.4 Zero shift: the constant-term formula (derived; checked)
Take α_j = β_j = z_j (any path to 0 gives the same limit: P(t; α, β) is holomorphic). Write
R = Π_j ζ(1+2z_j) · h(z), h(z) = Π_j c_j(2z_j) Π_{i<j} Λ_ij(1+z_i+z_j)² A(z) (on this slice α_i+β_j = z_i+z_j is symmetric).
The swap S sends z_j → −z_j (j ∈ S), so with X^{−(1−s)z} = X^{−z} X^{sz}:
  P(t; z) = Π_j X_j^{−z_j} Σ_{s ∈ {±1}^k} Π_j [ζ(1+2 s_j z_j) X_j^{s_j z_j}] h(s z).
The function Φ(z) = Π_j ζ(1+2z_j) X_j^{z_j} h(z) has only simple poles in each z_j; the sum over sign flips picks its
even part, and as z → 0 only the constant term survives (×2^k):

  P(t) = 2^k CT_z Φ = Σ_{m ⊆ {1..k}} h[m] Π_{j∉m} (L_j + 2γ),   L_j = ln(q_j t/2π),        (2.2)

with h[m] = coefficient of Π_{j∈m} z_j in h (2ζ(1+2z)X^z = 1/z + (L + 2γ) + O(z)). Only the multilinear coefficients
are needed, so h is computed exactly in the algebra R[z_1..z_k]/(z_1², …, z_k²) (tmlib.py: jets; log Λ_ij from PARI
lfun derivatives of order ≤ 3 at s = 1; log A from the per-prime series to p ≤ 10^5 plus the principal quartic tail
−Σ_{p>10^5} Σ_{i<j} p^{−2−2(z_i+z_j)} by the prime number theorem; 10^5 → 10^6 changes P by < 1e-8).
For k = 1 (2.2) is (φ(q)/q)(L + 2γ + 2Σ_{p|q} ln p/(p−1)), the classical second moment (2.9) (checked exactly).
P is a polynomial of degree k in ln t.

### 2.5 Comparators (derived)
* lead: h(0) Π_j L_j (Heap's leading term with ln(q_j t/2π)).
* diag (the AFE diagonal; the only finite "truncation" of the swap sum, see §5): the product L_1L_2L_3 has an
  approximate functional equation with two sums of length √C, C = Π_j X_j; the mean square of their diagonals,
  with a sharp cutoff, is D(t) = 2 Res_{s=0} ζ(1+s)^k g(s) C^{s/2}/s, g(s) = h(s/2, …, s/2) (Taylor to order k by a
  univariate jet; Stieltjes constants γ_0, γ_1, γ_2). Leading term h(0)(Σ L_j)³/24 against the recipe's h(0) Π L_j:
  ratio ≥ 27/24 = 1.125 at leading order (AM–GM; equality for equal conductors). For k = 1, D = P exactly; for k = 2
  the leading terms differ by h(0)(L_1−L_2)²/4, which vanishes for equal conductors.
  The sharp cutoff is a choice: a smooth AFE weight G(s) = 1 + s²/δ² + … moves D by ≈ 0.02/δ² relative (derived);
  so H_diag is defined only to ~2 %, against a 10 % separation.
* ind: Π_j P_1(χ_j; t) (three independent L-functions: no coupling factor).
* kirk ("pairwise superposition"): Π_{i<j} P_2(χ_i, χ_j; t) / Π_j P_1(χ_j; t): the triple built from the theorem-grade
  pair formulas; it has the factor Π L(1, χ_iχ_j)² but A_3 → Π A_ij and its own lower-order terms.
* no3: the recipe (2.2) with A_3(z) replaced by Π_{i<j} A_ij(z_i+z_j): isolates the three-body Euler factor.

## 3. Checks of the machinery (all checked; scripts and outputs in the scratch directory)

| check | script | result |
|---|---|---|
| per-prime series of A_p for pairs vs closed form 1 − x1²x2² (jets, p ≤ 10^5) | check_lib1.py | max diff 7e-14 |
| k = 1: (2.2) vs the classical formula | check_lib1.py | identical; D(t) = P(t) |
| k = 2: (2.2) vs Topacogullari's theorem polynomial (5576) (contour integral, mpmath + PARI), 9 pairs incl. ζ·L(χ_D) for D = 1001, −1003, 1068 | check_topa.py | relative difference ≤ 3e-12 |
| k = 2 leading coefficient vs his closed form (6/π²)|L(1,χ1χ2)|² φ(q1)φ(q2)/φ(q1q2) Π_{p|q1q2}(1 − 1/(p+1)) | by hand (§2.3) | identical |
| k = 3: (2.2) vs the literal sum of the 8 swap terms at shifts a = ε(1, 1.31, 1.73), every factor evaluated directly (mpmath ζ, PARI L at real s, direct Euler product), ε = 1e-3 … 1.25e-4, Richardson | check_recipe.py | 3e-7 … 3e-6 (the extrapolation error) for 3 triples × t = 100, 1000 |
| R_p series vs residue closed form (random distinct x) | check_recipe.py | 1e-15 |
| fourth moment of one L: Topacogullari's (5508) polynomial, leading coefficient vs (1/2π²)(φ(q)/q)² Π(1 − 2/(p+1)) | dev_controls.py | identical |
| values: PARI lfun vs an independent Hurwitz-zeta sum (mpmath) | t_mp.py | 1e-16 |
| values: realprecision 19 vs 38 digits, 18001 points | t_prec.py | 2e-16 relative |
| quadrature: trapezoid h = 0.05 vs h = 0.1 | dev_analysis.py | 5e-5 (so ≈ 1e-5 at h = 0.05); tapered weight 2e-11 |
| Euler product truncation 10^4 / 10^5 / 10^6 | dev_pred.py | P changes 8e-8 / 3e-9 |

## 4. What is theorem and what is conjecture (the project's filter)

* k = 1 (all terms), k = 2 (the full quadratic polynomial, including the mixed swaps S = {1}, {2}): **theorem as T → ∞**
  (classical; Topacogullari Thm 1.3, 1.5). At T ≤ 1000 with q ≈ 10³ only the polynomial is theorem-grade; the error is
  not controlled by the theorem. Readings on pairs and ζ·L are **controls**: they test nothing new.
* k = 3, order of magnitude T(log T)³: lower bound unconditional by the standard short-polynomial argument
  (Rudnick–Soundararajan; Radziwiłł–Soundararajan, "Continuous lower bounds for moments of zeta and L-functions"):
  ∫|L|² ≥ |∫ L Ā_N|²/∫|A_N|², A_N the first N = T^θ terms of L_1L_2L_3 (θ > 0 small enough that the twisted first
  moment is evaluable), gives M_3 ≥ (1+o(1)) T Σ_{n≤N} |α(n)|²/n ~ (θ³/6) h(0) T (ln T)³ (derived as a sketch; I did
  not find this exact statement for three distinct Dirichlet L-functions in print). Note: this lower bound already
  carries the full arithmetic factor h(0) = a_L(1), including A_3. Upper bound of the same order under GRH (Hagen 2024).
  A reading that only says "M_3 ≍ T(log T)³" or "M_3 grows with h(0)" over the 1000-fold range of h(0) tests (almost)
  nothing.
* k = 3, conjectural content:
  (a) the leading constant: M_3/(T Π L_j) → h(0) exactly (Heap Conj. 4), i.e. the universal factor is 1 (g_L = 1 for
      three distinct degree-1 L-functions, Π G(2)²/G(3) = 1); the provable lower bound only gives a factor θ³/6 in
      front of h(0)(ln T)³, and no upper bound with an explicit dependence on the characters is known;
  (b) the mixed swaps S with |S| = 1, 2 (the off-diagonal of the product AFE): they make the leading term h(0)ΠL_j
      instead of the AFE-diagonal h(0)(ΣL_j)³/24; at the test's heights the off-diagonal is −10.4 % of P
      (median; −8.0 … −15.4 % over the dev triples). For pairs the same terms are theorem; for ζ·L they are −8 … −26 %;
  (c) the three-body factor A_3/Π A_ij (0.28 … 1.70): conjectural; it is shared by the diagonal sum Σ|α(n)|²/n and by
      the provable lower bound, but no theorem fixes the moment's dependence on it beyond the slack of the lower bound (a factor 6/θ³ ≫ 1);
  (d) all lower-order terms (at T ≤ 1000 the leading term is only 35–160 % of P, median 75 %; sizes_out.txt).
* Individual swap classes are infinite at zero shift (each term ~ h(0) Π_j (±1/(2z_j))); only the full sum and the
  sum {∅, {1,2,3}} regularised by the AFE cutoff (= D) are finite. A "one-swap truncation" (S with |S| ≤ 1) is not
  symmetric under the functional equation and diverges as the shifts → 0: it is not a meaningful alternative at
  zero shift. The meaningful alternative is H_diag (the diagonal-only / Conrey–Gonek-type moment), which the data can
  reject (§6). Sizes "per swap class" are therefore given as (i) diagonal vs off-diagonal (above) and (ii) by the
  multilinear order |m| in (2.2) (median shares at t = 316: |m| = 0: 101 %, 1: +2 % (−184 … +47 %), 2: +4 %, 3: +0.3 %;
  sizes_out.txt).

## 5. Development set (8 characters, 1000 ≤ |D| ≤ 1300; now used)

D = 1001, −1003, 1005, −1048, 1068, −1195, 1253, −1267 (chosen for coverage: both parities; shared primes 7 (three
characters), 3, 5, 2, so several χ_iχ_j are imprimitive; 2 | D twice). Values: moment_values.py → dev_values/
(plus ζ, D = 1). Analysis: dev_analysis.py (dev_analysis_out.txt), dev_stats.py, dev_controls.py, dev_noise.py,
dev_noise2.py; dry run of the registered readings: moment_test.py on predictions_dev.json (dryrun_dev_out.txt).

Cost: lfuninit + 18,001 lfun values on [100, 1000]: 35–85 s per character on the busy machine (load 6–7).

Results (measured/recipe over [100, 1000]):
* singles (theorem): 0.979 … 1.016, mean 1.0013; fourth moments (theorem): 0.853 … 1.085, mean 0.976.
* ζ·L(χ) (theorem incl. mixed swaps): 0.968 … 1.005, mean 0.989; the AFE diagonal would read 0.80–0.92: the data pick
  the theorem's mixed-swap terms.
* pairs (28): mean 0.998 ± 0.030 (jackknife over characters), sd per pair 0.072.
* triples (56): mean 0.9945 ± 0.118, sd per triple 0.240; slope of ln M on ln P_full 0.986 ± 0.078 over a range of
  e^7 in P. Comparators (mean of M/P_c, jackknife): diag 0.900 ± 0.107, kirk 0.903 ± 0.210, no3 0.908 ± 0.202,
  ind 10.1 ± 8.7, lead 1.56 ± 0.51; sd of ln(M/P_c): full 0.242, diag 0.242, no3 0.393, kirk 0.408, lead 0.524, ind 1.48.
  Amplitudes γ_c (1 = recipe, 0 = comparator): no3 1.006 ± 0.241, kirk 1.018 ± 0.234, ind 0.972 ± 0.064, lead 1.05 ± 0.13.

Noise structure (the |L| values are exact to 1e-15; the "noise" is the recipe's error term on a window of length 900):
* per triple sd of ln(M/P) 0.24; split into character effects (σ_a ≈ 0.10) and the rest (σ_e ≈ 0.19).
  Character effects dominate the pooled error: Var(mean) ≈ 9σ_a²/n_char + σ_e²/N.
* the triple's deviation follows the pair deviations: corr(ln r_ijk, Σ ln r_pairs) = 0.81 (slope 1.47). Control
  variate with β = 1 (fixed, parameter-free: r3 − Σ_pairs (r_ij − 1)): σ_a drops to ≈ 0.022, σ_e to 0.15; on dev the
  CV mean is 1.0006 ± 0.045 (8 characters only).
* windows of length 100: sd of ln ratio 0.65 (triples), adjacent windows anti-correlated (−0.35): the error term has
  the rigidity expected of an oscillatory remainder; it is not white noise.
* heavy tails and correlation time (dev_corrtime.py): for the triple integrand f = |L_1L_2L_3|², Var f/(E f)² ≈ 199,
  max/mean ≈ 630, kurtosis ≈ 1100 (medians); correlation time (integral of the autocorrelation to its first zero)
  0.18. The "random process" estimate sqrt(Var f · 2τ · 900)/(900 E f) gives 0.28 (0.14 … 0.54) relative error per
  triple — close to the actual scatter 0.24. For one L-function the same estimate gives 0.087 against an actual 0.011,
  for pairs 0.16 against 0.072: lower degree means stronger rigidity of the remainder (the autocorrelation has negative
  tails); the triple product is nearly an uncorrelated heavy-tailed process on this window.
* tapered weight (sin² ramps of 100 at both ends) does not reduce the noise (sd 0.264 vs 0.240): rejected.

## 6. Reserved set and power

Frozen: sets.json (sha256 b23a8b41…), predictions_triple.json (sha256
b674db97f04aeadf2ee19cc62dfa1ddd918a1ca0e54ea2e46e42191af6cdd5a6, also in predictions_triple.sha256).
R = 56 fundamental discriminants, 100 ≤ |D| ≤ 1300, distinct |D|, none in the development set (seeded draw from 724
candidates, freeze.py): −1288 −1243 −1220 −1191 −1187 −1119 −1111 −1099 −1043 −1039 −1011 −983 −820 −807 −795 −763 −687
−635 −552 −532 −419 −415 −407 −383 −379 −339 −296 −287 −148 −136 −115 120 181 248 253 265 305 317 341 412 456 461 472 601
664 757 861 885 904 917 952 993 1013 1048 1149 1173. All 1540 pairs; 2000 of the 27720 triples (seeded; 80–128 per
character).
Over the test triples: P_diag/P_full 1.072 … 1.271 (mean 1.113); A_3/Π A_ij 0.284 … 2.450 (median 1.096);
P_lead/P_full 0.25 … 1.96 (mean 0.64); P_ind/P_full 0.011 … 14.6; ∫P_full 3.2·10^4 … 1.3·10^8 (h(0) 0.0086 … 187).
Pairs: P_diag/P_full 0.988 ± 0.007 (too close to read); ζ·L: 1.084 … 1.347 (mean 1.178).

Power (power.py → power_out.txt; noise model from the dev set only; 1000 simulations per hypothesis):
control-variate level θ: se 0.012 (0.024 if σ_a is doubled): H_diag (1.113) at 9 se (4.6); raw mean se 0.042 (H_diag
2.7 se — the raw reading alone cannot decide). γ_no3 and γ_kirk (CV): se 0.014 (0.020): H_no3/H_kirk at 70 se; raw 31 se.
γ_ind se 0.003, γ_lead se 0.013: > 40 se. Slope (CV) 1.000 ± 0.003 under H_full, 1.12 under H_kirk/H_no3, 0.15 under
H_ind. The delete-one-character jackknife reproduces the level se (0.0136 vs 0.0142 on 40 simulated sets) and
overstates the γ se (0.023 vs 0.012): conservative. Under H_diag, γ_no3 = 1.004: the three-body reading does not
separate the recipe from the diagonal (both carry a_L(1)); only the level does.

What varies between triples and is it trivial? The predicted moment varies by a factor ≈ 4000 between test triples,
almost all of it through h(0) = a_L(1) (the L(1, χ_iχ_j)² and the local factors). That variation is shared by every
"diagonal-type" model (recipe, AFE diagonal, the provable lower bound), so the slope of ln M on ln P across triples is
theorem-grade in bulk and gets no decision. The non-trivial between-triple contrast is the three-body factor (a factor
0.28 … 2.45, from the alignment of χ_j(2), χ_j(3), …), which separates the recipe (and the diagonal) from every pairwise
model; and the level, which separates the recipe from the diagonal.

## 7. Prior art (verdict NOT FOUND for a numerical test of the k = 3 conjecture)

* Heap 2013 (arXiv 1303.6119; MPCPS 2021): conjecture, no numerics. Heap's thesis (York 2013): second moment of ζ_K
  (quadratic), theory.
* CFKRS, "Integral moments of L-functions" (Proc. LMS 2005; math/0206018): the recipe; numerics for ζ (sixth moment
  on blocks of length 50,000 up to 2.35·10^6: ratios 0.995–1.004, "consistent with a remainder of size |D−C|^{1/2}D^ε")
  and for families in the q-aspect; no products of distinct L-functions.
* Conrey–Iwaniec–Soundararajan, "The sixth power moment of Dirichlet L-functions" (GAFA 2012; arXiv 0710.5176): the
  sixth moment averaged over all primitive χ mod q, over q ≤ Q and over a t-interval (asymptotic large sieve),
  confirming the recipe's constant 42; relevance: a theorem for a degree-3 average with a triple average, not for a
  fixed product of three distinct characters. Chandee–Li (arXiv 1303.4482, eighth moment, GRH) and
  Chandee–Li–Matomäki–Radziwiłł (arXiv 2307.13194, unconditional): same family-average setting.
* Topacogullari (Math. Z. 2021, arXiv 1909.11517): k = 2 theorem (our control); no numerics.
* Numerics of single-L moments: CFKRS 2005 (ζ, k = 3, 4), Rubinstein–Yamagishi ("Computing the moment polynomials of the
  zeta function", Math. Comp. 2015), Hiary–Odlyzko (Math. Comp. 2012: moments of ζ near height 10^23 against random-matrix models; not re-checked here); quadratic-twist families (Rubinstein,
  CFKRS). Hagen 2024 and Milinovich–Turnage-Butterbaugh: bounds for products, no numerics.
* Searched (six web searches and the arXiv sources of Heap, Topacogullari and CFKRS, November 2026): "numerical … mean square product of distinct Dirichlet L-functions",
  "moments of products of L-functions numerics", "Dedekind zeta mean square numerical": nothing for k ≥ 2 distinct
  L-functions, not even for ζ_K = ζ·L(χ) (k = 2). Verdict: NOT FOUND (the k = 3 conjecture has no numerics in print that
  I could find; the k = 2 theorem has none either).

## 8. Doubts

* The published (2021) version of Heap's paper may number or state the conjecture differently; only arXiv v1 was read.
* The "Heap modification" (no cross swaps between different L_j) is part of the conjecture; it is confirmed for k = 2 by
  Topacogullari's theorem (no extra term for q1 ≠ q2, and none for real characters with q1 = q2).
* The noise model is from 8 characters: σ_a is uncertain by a factor ~2; the power numbers carry that uncertainty.
* The control variate assumes the pair error terms have mean ~0 at these heights (dev: 0.998 ± 0.030); a finite-T bias
  of the pairs would pass into the CV reading (the raw reading is reported alongside).
* The window [100, 1000] is low: the recipe's error O(T^{1/2+ε}) is ~20 % per triple; the test is statistical over
  ~2000 triples and 56 characters, not a per-triple check.
