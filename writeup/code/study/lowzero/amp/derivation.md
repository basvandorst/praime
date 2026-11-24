# The 7 % shortfall of the off-diagonal amplitude in test B: five candidates, derived and tested

Folder: `writeup/code/study/lowzero/amp/` (all files below are
there; nothing in the repo was changed). Scripts and their outputs:

| script | output | what it does |
|---|---|---|
| `gridcache.py` | `grid.npz` | builds the repo's frozen `bklib.Grid` (copy `bklib.py`) once and caches it |
| `common.py` | | predictions per L-function (as `predict.py`), Σ² exactly as `r134.py`, γ as `r134.py` |
| `theory.py` | `theory_out.txt` | candidates (a), (b), (d): ratios conjecture against BK, exact phase, drift, unfolding |
| `retro.py` | `retro_out.txt` | γ in test B by split; γ in the rounds 117/120 sample at the same heights |
| `influence.py` | `influence_out.txt` | per-L-function contributions to γ − 1, zero-set completeness, γ − pooled null |
| `decomp.py` | `decomp_out.txt` | between-group (indicator) and within-group parts of the off pattern |
| `alt.py` | `alt_out.txt`, `alt_*.npz` | BK, AMP (0.932) and LP3 (L + 0.3) predictions for five samples |
| `table.py` | `table_out.txt` | the retrospective table and the joint reading |
| `chiexp.py` | `chiexp_out.txt` | exploratory: dependence on χ_D(2), χ_D(3) |
| `incomplete.py` | `incomplete_out.txt` | the one incomplete zero set of rounds 117/120 |
| `register.py`, `regfinal.py` | `register_out.txt`, `regfinal_out.txt`, `register.json`, `predictions_alt.json` | numbers for `round_draft.txt` |

Labels: **proved** (algebra or calculus from stated formulas), **classical** (in the literature), **computed**,
**retrospective** (data already known), **exploratory** (a search after looking at the data), **heuristic**.

Abbreviations: BK = the Bogomolny–Keating law as frozen in `lowzero/bk/`; AMP = BK with the bad-prime off-diagonal change
multiplied by 0.932; LP3 = BK with the bad-prime off-diagonal change evaluated at L + 0.3 (the diagonal at L; the review's
construction). R120 = the 181 real L-functions with 100 ≤ |D| ≤ 400 of rounds 117/120 (zeros `r117_part*.npz`,
`r120_part*.npz`), the sample of round 123 and the discriminants of the registered test C.

## Summary

- No candidate explains a 7 % shortfall. (a), (b) and (d) give corrections of at most 1·10⁻⁴ in γ, and (c) gives exactly
  zero. Each is killed by its own number and, independently, by the data (next point).
- The shortfall is **not a property of the off-diagonal term**. On R120 (the same heights 100–600, the same statistic,
  the same code; retrospective) γ = 1.001 ± 0.019. Inside test B itself the part of the per-L-function off pattern that
  the indicators do not carry has slope 1.042 ± 0.028. Only B's indicator contrasts are low (0.938 ± 0.021), and their
  difference from B's within part is z = −2.9.
- Joint reading of the four numbers (B and R120, between and within): a common amplitude of 0.983 ± 0.011. Against
  them, χ² on 4 degrees of freedom is 11.7 for BK (p = 0.019), 32.5 for AMP (p = 1.5·10⁻⁶) and 34.0 for LP3
  (p = 7.5·10⁻⁷). The single-amplitude model has χ² 9.3 on 3 (p = 0.026).
- The economical reading: BK is complete, and B's registered γ is a low excursion (p = 0.003 for that one reading). That
  is not proved, and a conductor-dependent effect is not excluded (descriptive trend −0.058 ± 0.035 per unit ln|D|).
- **No corrected formula.** The parameter-free prediction is BK's own: γ = between = within = pooled = 1. AMP and LP3
  predict about 0.93 for every reading in every sample. `round_draft.txt` registers tests that separate these.

## 0. What γ measures

γ = Σ res·off / Σ off² over the 96 cells (prime p ∈ {2, …, 13} × ℓ ∈ {0.5, …, 8}). Here off is the regression-matched
predicted bad-prime off-diagonal coefficient and res = fit − (predicted diagonal coefficient). **Computed**
(`retro_out.txt`): p = 2 carries 69 % of Σ off², p = 3 19 %, and p = 2 at ℓ > 5 alone carries 58 %. So γ is mainly the
amplitude of the p = 2 period-1 swing at ℓ = 5.5–8. There ε = 2πℓ/L ≈ 3.4–4.9, ε ln 2 ≈ π, and |F₂ − 1| = 4 sin²(ε ln 2/2)
is near its maximum 4.

γ − 1 is linear in the residuals S_i − full_i. In B it is spread over the L-functions: the 10 most negative contributions
sum to −0.039, all the others to −0.029, and the expected sd of the sum is 0.023. Dropping the 5 or 10 most extreme
contributors on both sides gives γ = 0.920 or 0.922 (`influence_out.txt`). So no outlier drives it.

## 1. Candidate (a): lower-order terms of the two-point function at finite height

### 1.1 The ratios conjecture for one Dirichlet L-function (classical recipe, written out here for χ mod k)

Let χ be primitive mod k with parity a (χ(−1) = (−1)^a). The functional equation is L(s, χ) = ε_χ X_χ(s) L(1 − s, χ̄),
with X_χ(s) = (π/k)^(s−1/2) Γ((1 − s + a)/2)/Γ((s + a)/2) and ε_χ ε_χ̄ = 1. Put s = 1/2 + it. The recipe of
Conrey–Farmer–Zirnbauer (CFZ 2008), as applied by Conrey–Snaith (Proc. LMS 94 (2007) 594, §7, for ζ), works in four steps:

1. Write the ratio L(s+α,χ)L(1−s+β,χ̄)/(L(s+γ,χ)L(1−s+δ,χ̄)) as a Dirichlet series.
2. Add the term with both numerator factors swapped by the functional equation. Its X-product X_χ(s+α)X_χ̄(1−s+β) does
   not oscillate in t. Terms with one swap oscillate like (kt/2πe)^(∓it) and drop out.
3. Keep the diagonal of the t-average. Every coefficient enters as χ(m)χ̄(n) with m = n, so only |χ(n)|² = 1_(n,k)=1 is
   left. Real and complex characters, and every exponent of p | k, therefore give the same terms.
4. Differentiate in α and β at γ = α, δ = β.

The Euler product of the diagonal is G(α,β,γ,δ) = Π_p G_p. Its polar part is
L(1+α+β,χ₀)L(1+γ+δ,χ₀)/(L(1+α+δ,χ₀)L(1+β+γ,χ₀)), where χ₀ is the principal character mod k. For p | k every coefficient
vanishes, so G_p = 1.

In the swapped term the factors 1/L(1−β+δ,χ₀) and 1/L(1−α+γ,χ₀) vanish at the evaluation point. Each contributes the
derivative of 1/L(1+z,χ₀) at z = 0. Since L(1+z,χ₀) = (φ(k)/k)/z + O(1), that derivative is k/φ(k). The result
(**proved** from the recipe; the recipe itself is a **conjecture**):

  ∫ (L′/L)(s+α,χ)(L′/L)(1−s+β,χ̄) dt = ∫ { (L′/L)′(1+α+β, χ₀) − B_k(α+β)
        + (k/φ(k))² X_χ(s+α) X_χ̄(1−s+β) L(1+α+β,χ₀) L(1−α−β,χ₀) A_k(α+β) } dt + O(T^(1/2+ε)),

  A_k(z) = Π_(p∤k) (1 − p^(−1−z))(1 − 2/p + p^(−1−z))/(1 − 1/p)²,   B_k(z) = Σ_(p∤k) (ln p)²/(p^(1+z) − 1)².

The two-point function follows by the contour argument of Conrey–Snaith §7. The diagonal is (1/4π²)·2 Re[(L′/L)′(1+iε,χ₀)
− B_k(iε)]. The off-diagonal is (1/4π²)·{(k/φ(k))² |L(1+iε,χ₀)|² A_k(−iε) · Ph + c.c.}. The phase Ph is the X-product on the
critical line, so it has modulus 1. With N̄(t) := θ_χ(t)/π and θ_χ(t) = (t/2)ln(k/π) + Im lnΓ((1/2 + a + it)/2):

  Ph = exp(2πi [N̄(t + ε/2) − N̄(t − ε/2)]).

### 1.2 Term by term against BK (proved; checked numerically)

- **Diagonal.** (L′/L)′(s,χ₀) = (ζ′/ζ)′(s) − Σ_(p|k) (ln p)² p^(−s)/(1 − p^(−s))². B_k drops the p | k terms of B. The
  change from p | k is therefore −(1/2π²) Σ_m (ln p)² cos(mε ln p)/p^m, which is BK's Ψ term (`bk/derivation.md` §4,
  step 2). For k = 1 the diagonal is BK's −(1/4π²)∂² ln[|ζ|²Φ_diag].
- **Off-diagonal, good p.** A_p(−iε) = 1 − (p^(iε) − 1)²/(p − 1)², which is BK's Φ_off factor.
- **Off-diagonal, bad p.** The ratios factor is (p/(p−1))² |1 − p^(−1−iε)|², times the p-part |1 − p^(−1−iε)|^(−2) of
  |ζ|²; the total is (p/(p−1))². In BK's product, write x = p^(iε). Then
  1 − (x−1)²/(p−1)² = (p − x)(p − 2 + x)/(p−1)² and F_p = (p − x̄)/(p − 2 + x), and

  |1 − x/p|^(−2) · (p − x)(p − 2 + x)/(p−1)² · (p − x̄)/(p − 2 + x) = p²/((p − x)(p − x̄)) · (p − x)(p − x̄)/(p−1)² = (p/(p−1))².

  The two are identical.
- **Numerically** (`theory_out.txt` (a1)): the ratios form, evaluated with direct Euler products over p < 10⁶ and
  mpmath ζ, against `bklib.C_of_r` (splines). For k = 1, 5, 8, 12, 105, 440, 664 and r = 0.25, 1.3, 3.7, 6.2, the
  largest difference is 1.0·10⁻¹². **Computed.**

So BK's law is the ratios conjecture's two-point function with all its lower-order terms. BK's paper states this in one
sentence for ζ ("the expression … that includes the lower order terms may also be obtained by a heuristic recipe"). The
per-prime identity for p | k, including that it does not depend on the exponent, is worked out here.

### 1.3 What the recipe leaves out at finite height, computed

The arithmetic factors are exact within the conjecture. The X-factor has modulus 1, so there is **no amplitude correction
at any order**. Everything at finite height sits in N̄. Stirling gives Im lnΓ(σ+iτ) = τ ln τ − τ + (σ − ½)π/2 −
(σ² − σ + 1/6)/(2τ) + O(τ⁻³), with σ = (½ + a)/2 and τ = t/2. Since σ² − σ + 1/6 = −1/48 for both σ = 1/4 and σ = 3/4,

  N̄(t) = (t/2π) ln(kt/(2πe)) + (2a − 1)/8 + 1/(48πt) + O(t⁻³).   (**proved**; checked: `theory_out.txt` (d) gives
  −0.124934 at t = 100 for a = 0, against −1/8 + 1/(4800π) = −0.124934.)

The first term is the test's u. The constant ±1/8 is invisible to window counts. By the reflection formula, the parity
enters the density only through terms O(e^(−πt)). Consequences:

- Against the test's unfolded phase 2πr = 2πΔu, the exact phase differs by 2πΔ[1/(48πt)] = −ε/(24t²). That is
  2.2·10⁻⁵ rad at t = 100 and ε = 5.5 (8 spacings). **Computed.**
- Against BK's literal εL(E) at the midpoint, the exact phase differs by −ε³/(24t²) − …: −7.2·10⁻⁴ rad at t = 100,
  ε = 5.5. This changes the p = 2 off-diagonal Σ² term by at most 1.8·10⁻⁶ (the term itself reaches 7.8·10⁻³). The
  projection ratio, which is what γ sees, is 1.000076. **Computed** (`theory_out.txt` (a2)).
- The phase difference between a = 0 and a = 1 is 0 to double precision.
- The recipe's error O(T^(1/2+ε)) is not a term of the recipe. Its size here is unknown (§6).

**Verdict (a): killed.** The largest finite-height effect on γ is 8·10⁻⁵, against the 0.07 needed. There is no 1/L or 1/t
term of the amplitude and no parity dependence. Independently, R120 has the same heights and a smaller L (9.39 against
10.23), so any 1/L or 1/t effect would be as large or larger there, and R120 shows γ = 1.001 ± 0.019.

## 2. Candidate (b): the drift of the density across a window

The test unfolds with u(t) = (t/2π) ln(kt/(2πe)). Its derivative is exactly (1/2π) ln(kt/2π), the BK density. So in u
the local density is 1 everywhere (to O(t⁻²), §1.3), and the pair phase is 2πΔu exactly. The drift of the raw density is
absorbed by the unfolding: it does not dephase anything. Two remainders are left (**computed**, `theory_out.txt` (b)):

- The pair normalisation d̄(t₁)d̄(t₂)/d̄(t)² − 1 is −9.4·10⁻⁵ at t = 100 for 8 spacings, and −1.5·10⁻⁶ at t = 600.
- The midpoint phase term, ε³/(24t²), is 7.1·10⁻⁴ rad at t = 100 (and does not even apply in u, §1.3).

What the brief's damping would be if the windows were unfolded at a fixed density (they are not): a pair at separation r,
whose midpoint sits δ away from the window centre, picks up the phase error εδ/t, with ε = r/d̄ and δ uniform on
±(ℓ − r)/(2d̄). The rms of that error peaks near r = ℓ/2 at (πℓ/L)²/(√3 t). For ℓ = 8, L = 9.1 and t = 100 that is
0.044 rad, an amplitude loss of 0.1 %. It grows like ℓ⁴/t². The data show neither growth: γ is 0.902 ± 0.052 for ℓ ≤ 4
and 0.941 ± 0.027 for ℓ > 4, and 0.942 ± 0.046 at heights 100–300 against 0.928 ± 0.045 at 300–600.

**Verdict (b): killed** (≤ 10⁻⁴ in the test's unfolding; 0.1 % even in a wrong one; wrong ℓ and height dependence).

## 3. Candidate (c): the bad-prime factor for k not squarefree

- In BK's derivation, S(k) = Π_(p|k) (p/(p−1))² comes from Ramanujan sums over squarefree divisors. It depends only on
  rad(k). The Gauss-sum identity Σ_r g(r) e^(2πirs/k) = |χ(s)|²|τ(χ)|² = k·1_((s,k)=1) holds for every primitive χ,
  including 4 ∥ k and 8 | k. So Ψ_off does not depend on the exponent.
- In the ratios recipe (§1.1), the local factor at p | k is 1 for every exponent, and the residue k/φ(k) depends on
  rad(k). The bad-prime factor is (p/(p−1))² in both cases (§1.2).

**Proved:** no exponent dependence in either heuristic. **Data** (review42): γ₂ = 0.936 ± 0.034 for 4 ∥ D and
0.943 ± 0.035 for 8 | D.

**Verdict (c): killed** (the correction is exactly 0; the data agree).

## 4. Candidate (d): statistics of the estimate and the unfolding

- **Unfolding.** θ_χ/π − u = (2a − 1)/8 + 1/(48πt) + O(t⁻³) (§1.3). Its non-constant part changes by 3.5·10⁻⁶ spacings
  across 8 spacings, a phase of 2.2·10⁻⁵ rad. A 7 % loss of a period-1 term by Gaussian phase jitter needs an rms phase of
  √(2 × 0.07) = 0.37 rad, or 0.060 spacings. The real error is 4·10⁴ times smaller. **Computed.**
- **Estimator.** The review found the bias of np.var to be 2·10⁻⁶ in the regression coefficients. The strongest check is
  empirical: the identical estimator, applied to R120 (about 745 zeros per set at the same heights), gives
  γ = 1.001 ± 0.019.
- **Zero completeness, recomputed here with the exact smooth count θ_χ/π** (`influence_out.txt`). In B, mean(N − θ/π)
  over 100–600 is 0.000 ± 0.004 (range −0.012 to +0.011), so all 184 sets are complete. The contributions to γ do not
  correlate with completeness (r = 0.10), drift (0.04) or the smallest gap (0.02). In R120 one set, D = −283, misses two
  zeros. Dropping it leaves γ = 1.0012 (`incomplete_out.txt`).

**Verdict (d): killed.**

## 5. Candidate (e): what the data themselves say

### 5.1 The same statistic on the other sample at the same heights (retrospective)

R120 has 181 sets, 100 ≤ |D| ≤ 400, heights 100–600 and mean L 9.386. The predictions come from the same frozen code; the
cached grid reproduces the frozen B and C predictions to 0.0e+00 (`alt.py` check). The fits reproduce `r123.json` to
1.5·10⁻¹⁶. Results:

- γ = 1.001 ± 0.019; in 100–300, 0.998 ± 0.028; in 300–600, 1.014 ± 0.030;
- p = 2: 0.989 ± 0.021; ℓ ≤ 4: 0.997; ℓ > 4: 1.002;
- D > 0: 1.043 ± 0.034; D < 0: 0.972 ± 0.034. This is the opposite order to B's, so the parity hint of B is noise.

This γ ≈ 1 was known before test B: `bk/derivation.md` §7.4 reports a residual-on-off slope of 0.995 on these data. It
counts as replication only in the weak sense that no parameter of BK was tuned on it.

### 5.2 Inside B: between-group and within-group parts (computed; the split is post hoc)

Write x_i(ℓ) = full_i − diag_i and y_i = S_i − diag_i. Let H be the projection on the design (constant and indicators of
2, 3, 5, 7, 11, 13). The **between** slope is that of Hy on Hx, both centred: this is what γ sees. The **within** slope is
that of (1−H)y on (1−H)x: it covers other primes dividing D, the non-additivity of the bad-prime product, and the spread
of L. Under any change of the law's amplitude, both parts move together. LP3 predicts 0.921 and 0.942 in B. Results
(`decomp_out.txt`, `table_out.txt`, null by resampling residual rows):

| sample | between | within | pooled (reading 6) | between − within (z) |
|---|---|---|---|---|
| B | 0.938 (null sd 0.021) | 1.042 (0.028) | 0.975 ± 0.014 | −0.104 (−2.9) |
| R120 | 0.996 (0.017) | 0.977 (0.025) | 0.990 ± 0.015 | +0.018 (+0.6) |

So the 7 % sits in B's indicator contrasts only. B's own registered reading (6), the pooled per-L-function slope, is
0.975 ± 0.014. That is 3.1σ above 0.932 and 1.8σ below 1. Under the null, γ − pooled has sd 0.015, and B's −0.043 sits at
z = −2.8 (`influence_out.txt`).

### 5.3 Joint reading (computed)

Four numbers (between and within, in B and in R120) with their null covariance (correlations −0.08 and +0.03):

| hypothesis | χ² | p |
|---|---|---|
| BK (all 1) | 11.7 on 4 | 0.019 |
| AMP (all 0.932) | 32.5 on 4 | 1.5·10⁻⁶ |
| LP3 (its own predictions) | 34.0 on 4 | 7.5·10⁻⁷ |
| one common amplitude 0.983 ± 0.011 | 9.3 on 3 | 0.026 |

Neither "plain amplitude 0.93" nor "L + 0.3" survives the data already in hand. BK survives at p = 0.019: B's between part
is its one bad number. No single amplitude fits all four numbers well. That is a sign of a fluctuation, or of an effect
specific to B's contrasts, rather than of a change of the law.

### 5.4 Things looked at and not found (exploratory)

- **Dependence on χ_D(2), χ_D(3) at good primes** (BK and the ratios conjecture: none; `chiexp_out.txt`). The projections
  on the p = 2 off pattern are +0.004 ± 0.013 and +0.003 ± 0.013 in B, and −0.015 ± 0.012 and −0.013 ± 0.011 in R120:
  nothing. One unrelated oddity: the χ(3) coefficient vector in B has χ² 36.1 over 16 lengths. Its null with an
  estimated covariance has mean about 17.7, so p ≈ 0.01 before allowing for the 4 looks. Not pursued.
- **Conductor dependence** (descriptive only, nothing derived). γ in the four half-samples is 1.001 ± 0.030
  (|D| ≤ 248), 1.017 ± 0.043 (> 248), 0.956 ± 0.041 (B, ≤ 551) and 0.917 ± 0.039 (> 551). A weighted line in mean ln|D|
  has slope −0.058 ± 0.035, which is 1.7σ, with χ² 1.34 on 2 about the line (`regfinal_out.txt`). I found no mechanism.
  The one k-specific approximation in BK is the HL average over the residue classes mod k. It needs many n per class, but
  n/k ≈ E/(2πR) does not depend on k. The k/E variant is contradicted by B's equal γ in its two height bands, where k/E
  differs by a factor of 2.3.

**Verdict (e).** No mechanism found. The 7 % is not an amplitude of the off-diagonal term. It is an excursion of B's
indicator contrasts: z = −2.9 against B's own within part, and −2.8 against BK for the registered reading alone. Chance is
the economical reading but is not proved.

## 6. Doubts

1. R120 is retrospective. BK was adopted after round 123 found the period-1 residual on these zeros, and γ ≈ 1 there was
   on record (`bk/derivation.md` §7.4, slope 0.995). No parameter was tuned, but it is not a blind replication.
2. The between/within split was chosen after seeing B (post hoc). The joint p-values carry that caveat.
3. The recipe's error term O(T^(1/2+ε)) has an unknown constant, and nothing is known about its uniformity in k at T ≈ k
   (here k/T is 0.7–7 in B). I found no term of the recipe that depends on k other than through L and the bad primes. A
   genuine effect for k ≳ T is not excluded by theory; the descriptive trend in §5.4 is 1.7σ.
4. The expected errors for the fresh tests come from a noise model: R120's residual rows at 100–600, scaled by
   √(zeros ratio). For B it gives se 0.025; B's own residual rows give 0.023; the bootstrap gives 0.024. At L ≈ 11
   (test Bp) the per-zero noise may differ.
5. AMP and LP3 differ by less than 0.01 in γ in every sample, so no test of feasible size separates them from each other.
   The registration separates both from BK.
6. Test A (complex characters, 0.943 ± 0.104) cannot tell any of these apart.
7. My band split of B uses the first and last zero in each band (as `r134.py`): 0.942 and 0.928. The review used the band
   edges: 0.938 and 0.931.

## 7. What is classical and what is derived here

- **Classical:** the ratios conjecture and its recipe (CFZ 2008), and the pair-correlation formula derived from it for ζ
  (Conrey–Snaith 2007 §7); BK's law (2013) and its HL derivation; Stirling's series; the reflection formula.
- **Derived here** (elementary consequences, not found in print, not searched beyond the two papers):
  - the recipe's two-point function written out for primitive χ mod k, with the bad-prime factor (k/φ(k))² from the
    residue of L(1+z, χ₀);
  - the per-prime identity with BK's F_p, and the proof that the exponent of p | k does not enter;
  - N̄(t) = u(t) + (2a−1)/8 + 1/(48πt) + O(t⁻³), with no parity dependence beyond the constant, and the resulting phase
    correction −ε/(24t²) in the test's units;
  - the bounds for (b) and (d);
  - the between/within decomposition as a check that a change of the law must pass and a fluctuation of group contrasts
    need not.
- **Computed / retrospective:** every γ, slope, χ² and p-value above.
- **Not derived:** any formula producing 0.93. None was found.
