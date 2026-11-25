# The maximum and the frozen partition function of |L(1/2+it, χ_D)| per prime dividing the conductor

Scratch directory: `writeup/code/study/lowzero/freeze/`
(every file named below is in it, unless it has a repo path). Labels: **theorem** (in print, cited), **derived**
(worked out here, elementary), **checked** (verified numerically by a separate script), **conjecture**,
**computed**, **dev** (measured on development characters only).

## 0. Result in four lines

* Per missing prime p, the frozen observables (the interval maximum M, ln Z_β for β > 1) shift by about
  (1/2)·ln R_p and β·ln R_p, where R_p = (1 − 1/p)(1 + 2 ln p/((p − 1)L′)) is the per-prime factor of the **second
  moment including its lower-order term** (theorem), not of the β-th moment (the "not frozen" value) and not of
  the leading Euler factor (1 − 1/p) only (the "independent Euler factor" value). On development data the
  not-frozen alternatives are already rejected at more than 25 errors, so freezing itself is settled at this height.
* The open question that can be measured is how the missing prime couples to the extreme values. There are two
  parameter-free closures. **Removal**: |L|·|E_p| behaves like a genuine L-function without p. **Addition**:
  |L_missing|/|E_p| with an independent phase behaves like a genuine L-function with p (this is the GHK splitting).
  For the mean square, removal holds exactly and addition fails by the lower-order term. This is a theorem
  (twisted second moment), and the development data reproduce it to 0.1–0.5 %. For the maximum and the frozen
  ln Z_β the development data sit halfway (γ = 0.50 ± 0.29 on a scale where removal is 0 and addition is 1). For
  the variance of M over intervals they follow removal (γ = −0.02 ± 0.12).
* A test with 240 fresh characters (1500 < |D| ≤ 2400) separates removal from addition by 13–14 errors for M and
  for ln Z_3, and the not-frozen alternatives by about 100 errors.
* Prior art: the arithmetic of the maximum is called "not completely clear" by Fyodorov–Keating (2014). No per-prime
  statement and no numerics for Dirichlet L-functions in t-aspect were found (§9).

## 1. Objects

χ = χ_D is real primitive, q = |D|, L = ln(qt/2π) (the test has mean L̄ = 11.89, range 11.64–12.11; the development
set has L̄ ≈ 11.5). The values are |L(1/2+it, χ)| on t_k = 100 + k/20 (moment_values.py, unchanged; sha256
a2f4c8630c1a33bc04439d0b785417f527f3e332a8fcf92caed9ca17172e9de5). There are 143 intervals I_j = [100 + 2πj,
100 + 2π(j+1)), each holding about L ≈ 12 zeros. This is FHK's matching: an interval of length 2π corresponds to
CUE(N = L).

Per interval:
* M_j = max_{I_j} ln|L|;
* Z_β = (1/2π)∫_{I_j}|L|^{2β}, for β = 0.5, 1, 1.5, 2, 3.

Per character:
* the means over j of M and of ln Z_β (quenched);
* slope32 = mean(ln Z_3 − ln Z_2);
* varM = the variance of M_j over j;
* lnEZ_β = ln of the mean over j of Z_β (annealed; for β = 1 this is the mean square).

Contrast: the OLS coefficient c_p of I(p | D) in a regression on [1, ln q − mean, I(p|D) for p = 2, 3, 5, 7, 11, 13].
It is read as "p missing minus p present".

**Discretisation (checked, refine_check.py, 386 intervals of 3 characters against PARI on a step of 0.0025).** The
plain grid maximum is low by 0.0016–0.0018 on average and by at most 0.0066. The parabola through the grid maximum
and its two neighbours leaves an error of −0.00002 ± 0.00006 (at most 0.0005). About 10 % of the maxima lie at an
interval end, which is not a grid point. There the value at the exact end is used (quadratic interpolation through
the three nearest grid points; error about 10⁻³). So the step 0.05 is fine once refined. The bias that remains is
10⁻⁴ of the smallest contrast.

## 2. Theorem-level facts (these test nothing; they are used as controls)

(T1) **Second moment with every lower-order term**, for q fixed and T → ∞ (classical):
  (1/T)∫|L(1/2+it,χ)|² = (φ(q)/q)(ln(qt/2π) + 2γ + 2Σ_{ℓ|q} ln ℓ/(ℓ−1)) + o(1).
  So the per-prime ratio is
  **R_p = (1 − 1/p)(1 + 2 ln p/((p−1)L′))**, with L′ = L + 2γ + 2Σ_{ℓ|q, ℓ≠p} ln ℓ/(ℓ−1).
  Development check: the regression coefficient of lnEZ_1 is 0.993 ± 0.015 times ln R_p. The leading factor alone
  would give 1.165 (control C1).

(T2) **Twisted second moment** (the same diagonal argument; for ζ, Balasubramanian–Conrey–Heath-Brown 1985):
  (1/T)∫|L|²p^{it} = χ(p)p^{−1/2}(φ/q)(L′ − ln p) + o(1). Derived from it:
  * The **removal** ratio is mean|L·E_p|²/mean|L|² = 1 + 1/p − (2/p)(L′ − ln p)/L′ = R_p **exactly** (identity
    checked in theory_check.py).
  * The **addition** ratio is mean|L_miss/E_p|²/mean|L_miss|² = 1/(1 − 1/p) **exactly**, because χ(p) = 0 kills
    every twist by a power of p.
  * Development check (C2, C4): removal/twisted = 1.0054, 1.0012, 0.9987, 1.0004 (± 0.002) for p = 2, 3, 5, 7;
    addition·(1 − 1/p) = 0.996, 1.003, 0.999, 0.999. The leading-only removal value would be 0.908–0.956.

  **So for the annealed second moment, removal = genuine and addition ≠ genuine. This is a theorem.**

(T3) **The Fourier coefficient of ln|L(1/2+it)| at frequency ln p is χ(p)p^{−1/2}/2.** This follows from RH and
  Littlewood's lemma on the rectangle ½ ≤ σ ≤ 2 (derived). Development check: 1.000 ± 0.004 times that value for
  every p (dev_closure.py).

(T4) **The mean of ln|L| over t is O(ln(qT)/T)** (Littlewood's lemma and GRH). So the per-prime shift of mean_ln
  is 0. Development: |c_p| ≤ 0.002 ± 0.001 (control C3).

(T5) **Fourth moment of an individual Dirichlet L-function in t-aspect**, with every lower-order term: Topacogullari,
  Math. Z. 2021 (the error is unproved at our sizes). The annealed Z_2 is therefore theorem-grade. Its per-prime
  terms were not computed here (no reading). Development: removal is close to genuine and addition far from it
  (γ = 0.10 ± 0.06).

(T6) **Freezing at leading order** (ζ, Arguin–Ouimet–Radziwiłł 2021): the moments over intervals of length (ln T)^θ
  are (ln T)^{f(β)+o(1)} with a frozen branch above β_c. **The maximum up to tightness** (ζ, Arguin–Bourgade–
  Radziwiłł 2020/2023; Harper 2019 for the upper bound in typical intervals). Nothing at the O(1) level, nothing
  per prime, and nothing for L(s, χ) in t-aspect.

## 3. The prime as a probe of the Gibbs measure (derived)

Let Y_p(t) = −ln|1 − χ(p)p^{−1/2−it}|, the Euler factor in log form. Write the field with p present as
"base + Y_p", where the base is the field without p. For any interval statistic,

  ln Z_β(base + Y_p) − ln Z_β(base) = ln ⟨e^{2βY_p}⟩_{μ_β},     μ_β = Gibbs measure e^{2β·base}/Z_β on the interval.   (3.1)

So **c_p(β) = −E ln⟨e^{2βY_p}⟩_{μ_β}** (missing minus present), if the base is the same with and without p and is
independent of p's phase (the addition closure). The value depends only on the step law of p and on how μ_β sits
relative to the wavelength 2π/ln p (9.1, 5.7, 3.9, 3.2 for p = 2, 3, 5, 7, against the interval 2π). Limits:

* **High temperature, Gibbs measure spread** (μ uniform on the interval): c^{spread}_p(β) = −E_φ ln (1/2π)∫e^{2βY_p}.
  For many cycles this becomes −ln E e^{2βY_p} = −ln ₂F₁(β, β; 1; 1/p), the per-prime factor of the 2β-th moment
  ("naive", quenched = annealed).
* **Frozen, Poisson–Dirichlet extremal process** with intensity ∝ e^{−2x} and positions spread over many
  wavelengths: marking every point by an independent Y_p is a shift of the process by (1/2)ln E e^{2Y_p}. Hence
  ΔM = (1/2)ln(1 − 1/p) and c_p(β) = β·ln(1 − 1/p) for β > 1. The free energy per prime is frozen at its β = 1
  value. "A shift of the location of M by (1/2)·ln of the critical factor" and "the tail prefactor multiplied by
  the critical factor" are the same statement, because the tail is exponential with rate 2.
* **Trunk** (the prime constant over the interval, θ ln p ≪ 1): M shifts by exactly −Y_p. The mean shift is
  E Y_p = 0, and Var M drops by Var Y_p = Li₂(1/p)/2 (0.291, 0.183, 0.106, 0.074). A trunk prime changes the shape
  of the distribution of M; a branching prime changes its location. With 0.7–1.9 cycles per interval,
  p = 2, …, 7 are mixtures.
* **Small β** (cumulant expansion): c_p(β) = −2β²[Var Y_p − Var(interval mean of Y_p)] + O(β³). At β = 0.25 this
  is −0.0318, −0.0227, −0.0128, −0.0093, against the exact spread values −0.0336, −0.0237, −0.0131, −0.0094
  (theory_check.py). This part is the variance of ln|L| over an interval: the explicit formula plus pair
  correlation (Selberg's mean shift 0 and variance Li₂(1/p)/2 per prime). It is a restatement, so β = 0.5 serves
  only as a control. Development: c(0.5) = −0.124, −0.082, −0.046, −0.034; the second-order values are −0.127,
  −0.091, −0.051, −0.037.

## 4. The closure question and H_crit (conjecture; the content of the test)

T1–T2 show that a genuine L_missing is not "L_present with E_p divided out and an independent base" at order 1/L.
The removal surgery reproduces the lower-order term exactly. The addition surgery does not.

**Removal closure, mechanism (derived, heuristic).** By the explicit formula, N_present(t) = N_missing(t) + S_p(t),
with S_p = −(1/π)arg(1 − χ(p)p^{−1/2−it}) (round 147's wave). The zeros are sparser where Y_p is high:
dS_p/dt ≈ −(ln p/π)Y_p. Near-field: ln|L_present(t)| ≈ base(t + δ(t)) + Y_p(t), with δ = (2π/L)S_p. The maximum
does not depend on the parametrisation. But the density of extremal candidates is the zero density, so the
extremal process of L_present has intensity ∝ (1 + δ′(u))e^{−2x} with marks Y_p(u). In the many-cycle limit:

  ΔM_wave = −(1/2) ln E_ψ[(1 + δ′)|1 − w|^{−2}] = −(1/2) ln[(p/(p−1))(1 − 2 ln p/((p−1)L))],

using E(r² − r cos ψ)|1 − w|^{−4} = −p/(p−1)² (identity checked, theory_check.py). To first order in 1/L this is
(1/2)·ln R_p. At L = 11.89: −0.2846/−0.2914, −0.1543/−0.1585, −0.0765/−0.0788, −0.0490/−0.0505
(wave-PPP / (1/2)ln R_p), against the leading value −0.3466, −0.2027, −0.1116, −0.0771.

**Hypotheses** (per missing prime, regression-matched to the test design, L̄ = 11.89; freeze_predict.py,
predictions_freeze.json sha256 83c24b307b40d398b6a6be82aa6866863145448e4a87b14a6ae3a40bd02abfdb):

| hypothesis | M | ln Z_1.5 | ln Z_2 | ln Z_3 | slope32 = Z_3 − Z_2 |
|---|---|---|---|---|---|
| crit (frozen; exact critical factor R_p) | −0.302 −0.167 −0.086 −0.057 | −0.904 −0.502 −0.258 −0.171 | −1.206 −0.669 −0.344 −0.227 | −1.809 −1.004 −0.517 −0.341 | −0.603 −0.335 −0.172 −0.114 |
| lead (frozen; leading factor 1 − 1/p) | −0.347 −0.203 −0.112 −0.077 | −1.040 −0.608 −0.335 −0.231 | −1.386 −0.811 −0.446 −0.308 | −2.079 −1.216 −0.669 −0.463 | −0.693 −0.406 −0.223 −0.154 |
| naive (not frozen: 2β-th moment factor) | no finite value | −1.508 −0.893 −0.496 −0.344 | −2.485 −1.504 −0.852 −0.596 | −4.644 −2.921 −1.726 −1.236 | −2.160 −1.417 −0.874 −0.640 |
| spread (not frozen: uniform Gibbs measure, finite θ) | −1.071 −0.861 −0.593 −0.475 | −1.220 −0.886 −0.481 −0.343 | −2.009 −1.495 −0.830 −0.596 | −3.791 −2.908 −1.694 −1.235 | −1.782 −1.413 −0.864 −0.639 |
| zero (no dependence) | 0 | 0 | 0 | 0 | 0 |

Relational hypotheses, needing no numbers:
* **H_W (removal):** c_p = r_p, the mean removal contrast over the test characters with p ∤ D. Exact for M under the
  wave closure. For ln Z_β the warp Jacobian enters at second order (|.| ≲ 0.005).
* **H_G (addition, GHK):** c_p = a_p, the mean addition contrast over the characters with p | D.

Reading: γ = Σ_p w_p (c_p − r_p)(a_p − r_p) / Σ_p w_p (a_p − r_p)². H_W gives 0, H_G gives 1.

What is conjecture here:
* the FHK-type statement that the extremal process carries the critical (β = 1) factor;
* that the exact finite-L factor (with its lower-order term) is the one that matters;
* the closure. The removal closure is proved only for the annealed second moment.

What is not: the numbers R_p (theorem) and the algebra (checked).

**Theory's own uncertainty (not in the measurement bands):**
* crit is a many-cycle, exponential-tail (PPP) formula. With 0.7 cycles of p = 2 per interval and N ≈ 12 the
  finite-size factors are visible. On development data the removal amplitudes on the crit pattern are 0.94 (M),
  0.79 (ln Z_1.5), 0.84 (ln Z_2), 0.89 (ln Z_3) and 0.99 (slope32). The frozen slope converges fastest.
* In both CUE-based synthetic worlds the GHK world gives amplitude 1.00 for M, because the finite-interval Jensen
  reduction mimics the lower-order term. **So the absolute reading cannot tell the mechanisms apart. Only the
  relational reading γ can.**
* Stated bands: crit M [0.92, 1.08], crit ln Z_3 [0.85, 1.10], frozen slope32 [0.90, 1.20], H_W |γ| ≤ 0.15.

## 5. Moments of moments: the arithmetic factor per bad prime (derived; conjecture except where marked)

MoM(k, β) = (1/T)∫ Z_β(t)^k dt, with Z_β = (1/2π)∫_{I(t)}|L|^{2β}. For integer k and β this is a shifted 2kβ-th
moment integrated over the k positions. Its local factor at a bad prime is 1. At a good prime it depends on the
shifts through p^{i(h_a − h_b)}. Bailey–Keating (ζ, from CFKRS) give degree k²β² − k + 1 for kβ² > 1 and kβ² for
kβ² < 1.
* For kβ² > 1 the leading order comes from coalesced configurations (all shifts ≲ 1/L). There p's factor is that of
  the 2kβ-th moment, so the per-prime factor of MoM is 1/₂F₁(kβ, kβ; 1; 1/p) at leading order, with O(1/L)
  corrections from non-coalesced shifts.
* For kβ² < 1 the k points decorrelate. A branching prime then gives [1/₂F₁(β, β; 1; 1/p)]^k, and a trunk prime
  gives 1/₂F₁(kβ, kβ; 1; 1/p), because it is common to all k copies.
* I could not open the full text of Bailey–Keating, so their α_{k,β} was not compared with this.

Status:
* k = 1, β = 1: theorem (T1).
* k = 1, β = 2: theorem for fixed q (T5).
* k = 2, β = 1: a shifted fourth moment. Theorem for ζ. For a fixed Dirichlet L in t-aspect it is expected from
  Topacogullari's method, **not checked**.
* All else: conjecture (CFKRS / Bailey–Keating; Curran 2023 proves the k = 2 freezing at β = 1/√2 for ζ).

The MoM are annealed and dominated by rare spikes (dev: the residual sd of lnEZ_2 is 0.10 per character), and they
are theorem-grade where computable, so **no MoM reading is registered** apart from the controls C1, C2, C4.

## 6. Random models (computed; used only as dry runs, not for predictions)

modellib.py, model_first.py:
* hyb_k0: CUE(N = L) on the circle mapped onto the interval, Fourier modes k < k0 replaced by explicit Euler factors
  of the primes with ln p < k0 − ½;
* cue: probes on top of CUE;
* ep: randomized Euler product up to X = exp(e^{−γ}L).

At N = 11 the hyb models give ΔM = −0.249, −0.145, −0.093, −0.068 and c(ln Z_3) = −1.43, −0.86, −0.53, −0.39
(hyb2/3/4/ep agree within 5 %; cue is lower). But their levels are off: M 2.36 against 2.19 real, var M 0.34–0.37
against 0.20.

synth.py builds |L| on the data grid in two worlds: CUE fine structure, primes ≤ 31 exact, and the CUE zeros either
warped by the primes' counting waves ("wave") or not ("ghk"). Both worlds fail the theorem control C2
(removal/twisted 0.91–0.96). The CUE+Euler-factor models lack the lower-order term of the mean square. In both
worlds c_p = r_p within 0.01 (M), so freeze_test.py recovers the closure that holds there. |a_p − r_p| there is
only 0.005–0.02, against 0.03–0.06 in the real data.

## 7. Development results (32 characters: round 149's 8 dev and 24 new with 1300 < |D| ≤ 1500; dryrun_dev_out.txt)

| statistic | c_p (p = 2, 3, 5, 7) | removal r_p | addition a_p | γ |
|---|---|---|---|---|
| M | −0.305(19) −0.177(22) −0.120(26) −0.093(22) | −0.284 −0.160 −0.086 −0.054 | −0.346 −0.213 −0.124 −0.085 | 0.50 ± 0.29 |
| ln Z_3 | −1.70(10) −1.03(11) −0.66(13) −0.52(11) | −1.590 −0.946 −0.498 −0.333 | −1.878 −1.186 −0.666 −0.470 | 0.53 ± 0.33 |
| var M | −0.119(6) −0.048(6) −0.010(7) −0.009(8) | −0.121 −0.044 −0.015 −0.009 | −0.172 −0.078 −0.040 −0.024 | −0.02 ± 0.12 |
| lnEZ_2 | −1.970(32) −1.082(36) −0.579(38) −0.377(42) | −1.919 −1.040 −0.539 −0.306 | −2.472 −1.539 −0.848 −0.590 | 0.10 ± 0.06 |

Amplitudes on the crit pattern: M 1.04 ± 0.08, ln Z_3 0.97 ± 0.06, slope32 1.10 ± 0.09.
* naive (3.78) and spread (3.30) on slope32 are 30 and 25 errors away: **freezing is settled already on development
  data.**
* lead (1.165) is 1.6 errors away for M.

Noise per character: residual sd 0.042 (M), 0.22 (ln Z_3), 0.014 (var M), 0.035 (lnEZ_1). The within-character
error of a mean over the 143 intervals is 0.025 (M); neighbouring intervals are anti-correlated (lag-1 −0.2 to
−0.6, from p = 2: half a turn of 2^{−it} per interval).

## 8. Power (power.py; resampling the development characters)

At n = 240:
* se(γ) = 0.070 (M), 0.078 (ln Z_3), 0.034 (var M). H_W and H_G are 14, 13 and 29 errors apart.
* se(amplitude) = 0.022 (M), 0.018 (ln Z_3), 0.024 (slope32).
* lead against crit: 7.7 and 9.3 errors, but only 4.0 and 4.8 errors from the edge of the crit band. So the absolute
  reading decides crit against lead only marginally.

To separate H_W from H_G at 5 errors needs n ≈ 80–100. To exclude both at 5 errors if the truth is halfway
(γ = 0.5, as on dev) needs n ≈ 310–380. 240 gives about 7 errors from each in that case.

## 9. Prior art (mathematical searches; abstract pages opened as stated)

* Fyodorov–Keating, Phil. Trans. R. Soc. A 372 (2014), arXiv 1211.6063 (opened, ar5iv text):
  * conjectures freezing for Z_T(β) on [T, T+2π] with N = ln(T/2π);
  * the arithmetic factor a(x) = Π_p(1 − 1/p)^{x²}Σ_m (Γ(x+m)/(m!Γ(x)))²p^{−m} appears in the high-point measure;
  * for the maximum they say it is "not at this stage completely clear" how arithmetic modifies the formulas.
  * KNOWN IN SUBSTANCE for H_lead as a ζ statement (the a(1) factor at the critical level). Not stated per prime
    or for L(s, χ).
* Fyodorov–Hiary–Keating, PRL 108 (2012): not opened, known from FK 2014's account.
* Arguin–Belius–Bourgade–Radziwiłł–Soundararajan (CPAM 2019, arXiv 1612.08575; search result only): leading order
  for ζ.
* Arguin–Bourgade–Radziwiłł I/II (arXiv 2007.00988, 2307.00982; search results): upper bound with the y·e^{−2y}
  tail and tightness, for ζ.
* Harper, arXiv 1906.05783 (abstract opened): the critical partition function and the upper bound for the maximum
  in typical intervals, for ζ.
* Arguin–Ouimet–Radziwiłł, Ann. Probab. 49 (2021) (search result): freezing for moments on intervals (ln T)^θ.
* Bailey–Keating, J. Number Theory (arXiv 2006.04503, abstract opened): MoM for integer k, β from CFKRS, with an
  arithmetic factor α_{k,β}. The full text was not opened.
* Curran, arXiv 2301.10634 (abstract opened): MoM(2, β) freezes at β = 1/√2.
* Arguin–Dubach–Hartung, arXiv 2103.04817 (abstract opened): random model, intervals (ln T)^θ.
* Amzallag–Arguin–Bailey–Hui–Rao, arXiv 2104.07403 (title seen): numerics of random-matrix corrections in Selberg
  large deviations, ζ.
* Atherfold, arXiv 2608.09906 (Aug 2026, abstract opened): β = 2 partition function and typical maxima of L(1/2+ih,
  χ) in q-aspect for prime q. There are no bad primes there.
* Searches for FHK statements or numerics for Dirichlet L-functions in t-aspect, or per prime dividing the
  conductor: none found. Arguin–Tai ("1-RSB", 2018) not opened.

**Verdict.**
* Per-prime frozen shift with the exact critical factor (H_crit), and the removal-versus-addition closure for
  extreme values: NOT FOUND.
* The leading factor at the critical level (H_lead): KNOWN IN SUBSTANCE for ζ (FK 2014 high points), not per prime.
* Any test: NOT FOUND.

## 10. Doubts

1. **The relational reading mixes between-character and within-character quantities.** The removal and addition
   groups differ in their other bad primes. In the GHK synthetic world this biases a_p − r_p by up to +0.02
   (M, p = 2), against −0.06 in the real data. c_p − r_p is unbiased in both synthetic worlds (≤ 0.008).
2. **The absolute crit reading is not diagnostic of the mechanism** (§4). It is secondary.
3. **The not-frozen alternatives are straw men at N ≈ 12.** Reading T4 confirms on fresh data what development data
   already show.
4. **Height.** The test L̄ is 11.9, not the 9–11 asked. Lower heights need |D| < 1300, where round 149's reserved
   characters and the zero sets of rounds 134–147 live. 9–11 would need t ∈ [10, 400] with q ≈ 300–1300: a
   different grid and a new value script.
5. **Discretisation and the anti-correlated intervals** are handled. The bootstrap is over characters.
6. **The mechanism of the removal closure for the max (candidate density modulated by the zero wave) is heuristic.**
   The synthetic wave world reproduces only part of the real removal/addition gap. The real gap also contains the
   approximate-functional-equation cutoff effect that produces the lower-order term of T1.
7. The 24 new development characters (1300 < |D| ≤ 1500) are now used. The test set is disjoint from them, from
   round 149's reserved 56 and from its development 8 (checked).
