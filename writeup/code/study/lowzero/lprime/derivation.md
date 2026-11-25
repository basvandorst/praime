# L′(ρ, χ_D) at the zeros: what the conductor's primes predict

Status labels: **proved** (a theorem, or an identity checked to machine precision), **conjecture** (the CFKRS recipe, the
ratios conjecture, the HKO/hybrid model, or a published conjecture), **heuristic** (derived here with a model step that
is not a theorem), **computed** (numbers from the scripts here), **measured** (data, dev set B only, noise only).

Setting: χ = χ_D, the real primitive character of a fundamental discriminant D, q = |D|. Zeros ρ = 1/2 + iγ with
100 ≤ γ ≤ 600, L = L(γ) = ln(qγ/2π) (7.4 to 11.7 across the sets; zero-weighted mean 1/L = 0.107 for R120, 0.098 for
B, 0.094 for Dn, 0.091 for Fn). Y = log|L′(ρ)|. "Bad prime" = p | q (χ(p) = 0: a missing Euler factor, a missing orbit).

Scripts (all in this folder): `dev/theory.py` (sections 1 and 2), `dev/cue_sim.py` (CUE check), `dev/recipe_k1b.py`
(Cauchy-integral recipe), `lprime_predict.py` (every frozen number), `lprime_values.py`, `lprime_test.py`,
`check_deriv.py`.

## 1. Random-matrix side (proved)

For A in CUE(N) with characteristic polynomial Λ(z) = ∏(1 − z e^{−iθ_j}), |Λ′(θ_n)| = ∏_{j≠n}|e^{iθ_n} − e^{iθ_j}|. By
rotation invariance put θ_n = 0; the other N − 1 eigenvalues then have the CUE(N − 1) density times ∏_j|1 − e^{iφ_j}|²
(Palm measure). Hence

  E|Λ′(θ_n)|^{2k} = M_{N−1}(1 + k) / M_{N−1}(1),   M_n(s) = E_{CUE(n)}|Λ(1)|^{2s} = ∏_{j=1}^{n} Γ(j)Γ(j+2s)/Γ(j+s)²,

with M_n(1) = n + 1, so E|Λ′|^{2k} = M_{N−1}(1+k)/N (this is the exact finite-N form of Hughes–Keating–O'Connell; at
large N it is [G²(k+2)/G(2k+3)]·N^{k(k+2)}). The cumulant generating function of Y_N = log|Λ′(θ_n)| is
K_N(t) = ln M_{N−1}(1 + t/2) − ln M_{N−1}(1). Differentiating (ψ = digamma):

* mean: K′(0) = Σ_{j=1}^{N−1}[ψ(j+2) − ψ(j+1)] = H_N − 1 = ln N + γ_E − 1 + 1/(2N) − … (exact);
* variance: K″(0) = Σ_{j=1}^{N−1}[ψ′(j+2) − ½ψ′(j+1)] = ½ ln N + (γ_E + 1)/2 + 1 − π²/4 + O(1/N)
  = ½ ln N − 0.67879 + O(1/N);
* κ₃ → −0.2583, κ₄ → +0.1347 (N → ∞); left-skewed;
* E|Λ′|^{−2} = M_{N−1}(0)/N = 1/N exactly;
* left tail: E|Λ′|^{2k} has a simple pole at k = −3/2, so P(|Λ′| < y) ~ C_N y³ with C_N = (2/3)·Res/N
  (C_N = 0.0242, 0.0207, 0.0182 at N = 8, 10, 12).

| N | mean | var | κ₃ | κ₄ | E\|Λ′\|^{−2} | E\|Λ′\|^{−1} | E\|Λ′\| | E\|Λ′\|² |
|---|---|---|---|---|---|---|---|---|
| 8 | 1.7179 | 0.4791 | −0.1805 | 0.1243 | 0.12500 | 0.23655 | 6.901 | 67.5 |
| 10 | 1.9290 | 0.5681 | −0.1944 | 0.1277 | 0.10000 | 0.20074 | 8.891 | 121.0 |
| 12 | 2.1032 | 0.6439 | −0.2042 | 0.1296 | 0.08333 | 0.17546 | 10.975 | 197.2 |

CUE simulation (`dev/cue_sim.py`, Haar via QR, 30,000 matrices N = 10 and 12,000 N = 20): mean 1.9310 / 2.5972 (exact
1.9290 / 2.5977), variance 0.5656 / 0.8690 (0.5681 / 0.8679), κ₃ −0.192 / −0.223 (−0.194 / −0.225), κ₄ 0.122 / 0.129
(0.128 / 0.133), E|Λ′|^{2k} for k = −½, ½, 1 within 0.5 %; k = −1: 0.094 / 0.049 against 0.100 / 0.050 (the infinite
variance of |Λ′|^{−2} makes the sample mean low and erratic: the same heavy tail as in the data). **Computed; agrees.**

**N_eff.** The identification N ↔ L is only a leading-order statement (HKO: E over zeros of |ζ′(ρ)|^{2k} ~
a_k·[G²(k+2)/G(2k+3)]·L^{k(k+2)}). The constants of the mean and variance at finite L are *not* fixed by it: there are
1/L terms from every prime (section 3) whose sum over all primes is not convergent and merges with the random-matrix
constant. Two data points for orientation (not used in any prediction): for ζ (813 zeros below 1200, ⟨1/L⟩ = 0.21)
mean(Y − ln L) = −0.262 against γ_E − 1 = −0.423; Hiary–Odlyzko found at height 10²³ mean 3.49 and sd 1.10 against
the CUE values 3.47 and 1.12. **The absolute level of the mean and variance at L ≈ 10 is therefore not predicted**;
every reading below is a per-prime difference with the level as a nuisance, except the k = 1 moment (section 4), whose
recipe gives all lower-order terms.

## 2. Arithmetic at leading order (conjecture: CFKRS / HKO recipe, hybrid Euler–Hadamard model)

The recipe for discrete moments gives E_zeros|L′(ρ)|^{2k} ≈ a_k(q)·E_CUE(N=L)|Λ′|^{2k} with
a_k(q) = ∏_p (1 − 1/p)^{k²} Σ_m d_k(p^m)²|χ(p)|^{2m} p^{−m}. For p ∤ q the sum is ₂F₁(k, k; 1; 1/p); for p | q it is 1. A
bad prime therefore multiplies every moment by

  r_p(k) = 1/₂F₁(k, k; 1; 1/p):   k = −1: p/(p+1) (Pearce-Crump's A_q);  k = 1: 1 − 1/p;  k = −3/2: 0.4628 (p = 2).

Equivalently (hybrid model, Gonek–Hughes–Keating; Bui–Gonek–Milinovich for zeros): Y = Σ_p A_p + U with
A_p = −ln|1 − χ(p)p^{−1/2−iγ}| (one Euler factor) and U the zeros part, at leading order independent, the phase γ ln p
equidistributed. ln E e^{tA_p} = ln ₂F₁(t/2, t/2; 1; 1/p) (proved, Parseval), so a bad prime removes the cumulants of
A_p (proved for the uniform phase; their use for L′ is the conjecture):

| p | κ₁ | κ₂ = ½Li₂(1/p) | κ₃ | κ₄ |
|---|---|---|---|---|
| 2 | 0 | 0.29112 | +0.14213 | −0.03747 |
| 3 | 0 | 0.18311 | +0.05374 | −0.03014 |
| 5 | 0 | 0.10550 | +0.01732 | −0.01309 |
| 7 | 0 | 0.07416 | +0.00846 | −0.00703 |

(κ₃ = (3/2)Σ_{m≥2} H_{m−1}p^{−m}/m², from ((k)_m/m!)² = (k²/m²)(1 + 2kH_{m−1} + …); checked numerically.)

## 3. The first order in 1/L: what the zeros know about γ ln p

### 3.1 Landau–Gonek (theorem, Gonek's uniform form)

Σ_{0<γ≤T} x^{iγ} = −(T/2π)χ(x)Λ(x)/√x + small, so over the zeros at height t:
E[p^{−imγ}] = −χ(p)^m ln p/(p^{m/2}L). The phase of γ ln p is not uniform at a zero; its density is
f(φ) = (1/2π)(1 + ε_p(φ)) with

  ε_p(φ) = −(2 ln p/L)·Re[w/(1 − w)],  w = χ(p)p^{−1/2}e^{−iφ}.

The same ε_p is exactly the p-part of the smoothed density of zeros, S′(t) = (1/π)Re(L′/L)(1/2+it) divided by L/2π
(identity). So: the zeros are locally denser by the factor 1 + ε_p where the Euler factor is small.

### 3.2 The mean: a per-prime shift (derived here; conjecture-level input)

Hadamard: Y = Σ_{ρ′≠ρ} ln|γ − γ′| + (terms smooth in γ), so E Y is linear in the pair density R₂. The ratios-conjecture
pair correlation of L(s, χ) (Bogomolny–Keating / Conrey–Snaith with p | q deleted; confirmed in rounds 134/136) has the
diagonal part (1/2π²)Σ_n Λ(n)²|χ(n)|² n^{−1} cos(r ln n). Removing p deletes
(1/2π²)Σ_m (ln²p/p^m) cos(m r ln p). With ∫ ln|r| cos(ωr) dr = −π/|ω| and density L/2π:

  Δμ_p = (2π/L)·(1/2π²)·Σ_m (ln²p/p^m)·π/(m ln p) = + ln p · ln(p/(p−1)) / L.

The same number follows from the Euler part alone, E_f[A_p] = −ln p ln(p/(p−1))/L for a good prime (Landau–Gonek),
and the GHK zeros part Z_X is blind to p < X (∫Ci(|y|)cos(ωy)dy = 0 for ω < 1): the two routes agree. A bad prime
*raises* the mean of log|L′(ρ)| by

  +0.4805/L (p = 2), +0.4454/L (3), +0.3591/L (5), +0.3000/L (7), +0.2285/L (11), +0.2053/L (13),

i.e. +0.048, +0.044, +0.035, +0.029, +0.022, +0.020 at L ≈ 10. Order: ln p·ln(p/(p−1))/L ≈ ln p/(pL), not
ln p/(√p L): the p^{−1/2} of the Euler factor multiplies the p^{−1/2} of the Landau–Gonek correlation. The shift does not
depend on χ(p) = ±1 for good primes (it uses χ(p)² only). **Status: conjecture-level (ratios R₂), derivation ours,
first order.** Second order: the off-diagonal (e^{−irL}) part of ΔR₂ has a regular envelope; its leading contribution
is −ln²p/((p−1)L²), i.e. −10 % (p = 2), −13.5 % (3), −18 % (5) of the first-order shift at L = 10 (heuristic
estimate). Theory allowance on the amplitude: 0.2.

### 3.3 The local-density-modulation (LDM) model (heuristic, ours)

At a zero with phase φ the Euler factor contributes A_p(φ), and the zeros part U sees a local density L(1 + ε_p(φ)).
Take U's moments to scale like the HKO N^{k(k+2)} with the local density: E[e^{tU} | φ] ∝ e^{(t + t²/4)ε_p(φ)}. With
the zero weight (1 + ε_p) the per-prime cumulant generating function becomes, to first order in 1/L,

  K_p(t) = ln ₂F₁(t/2, t/2; 1; 1/p) − (2 ln p/L)·(1 + t/2)²·G_p(t),  G_p(t) = E_u[|1−w|^{−t} Re(w/(1−w))]/E_u[|1−w|^{−t}],

and a bad prime removes K_p. Three exact checks at first order:

1. the mean: K′_p(0) = −(2 ln p/L)E_u[A_p Re(w/(1−w))] = −ln p ln(p/(p−1))/L, the Hadamard result of 3.2;
2. k = 1: the bad-prime ratio becomes (1 − 1/p)(1 + 8 ln p/((p−1)L)) since G_p(2) = 1/(p−1). The exact recipe
   (section 4) has (ratio − (1−1/p))·L/(ln p (1−1/p)) = 7.81 at L = 50 → 8 (p = 2) and 3.84 → 4 (p = 3): **the 1/L
   coefficient agrees exactly**;
3. k = −1: (1 + k)² = 0, no first-order correction: the bad-prime factor of Σ1/|L′(ρ)|² is p/(p+1) at
   O(1/L) (consistent with Pearce-Crump's constant and with E_CUE|Λ′|^{−2} = 1/N exactly).

Second order is not controlled. At L = 10 the first-order term overshoots the exact k = 1 correction by 15 % (p = 2) and
25 % (p = 3) (recipe 0.7415 / 0.9013 against LDM 0.7773 / 0.9596). That sets the theory allowance of every LDM 1/L
term: ±0.35 of the term.

Per bad prime, per unit (ln p)/L, the LDM terms of the cumulants are (removed with the prime; `dev/theory.py`):

| p | Δκ₁ | Δκ₂ | Δκ₃ | Δκ₄ |
|---|---|---|---|---|
| 2 | −0.6931 | −1.8667 | −2.5298 | −1.0791 |
| 3 | −0.4055 | −0.9753 | −1.0606 | −0.1826 |
| 5 | −0.2231 | −0.4961 | −0.4606 | −0.0268 |

so the change of the variance by a bad prime 2 is −0.2911 + 1.8667·0.693/L = −0.291 + 0.129 (L = 10): the leading
order is cut almost in half at these heights. Moment factors (bad/good) r_p(k)[1 + c_p(k)/L]:

| p | k = −3/2 | k = −1 | k = −½ | k = ½ | k = 1 |
|---|---|---|---|---|---|
| 2 | 0.4628 (c −0.128) | 0.6667 (0) | 0.8853 (−0.083) | 0.8472 (+1.425) | 0.5000 (+5.545) |
| 3 | 0.5663 (−0.162) | 0.7500 (0) | 0.9215 (−0.088) | 0.9059 (+1.173) | 0.6667 (+4.394) |
| 5 | 0.6870 (−0.170) | 0.8333 (0) | 0.9518 (−0.079) | 0.9465 (+0.880) | 0.8000 (+3.219) |

### 3.4 Within one L-function (heuristic, ours; follows from 3.3)

E[Y | γ] = const + Σ_{p∤q}[A_p(γ) + ε_p(γ)] + (no term for p | q). So in a pooled regression of Y − ln L on A_p and
E_p := ε_p over the zeros of each L-function: coefficient of A_p = 1, of E_p = 1 (0 if the zeros part were independent
of the phase, "P-only"), and for p | q no dependence on the phase γ ln p at all (cos, sin coefficients 0; if the orbit
were kept, the squared amplitude at frequency ln p would equal that of a good prime, p^{−1}(1 − 2 ln p/L)²).

## 4. The k = 1 moment: all lower-order terms (conjecture: ratios / recipe, Conrey–Snaith §7.1; adaptation ours)

Conrey–Snaith derive Σ_{γ<T}|ζ′(ρ)|² = ∫₀^T m₁(t)dt with all lower-order terms (Milinovich verified them by Gonek's
method). For L(s, χ) the arithmetic sums run over (hmn, q) = 1, so ζ(1+s) → ζ_q(1+s) = ζ(1+s)∏_{p|q}(1 − p^{−1−s}) and
t/2π → qt/2π. Writing ζ_q(1+s) = c_q(1/s + c₀ + c₁s + c₂s² + c₃s³ + …), c_q = ∏(1 − 1/p):

  2π m₁/c_q = L⁴/12 + (2/3)c₀L³ + (c₀² + 2c₁)L² + (−2c₀³ + 10c₀c₁ − 2c₂)L + 2c₀⁴ − 12c₀²c₁ + 14c₁² + 16c₀c₂ − 20c₃.

Checks (computed): (i) q = 1 reproduces Conrey–Snaith's polynomial (17.880885063 at L = 5); (ii) an independent
evaluation of the three-term recipe by Cauchy integrals on the torus |β| = r, |γ| = 2r, |δ| = 3r with the true ζ_q
(`dev/recipe_k1b.py`) gives the same bad-prime ratios (p = 2: 0.741501, p = 3: 0.901260 at L = 10); (iii) ζ data: 813
zeros below 1200 give Σ|ζ′(ρ)|²/∫m₁ = 0.9989 on (100, 1200] (1.014, 0.997, 0.996 on the three thirds). A bad prime
changes the k = 1 moment far from its leading factor 1 − 1/p at these heights: ratio 0.74 for p = 2 (leading 0.50),
0.90 for 3 (0.67), and > 1 for p ≥ 13 (D = −403 = −13·31: 1.026), because the 1/L term 8 ln p/((p−1)L) beats 1/p.

## 5. Answers to the four questions

**(a) Variance.** V = ½⟨ln L⟩ + c(L) − Σ_{p|D}[½Li₂(1/p) + Δκ₂(p)·(ln p)⟨1/L⟩] with Δκ₂(p) the LDM coefficients of
3.3 (−1.8667, −0.9753, −0.4961 for p = 2, 3, 5). The per-prime difference at leading order is ½Li₂(1/p) (conjecture); at L ≈ 10 the LDM first-order term reduces it
to 0.165 (p = 2), 0.078 (3), 0.027 (5), 0.011 (7) (heuristic, theory allowance ±35 % of the 1/L term). c(L) =
−0.6788 (CUE) + c_A (Σ_p ½[Li₂(1/p) + ln(1 − 1/p)] = −0.088) + good-prime 1/L terms whose sum does not converge: **not
predicted** at L ≈ 10. N_eff: not determined by the theory at this order (section 1).

**(b) Mean.** Yes: a missing prime shifts the mean of log|L′(ρ)| up by ln p·ln(p/(p−1))/L (+0.048 for p = 2 at L = 10),
via the Ford–Zaharescu/Landau–Gonek correlation between γ ln p and being a zero; its size is ln p/(pL), not
ln p/(√p L). New prime-by-prime signature, first order exact given the ratios two-point function.

**(c) Moments.** Bad-prime factor r_p(k) = 1/₂F₁(k,k;1;1/p) (verified: p/(p+1) at k = −1, 1 − 1/p at k = 1) times
[1 + 2 ln p (1+k)² G_p(2k)/L] (LDM). k = −1: no first-order correction, so Pearce-Crump's factor p/(p+1) is predicted
to hold already at L ≈ 10. k = 1: the exact recipe replaces the LDM. Level: predicted only for k = 1 (recipe);
for k = −1 the conjecture's constant at leading order (theory allowance unknown); for k = ±½ not predicted.
Estimator: per-L-function sums, quasi-Poisson ratio regression E S_i = P_i e^{X_iβ}, bootstrap over L-functions. The
k = −1 tail is not reached in practice: on B (149,523 zeros) no |L′(ρ)| below 0.2 (the y³ law sets in further out),
so a truncation at 0.2 is registered but inactive; m-out-of-n considerations are moot at this size, and k = −½ (finite
variance) is read alongside.

**(d) Shape.** A bad prime changes κ₄ by +0.037 (p = 2, leading; the Euler factor has κ₄ < 0) and +0.111 with the LDM
term; κ₃ by −0.142 (leading) → +0.032 (LDM): the 1/L term reverses the sign of the κ₃ change at L ≈ 10. The left tail
P(|L′| < y) ~ C_q y³ carries the arithmetic factor a_{−3/2}(q): a bad prime 2 multiplies its amplitude by 0.463 at
leading order. So the close-pair statistic *is* affected by bad primes, not only through R₂: through the Euler factor
multiplying the close pair's product. The y³ regime is below |L′| = 0.2 and not reached by these sets: no reading.

## 6. What is firm and what is not

Firm (parameter-free, per-prime differences or exact recipe): the k = 1 moment with all lower-order terms; the
Pearce-Crump factor p/(p+1); the per-prime mean shift at first order; the presence/absence of the Euler factor within
each L-function (A_p coefficient 1, bad-prime phase coefficient 0). Model-dependent at the 1/L level (allowance ±35 %):
the per-prime variance, κ₃, κ₄, k = ±½ and the E_p coefficient. Not predicted: absolute levels of mean and variance,
N_eff, the k = −1 level beyond leading order.

## 7. Doubts

* The LDM model is a heuristic; it passes three first-order checks but its second order is wrong by 15–25 % at k = 1,
  L = 10; the variance reading may land between the LDM (γ₁ = 1) and 0.65.
* The mean-shift second order (−10 to −18 %) is an asymptotic estimate of the off-diagonal term, not a computation.
* The recipe's error term (O(t^{−1/2+ε}) in Conrey–Snaith) is unquantified at t = 100–600; the ζ check bounds it to
  about 1 % over (0, 400] and 0.1 % over (100, 1200].
* E_p is 91 % correlated with A_p for p = 2 (E_u[A R] / √(E A² E R²)); the within-function coefficients are read jointly.
* The k = −1 bootstrap is formally unsafe for an infinite-variance variable; in these sample sizes the tail is not
  reached (no |L′| < 0.2 on B), but one extreme close pair in a fresh set could move a per-prime coefficient by
  several se: the reading reports the largest single term.

## 8. Noise and power (measured on B, noise only; `noise_B_out.txt`)

Bootstrap se of the per-prime coefficients on B (184 L-functions; p = 2, 3, 5): mean 0.0013, 0.0013, 0.0015; variance
0.0027, 0.0026, 0.0031; k = 1 log-sum 0.0012, 0.0015, 0.0022 (level 0.0007); k = −1 log-sum 0.034, 0.026, 0.035
(level 0.015); within-function A_2 0.0044, E_2 0.025, bad-prime cos/sin 0.002. The per-L-function scatter of the k = 1
sum around the recipe is only 0.9 % (Σ|L′(ρ)|² over a band is as rigid as a mean value), and the mean of log|L′| is
three times less noisy than independent zeros would give (synthetic check with iid noise: se 0.0038 instead of 0.0013).
Scaled to the 548 fresh L-functions (×0.58), every primary reading separates its hypothesis from "no conductor
dependence" by 26 (k = −1) to 420 (k = 1) expected se, and from "half the size" by 13 to 210 (round_draft.txt).
Validation: PARI's L′ agrees with an mpmath Hurwitz-zeta computation to 4·10⁻¹⁴ in log|L′| (`check_deriv_out.txt`);
|L(ρ)|/|L′(ρ)| ≤ 6·10⁻¹⁴ at every zero of B; synthetic values with planted LDM effects on B's zero heights are recovered
by `lprime_test.py` (`synth_B_out.txt`: within-function A_2 1.03 ± 0.01, E_2 1.14 ± 0.06, γ_1 0.85 ± 0.02, where the
synthetic truth contains second-order terms the first-order formula lacks).
