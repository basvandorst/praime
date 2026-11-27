# Entire L-functions of real archimedean type: a zero with |γ| < 9.225 under GRH and Ramanujan at p = 2 (ball-arithmetic certificate, 27 November 2026)

**In short.** Round 175 registered the claim: under GRH, and Ramanujan at the primes p ≤ 31, every entire L-function of real archimedean type has a zero 1/2 + iγ with |γ| < 9.23, possibly < 9.225. The Ramanujan Δ function has its first zero at 9.2224. The float dual certificates are turned into one explicit test function G each, with exact rational data, and checked in ball arithmetic.

1. **Proved (Lemma 1).** The certificate is a function G with Ĝ(t) ≥ 0 for |t| ≥ T. This holds by construction, from 9 + 7 (T = 9.225) or 8 + 6 (T = 9.23) rational vectors.
2. **Cited, hypotheses checked (Lemma 2).** Weil's explicit formula holds for G. The archimedean term in the u-form is proved here. Our normalisation is checked against 982 zeros of L(s, Δ) and 630 zeros of L(s, χ₋₃): it agrees to 4·10⁻¹¹ or better.
3. **Proved (Lemmas 0, 3).** Assume no zero has |γ| < T. Then 0 ≤ W(G)/d ≤ B := sup_κ F(κ) + Σ_{p≤31} max P_p. This uses c₀ = G(0) < 0, N ≥ 1, the averaging over the degree, and the maximum principle on the disk |α| ≤ R_p.
4. **Proved (Lemma 4).** The sup over all κ ≥ 0 reduces to [0, 12]. Beyond 11.59 (T = 9.225) or 11.66 (T = 9.23), F′ < 0. On [0, 12] a grid plus a bound on F″ suffices.
5. **Computer-assisted (python-flint arb, one code; 80 checks, 0 failed).** For T = 9.225: c₀ ≤ −0.027618, sup F ≤ −2.1015655e−3, and B ≤ **−2.0525538e−4**. For T = 9.23: B ≤ −9.2840941e−4.

**Result.** The theorem holds with T₀ = 9.225. For T₀ = 9.225, Ramanujan is used only at p = 2: |α_j(2)| ≤ 1 (1.0365 suffices if Ramanujan holds at the other p ≤ 31). At 3 ≤ p ≤ 31 only the bound |α_j(p)| ≤ p is used, and every L-function of the class has it. The statement also holds for real κ_j > −1/2.

## 1. The class and the theorem

The degree d ≥ 1 is an integer. The axioms:

- **(A1) Euler product.** L(s) = Σ a(n) n^{−s} = ∏_p ∏_{j=1}^{d} (1 − α_j(p) p^{−s})^{−1}, absolutely convergent for Re s > 1. So −L′/L(s) = Σ c(n) n^{−s}, with c(p^m) = log p · Σ_j α_j(p)^m and c(n) = 0 when n is not a prime power. Some α_j(p) may be 0.
- **(A2) Gamma factor.** Λ(s) = N^{s/2} ∏_{j=1}^{d} Γ_R(s + κ_j) L(s), with Γ_R(s) = π^{−s/2} Γ(s/2), N ≥ 1 real, and κ_j real, κ_j ≥ 0. Section 6 extends this to κ_j > −1/2.
- **(A3) Functional equation.** Λ is entire of order 1, and Λ(s) = ε \overline{Λ(1 − s̄)} with |ε| = 1.
- **(A4) Growth.** L has polynomial growth in vertical strips. This follows from (A1)–(A3) by Phragmén–Lindelöf, but it is listed so that the cited explicit formula applies verbatim.
- **(GRH)** Every zero of Λ has real part 1/2.

**Theorem (computer-assisted).** Let L satisfy (A1)–(A4) and GRH, and let |α_j(2)| ≤ 1 for all j. Then Λ(1/2 + iγ) = 0 for some real γ with |γ| < 9.225.

- If Ramanujan holds at every other p ≤ 31, the bound at 2 can be weakened to |α_j(2)| ≤ 1.0365 (computed, certified at 1.036574 with all other radii 1).
- With T₀ = 9.23 the theorem holds under GRH together with either of two hypotheses:
  - Ramanujan at p = 2 and 5;
  - Ramanujan at p = 2, and |α_j(p)| ≤ √p at 3 ≤ p ≤ 31.
- Both statements hold with κ_j real and κ_j > −1/2.

**Hypotheses, and where each is used.**

| hypothesis | used in |
|---|---|
| (A1) Euler product with exactly d roots at each p (zeros allowed) | Lemma 3: averaging over j at p = 2, …, 31; Lemma 0: \|α_j(p)\| ≤ p |
| (A2) the number of Γ_R factors equals d; κ_j real | Lemma 3: averaging (1/d) Σ_j F(κ_j) ≤ sup F; complex κ would need the sup over a half-plane (class C, not covered) |
| (A2) κ_j > −1/2 (≥ 0 in the main statement) | Lemma 2: Gauss's integral and the contour shift; Lemma 4: the range of the sup |
| (A2) N ≥ 1 | Lemma 3: (log N/d) c₀ ≤ 0, because c₀ < 0 |
| (A3), (A4), entire (no pole) | Lemma 2 |
| GRH | proof of the Theorem: every γ is real, so Ĝ(γ) ≥ 0 when \|γ\| ≥ T |
| \|α_j(2)\| ≤ 1 | Lemma 3 at p = 2; the certificate fails with the trivial \|α_j(2)\| ≤ 2 (section 5) |

**Not used:**

- Ramanujan at any p ≥ 3 (for T₀ = 9.225);
- any coefficient c(n) with n ≥ 32;
- the root number;
- self-duality: L may have complex coefficients, and only Re c(n) enters;
- whether L is automorphic.

**Lemma 0 (proved).** (A1) implies |α_j(p)| ≤ p.

*Proof.* Σ_m a(p^m) X^m converges absolutely for |X| < 1/p, by absolute convergence of the Dirichlet series at real s > 1. It is the Taylor series of the rational function ∏_j (1 − α_j(p) X)^{−1}. This function has a pole at every 1/α_j(p), because there is no numerator to cancel it. The radius of convergence is the distance to the nearest pole, so |1/α_j(p)| ≥ 1/p. □

## 2. The test function (Lemma 1, proved)

**Set-up.**

- a = (log 32)/2, so 2a = log 32 = log X with X = 32.
- ω_k = kπ/a for |k| ≤ 23; a vector has 47 entries.
- T is the exact rational 9225/1000 or 923/100.
- The certificate of record consists of rational vectors s_i, r_j ∈ ℚ⁴⁷ with Σ_k (−1)^k r_{jk} = 0 exactly, stored in gap_cert_rational_T*.json:
  - 9 vectors s_i and 7 vectors r_j for T = 9.225;
  - 8 vectors s_i and 6 vectors r_j for T = 9.23.
- How they were made: s_i and r_j are the float eigen-factors √λ v of P_H and of Z P_G Zᵀ, converted exactly (dyadic). Each r_j is then projected exactly onto Σ_k (−1)^k r_k = 0. The matrices Σ s_i s_iᵀ and Σ r_j r_jᵀ reproduce the npz matrices to 1.6·10⁻¹⁷.
- Once the vectors exist, the rounding plays no further role: any real vectors with the exact side condition define an admissible G.

**Definitions.**

- f_i(u) = Σ_k s_{ik} e^{iω_k u} on [−a, a], and 0 outside.
- ψ_j is defined the same way from r_j. ψ_j′ is its derivative on (−a, a), extended by 0.
- F ⋆ F̃ (u) = ∫ F(x + u) \overline{F(x)} dx.
- g = Σ_i f_i ⋆ f̃_i + Σ_j (ψ_j′ ⋆ ψ̃_j′ − T² ψ_j ⋆ ψ̃_j), and G = Re g.

**Lemma 1.**

- (a) G is real, even, continuous and Lipschitz, with support in [−2a, 2a] and G(±2a) = 0. On [0, 2a],
  `G(u) = (2a − u) Σ_{j=0}^{23} α_j cos ω_j u + Σ_{k=1}^{23} b_k sin ω_k u.`
  Here D = Σ s_i s_iᵀ + Ω R Ω − T² R, with R = Σ r_j r_jᵀ and Ω = diag(ω_k). The coefficients are α_0 = D_00, α_j = D_jj + D_{−j,−j}, and b_k = β_k − β_{−k}, with β_k = (2a/π) Σ_{l≠k} D_lk (−1)^{l−k}/(l − k). On each of [−2a, 0] and [0, 2a], G is the restriction of an entire function.
- (b) Ĝ(t) := ∫ G(u) e^{itu} du = (ĝ(t) + ĝ(−t))/2, where ĝ(t) = Σ_i f̂_i(t)² + (t² − T²) Σ_j ψ̂_j(t)². Here f̂(t) = 2 sin(at) Σ_k (−1)^k s_k/(t + ω_k) is real. Hence **Ĝ(t) ≥ 0 for |t| ≥ T.**
- (c) Ĝ is entire, and |Ĝ(t + iy)| ≪ e^{2a|y|} min(1, |t|^{−2}).

*Proof.*

- (a) For 0 ≤ u ≤ 2a, ∫ e^{iω_l(x+u)} e^{−iω_k x} dx over the overlap [−a, a − u] equals:
  - (2a − u) e^{iω_k u} when l = k;
  - a(−1)^{l−k}(e^{iω_k u} − e^{iω_l u})/(iπ(l − k)) when l ≠ k.

  The derivative ψ_j′ has coefficients iω_k r_{jk}, and the factor i cancels in the product. So g(u) = Σ_{l,k} D_lk (the integral above). The real part, with D symmetric, gives the closed form.
  - G(2a) = 0·A + Σ b_k sin 2πk = 0.
  - G(−u) = Re g(−u) = Re \overline{g(u)} = G(u), so G is even.
  - Lipschitz: an autocorrelation of bounded functions of bounded variation.
- (b) ĝ(t) = |f̂(t)|²-type terms by Fubini. Since f_i(−u) = \overline{f_i(u)}, f̂_i is real. The ψ′ terms use integration by parts: ∫_{−a}^{a} ψ′ e^{itu} = [ψ e^{itu}]_{−a}^{a} − it ψ̂(t). The boundary term vanishes because ψ_j(±a) = Σ_k (−1)^k r_{jk} = 0 exactly; this is the only place the side condition is used. Taking the real part of g symmetrises ĝ.
- (c) Integrate by parts twice on [−2a, 0] and on [0, 2a]. The first boundary terms cancel, by continuity at 0 and G(±2a) = 0. The second ones and the remaining integral are O(e^{2a|y|}/|z|²). □

## 3. The explicit formula (Lemma 2)

**Lemma 2.** Let L satisfy (A1)–(A4), with real κ_j > −1/2. No GRH is needed. Then

`Σ_ρ Ĝ(γ_ρ) = W(G) := G(0) log N + Σ_j A_{κ_j}(G) − 2 Σ_{n<32} Re c(n) n^{−1/2} G(log n),`

where ρ = 1/2 + iγ_ρ runs over the zeros of Λ with multiplicity, and the sum converges absolutely. The archimedean term is

`A_κ(G) = (1/2π) ∫ Ĝ(t) [−log π + Re ψ((1/2 + κ + it)/2)] dt`
`       = −(γ_E + log π) G(0) − G(0) log(1 − 2^{−10}) + ∫_0^{2a} [2e^{−2u} G(0) − 2e^{−(1/2+κ)u} G(u)] / (1 − e^{−2u}) du.`

**Proof.**

- *First line: cited.* This is Weil's explicit formula.
  - Source 1: Bober, Conrey, Farmer, Fujii, Koutsoliotas, Lemurell, Rubinstein, Yoshida, *The highest lowest zero of general L-functions*, J. Number Theory 147 (2015), Lemma 4.1 (arXiv 1211.5996). It requires a test function holomorphic in |Im| < 1/2 + δ with f ≪ min(1, |t|^{−1−ε}), L entire with the functional equation and polynomial growth, and absolute convergence of the n-sum. Lemma 1(c) gives the first hypothesis, (A3)–(A4) the second, and the n-sum is finite because G(log n) = 0 for n ≥ 32.
  - Source 2: Iwaniec–Kowalski, Thm 5.12, same hypotheses; not reopened.
  - Neither source was read verbatim (see section 8). The normalisation used here is therefore not taken from print; it is checked numerically below.
- *Reduction to smooth test functions (proved).* The cited input can be limited to even real C_c^∞ test functions, the standard case.
  - Let η ∈ C_c^∞[−1, 1] be even, η ≥ 0, ∫η = 1, and set G_ε = G ∗ η_ε with η_ε(x) = η(x/ε)/ε. Then G_ε ∈ C_c^∞ is even and real, and Ĝ_ε(z) = Ĝ(z) η̂(εz) with |η̂(εz)| ≤ e^{ε|Im z|}.
  - Let ε → 0 in the formula for G_ε, with ε < log(33/32) so that only n ≤ 32 enter:
    - zero side: dominated convergence, using Σ_ρ |Ĝ(γ_ρ)| < ∞. This holds by Lemma 1(c), |Im γ_ρ| ≤ 1/2, and O(log T) zeros per unit height (classical, e.g. Iwaniec–Kowalski Prop. 5.7; recalled).
    - G_ε(0) → G(0) and G_ε(log n) → G(log n), by continuity; G(log 32) = 0.
    - archimedean integral: dominated convergence, since |η̂| ≤ 1 on ℝ and Ĝ(t) log(2 + |t|) ∈ L¹.
  - The n-sum is Σ n^{−1/2}[c(n) G(log n) + \overline{c(n)} G(−log n)] = 2 Σ Re c(n) n^{−1/2} G(log n), because G is even and real.
- *Second line: proved here.*
  - For Re z > 0, Gauss: ψ(z) = −γ_E + ∫_0^∞ (e^{−v} − e^{−zv})/(1 − e^{−v}) dv (DLMF 5.9.16).
  - Put z = (1/2 + κ + it)/2 and v = 2u. Re z > 0 is exactly κ > −1/2.
  - On u ≥ δ the double integral converges absolutely, because Ĝ ∈ L¹ (Lemma 1(c)). Exchange the integrals there and use (1/2π)∫Ĝ = G(0) and (1/2π)∫Ĝ(t) cos(tu) dt = G(u) (Fourier inversion: G continuous, Ĝ ∈ L¹).
  - The piece u < δ is bounded by C(1 + log(2 + |t|)) uniformly in δ ≤ 1 and tends to 0 pointwise. Dominated convergence applies because Ĝ(t) log(2 + |t|) ∈ L¹.
  - On the u-side the integrand is bounded at u = 0 (G is Lipschitz), so δ → 0 gives the full integral.
  - Since G = 0 for u ≥ 2a, the rest of ∫ 2e^{−2u}G(0)/(1 − e^{−2u}) is −G(0) log(1 − e^{−4a}), and e^{−4a} = 2^{−10}.

**Where each constant comes from.**

- −log π: the factor π^{−s/2} in Γ_R.
- −γ_E: Gauss's integral.
- log(1 − 2^{−10}): the part of the g(0)-integral beyond the support.
- log N: the factor N^{s/2}.
- n^{−1/2}: evaluation on the critical line.
- The factor 2 in the n-sum: n and its mirror −log n.

**Numerical check of the normalisation (computed; section H of the script).** The zeros were computed by PARI lfunzeros up to height 800:

| | zeros | default vs divz = 32 | Riemann–von Mangoldt main term | zero side incl. tail | W(G) from c₀, F, primes | difference |
|---|---|---|---|---|---|---|
| L(s, Δ) (κ = 11/2, 13/2) | 982 | same | 982.3 | −2.2143054222e−3 | −2.2143054076e−3 | 1.5e−11 |
| L(s, χ₋₃) (κ = 1, N = 3) | 630 | same | 629.8 | −7.0747789603e−2 | −7.0747789647e−2 | 4.3e−11 |

The table is for the T = 9.23 function. With the T = 9.225 function the differences are 1.5e−12 (Δ) and 4.5e−12 (χ₋₃).

The tail above 800 is estimated, so this check is not exact:

- the smooth density part is 6.7e−8 (Δ, T = 9.23);
- the resonance of the prime power 32 is −1.4e−9. It is resonant because log 32 = 2a matches the factor sin²(at) in Ĝ = 4 sin²(at) Q(t). Without this term the agreement is worse.

## 4. The reduction (Lemma 3, proved)

**Lemma 3.** Assume (A1)–(A4) with real κ_j > −1/2, and |α_j(p)| ≤ R_p for p ≤ 31. Then

`W(G)/d ≤ (log N/d) c₀ + sup_{κ>−1/2} F(κ) + Σ_{p≤31} M_p(R_p),`

with the notation:

- c₀ = G(0) and F(κ) = A_κ(G);
- M_p(R) = max_θ Σ_m w_{p,m} R^m cos mθ;
- w_{p,m} = −2 log p · p^{−m/2} G(m log p), for p^m < 32.

*Proof.*

- The n-sum equals Σ_{p≤31} Σ_j Re P̃_p(α_j(p)), where P̃_p(z) = Σ_m w_{p,m} z^m has real coefficients.
- Re P̃_p is harmonic. By the maximum principle, its maximum over |z| ≤ R lies on |z| = R, where Re P̃_p(Re^{iθ}) = Σ_m w_{p,m} R^m cos mθ. So each of the d terms is at most M_p(R_p). A root α_j(p) = 0 contributes 0 ≤ M_p, because the mean over the circle is 0.
- (1/d) Σ_j F(κ_j) ≤ sup F.
- Divide by d. □

**Degree-free.** d never enters: the bound is on W/d.

**Proof of the Theorem.** Suppose Λ has no zero with |γ| < T. By GRH every γ_ρ is real with |γ_ρ| ≥ T, so Ĝ(γ_ρ) ≥ 0 by Lemma 1(b). Then 0 ≤ Σ_ρ Ĝ(γ_ρ) = W(G) by Lemma 2. Lemma 3 with c₀ < 0 and log N ≥ 0 gives W(G)/d ≤ B < 0, a contradiction. □

## 5. The certified numbers (computer-assisted)

The code is [validate_gap_cert.py](validate_gap_cert.py) (python-flint 0.9.0, arb, 128 bits), with output in [validate_gap_cert_out.txt](validate_gap_cert_out.txt).

**What is computed, and how.**

- α_j and b_k: exact rationals combined with π and log 2 as balls.
- G(u): evaluated from the closed form. Two independent float routes agree: direct autocorrelation at u = 0, log 2, log 3, log 5, 2.5 (to 10⁻¹²), and Fourier inversion at log 2, log 3.
- **c₀** = 2a Σ_j α_j, a ball.
- **max P_p**: P_p is converted to a polynomial in x = cos θ with exact Chebyshev coefficients, and branch and bound on [−1, 1] uses the mean-value form P(c) + P′(X)(x − c) on dyadic boxes.
- **F(κ):** acb.integral of an integrand with no removable singularity: (1 − e^{−2u}) = 2u e^{−u} sinh(u)/u, and sinh(x)/x = sinc(ix). The integrand is analytic in |Im u| < π.
- **Lemma 4 (proved).**
  - F′(κ) = ∫_0^{2a} e^{−κu} h(u) du and F″(κ) = −∫ u e^{−κu} h(u) du, where h(u) = 2u e^{−u/2} G(u)/(1 − e^{−2u}) and h(0) = c₀.
  - (a) Tail. If h ≤ −η on [0, u₁] and |h| ≤ M on [u₁, 2a], then F′ < 0 for κ > log(1 + M/η)/u₁.
  - (b) Grid. On [κ_l, κ_r], F ≤ max(F(κ_l), F(κ_r)) + M₂(κ_l)(κ_r − κ_l)²/8, with M₂(κ_l) = ∫ u e^{−κ_l u} |h|.
  - Since F ≈ c₀ log κ, F → −∞ as κ → ∞.
- Bounds on h come from arb balls on 4096 pieces of [0, 2a]. The Riemann sums are added in floats with a relative margin of 10⁻⁹.
- The grid has step 1/32 on [0, 12], with bisection where the bound comes within 10⁻¹⁰ of the best value.

| | T = 9.23 | T = 9.225 |
|---|---|---|
| c₀ = G(0) | [−0.0287445310727, −0.0287445310726] | [−0.0276185373065, −0.0276185373064] |
| max P₂ | ≤ 1.7319150e−3 | ≤ 1.8952439e−3 |
| max P₃ | ≤ 1.2124706e−5 | ≤ 1.0662782e−6 |
| max P₅ | ≤ 5.2885328e−5 | ≤ 1.91e−16 |
| max P_p, 7 ≤ p ≤ 29 | each ≤ 3.8e−14 | each ≤ 5.1e−16 |
| max P₃₁ | ≤ 4.7190452e−8 | ≤ 6.1925120e−13 |
| tail: F decreasing from | 11.660 | 11.589 |
| sup_{κ≥0} F | [−2.725381640e−3, −2.725381577e−3] at κ ≈ 5.8738 | [−2.101565611e−3, −2.101565549e−3] at κ ≈ 5.8738 |
| sup_{−1/2≤κ≤0} F | ≤ −0.04727967 (at κ = 0) | ≤ −0.04537750 (at κ = 0) |
| **B** | **[−9.2840948e−4, −9.2840941e−4]** | **[−2.0525545e−4, −2.0525538e−4]** |

**Float cross-check (cert.py, cert_real_out.txt).**

- c₀ agrees to 10⁻¹².
- sup F: cert.py's grid value −2.7253827e−3 at κ = 5.875 lies 1.0e−9 below the certified sup, which is attained at 5.8738.
- The per-prime maxima agree to all printed digits for p = 2, 3, 5. The values below 10⁻¹² (p ≥ 7) differ in the second or third digit, from float cancellation in cert.py; this is irrelevant at the size of B.
- B: −9.284105e−4 and −2.052565e−4, against the certified [−9.2840948e−4, −9.2840941e−4] and [−2.0525545e−4, −2.0525538e−4].

**Second route for F (float).** The Fourier side (1/π)∫_0^∞ Ĝ(t)[−log π + Re ψ((1/2+κ+it)/2)] dt, with the tail beyond t = 1000 taken from the exact rational form Ĝ = 4 sin²(at) Q(t), agrees with the certified u-integral at κ = 0, 1, 5.5, 5.8738, 6.5 and 20 to better than 10⁻¹¹ (observed ≤ 2·10⁻¹⁵). c₀ and G(log 2), G(log 3) by Fourier inversion agree to 10⁻¹⁵.

**Which local hypotheses are needed (certified; B with radius R_p).**

| Ramanujan (R = 1) at | other p ≤ 31 | T = 9.23 | T = 9.225 |
|---|---|---|---|
| all p ≤ 31 | — | −9.284e−4 | −2.053e−4 |
| 2, 5 | R = p (Lemma 0) | −6.282e−4 | −1.813e−4 |
| 2 | R = p | **+4.06e−4 (fails)** | **−1.813e−4** |
| 2 | R = √p | −7.063e−4 | −2.012e−4 |
| none | R = p | +1.28e−2 (fails) | +1.29e−2 (fails) |

- Uniform relaxation |α_j(p)| ≤ R at every p ≤ 31: enough for R ≤ 1.155999 (T = 9.23) and R ≤ 1.036555 (T = 9.225).
- Relaxing p = 2 alone: |α_j(2)| ≤ 1.159093 and ≤ 1.036574.
- Relaxing p = 3 alone: R = 3 is enough. Relaxing p = 5 alone: 4.747516 (T = 9.23) and 5 (T = 9.225).
- So the prime 2 is where arithmetic enters. The certificate fails with the trivial bound there.

## 6. Extending the parameters, and what enters

- **κ below 0 (proved and certified).** On [−1/2, 0], F stays below its sup over [0, ∞): the maximum there is at κ = 0, and F(−1/2) = −0.05225 (T = 9.23). So the Theorem holds for real κ_j > −1/2. This range includes exceptional-eigenvalue parameters, which are real and lie in (−1/2, 1/2).
  - At κ = −1/2 something new enters. The pole of Γ_R(s + κ) at s = −κ reaches the critical line, and Gauss's integral needs Re z > 0, which fails.
  - For κ < −1/2, (A3) forces L(−κ) = 0 at a point right of the line, and the explicit formula gains a residue term.
- **Ramanujan radius.** See section 5. At p = 2 the radius can grow to 1.0365 (T = 9.225, radius 1 at the other primes). Beyond that this certificate fails: the local factor at 2 is where arithmetic enters.
- **T.** Only the two given certificates were checked. Between 9.2224 and 9.225 nothing is certified, and other T need the LP again. Under GRH for L(s, Δ), Δ belongs to the class: its local roots satisfy |α| = 1 by Deligne, κ = 11/2, 13/2, and N = 1. Then no T₀ ≤ 9.22237… is possible, and 9.225 is within 0.0027 of the optimum. This is conditional on GRH for Δ, which is not proved.
- **Conductor.** B < 0 leaves slack. The bound needs only log N/d > B/|c₀| = −0.0074 (T = 9.225), so N ≥ 1 is not tight.
- **Complex κ (class C).** Not covered: Lemma 3 needs sup F over complex κ, and the round's own pilot puts that bound near 17.

## 7. Consistency with Δ (computed)

- For L(s, Δ) and the T = 9.23 function: W(G)/2 = −1.1072e−3 ≤ B. The zero side is dominated by 2Ĝ(9.22238) = −2.2165e−3: Δ's first zero lies inside (−T, T), where Ĝ < 0. This is the reason Δ does not contradict the certificate.
- Ĝ ≥ 0 on [T, 400], sampled; this holds by construction.
- The minimum of Ĝ on [0, T] is negative.

## 8. What is not proved, sources, and what a referee should attack

**Labels.**

- Proved here: Lemma 0, Lemma 1, the second line of Lemma 2, Lemma 3, Lemma 4, the assembly.
- Cited: the first line of Lemma 2, i.e. Weil's explicit formula for this G.
- Computer-assisted (one code): every number in section 5.

**Not proved or not checked:**

- **Explicit formula.**
  - Bober et al., Lemma 4.1: seen only through an automated summary of the ar5iv HTML page, not read verbatim. The hypotheses quoted above come from that summary, and the exact form of their archimedean term could not be confirmed.
  - Iwaniec–Kowalski Thm 5.12: recalled, not reopened.
  - Our normalisation rests on the derivation in section 3 and on the zero check, to 10⁻¹⁰ at two L-functions.
- **Software.** arb and python-flint are trusted.
  - Float steps with margins: the Riemann sums for M₂, and the tail constants η, M, K₁. Each is a float sum or float max of arb upper bounds, with relative margin 10⁻⁹ or 10⁻¹² and an additive 10⁻¹⁵.
  - The branch and bound accepts a box when its bound is within 10⁻¹³ of the best value. It is rigorous, but not a second code.
- **No second rigorous code.** The float routes are independent, but they are not proofs.
- **Sharpness.** It is conditional on GRH for Δ, and T₀ ∈ (9.2224, 9.225) is open.
- **Literature.** Not found in print: a bound below 14.13 for this class with local hypotheses at p = 2.
  - Searched: two web queries ("highest lowest zero of general L-functions" explicit formula; lowest zero bound with Ramanujan/Δ/9.2224/real archimedean type), and the registration's own check.
  - Miller (Duke 2002, arXiv math/0112196): search snippets only, paper not read. Bober et al.: two automated summaries of the ar5iv HTML; the arXiv PDF could not be read here (no PDF text tool). From those summaries and the registration: Miller's 14.1347 assumes GRH and real archimedean type, and its test functions have support in [−log 2, log 2], so no prime terms enter.

**Attack first:**

1. Lemma 2. The cited formula for C_c^∞ test functions is used in our normalisation, which is checked only numerically. Also check the smoothing limit for the Lipschitz G, which has kinks at 0 and ±2a.
2. Lemma 3's use of exactly d local roots at p = 2, equal to the number of Γ_R factors. The theorem is false without such an axiom.
3. The float-with-margin steps in section 5 (M₂, η, M).
4. The regularised integrand for F. It is checked against the difference quotient at three points and against the Fourier route to 10⁻¹¹.
