# Second moment of the one-point sum over the Petersson family: experiment, mechanism, and where a proof stops

Folder: `scratchpad/second/` (nothing in the repository was changed). Builds on `writeup/code/study/lowzero/ecross/family/family_proof.md`
("FP": Theorems A–D, Lemmas 1–6, §7, §9a) and facts 2052–2055. Code: `smlib.py` (coefficients, Petersson/Kloosterman machinery),
`moment2.py` (M₁₁, M₂ per modulus), `mixed.py` (X_q, C₂, V), `forms_side.py` (eigenform side at small level), `xprofile.py`
(where inside a Kloosterman sum the off-diagonal lives), `analysis.py` (tables), `check_second.py` (assertions, sections a, b, c).
Labels: **proved**, **proved (cited)**, **sketch**, **GAP**, **computed**, **heuristic**.

## 0. Verdict in short

1. **Experiment (computed, exact Petersson sums, no eigenform data; validated against eigenforms at N = 61: M₁₁ to 6 digits, X_q to 10⁻⁸).**
   V := Σ^h|Z_f + ε_f c_f^{(10)}(1/2π)∫H|² at N = 10007: V/T² = 0.095, 0.060, 0.040, 0.027 at T = 6, 12, 20, 30 (k = 4 at T = 6;
   k = 14 at T ≥ 12; k = 4 at T = 12 gives 0.059). V is small and decreasing, **but V/D₁₁ ≈ 0.79–0.84 throughout**: the residual
   has the size of the diagonal D₁₁ ≍ T log(NT). At these heights V/T² → 0 is the trivial decay of the diagonal; what decides is
   item 2. Second level N = 5003: see Table 1 (same picture).
2. **Where the T² term lives (computed + heuristic derivation; the arithmetic identity proved).** The off-diagonal of
   Σ^h|Z_f|² is carried by the moduli c = N·m with 1 ≤ m ≲ T/2π (one modulus per "lag" t − t′ = 2πm of the heights), and inside
   each Kloosterman sum S(p, p′; Nm) = Σ_x e((px + p′x̄)/Nm) by a handful of units only, the **arcs** x ≡ a·(Nm/q) ± 1,
   x̄ ≡ a′·(Nm/q) ± 1 (q | m). The x = ±1 part equals, to 3–4 digits, the same modulus of the *integer* model (Λ replaced by 1),
   whose total is a theorem; the q ≥ 2 arcs multiply it by w(m) = ∏_{p∥m}(1 + (p−2)/(p−1)²)·∏_{p²|m}(1 + 1/(p−1)), and **the mean
   of w(m) is ∏_p(1 + p⁻²) = ζ(2)/ζ(4) = lim Σ^h c_f²** (proved). The q-th arc family carries exactly the q-th term μ²(q)/q² of
   Σ^h c_f² = Σ_q μ²(q)/q²; arcs with different q never pair (proved). Generic pairs (p, p′) do not contribute: the mod-N factor
   S(1, pp′m̄²; N) is a genuine Kloosterman sum and averages out; the main term is two single exponentials of it.
   For Σ^h Z_f² the same arcs carry the T² term at the moduli m ≈ (t + t′)/2π ∈ [T/π, 2T/π].
3. **Proof.** Proved: the exact quadratic-form identities, the one-modulus large sieve (sharp for general coefficients), the
   diagonal, the integer-model second moment, the arithmetic of the arcs (w(m), its mean, orthogonality of arcs). Sketch (GRH):
   the arc main terms. **GAP (the route is closed with known tools): the minor arcs**, i.e. Σ_m |O_m − O_m^{arc}| = o(T²) for the
   bilinear forms Σ_{p,p′}β_pβ̄_{p′}S(p, p′; Nm)J_{k−1}(4π√(pp′)/Nm), p, p′ ≍ NT², 1 ≤ m ≲ T². The best per-modulus bound we can
   prove is the Parseval (large-sieve) bound O((1 + T²/m)·T log(NT)), sharp for general coefficients; summed over m it gives
   O(T³log(NT)log²T), short of o(T²) by a factor ≍ T·log(NT)log²T. The bilinear Kloosterman bounds of
   Fouvry–Kowalski–Michel (Duke 2014, Thm 1.17) and Kowalski–Michel–Sawin (Annals 2017, Thm 1.1) are stated against the trivial
   bound ‖α‖‖β‖(MN)^{1/2} and, in our range (both variables ≍ NT² > N), are weaker than Parseval by N^{1/4}; KMS needs
   MN < q^{5/4} and the second variable in an interval of length < q. GRH for Dirichlet L-functions controls the arcs, not the
   minor arcs (≍ Nm characters per modulus). Observed: the minor arcs of Σ^h|Z_f|² stay O(1) in total (−0.27, +0.23, −0.33 at T = 12, 20, 30;
   0.1–0.4 per small modulus against a large-sieve bound ≍ 10²), consistent with (MA); those of Σ^h Z_f² are ≈ −6 % of its arcs
   at T = 20, 30 (not decided whether O(T) or a small multiple of T²; V does not involve them). Not provable by this route now.
4. **Theorem we can state** (§4.6): the conditional reduction "(MA) + GRH (arcs) ⇒ V = o(T²) ⇒ almost all f satisfy the law",
   with (MA) the explicit minor-arc hypothesis; plus the bounds L²(1 + o(1)) ≤ Σ^h|Z_f|² ≪ T³ log(NT) log²T (lower bound in the ranges of FP Theorems B/C).

## 1. Setting

N prime, k ∈ K = {2, 4, 6, 8, 10, 14} (all forms new), harmonic weights ω_f, Σ^h := Σ_{f∈H_k(N)} ω_f; ε_f = −i^kλ_f(N)√N;
H = W(·/T) ∗ g₁ (W the bump exp(1 − 1/(1 − (2u−3)²)) on (1, 2)); L := (1/2π)∫H; θ, I_u as FP §1.
**Coefficients.** w_{p^j} := −(1/2π)(p^j)^{−½}I_{p^j}ln p; β_{p^j} := w_{p^j} − w_{p^{j+2}} (p ≠ N), β_{N^j} := w_{N^j},
β₁ := −Σ_{p≠N}w_{p²}. Since Λ_f(p^j) = (λ_f(p^j) − λ_f(p^{j−2}))ln p (p ≠ N) and λ_f(N)^j ln N, Proposition 1 of the note gives
**Z_f(H) = Σ_n β_nλ_f(n) + E₀(f)**, |E₀| ≪ N^{−3/4}log(NT) (proved; FP Theorem A(1)). Numerically n runs over prime powers with
t_n := 2π√(n/N) ∈ [0.3T, 3.5T] (k = 4) or [0.3T, 3.0T] (k = 14); the dropped tails change M₁₁ by ≤ 0.05 at T = 12 (§6, doubt 2).
**Quantities.** M₁₁ := Σ^h|Z_f|², M₂ := Σ^h Z_f² (= Σ^h ε_f²Z_f², ε_f² = 1), X_q := Σ^h ε_fλ_f(q²)Z_f,
c_f^{(Q)} := Σ_{q≤Q,(q,N)=1}μ(q)λ_f(q²)/q, C₂(Q) := Σ^h(c_f^{(Q)})², X_Q := Σ_{q≤Q}μ(q)q⁻¹X_q and
**V(Q) := Σ^h|Z_f + ε_f c_f^{(Q)}L|² = M₁₁ + 2L·Re X_Q + L²C₂(Q)** (ε_f, c_f, L real). The conjecture (fact 2055) says
Z_f = −ε_f c_f L + o(T) per form, hence predicts V = o(T²) once Q → ∞ (c_f^{(Q)} → c_f in mean square: Σ^h(c_f − c_f^{(Q)})² ≈
Σ_{q>Q}μ²(q)/q² ≍ Q⁻¹ by Petersson at small indices).
**Petersson (FP Lemma 1).** Σ^hλ_f(n)λ_f(n′) = Δ(n, n′) = δ_{nn′} + Σ_{m≥1}(2πi^{−k}/Nm)S(n, n′; Nm)J_{k−1}(4π√(nn′)/Nm).
Hence (proved, exact up to E₀):

    M₁₁ = D₁₁ + Σ_{m≥1} O_m + (n = 1 terms),   D₁₁ := Σ_{n≥2}|β_n|²,
    O_m := (2πi^{−k}/Nm) Σ_{n,n′≥2} β_n β̄_{n′} S(n, n′; Nm) J_{k−1}(4π√(nn′)/Nm),                     (1.1)

and the same for M₂ with β̄_{n′} → β_{n′} (O_m^{(2)}). The n = 1 terms (β₁ ≈ 0.05) contribute ≤ 5·10⁻³. X_q is linear in the
root-number Petersson sums (FP Lemma 2): Σ^hε_fλ_f(ℓ) = 2πN^{−½}Σ_{(m,N)=1}S(ℓ, N̄; m)m⁻¹J_{k−1}(4π√(ℓ/N)/m), with
λ_f(q²)λ_f(n) = Σ_{d|(q²,n)}λ_f(q²n/d²). C₂(Q) needs only Δ(q², q′²) = δ + O(N^{−3/2}).

## 2. The experiment (task 1)

### 2.1 Method (all exact identities; checked)
- **Bilinear Kloosterman forms by FFT (proved, elementary).** Σ_{n,n′}x_ny_{n′}S(n, n′; c) = Σ_{x∈(ℤ/c)^×}x̂(x)ŷ(x̄), x̂(x) :=
  Σ_n x_n e(nx/c): one FFT of length c per factor after binning n mod c; inverses mod Nm by CRT from tables mod N and mod m
  (check a2), moduli N²m′ handled directly at small N.
- **The Bessel kernel is separated by a low-rank factorisation.** With s = t_n, J_{k−1}(4π√(nn′)/Nm) = J_{k−1}(s s′/(πm)); on
  [s_min, s_max]² it is interpolated at Chebyshev nodes (0.6(s_max − s_min)s_max/(πm) + 40 nodes) and eigen-decomposed (eigenvalues
  kept down to 10⁻¹³ of the largest): max interpolation error ≤ 2·10⁻⁸ (t2_kernel.py); rank 830 (m = 1, T = 30) down to 6.
  Check a3: FFT + low-rank = direct double sum to 10⁻⁸ relative.
- **Convergence in m.** For k = 4 the terms decay slowly (J₃(y) ~ y³/48 below the transition y ≈ 1, i.e. m ≳ 0.7T²); at T = 20 the
  k = 4 tail m ∈ [31, 110] still moved M₁₁ by +0.37 with steps ±0.05 (not converged). **For T ≥ 12 we use k = 14**, where J₁₃ is
  negligible below y ≈ 5: terms are < 10⁻⁴ beyond m ≈ 50 (T = 20) and m ≈ 100 (T = 30); runs go to m = 60, 150, 140.
- **I_u by trapezoid quadrature** on a step-0.05 grid (agrees with step 0.02 to 4·10⁻¹³).
- **Validation against eigenforms (check b; N = 61, k = 4, T = 6, all 15 newforms, mfpetersson weights, Z_f summed form by form):**
  M₁₁ 1.96555 (Petersson, c ≤ 61·4000, moduli N²m′ included) vs 1.965554 (forms); M₂ 0.03650 + 0.01840i vs 0.036503 + 0.018396i;
  X_q (q ≤ 10) to 1.3·10⁻⁸; V(10) 1.84651 vs 1.846613 (the 10⁻⁴ is C₂(10): Δ(q², q′²) truncated at 6 moduli on the Petersson side).
- **Validation against a theorem (check c1): the integer model** Y_f := −(1/2π)Σ_{n≥1}λ_f(n)n^{−½}I_n satisfies (§4.3)
  Σ^h|Y_f|² = (1/4π²)Σ_ℓ|Ĥ(ln ℓ)|²/ℓ + O(N^{−1}). Petersson pipeline at N = 10007 (n ∈ [0.3T, 3.5T] in t): 0.40141 vs 0.40445
  (T = 6), 1.33876 vs 1.34265 (T = 12): 0.3–0.8 %, the size of the truncation of the n-range.

### 2.2 Table 1 — the deliverable (analysis.py; Q = 10, C₂(10) = 1.4593 at N = 5003, 10007)

D₁₁ = diagonal Σ|β_n|² (in brackets its asymptotic (1/2π)∫H²ln(Nt²/4π²), Prop. 4.4); "off" = Σ_mO_m of (1.1) = arcs + minor
(arcs: units x ≡ a·Nm/q ± 1, Lemma 3.1; in brackets the x = ±1 part alone); L² − D_int := L² − (1/2π)∫H² is the integer-model
off-diagonal (theorem 4.3); D_int = 0.903, 1.542, 2.332 at T = 12, 20, 30.

| N | k | T | D₁₁ [asym.] | off = arcs + minor [x = ±1] | L² − D_int | M₁₁ | M₂ | L² | X₁₀ (−LC₂) | V(10) | V/T² | V/D₁₁ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 10007 | 4 | 6 | 4.063 [4.070] | −0.112 (m ≤ 60) | −0.078 | 3.953 | 0.440+0.156i | 0.332 | −0.884−0.034i (−0.841) | 3.419 | 0.0950 | 0.841 |
| 10007 | 4 | 12 | 10.202 [10.210] | +0.107 (m ≤ 100) | 0.425 | 10.314 | 1.433−0.246i | 1.328 | −1.614+0.008i (−1.682) | 8.532 | 0.0593 | 0.836 |
| 10007 | 14 | 12 | 10.298 [10.210] | 0.159 = 0.431 − 0.272 [0.430] | 0.425 | 10.460 | 1.778+0.205i | 1.328 | −1.651+0.040i (−1.682) | 8.592 | 0.0597 | 0.834 |
| 10007 | 14 | 20 | 19.064 [19.008] | 2.346 = 2.113 + 0.233 [2.135] | 2.147 | 21.413 | 5.067+0.125i | 3.690 | −2.819−0.008i (−2.803) | 15.968 | 0.0399 | 0.838 |
| 10007 | 14 | 30 | 30.670 [30.632] | 5.760 = 6.094 − 0.334 [5.979] | 5.970 | 36.433 | 11.126−0.268i | 8.302 | −4.202+0.054i (−4.205) | 24.332 | 0.0270 | 0.793 |
| 5003 | 14 | 12 | 9.670 [9.584] | 0.059 = 0.408 − 0.350 [0.408] | 0.425 | 9.745 | 1.667+0.284i | 1.328 | −1.631+0.026i (−1.682) | 7.925 | 0.0550 | 0.819 |
| 5003 | 14 | 20 | 17.991 [17.938] | 2.319 = 2.127 + 0.192 [2.144] | 2.147 | 20.321 | 5.149+0.507i | 3.690 | −2.812−0.001i (−2.803) | 14.902 | 0.0373 | 0.828 |
| 5003 | 14 | 30 | 29.045 [29.016] | 5.718 = 6.015 − 0.297 [5.937] | 5.970 | 34.766 | 11.435+0.124i | 8.302 | −4.190−0.007i (−4.205) | 22.734 | 0.0253 | 0.783 |
| 61 | 4 | 6 | 1.920 [1.977] | 0.053 (m ≤ 4000, all c) | −0.078 | 1.966 | 0.037+0.018i | 0.332 | −0.529−0.034i (−0.852) | 1.847 | 0.0513 | 0.962 |

Reading. (a) **V/T² is small and decreasing** at both levels (0.060 → 0.040 → 0.027 at N = 10007; 0.055 → 0.037 → 0.025 at
N = 5003), **but V ≈ 0.78–0.84·D₁₁**: the residual R_f = Z_f + ε_fc_f^{(10)}L has the size of the diagonal (random-phase variance
≍ T log(NT)); V/T² falls like log(NT)/T. (b) **The off-diagonal of M₁₁ carries a T² term**: 0.16 → 2.35 → 5.76 (T = 12, 20, 30),
and its x = ±1 part alone equals the integer model's L² − D_int (0.430/0.425, 2.135/2.147, 5.979/5.970): the T² term of the prime
second moment is, modulus by modulus, the integer model's, times w(m). (c) The minor arcs of M₁₁ do not grow: −0.27, +0.23, −0.33
(N = 10007), −0.35, +0.19, −0.30 (N = 5003). (d) X_q follows the uniform line sums of Theorem D within 0.035–0.068 (q ≤ 10; L =
1.15–2.88) also at T = 20, 30, beyond the proved range T ≤ N^{1/3}. (e) **If the off-diagonal had no T² term**, V would be
D₁₁ + 2L Re X_Q + L²C₂ ≈ D₁₁ − C₂L², i.e. 18.6 instead of 24.3 at T = 30, and negative from T ≈ 90 on (D₁₁ ≍ 1.02T·log-growth against C₂L² = 0.0135T²;
impossible, V ≥ 0); the
arcs are exactly what keeps V ≥ 0. At finite T the T² parts of V cancel only up to (C₂ − w_eff)L² with w_eff = 1.02 at T = 30
(the lags 2πm ≤ T reach m = 1, 2, 3 only, where w = 1, 1, 1.25): V = D₁₁ − w_eff D_int − (C₂ − w_eff)L² + minor + (X_Q − line terms)
= 30.67 − 2.38 − 3.64 − 0.33 + 0.02 = 24.33 (N = 10007, T = 30).

### 2.3 Table 2 — which moduli (M₁₁: x = ±1 part against the lag formula (1/π)∫H(t)H(t − 2πm)dt of (3.3))

| N, T | m = 1: total [arcs, x=±1, lag] | m = 2 | m = 3 | m ≥ 4 total (arcs) |
|---|---|---|---|---|
| 10007, 12 | 0.240 [0.488, 0.488, 0.434] | −0.072 [−0.062, −0.062, 0.000] | −0.006 [0.010, 0.008, 0.000] | −0.002 (−0.004) |
| 10007, 20 | 1.975 [1.844, 1.844, 1.842] | 0.338 [0.369, 0.369, 0.308] | −0.088 [−0.061, −0.048, 0.000] | 0.121 (−0.039) |
| 10007, 30 | 3.727 [3.670, 3.670, 3.669] | 1.822 [1.845, 1.845, 1.853] | 0.596 [0.615, 0.493, 0.441] | −0.385 (−0.036) |
| 5003, 30 | 3.647 [3.670, 3.670, 3.669] | 1.791 [1.844, 1.844, 1.853] | 0.554 [0.615, 0.492, 0.441] | −0.274 (−0.114) |

The x = ±1 values do not depend on N (3.6701/3.6702) — they are the integer model's — and at m = 3 the q = 3 arcs add
0.25 × 0.492 (w(3) = 1.25). **M₂**: its T² term sits at m ∈ [T/π, 2T/π] (T = 30: the arcs at m = 10–19 sum to 12.68, all arcs to 12.22;
T = 20: m = 6–14), and there arcs = Σ_m w(m)·(x = ±1)_m to 0.002 (12.223 + 0.178i against 12.221 + 0.176i at T = 30, N = 10007);
the x = ±1 parts sum to 8.307 ≈ L² = 8.302 (integer model: Σ^hY_f² ≈ L²). The M₂ minor arcs are larger than M₁₁'s and negative,
concentrated in the same moduli: −0.17, −0.41 − 0.06i (m 6–14), −0.72 − 0.07i (m 10–19) at T = 12, 20, 30 (N = 10007), about
−6 % of the arcs there; with ≈ T/π main moduli this is ≈ −0.05 to −0.07 per modulus (heuristic reading: isolated, non-degenerate
stationary points of "near-arc" units x ≡ a·c/q + α, x̄ ≡ a′·c/q′ + α′ with αα′ ≠ 1, each O(1)); not decided by the data whether
this total is O(T) or a small multiple of T².

**Checks (check_second.py; outputs check_ab_out.txt, check_c_out.txt):** a 7/7, b 4/4, c 9/9 PASS (+ 3 INFO lines, §3.3).
c4: on all 48 moduli (8 runs, M₁₁ and M₂) where the x = ±1 part exceeds 0.3, arcs/(x = ±1 part) = w(m) to 0.004.

## 3. Where the T² term comes from (task 2a)

### 3.1 Not from generic pairs
For (m, N) = 1 and N ∤ nn′, twisted multiplicativity (check a1) gives

    S(n, n′; Nm) = S(n m̄, n′ m̄; N)·S(n N̄, n′ N̄; m) = √N·Kl₂(nn′m̄²; N)·S(nN̄, n′N̄; m),                   (3.1)

so the mod-N factor is a genuine normalised Kloosterman sum of the product nn′. Averaged over residue classes it is
Σ*_{a,b}S(a, b; c)/φ(c)² = μ²(c)/φ(c): if the primes p, p′ were replaced by their averages over classes mod Nm (the "PNT in
both variables" approximation), the m-th term would be (2π/Nm)·μ²(Nm)/φ(Nm)·∫∫(…) and, by the Bessel–Mellin identity (FP (★))
in both variables, equal for Σ^hZ_f² to (i^{−k}/2π)(μ²(m)/((N−1)φ(m)))·∫H(t)²e^{i(2θ(t) − t ln(Nm²))}dt ≪ T^{½}m^{½}/N, and for
Σ^h|Z_f|² to O_A(T^{−A}) (the v-integral produces δ(t + s), off the support of H) (heuristic derivation, not separately checked):
**generic pairs give nothing of order T²**. (Contrast with the first moment, where
S(Np, 1; N) = μ(N) is a Ramanujan sum for every p: all N − 1 terms of the Kloosterman sum conspire. Here no such degeneration
occurs, since N ∤ pp′.)

### 3.2 The arcs (computed: xprofile.py at T = 12; moment2.py records them for every run)
Opening the Kloosterman sum, O_m = (2πi^{−k}/c)Σ_{x∈(ℤ/c)^×}𝔅_m(x), 𝔅_m(x) := Σ_{n,n′}β_nβ̄_{n′}e((nx + n′x̄)/c)J(…), c = Nm.
At N = 10007, k = 4, T = 12, m = 1 (c = 10007): the two units x = ±1 (x̄ = x) give **0.48540**, the other 10004 units give
−0.278 (sum of their absolute values 3.39, r.m.s. 4·10⁻³ per unit: random-walk size √φ(c)·rms = 0.40). The integer model at
the same modulus gives **0.48516** in total (Poisson: only x = ±1 survive). At m = 4, Σ^hZ_f²: x = ±1 give 0.10939 + 0.23861i
(integer model at m = 4: 0.10923 + 0.23837i), all arcs give 0.21873 + 0.47725i = **2.000×** (w(4) = 2); at m = 6:
0.44842 − 0.29100i (integer 0.44851 − 0.29093i) and all arcs 0.56057 − 0.36334i = **1.250×** (w(6) = 1.25).

**Definition (arcs).** For q | m, d := c/q, a ∈ (ℤ/q)^×, j ∈ {±1}: the units x ≡ ad + j (mod c) whose inverse is again of the
form x̄ ≡ a′d + j (mod c).
**Lemma 3.1 (proved).** (i) For such a pair the signs agree (j′ = j) and a′ is determined by a: aa′d + a + a′ ≡ 0 (mod q) for
j = +1 (sign change for j = −1). (ii) The number of admissible a is A_q(m) = ∏_{p|q}A_p(m), A_p = p − 2 if p ∥ m, A_p = p − 1
if p² | m. (iii) Arcs of different q never pair: there are no units x ≡ ad + j, x̄ ≡ a′d′ + j′ with d = c/q ≠ d′ = c/q′.
*Proof.* xx̄ ≡ 1 mod gcd(d, d′) forces jj′ ≡ 1 modulo a number ≥ N > 2, so j = j′. (i): (ad + j)(a′d + j) = aa′d² + j(a + a′)d + 1
≡ 1 (mod qd) iff aa′d + j(a + a′) ≡ 0 (mod q). If p ∥ m (so p ∥ q·(stuff) with (d, p) = 1), a′(ad + j) ≡ −ja (mod p) has a
solution iff ad + j is a unit mod p: p − 2 values of a mod p; if p² | m then p | d and a′ ≡ −a: p − 1 values; CRT gives (ii).
(iii): take a prime ℓ with v_ℓ(q) ≠ v_ℓ(q′), say v_ℓ(q) > v_ℓ(q′), v := v_ℓ(c). In (ad + j)(a′d′ + j) − 1 = aa′dd′ + j(ad + a′d′)
(j² = 1), v_ℓ(ad) = v − v_ℓ(q) < v ≤ v_ℓ(a′d′), v_ℓ(aa′dd′) ≥ v, so the left side is not divisible by ℓ^v. ∎

**Weights.** Replacing Λ on each arc by its mean in progressions mod q, Σ_p ln p e(pa/q)Φ(p) ≈ (μ(q)/φ(q))∫Φ (heuristic; the
phase e(pj/c) of the arc is absorbed into the weight, which is legitimate because on the stationary set the total phase is slowly
varying, see 3.3), the x = ad ± 1 arc contributes (μ(q)/φ(q))² times the x = ±1 term of the integer model. Summing over the
admissible a:

    O_m^{arc} ≈ w(m)·O_m^{int},   w(m) := Σ_{q|m} μ²(q)A_q(m)/φ(q)² = ∏_{p∥m}(1 + (p−2)/(p−1)²)·∏_{p²|m}(1 + 1/(p−1)).   (3.2)

w(1) = w(2) = 1, w(3) = 1.25, w(4) = 2, w(5) = 1.1875, w(6) = 1.25, w(8) = 2, w(9) = 1.5.
**Proposition 3.2 (proved; check a6, a7).** w is multiplicative, (1 − 1/p)Σ_e w(p^e)p^{−e} = 1 + p⁻², hence
lim_{M→∞}M⁻¹Σ_{m≤M}w(m) = ∏_p(1 + p⁻²) = ζ(2)/ζ(4) = 15/π² (computed mean over m ≤ 2·10⁵: 1.51979). The q-th summand has mean
density μ²(q)A_q(m)1_{q|m}/φ(q)² → μ²(q)/q² (p-locally ((p−2)(1/p − 1/p²) + (p−1)/p²)/(p−1)² = 1/p²).
So **the q-th arc family carries the q-th term μ²(q)/q² of Σ^hc_f² = Σ_{q,q′}μ(q)μ(q′)(qq′)⁻¹Σ^hλ_f(q²)λ_f(q′²) → Σ_qμ²(q)/q²**,
and Lemma 3.1(iii) (no pairing of arcs with q ≠ q′) is the Petersson-side image of the diagonal q = q′ on the constant side.
Hecke products with d > 1 (p = p′) are the diagonal D₁₁ only; prime powers contribute O(T^{½}) per form and ≤ 10⁻² to M₁₁.

### 3.3 Which moduli (stationary phase; heuristic, computed)
Put u = √(n/N), v = √(n′/N) (so t_n = 2πu). By stationary phase β_n ≈ −(1/2π)n^{−½}√(πt_n)H(t_n)e(−2u)e^{i((k−1)π/2+π/4)}
(FP Lemma 5, fact 2052), and J_{k−1}(z) ≈ (2πz)^{−½}Σ_±e^{±i(z − (2k−1)π/4)}, z = 4πuv/m. On the arc x = x̄ = 1 the total phase
of 𝔅_m(1) (in units of 2π) is
- for M₁₁ (β̄_{n′}, minus branch): −2u + 2v + (u² + v²)/m − 2uv/m = ((u − v − m)² − m²)/m: **stationary on the whole line
  u − v = m, i.e. t − t′ = 2πm**, with e(−m) = 1;
- for M₂ (β_{n′}, plus branch): ((u + v − m)² − m²)/m: **stationary on t + t′ = 2πm**.
The Fresnel integral across the line (∫e((w − m)²/m)dw = √(m/2)e^{iπ/4}) and the amplitudes combine, for M₁₁, to the constant
1 (all phases cancel: i^{−k}e^{i(2k−1)π/4}e^{iπ/4} = 1), so (x = 1 and its mirror x = −1)

    O_m^{int} ≈ (1/π)∫H(t)H(t − 2πm)dt   (M₁₁, m ≥ 1),        O_m^{int,(2)} ≈ (1/2π)∫H(t)H(2πm − t)dt·e^{iψ_m}  (M₂).   (3.3)

By Poisson in the lag, Σ_{m≥1}(1/π)∫H(t)H(t − 2πm)dt = L² − (1/2π)∫H² + O_A(T^{−A}), the exact integer-model total (§4.3), and
(1/2π)Σ_m(H∗H)(2πm) = L² for M₂. So **M₁₁'s T² term sits at 1 ≤ m ≤ T/2π (lags t − t′ = 2πm inside the window), M₂'s at
m ∈ [T/π, 2T/π]**. Check c3: the computed x = ±1 parts agree with (1/π)∫H(t)H(t − 2πm)dt to 0.1 % at m = 1 (T = 20: 1.8437 vs 1.8418; T = 30:
3.6702 vs 3.6691) and 0.5 % at m = 2, T = 30 (1.8447 vs 1.8530); at lags close to the window length the leading order is poor
(T = 20, m = 2: 0.369 vs 0.308; T = 30, m = 3: 0.493 vs 0.441; T = 20, m = 3: −0.049 vs 0.000). (A first version of c3 asserted
10 % for all m ≤ 3 and failed at (T, m) = (20, 2); the check now asserts only the two claims above and prints the rest.) The phases
e^{iψ_m} of the M₂ terms are not constant at these T (the Fresnel zone 2π√m in t is as wide as the window when m ≈ T/π); only the
sum over m is robust there.
**Consequence (heuristic, given (3.2), (3.3)).** Σ_m O_m^{arc} ≈ Σ_m w(m)(1/π)∫H(t)H(t − 2πm)dt → (ζ(2)/ζ(4))(L² − (1/2π)∫H²)
+ O(T log T) (the error from the discrepancy of w: Σ_{m}(w(m) − w̄)g(m/M) = −(g(0)/2)Σ_{d≤M}h(d) + …, h = μ ∗ w, Σ_{d≤M}h(d) ≍ log M).
Hence, if the minor arcs are o(T²) (they are numerically O(1) per modulus), M₁₁ = (ζ(2)/ζ(4))L² + D₁₁ − (ζ(2)/ζ(4))(1/2π)∫H² + O(T log T)
+ (minor), X_Q = −LC₂(Q)(1 + o(1)) (FP Theorem D) and V = D₁₁ − O(T log T) + o(T²) = O(T log(NT)): the observed V/D₁₁ ≈ 0.8.
At finite T the effective weight Σ_m w(m)g(m)/Σ_m g(m) is a short average of w (at T = 20 only m = 1, 2 matter, where w = 1), so
the cancellation of the T² terms in V is incomplete by L²(C₂ − w_eff) ≍ T²·O(1/T): this is part of V at the sizes computed.

## 4. What is proved (task 2c)

**4.1 Lemma (exact quadratic forms; proved).** (1.1) holds with an error ≪ N^{−3/4}log(NT)·M₁₁^{½} from E₀ (Cauchy–Schwarz,
Σ^h1 = 1 + O(N^{−3/2+ε})), and V(Q) = M₁₁ + 2L Re X_Q + L²C₂(Q) exactly (expand the square; ε_f, c_f^{(Q)}, L real). *Proof:*
FP Theorem A(1), Lemma 1, Hecke; the double series converge absolutely (Weil, |J_ν(y)| ≤ (y/2)^ν/ν!). ∎

**4.2 Proposition (one-modulus large sieve; proved).** Let c ≥ 1, x, y supported on n ∈ [1, X], A_x(b) := Σ_{n≡b (c)}x_n. Then
(i) Σ_{n,n′}x_ny_{n′}S(n, n′; c) = Σ_{x∈(ℤ/c)^×}x̂(x)ŷ(x̄) and its modulus is ≤ c‖A_x‖₂‖A_y‖₂ ≤ (c + X)‖x‖₂‖y‖₂;
(ii) the matrix (S(a, b; c))_{a,b mod c} has operator norm exactly c (it is U P Uᵀ with U = (e(ax/c))_{a,x∈(ℤ/c)^×}, UᵀŪ = c·I on
the units, P the permutation x ↦ x̄), so (i) is sharp for general coefficients;
(iii) with the Bessel weight, for every σ ∈ (−(k−1)/2, 0):
|Σ_{n,n′}x_ny_{n′}S(n, n′; c)J_{k−1}(4π√(nn′)/c)| ≤ C_k|σ|⁻¹(c + X)‖x‖‖y‖(2πX/c)^{2|σ|}.
*Proof.* (i) open S and use Parseval on (ℤ/c), Cauchy–Schwarz within residue classes (≤ X/c + 1 terms each). (ii) the columns
e(·x/c), x a unit, are orthogonal of norm √c. (iii) Mellin–Barnes: J_ν(z) = (1/2πi)∫_{(σ)}Γ(ν/2 + s)Γ(ν/2 + 1 − s)⁻¹(z/2)^{−2s}ds,
absolutely convergent for −ν/2 < σ < 0 (Gradshteyn–Ryzhik 6.561.14 inverted, as FP Lemma 4; |Γ-ratio| ≍ |t|^{2σ−1}, integral ≍ |σ|⁻¹);
for each s the kernel (2π√(nn′)/c)^{−2s} factors and (i) applies to x_nn^{−s}, y_{n′}n′^{−s}. ∎
**Corollary 4.2′ (proved).** With x = y = β (‖β‖² = D₁₁ ≍ T log(NT), X ≍ NT²): |O_m| ≪ (1 + T²/m)D₁₁·log(2 + T²/m) for m ≤ T²
(σ = −1/log(2 + T²/m)) and ≪ (1 + T²/m)D₁₁(T²/m)^{k−1−ε}·… for m > T² (σ near −(k−1)/2): Σ_m|O_m| ≪ T²log²T·D₁₁ ≍
**T³log(NT)log²T**. Together with Cauchy–Schwarz and FP Theorems B/C: L²(1 + o(1)) ≤ M₁₁ ≪ T³log(NT)log²T. Useless for "almost
all" (needs o(T²)); and no per-modulus bound valid for general coefficients can do better than (i) (by (ii)).

**4.3 Theorem (integer model; proved).** Y_f := −(1/2π)Σ_{n≥1}λ_f(n)n^{−½}I_n = −(1/2π)ε_fΣ_{ℓ≥1}λ_f(ℓ)ℓ^{−½}Ĥ(−ln ℓ), and
Σ^h|Y_f|² = (1/4π²)Σ_ℓ|Ĥ(ln ℓ)|²/ℓ + O(N^{−1}) = L²(1 + O_A(T^{−A})) + O(N^{−1}).
*Proof.* First identity: the weight-k, level-N form of fact 2053 (r = 1): write the left side as (1/2πi)∫_{(1+δ)}h(s)X(s)⁻¹L(s, f)ds,
move to Re s = −δ (X⁻¹ is holomorphic for −(k−1)/2 < σ < (k+1)/2, FP §1; L(s, f) entire; horizontal edges vanish by the Gaussian
decay of h), use L(s, f) = ε_fX(s)L(1 − s, f) and expand L(1 − s, f) — the proof of fact 2053 verbatim with the weight-k X.
Second: Petersson for (ℓ, ℓ′) with Ĥ(ln ℓ) ≪ T(T ln ℓ)^{−A}e^{−(ln ℓ)²/4}: the off-diagonal Δ(ℓ, ℓ′) − δ ≪ N^{−1/2−(k−1)+ε}(ℓℓ′)^{(k−1)/2} (Weil and |J_ν(y)| ≤ (y/2)^ν/ν!, c = Nm)
summed against these weights is O(N^{−1}). ∎ (Computed: §2.1.) Hence the integer-model off-diagonal Σ_mO_m^{int} equals
L² − (1/2π)∫H²(1 + O(1/T)) (its diagonal is (1/2π)∫H², stationary phase): **in the integer model the T² term of the off-diagonal
is a theorem**, proved through the spectral side (fact 2053), not by evaluating Kloosterman sums.

**4.4 Proposition (diagonal; proved (cited) FP §7).** D₁₁ = (1/2π)∫H(t)²ln(Nt²/4π²)dt·(1 + O(1/T)) + O(T^{½}) (stationary phase
|I_p|² = πt_pH(t_p)²(1 + O(1/T)) and the PNT). Computed (analysis.py): D₁₁ against the formula at every run (Table 1, columns 4–5).

**4.5 Arc main terms (sketch, GRH).** For each m ≤ M₁ and each arc pair (x, x̄) of Lemma 3.1, separate the Bessel kernel by
Mellin–Barnes and evaluate Σ_p ln p e(pa/q)·[e(pj/c)×smooth] by the explicit formula for L(s, χ), χ mod q, under GRH:
main term (μ(q)/φ(q))∫, error ≪ √X log²(qX)·Var(weight), Var ≍ sup·(X/c + T). Bookkeeping (not written in full) gives a relative
error ≍ T/(m√N)·log², so Σ_m O_m^{arc} = Σ_m w(m)O_m^{int}(1 + O(T N^{−½}log²)) for T ≤ N^{½−ε}; with (3.3) and Prop. 3.2 the
total is (ζ(2)/ζ(4))(L² − (1/2π)∫H²) + O(T log T). **Status: sketch** (the per-arc error analysis and (3.3) beyond stationary
phase are not written; the numerics of §2–3 support every step).

**4.6 Theorem (conditional reduction; proved, trivial given its hypotheses).** Assume GRH for Dirichlet L-functions mod q ≤ T,
the arc evaluation 4.5, and
(MA) Σ_{m≥1}|O_m − O_m^{arc}| = o(T²) (T → ∞, T ≤ N^{δ}).
Then, with Q = Q(T) → ∞ slowly, V(Q) = O(T log(NT)) + o(T²), and for every η > 0
Σ^h{f : |Z_f(H) + ε_f c_f (1/2π)∫H| > ηT} → 0 if T/log N → ∞ (harmonic measure; c_f^{(Q)} → c_f in Σ^h-mean square, FP Prop. 5
argument, GRH(Sym²f)).
*Proof.* M₁₁ = D₁₁ + Σ_m O_m^{arc} + Σ_m(O_m − O_m^{arc}) = (ζ(2)/ζ(4))L² + O(T log(NT)) + o(T²) by 4.4, 4.5, (MA); X_Q =
−L·Σ_{q≤Q}μ²(q)q⁻²(1 + o(1)) by FP Theorem D (uniform in q ≤ Q for fixed Q; let Q → ∞ slowly); C₂(Q) = Σ_{q≤Q}μ²/q² + O(N^{−1}).
So V(Q) = (ζ(2)/ζ(4) − Σ_{q≤Q}μ²/q²)L² + O(T log(NT)) + o(T²) = o(T²). Chebyshev. ∎
The hypothesis (MA) is the whole difficulty; everything else is proved or sketched under GRH.

## 5. The obstruction (task 2b)

**What is needed.** (MA): after removing the arcs, Σ_m|O_m^{min}| = o(T²), where by (1.1), (3.1)

    O_m^{min} = (2πi^{−k}/(√N m)) Σ_{p,p′≍NT²} β_pβ̄_{p′} Kl₂(pp′m̄²; N) S(pN̄, p′N̄; m) J_{k−1}(4π√(pp′)/Nm)  − (arcs),     (5.1)

for 1 ≤ m ≲ T² (k fixed; for m ≫ T² the Bessel factor decays like (T²/m)^{k−1}). Normalisation: |β_p| ≍ T^{½}X^{−½}ln p,
X ≍ NT², ‖β‖₂² ≍ T log(NT), ‖β‖₁ ≍ (XT)^{½}.
**What is known per modulus.**
- Trivial (Weil, |Kl₂| ≤ 2, |S(·;m)| ≤ τ(m)√m): |O_m| ≪ N^{−½}m^{−½}‖β‖₁² ≍ √N·T³/√m.
- Parseval/large sieve (Prop. 4.2): |O_m| ≪ (1 + T²/m)T log(NT) — better than Weil by ≈ √N, and **optimal for general
  coefficients** (Prop. 4.2(ii)).
- Fouvry–Kowalski–Michel, *Algebraic trace functions over the primes*, Duke Math. J. 163 (2014), Theorem 1.17(1) (opened,
  arXiv:1211.6043v3): for an isotypic non-exceptional trace weight K mod p and α, β supported on [M/2, 2M], [N/2, 2N],
  Σ_{m,n}α_mβ_nK(mn) ≪ ‖α‖‖β‖(MN)^{½}(p^{−1/4} + M^{−½} + p^{1/4}log^{½}p·N^{−½}). Here p ↦ N (the level), M = N′ ≍ NT²: the saving
  over the trivial bound is N^{−1/4} — i.e. |O_m| ≪ N^{1/4}T³m^{−½}: **worse than Parseval by N^{1/4}**. (Remark 1.2 of KMS quotes
  the same bound as the Pólya–Vinogradov-type estimate.)
- FKM, same paper, Theorem 1.5: Σ_{q prime}K(q)V(q/X) ≪ X(1 + p/X)^{1/6}p^{−η}, η < 1/24: a linear sum over primes; for X ≥ p it
  saves p^{−1/24}; our sums are bilinear with a second prime variable and the coupling Bessel kernel.
- Kowalski–Michel–Sawin, *Bilinear forms with Kloosterman sums and applications*, Ann. of Math. 186 (2017), Theorem 1.1 (opened,
  arXiv:1511.01636v5): for q prime, 1 ≤ M ≤ Nq^{1/4}, q^{1/4} < MN < q^{5/4}, the n-variable in an interval of length N ⊂ [1, q−1]:
  B(Kl_k, α, β) ≪ q^ε‖α‖‖β‖(MN)^{½}(M^{−½} + (MN)^{−3/16}q^{11/64}). Designed for both variables **below** the Pólya–Vinogradov range;
  our M = N′ ≍ qT² violate MN < q^{5/4} and N < q.
- Blomer–Fouvry–Kowalski–Michel–Milićević, *On moments of twisted L-functions*, Amer. J. Math. 139 (2017), Theorem 5.1 (opened,
  arXiv:1411.4467v3): for M, N ≤ q, Σ_{m≤M, n∈𝒩}α_mβ_nKl₂(amn; q) ≪ (qMN)^ε(MN)^{½}‖α‖‖β‖(M^{−½} + q^{1/4}N^{−½}); again M, N ≤ q,
  and in the normalisation of the trivial bound.
- Fouvry–Michel 1998 (Ann. Sci. ÉNS, "Sur certaines sommes d'exponentielles sur les nombres premiers"): **not opened**; cited in
  KMS Remark 1.4(2) as the starting point of the Type I method for Kloosterman sums; same range issue.
**Why GRH does not help.** GRH for L(s, χ) gives the arcs (χ mod q | m, q small). The minor arcs need the characters mod Nm of
the expansion Kl₂(ab) = Σ_χĉ_χχ(a)χ(b), |ĉ_χ| ≈ 1 for ≈ Nm characters; GRH bounds each |Σ_pβ_pχ(p)| by ≍ T, and the large sieve over
χ mod Nm is again Prop. 4.2(i). Equivalently (5.1) is a **Kloosterman-weighted variance of primes in progressions to one modulus
c = Nm < x = NT²**: Σ_{a,b mod c}S(a, b; c)E(a)Ē(b), E(a) the deviation of the β-weighted prime count in the class a from the arcs'
prediction. The operator norm c (Prop. 4.2(ii)) turns even an optimal (Hooley-type) variance asymptotic Σ_a|E(a)|² ~ (x/c)·log c
into a bound of the Parseval size; what is needed is that the E(a) do not correlate with the Kloosterman kernel. Already one
modulus m ≍ 1 has Parseval bound ≍ T³log(NT) against the o(T²) needed (a saving ≫ T log(NT)); in total the saving needed over
Parseval is a factor ≫ T·log(NT)log²T.
**How far off.** Parseval total ≍ T³log(NT)log²T; needed o(T²): a factor ≫ T·log(NT)log²T. The observed minor arcs (Table 2)
are 0.1–0.4 per modulus at N = 10007, i.e. ≍ D₁₁/√(Nm)-size (square-root cancellation over the φ(c) units: at m = 1, T = 12 the
non-arc units sum to −0.278 against the random-walk size 0.40).
**Relation to known ranges.** In the language of low-lying zeros (Hughes–Miller, Duke Math. J. 136 (2007); abstract opened,
paper not): the n-th centred moments of the one-level density of this family are Gaussian for Fourier support in (−1/n, 1/n)
with the sign split; our test function e^{2iθ(γ)}H(γ) has its Fourier mass at exactly the mean density (normalised frequency 1),
so the second moment needs support ≈ 2 — outside every range where the Petersson off-diagonal is controlled by Weil bounds.

## 6. GAPs and doubts

**GAPs.** (1) **(MA), the minor arcs** (§5): the route is closed here. (2) The arc evaluation 4.5 under GRH (error bookkeeping
not written). (3) The lag formula (3.3) beyond leading-order stationary phase (two-variable degenerate stationary set; the
Fresnel zone 2π√m is comparable to the window when m ≍ T/π, which is why the M₂ phases e^{iψ_m} are not constant at T ≤ 30).
(4) The derivation of the weight w(m) replaces Λ on each arc by its mean in progressions mod q with a weight that contains
e(pj/c); this is the standard major-arc step but the stationary-set argument that justifies it is a sketch. (5) The T log T
term of Σ_m(w(m) − w̄)·lag_m (§3.3) is computed heuristically, not proved.
**Doubts.** (1) Range: FP Theorem B is proved for T ≤ N^{1/3}(log N)^{−3} (N^{1/3} = 21.5 at N = 10007, 17.1 at N = 5003); the runs at
T = 20, 30 are outside the proved range of the first moment — the identities used are exact at every T, only the asymptotic
interpretation is not covered; X_q still match the uniform line sums of Theorem D there (Table 1, last column).
(2) Truncation of the n-range: t_n ∈ [0.3T, 3.0T] (k = 14) or [0.3T, 3.5T] (k = 4). At T = 12 (k = 4) the cut at 2.8T instead of
3.5T moved D₁₁ by −0.018 and the m ≤ 6 off-diagonal by +0.04: a systematic error ≈ 0.5 % of M₁₁; the integer model, which is far
more sensitive (dense n), agrees with its theorem to 0.3–0.8 %. (3) k = 4 at T = 20: the m-sum is not converged at m = 110 (tail
steps ±0.05); those numbers are not used in Table 1. k = 14 converges (terms < 10⁻⁴ by m = 50–100). (4) V at fixed Q = 10:
C₂(10) = 1.4593 against ζ(2)/ζ(4) = 1.5198; V(Q) varies with Q by ≈ 0.04·L² per doubling of Q (Table in mixed outputs);
c_f^{(Q)} → c_f only in mean square. (5) The low-rank kernel tolerance 10⁻¹³ and interpolation error ≤ 2·10⁻⁸ are checked on
random points, not proved; the N = 61 agreement (4·10⁻⁶) is the end-to-end check. (6) "Square-root cancellation of the minor arcs"
is an observation at two levels and T ≤ 30, not a statement about N → ∞.

## 7. Prior art

- **Petersson / Kloosterman machinery, root number via λ_f(N):** Iwaniec–Luo–Sarnak (FP §10, opened there). Known.
- **Second moments/variances of linear statistics of zeros over this family:** Hughes–Miller, Duke Math. J. 136 (2007) (abstract
  opened): centred moments of the one-level density, Gaussian for support in (−1/n, 1/n) with the sign split, non-Gaussian
  corrections up to (−1/(n−1), 1/(n−1)). Hughes–Rudnick (mock-Gaussian behaviour, not opened). Our statistic sits at normalised
  frequency exactly 1 at height T (support ≈ 2 for the second moment): **not covered; no source found** that computes
  Σ^h|Σ_γH(γ)e^{2iθ_f(γ)}|² or locates its off-diagonal main term.
- **Bilinear forms in Kloosterman sums:** FKM 2014 (Thm 1.5, 1.17), KMS 2017 (Thm 1.1, 1.3, Remark 1.2), BFKMM 2017 (Thm 5.1, Prop. 5.5):
  **opened** (arXiv texts in `dl/`), statements quoted in §5. FKM *Algebraic twists of modular forms and Hecke orbits* (GAFA 2015,
  arXiv:1207.0617): downloaded and searched, not used (it bounds Σλ_f(n)K(n)V(n/X) for one form, Theorem 1.2 there; not bilinear).
  Fouvry–Michel 1998: **not opened**.
- **The arc mechanism** (the T² term from the units x ≡ a·Nm/q ± 1 of S(p, p′; Nm), weights w(m) with mean ζ(2)/ζ(4)): **not found**.
  It is the Petersson-side image of the Hardy–Littlewood singular series of FP §9a (constant side Σ^hc_f² = Σμ²(q)/q²).
- **Integer model (5.3):** a direct consequence of fact 2053 (Ren–Ye for level 1, per FACTS) + Petersson; **known in substance**.
- **Hooley/Montgomery variance of primes in progressions** (the analogy of §5): not opened; quoted from memory as context only.

## 8. What a referee should check first

1. Lemma 3.1 (arcs: signs agree, count A_q(m), no pairing of different q) and Prop. 3.2 (mean of w = ζ(2)/ζ(4)); check a6, a7.
2. The claim that the x = ±1 part of O_m for primes equals the integer model's O_m (0.48540 vs 0.48516 at m = 1; M₂ at m = 4, 6), and
   that arcs/(x = ±1) = w(m) at every modulus with a sizeable term (check c2, c4; Table 2).
3. Prop. 4.2(ii): the operator norm of (S(a, b; c)) is c — this is why no per-modulus bound for general coefficients beats Parseval,
   and why the FKM/KMS bounds (normalised against ‖α‖‖β‖(MN)^{½}) are weaker here by N^{1/4}.
4. The N = 61 end-to-end check (check b): Petersson double sums including the moduli N²m′ against eigenform data.
5. §5: that (MA) is genuinely needed — i.e. that no averaging over m (Kuznetsov / spectral large sieve) is available that is not
   circular: the spectral side of Σ_{c≡0(N)}S(n, n′; c)c⁻¹J(4π√(nn′)/c) is the family itself.
6. Table 1: V/T² and V/D₁₁ from the saved runs (analysis.py), and the statement that V/D₁₁ ≈ 0.8 means the decrease of V/T² is the
   diagonal's, not a test of the T² cancellation by itself.

## 9. Files and commands

All in ``:
`smlib.py`, `moment2.py` (usage `moment2.py N k T m0:m1 primes lo hi`), `moment.py` (first version, k = 4 runs at N = 10007,
T = 6, 12 and the integer model), `mixed.py` (`mixed.py N k T Q MG primes lo hi momentfiles`), `forms_side.py` (`N k T`),
`xprofile.py`, `analysis.py`, `check_second.py` (`a`, `b`, `c`), `run_5003.sh`; outputs `out_*.txt`, `mix_*.txt`, `mom_*.npz`,
`mixed_*.npz`, `forms_61_4_6.npz`; papers in `dl/` (kms, fkm_primes, fkm_twists, bfkmm: pdf + extracted txt).
Every job ran as `(ulimit -v 1500000; timeout 600 nice .venv/bin/python -I …)`.
