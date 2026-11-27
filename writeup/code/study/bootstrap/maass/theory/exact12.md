# Level 1, one form L(f): no admissible weight below 12 − 3.7·10⁻⁸, from Ξ_Δ²·B, unconditionally (27 November 2026)

**In short.** This is the problem of [../maass.md](../maass.md), sections 2 and 4, with one form only. The Gamma factor is Γ_C(s + (κ−1)/2), κ real; the conductor is 1; there is no pole. The unknown is one probability measure μ_p per prime on θ ∈ [0, π], and the coefficient at p^m is ⟨2 cos mθ⟩ log p. Take ĝ = Ξ_Δ(t)²·B(t) with Ξ_Δ(t) = Λ(1/2 + it, Δ) and B ≥ 0 an even polynomial.

1. **Proved (Lemma 1).** A^{(12)}(g) = 2 Σ_p P_p(x_p), with x_p = τ(p)/(2p^{11/2}). This is the classical explicit formula of L(s, Δ) applied to this kernel; the zero side is identically 0 and no hypothesis on the zeros is used. Deligne is not used either.
2. **Proved (Lemma 5, Theorem 6).** Re ψ′(x + iy) ≥ 0 for x ≥ 2 (Pólya's argument), so A^{(κ)}(g) does not decrease in κ on [4, 12]. A family that is feasible at weight κ must satisfy A^{(12)}(g) − A^{(κ)}(g) ≤ 2 Defect(g). With B = 1 this works on all of (0, 12), because there g > 0.
3. **Computer-assisted, one code.** No family is feasible at any κ ∈ (0, 12 − ε):
   - B = 1: ε = 4.09;
   - degree 4: ε = 0.724;
   - degree 8: ε = 8.8·10⁻³;
   - degree 12: ε = 3.5·10⁻⁵;
   - degree 18: **ε = 3.7·10⁻⁸**.
   The A′ integral is bounded below in closed form (Hurwitz ζ values times exact moments), so no quadrature is needed.
4. **Caveats.**
   - Without GRH for L(s, Δ), and the real angles that Deligne's theorem provides, no feasible family is known at κ = 12. So the statement may be vacuous.
   - Exactness (κ*_∞ = 12) needs Defect → 0, which is unproved. Computed: it falls steadily, 0.36 → 2.2·10⁻⁹ over degrees 0-18.
   - Only ε = 3.7·10⁻⁸ is shown, not every ε > 0.
   - The repo's 11.9999983 (support 8) and 12.0000000 (support 16) are finite-support values of a problem with four forms; nothing here bounds them.

## 1. Set-up (maass.md section 2, κ₁ = (κ−1)/2, κ₂ = (κ+1)/2)

Q_{μ,κ}(f) = A^{(κ)}(g) − 2 Re Σ_p Σ_m ⟨2 cos mθ⟩_{μ_p} log p · p^{−m/2} g(m log p), with g = f*f̃, f ∈ C_c^∞(ℝ). The archimedean term is

`A^{(κ)}(g) = A_{κ₁}(g) + A_{κ₂}(g) = −2(γ + log π) g(0) + ∫₀^∞ [4e^{−2u} g(0)/(1 − e^{−2u}) − 2 e^{−κu/2} g(u)/(1 − e^{−u})] du`

(even g). This uses (e^{−κu/2} + e^{−(κ/2+1)u})/(1 − e^{−2u}) = e^{−κu/2}/(1 − e^{−u}). Equivalently A^{(κ)}(g) = (1/π) ∫ ĝ(t) [−log 2π + Re ψ(κ/2 + it)] dt. This follows from Γ_R(s)Γ_R(s+1) = Γ_C(s) = 2(2π)^{−s}Γ(s) and Gauss's integral; the two constants agree because ∫₀^∞ e^{−u}/(1 + e^{−u}) du = log 2. Since κ > 0 is assumed, A^{(κ)} is defined. maass.md section 4 imposes four forms (L, Sym², Sym³, Sym⁴). Positivity of L(f) alone is a weaker requirement, so every bound here also holds for the four-form problem.

## 2. The objects (derived here; checked numerically, section 6)

- **Normalisation.** L(s, Δ) = Σ τ(n) n^{−11/2−s} (analytic normalisation) and Λ(s) = Γ_C(s + 11/2) L(s, Δ). Classical inputs:
  - (C1Δ) Δ(i/y) = y^{12} Δ(iy);
  - (C2Δ) the Euler product L(s, Δ) = ∏_p (1 − τ(p)p^{−11/2}p^{−s} + p^{−2s})^{−1} (Mordell), together with Hecke's bound τ(n) = O(n^6), so that everything converges absolutely for Re s > 3/2;
  - (C3) the product formula Δ = q ∏(1 − q^n)^{24}, so that Δ(iy) > 0.
- **Mellin (proved).** Λ(s) = 2 ∫₀^∞ Δ(iy) y^{s+11/2} dy/y for every s: termwise, ∫ e^{−2πny} y^{w−1} dy = (2π)^{−w} Γ(w) n^{−w}, and Γ_C(w) = 2(2π)^{−w}Γ(w). Δ(iy) decays like e^{−2πy} at ∞, and by (C1Δ) doubly exponentially at 0, so Λ is entire. Substituting y ↦ 1/y gives Λ(1 − s) = Λ(s), i.e. root number +1 (i^{12} = 1).
- **φ and Ξ.** With y = e^u, Ξ(t) := Λ(1/2 + it) = ∫ φ(u) e^{itu} du, where φ(u) = 2e^{6u} Δ(ie^u). φ is even by (C1Δ) and positive by (C3). φ and all its derivatives are O(exp(−e^{|u|})), so Ξ is real, even and entire, and decays faster than any power in horizontal strips.
- **Closed form (Lemma 2, proved).** (φ*φ)(w) = 8 e^{6w} Σ_{N≥1} (τ*τ)(N) K₀(4π√N e^{w/2}), with (τ*τ)(N) = Σ_{nm=N} τ(n)τ(m). Proof: write φ(v+w) with Δ(ie^{v+w}) and φ(v) = φ(−v) with Δ(ie^{−v}). Each pair gives ∫ exp(−2πn e^{w} e^{v} − 2πm e^{−v}) dv = 2K₀(4π√(nm e^w)). Tonelli applies via Hecke's bound.
  - Derivatives: z = 4π√N e^{w/2}, d/dw = (1/2) z d/dz, D^j[z^{12}K₀] = z^{12}(p_j K₀ + q_j K₁), with p_{j+1} = 12p_j + z p_j′ − z q_j and q_{j+1} = 11q_j + z q_j′ − z p_j. Hence (d/dw)^j(φ*φ)(log n) = 8n^6 Σ_N (τ*τ)(N) 2^{−j}[p_j K₀ + q_j K₁](4π√(Nn)), and g_B = Σ_k b_k(−1)^k(φ*φ)^{(2k)}.
  - Constants: 8 = 2·2·2 (twice the 2 of φ, and 2K₀); 4π√(nm) is twice the geometric mean of 2πn and 2πm; the exponent 12 is the weight; τ*τ appears because τ is multiplicative but not completely multiplicative.
- **Decay.** g_B(log n) ~ exp(−4π√n): much slower than for the characters (e^{−2πn/3}, e^{−2πn/5}), so many more primes matter.
- **Local polynomials.** x = cos θ ∈ [−1, 1] and 2 cos mθ = 2T_m(x), so P_p(x) = Σ_m 2T_m(x) log p · p^{−m/2} g(m log p). Defect = Σ_p [P_p(x_p) − min_{[−1,1]} P_p], with x_p = τ(p)/(2p^{11/2}).

## 3. Lemma 1 (proved, unconditional): A^{(12)}(g_B) = 2 Σ_p P_p(x_p)

This is route 1 of [exact3.md](exact3.md), with these changes:
- Λ′/Λ(s) = ψ_C(s + 11/2) + L′/L(s, Δ), with ψ_C(s) = −log 2π + ψ(s).
- On Re s = 2: −L′/L(s, Δ) = Σ_p Σ_m (α_p^m + β_p^m) log p · p^{−ms}, where α_p + β_p = τ(p)p^{−11/2} and α_pβ_p = 1. This converges absolutely for Re s > 3/2 by (C2Δ).
- α^m + α^{−m} = 2T_m((α + α^{−1})/2) = 2T_m(x_p) is an identity of polynomials, so P_p(x_p) is meant algebraically. Deligne is not needed: x_p ∈ [−1, 1] is not used.
- The poles of ψ_C(s + 11/2) lie at s = −11/2, −13/2, … (the trivial zeros, absent from Λ), outside 1/2 ≤ Re s ≤ 2.
- On Re s = 1/2 the Gamma piece becomes (1/π)∫ĝ Re ψ_C(6 + it) dt = A^{(12)}(g).

F = ΛΛ′B̃ is entire, F(1 − s) = −F(s), and F decays in strips (section 2), so ∫_{(2)} F = 0, exactly as in exact3.md. Route 2: the classical explicit formula for L(s, Δ), whose zero side Σ_ρ Ξ(γ_ρ)²B(γ_ρ) is 0.

**Where Deligne enters (not in Lemmas 1-6):**
- (i) x_p ∈ [−1, 1], so that Defect_p ≥ 0 and the Satake data of Δ define a competitor;
- (ii) the tail bound |(τ*τ)(N)| ≤ 8N^7 in the certificate (from |τ(n)| ≤ d(n)n^{11/2}, d₄ ≤ d³ ≤ 8N^{3/2}); Hecke's bound with an explicit constant would serve as well;
- (iii) the upper side: Δ is feasible at κ = 12 only if its angles are real (Deligne) and GRH for L(s, Δ) holds.

## 4. Lemma 3 and the monotonicity lemma

**Lemma 3 (as in exact3.md, proved).** Feasibility for all f ∈ C_c^∞ gives Q_{μ,κ}(g_B) ≥ 0. F = b(i d/du)φ decays super-exponentially, Ω_κ(t) = −log 2π + Re ψ(κ/2 + it) = O(log(2 + |t|)) for κ > 0, and |⟨2 cos mθ⟩| ≤ 2.

**Lemma 5 (proved).** Re ψ′(x + iy) ≥ 0 for x ≥ 2 and every real y.
- Write ψ′(z) = ∫₀^∞ h(t) e^{−iyt} dt with h(t) = t e^{−xt}/(1 − e^{−t}), so Re ψ′ = ∫ h(t) cos(yt) dt.
- With k(t) = t/(1 − e^{−t}) = (t/2)coth(t/2) + t/2: k ≥ 1, 1/2 ≤ k′ < 1 and k″ ≥ 0, because (u coth u)″ = 2(u cosh u − sinh u)/sinh³u ≥ 0.
- Then h′ = e^{−xt}(k′ − xk) < 0 and h″ = e^{−xt}(k″ − 2xk′ + x²k) ≥ e^{−xt}x(x − 2) ≥ 0. So h is positive, decreasing, convex and tends to 0.
- Integrating by parts, ∫h cos(yt) dt = (1/|y|) ∫(−h′) sin(|y|t) dt ≥ 0, since −h′ ≥ 0 is decreasing.
- The restriction matters: Re ψ′(x + iy) ≈ (x − 1/2)/y² for large y, and Re ψ′(0.3 + 50i) < 0 (checked).

Consequence: Δ_B(κ) := A^{(12)}(g) − A^{(κ)}(g) = (1/π)∫ĝ [Re ψ(6+it) − Re ψ(κ/2+it)] dt does not increase in κ on [4, 12]. Equivalently, Δ_B(κ) = ∫_κ^{12} A′(s) ds, with A′(s) = ∫₀^∞ u e^{−su/2} g(u)/(1 − e^{−u}) du = (1/2π)∫ĝ Re ψ′(s/2 + it) dt.

## 5. Theorem 6 and the certificate

**Theorem 6 (proved).** If a family is feasible at weight κ ∈ (0, 12], then Δ_B(κ) ≤ 2 Defect(g_B) for every B.

*Proof.* Lemma 3 gives A^{(κ)}(g) ≥ 2Σ_p ∫P_p dμ_p ≥ 2Σ_p min P_p, and Lemma 1 gives A^{(12)}(g) = 2Σ P_p(x_p). Subtract. □

**Lower bound for Δ_B (proved, closed form).** A′(s) = ∫ h_s g with h_s ≥ 0 and g(u) = (1/2π)∫ĝ cos(tu) dt, ĝ ≥ 0. From cos x ≥ Σ_{j≤2m+1}(−1)^j x^{2j}/(2j)! and ∫₀^∞ h_s(u) u^{2j} du = (2j+1)! ζ(2j+2, s/2):

A′(s) ≥ Σ_{j≤2m+1} (−1)^j (2j+1) ζ(2j+2, s/2) M_{2j}, with M_{2j} = (1/2π)∫ĝ t^{2j} dt (exact, from the closed form at n = 1).

Integrating over s ∈ [κ, 12], using ∂_a ζ(n, a) = −nζ(n+1, a) and ζ(2, a) = ψ′(a):

`Δ_B(κ) ≥ 2M₀[ψ(6) − ψ(κ/2)] + 2 Σ_{j=1}^{2m+1} (−1)^j M_{2j} [ζ(2j+1, κ/2) − ζ(2j+1, 6)]`

for every 2m+1; the best truncation is taken.

**Assembly.**
- B = 1: g = φ*φ > 0, so Δ_1 decreases in κ on all of (0, 12) (u-side form). Hence no κ < 12 − ε₁ is feasible, where ε₁ is certified by the bound above at κ = 12 − ε₁.
- General B: for κ ∈ [4, 12 − ε], Lemma 5 gives Δ_B(κ) ≥ Δ_B(12 − ε) > 2 Defect.
- Since 12 − ε₁ = 7.91 > 4, the union covers (0, 12 − ε).

**Certificate ([validate_exact12.py](validate_exact12.py), arb at 600 bits).**
- B > 0 by Sturm.
- g(log n) for prime powers n ≤ 700, with the N-tail bounded by an incomplete-Gamma integral, and τ from PARI (checked against the q-expansion).
- Branch and bound for min P_p on [−1, 1], with P_p converted from the Chebyshev basis with exact integer coefficients, and P_p(x_p) as a ball.
- Everything above 700 is bounded by 4 Σ_{n>700} log n · n^{−1/2}|g(log n)| (factor 4 since |2T_m| ≤ 2), below 10⁻⁷⁴ at every degree.
- The moments M_{2j} and Hurwitz ζ values come from arb; ε is found by bisection and verified in ball arithmetic.

| deg B | Defect/g(0) | A′(12)/g(0) ≥ (rigorous) | no feasible κ < 12 − ε, ε = | largest share |
|---|---|---|---|---|
| 0 | [0.3583, 0.3584] | 0.1600 | 4.086 | p = 2 |
| 2 | [0.1827, 0.1828] | 0.1526 | 2.317 | 2 |
| 4 | [0.04739, 0.04740] | 0.1324 | 0.7236 | 5 |
| 6 | [8.043, 8.044]·10⁻³ | 0.1244 | 0.1293 | 7 |
| 8 | [5.402, 5.403]·10⁻⁴ | 0.1230 | 8.785·10⁻³ | 11 |
| 10 | [3.847, 3.848]·10⁻⁵ | 0.1224 | 6.285·10⁻⁴ | 13 |
| 12 | [2.111, 2.112]·10⁻⁶ | 0.1221 | 3.459·10⁻⁵ | 19 |
| 14 | [1.331, 1.332]·10⁻⁷ | 0.1219 | 2.185·10⁻⁶ | 19 |
| 16 | [2.323, 2.324]·10⁻⁸ | 0.1218 | 3.813·10⁻⁷ | 19 |
| 18 | [2.248, 2.249]·10⁻⁹ | 0.1217 | 3.694·10⁻⁸ | 5 |

All b_k > 0 at every degree. The polynomials come from [c12refine.py](c12refine.py) (SVD-coordinate LP; logs c12refine_d0_8.log, c12refine_d10_18.log), stored as cert12_d*.json.

**Does the construction fail because only one form is used? No (computed).**
- The worry was that min_θ P_p is not at θ_p for small p. But P_p contains p, p², p³, … tied to one measure. At degree 8, P_2 has its minimum on [−1, 1] at x = −0.268, next to x_2 = −0.2652, and Defect_2/g(0) = 1.2·10⁻⁷: the search tunes g(log 2)/g(log 4) so that the quadratic part of P_2 has its vertex at x_2.
- The defect is carried by the first primes whose powers are not yet controlled. It falls with every degree tried, more slowly than for the characters because g decays only like e^{−4π√n}.
- So the two-form version (L and Sym²) is not needed for the bound. **Sketch only:** add the Lemma-1 identity for Sym²Δ (conductor 1, Gamma factor Γ_R(s+1)Γ_C(s+11), entire by Shimura) with its own kernel Ξ_{Sym²}²B′. The local polynomial becomes P_p + P′_p, with P′_p built from 1 + 2cos 2mθ; this gives a second lever at each prime. Not done.

## 6. Checks (computed; validate_exact12.py section A)

- τ(n) from PARI equals the own q-expansion for n ≤ 150; Δ(i/y) = y^{12}Δ(iy) at y = 1.3, to 10⁻²⁵.
- Mellin Λ(8) against Γ_C(8 + 11/2)L(8, Δ), to 10⁻²⁰.
- PARI's own functional-equation check for the L-data (lfuncheckfeq ≤ −100).
- Ξ(t) from the own Mellin integral against PARI's lfunlambda(Δ, 6 + it) (arithmetic normalisation, Ξ(t) = Λ_PARI(6 + it)) at t = 0 and 5.5: 0.0030897587207900544 and 0.00050388323039607441, agreeing to about 19 digits.
- Closed form of φ*φ against quadrature at log 1, log 2, to 20 digits.
- Lemma 1 with A^{(12)} and g(0) on the Ξ side (PARI) against 2Σ_p P_p(x_p) from the closed form. Residuals −1.9·10⁻²⁶ (B = 1, A^{(12)} = −1.15·10⁻⁶) and −1.8·10⁻²⁶ (B = 1 + t²/50).
- Lemma 5 sampled at x = 2, 3, 6, y ≤ 199.5 (minimum 3.8·10⁻⁵ ≥ 0), and its failure at x = 0.3: Re ψ′(0.3 + 50i) = −8.0·10⁻⁵.

Not checked by a second code: the branch and bound, the tails and the moment bound. Result: 88 checks, 0 failed ([validate_exact12_out.txt](validate_exact12_out.txt)).

## 7. What is missing

- Exactness: inf_B Defect/g(0) = 0 is unproved (computed decay through degree 18, no stall seen), and so is the upper side (GRH and Deligne).
- A finite-support version (Lemma 3 with explicit constants).
- A second code for the certificate.
- Literature: no search was made for this weight version. The positivity method for weights is in print (Chenevier–Lannes, maass.md section 1).
