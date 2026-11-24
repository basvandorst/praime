# The mean variance of primes in (x, (1 + δ)x]: theorem and proof

> **In short:** The repo states that the mean over ln X of R(X, δ), the normalised variance of ψ over (x, (1 + δ)x], is a closed form 𝔇(δ). Proof: (1) under RH the mean exists and equals (1/δ) Σ_λ m(λ)² |c_δ(λ)|², over the distinct zeros with squared multiplicities; (2) with no hypothesis, the same sum with multiplicity m equals δ·𝔇(δ), by the explicit formula for Σ_{n≤y} Λ(n)(y/n − 1). Two corrections follow: under RH alone the mean equals 𝔇(δ) exactly when all zeros are simple, and for δ > 1 prime powers n ≤ 1 + δ enter 𝔇.

Proof note for RIEMANN.md §23 (rounds 95 and 102, where the formula has status "computed"), revised after a referee report. Every number in this note is recomputed and asserted by one script, [check_proof.py](check_proof.py), with its output in [check_proof_output.txt](check_proof_output.txt). The script imports nothing from the repo; its only repo input is the zero table `writeup/data/riemann/zeros_40000_flint.txt`. Each step carries a label: **proved** (a complete argument from the cited inputs), **computed** (numerical, with the error stated) or **heuristic**.

---

## 1. Statement

**Notation.** γ is Euler's constant and γ₁ the first Stieltjes constant: ζ(s) = 1/(s − 1) + γ − γ₁(s − 1) + …, γ₁ = −0.0728158454836767. Λ(n) is von Mangoldt's function. A nontrivial zero is ρ = ½ + iγ_ρ (γ_ρ is complex if ρ is off the line), counted with multiplicity. λ runs over the *distinct* ordinates, with multiplicity m(λ). A height is written T; the function T(y), always written with its argument, is the trivial-zero sum of §4.1. Two constants recur:

```
C = 2 + γ − ln 4π = 0.0461914179322421,        B = 1 − γ − ln 2π = −1.41509273131.
```

**The objects.** Put y = 1 + δ and

```
Δ_δ(x) = ψ((1 + δ)x) − ψ(x) − δx,          R(X, δ) = (2/(δX²)) ∫₁^X Δ_δ(x)² dx,
M(δ)   = lim_{L→∞} (1/L) ∫_a^{a+L} R(e^U, δ) dU                 (the mean over ln X),
c_δ(λ) = ((1 + δ)^(½+iλ) − 1)/(½ + iλ),     H_δ(s) = (y^s − 1)(y^(1−s) − 1)/(s(1 − s)),
𝔇(δ)   = ln(1/δ) − γ − ln 2π + 2 ln(1 + δ) + ((2 + δ)/δ) ln(1 + δ/2) − (2/δ) Σ_{n ≤ 1+δ} Λ(n) ((1 + δ)/n − 1).
```

H_δ is entire. The prime sum in 𝔇 is empty for 0 < δ < 1. At δ = 1 it holds n = 2 with weight (1 + δ)/2 − 1 = 0, so it contributes nothing. For 1 < δ ≤ 2 it adds −((δ − 1)/δ) ln 2.

**Theorem A (the mean exists; RH).** Assume RH and let δ > 0. Then R(X, δ) is bounded for X ≥ 1, the limit M(δ) exists uniformly in a ≥ 0, and

```
M(δ) = (1/δ) Σ_λ m(λ)² |c_δ(λ)|²        (finite; summed over the distinct ordinates, both signs).
```

**Theorem B (the closed form; no hypothesis).** For every δ > 0, summing over all nontrivial zeros with multiplicity (the sum converges absolutely),

```
Σ_ρ H_δ(ρ) = δ·𝔇(δ).
```

Under RH, 1 − ρ = ρ̄, so H_δ(ρ) = |c_δ(γ_ρ)|², and Theorem B says (1/δ) Σ_λ m(λ) |c_δ(λ)|² = 𝔇(δ).

**Corollary 1 (the role of simple zeros; RH).** Theorem A minus Theorem B:

```
M(δ) − 𝔇(δ) = (1/δ) Σ_λ m(λ)(m(λ) − 1) |c_δ(λ)|² ≥ 0.
```

Since c_δ(λ) ≠ 0 for every λ, equality at one δ > 0 holds if and only if every zero is simple, and then it holds at every δ. Whether all zeros are simple is not known, so "the mean is 𝔇(δ)" is conditional on it.

**The same theorem in the notation of the tube (fact 2048).** Let G(y) = Σ_{n≤y} Λ(n)(y/n − 1) − y ln y + (1 + γ)y and T(y) = Σ_{k≥1} y^(−2k)/(2k(2k + 1)), as in fact 2048. Then Theorem B reads (§4.1, Step 3)

```
δ·𝔇(δ) = (2 + δ)(2 + γ − ln 4π) − 2 [G(1 + δ) − ln 2π + T(1 + δ)].
```

For δ ≤ 1 this is the closed form of round 95. For δ > 1 the round-95 formula, which RIEMANN.md states for δ ≤ 1, does not extend: at δ = 2 it gives 0.4753, while the zeros (§7) and the formula above give 0.1287.

**Characters (§5).** Let χ be a primitive character mod q > 1 (so q is its conductor) with χ(−1) = (−1)^a. Theorem B_χ (Σ_ρ H_δ(ρ) over the zeros of L(s, χ) equals δ·𝔇_χ(δ), with 𝔇_χ given in §5) needs no hypothesis. Theorem A_χ (the mean of R_χ) needs GRH for L(s, χ), and so does reading B_χ as a sum of |c|².

**Where RH is used.** These are all the places:

- (R1) Lemma 3: |x^ρ| = √x, so Σ c_ρ x^ρ = √x Σ_λ m(λ) c_δ(λ) e^{iλ ln x} is a trigonometric series in ln x.
- (R2) Lemmas 1 and 2: the frequencies λ are real. This gives |∫_I e^{iθy} dy| ≤ min(|I|, 2/|θ|).
- (R3) The bounds (√(1 + δ) − 1)/|½ + iλ| ≤ |c_δ(λ)| ≤ K_δ/(1 + |λ|), K_δ = 2√2(1 + √(1 + δ)). Both use |(1 + δ)^ρ| = √(1 + δ). The upper bound drives Lemmas 1 to 3; the lower bound gives c_δ(λ) ≠ 0 in Corollary 1.
- (R4) H_δ(ρ) = |c_δ(γ_ρ)|². It turns Theorem B into a statement about Σ m|c|² (Corollary 1, C3, and the zero sums of §7).
- (R5) The bound |J(y)| ≤ C√y in C7, and the mean squares of C8.

Theorems B and B_χ, the converse in C7, and C2, C4, C5, C6 use no hypothesis. Not used anywhere: simple zeros (a condition in Corollary 1 only), linear independence of the ordinates, Montgomery's pair-correlation conjecture, and any bound on multiplicities beyond m(λ) ≪ ln|λ|, which follows from (I2). RH and GRH for Dirichlet L-functions have a formal proof (fact 2047, not yet refereed); the hypotheses are kept explicit here.

**Cited inputs.** H. Davenport, *Multiplicative Number Theory*, ch. 15–17 and 19.

- (I1) *Explicit formula* (Davenport ch. 17). For x > 1, ψ₀(x) = x − lim_{T→∞} Σ_{|γ_ρ|≤T} x^ρ/ρ − ln 2π − ½ ln(1 − x⁻²), zeros with multiplicity. For χ primitive mod q > 1 (ch. 19): ψ₀(x, χ) = −lim_T Σ_{|γ_ρ|≤T} x^ρ/ρ + A_χ(x), A_χ(x) = −(1 − a) ln x − b(χ) + Σ_{k≥1} x^(a−2k)/(2k − a), where b(χ) = L′/L(0, χ) for odd χ and b(χ) = B₀ (defined in §5) for even χ.
- (I2) *Zero counting* (ch. 15, 16). A unit interval of heights [t, t + 1] holds ≪ ln(|t| + 2) zeros (≪ ln(q(|t| + 2)) for L(s, χ)), and N(T) ≪ T ln T. Hence Σ_ρ 1/|ρ(1 − ρ)| < ∞.
- (I3) *Growth of ζ′/ζ* (ch. 17; ch. 19 for L′/L). There are heights T_j ∈ [j, j + 1] with ζ′/ζ(σ ± iT_j) ≪ ln² T_j for −1 ≤ σ ≤ 2. For σ ≤ −1 and |s + 2k| ≥ ½ (k = 0, 1, 2, …), ζ′/ζ(s) ≪ ln(|s| + 2).
- (I4) Only for the check route of §4.2: Weil's explicit formula (A. Weil 1952). For compactly supported continuous test functions that are only piecewise smooth, the reference is J. Barner, J. reine angew. Math. 323 (1981) 139–152; its exact conditions were not checked here, so §4.2 also gives a mollification argument.

---

## 2. Normalisation (read from the repo's code)

- RIEMANN.md §23 and `writeup/code/study/r91.py`, `r93.py`, `r96.py` use Δ_δ(x) = ψ((1 + δ)x) − ψ(x) − δx, with prime powers in ψ, I(X, δ) = ∫₁^X Δ_δ(x)² dx and R(X, δ) = I(X, δ)/(½δX²).
- The zero form of r93.py is R(X, δ) = (2/(δX²)) Re Σ_{ρ,ρ′} c_ρ c̄_ρ′ (X^(2 + i(γ_ρ − γ_ρ′)) − 1)/(2 + i(γ_ρ − γ_ρ′)), c_ρ = ((1 + δ)^ρ − 1)/ρ. Its diagonal ρ = ρ′ is (1/δ) Σ_ρ |c_ρ|² (1 − X⁻²).
- "Mean over ln X": r96.py averages R at 200 log-spaced X over three decades. That is a Riemann sum for (1/L) ∫_a^{a+L} R(e^U) dU with a finite L. Theorem A is about L → ∞.
- r95.py computes D(δ) = (2/δ) Σ_{γ>0} h_δ(γ), with h_δ(t) = |(1 + δ)^(½+it) − 1|²/(¼ + t²) = ((s − 1)² + 4s sin²(ut/2))/(¼ + t²), s = √(1 + δ), u = ln(1 + δ). Under RH this is (1/δ) Σ_ρ |c_ρ|², each zero counted with its multiplicity: it is m·|c|², not m²·|c|².
- For characters (`writeup/code/riemann/race/r102.py`) R_χ = ∫(ψ((1 + δ)x, χ) − ψ(x, χ))² dx/(½δX²), with no main term. For complex χ the square is |·|².

The diagonal of R is its mean only if the off-diagonal terms average out and only if the multiplicities are counted correctly. Theorem A settles both.

---

## 3. Proof of Theorem A

Assume RH. Let Λ* be the set of distinct ordinates (real, symmetric about 0). Write β(λ) = m(λ) c_δ(λ), w(λ) = m(λ)/(1 + |λ|) and κ(θ) = min(1, 2/|θ|).

**The bound on c_δ (proved; uses R3).** |(1 + δ)^(½+iλ) − 1| ≤ √(1 + δ) + 1, so |c_δ(λ)| ≤ (√(1 + δ) + 1)/|½ + iλ|. Also (1 + |λ|)² ≤ 8(¼ + λ²), because the difference 7λ² − 2|λ| + 1 has discriminant 4 − 28 < 0. So |c_δ(λ)| ≤ K_δ/(1 + |λ|) with K_δ = 2√2(1 + √(1 + δ)), and |β(λ)| ≤ K_δ w(λ). In the same way |(1 + δ)^(½+iλ) − 1| ≥ √(1 + δ) − 1 > 0 gives the lower bound in (R3).

**Lemma 1 (the pair sum converges; proved; uses R2).** W := Σ_{λ,λ′ ∈ Λ*} w(λ) w(λ′) κ(λ − λ′) < ∞.

*Proof.* By symmetry, count the pairs with |λ′| ≤ |λ| twice.
(i) Pairs with |λ − λ′| ≥ (1 + |λ|)/4. Then κ ≤ 8/(1 + |λ|), and these pairs contribute at most 16 Σ_λ m(λ)/(1 + |λ|)² · Σ_{|λ′|≤|λ|} m(λ′)/(1 + |λ′|). By (I2) and partial summation, Σ_{|λ′|≤T} m(λ′)/(1 + |λ′|) ≪ ln²(T + 2). The total is ≪ Σ_λ m(λ) ln²(|λ| + 2)/(1 + |λ|)², which converges by partial summation from N(T) ≪ T ln T.
(ii) Pairs with |λ − λ′| < (1 + |λ|)/4 and |λ′| ≤ |λ|. Then |λ′| > |λ| − (1 + |λ|)/4, so 1 + |λ′| > ¾(1 + |λ|), and w(λ)w(λ′) ≤ (4/3) m(λ) m(λ′)/(1 + |λ|)². Group the λ′ by j = ⌊|λ − λ′|⌋. For each j the λ′ lie in two unit intervals at height ≤ 2|λ| + 1, which hold ≪ ln(|λ| + 2) zeros (I2), and κ ≤ 2/j for j ≥ 1. So Σ_{λ′} m(λ′) κ(λ − λ′) ≪ ln(|λ| + 2)·(1 + Σ_{1≤j≤|λ|} 2/j) ≪ ln²(|λ| + 2). The contribution is ≪ Σ_λ m(λ) ln²(|λ| + 2)/(1 + |λ|)² < ∞. ∎

**Lemma 2 (mean square of the zero series; proved; uses R2).** Let Q_T(y) = Σ_{|λ|≤T} β(λ) e^{iλy}.
(a) For an interval I of length L ≥ 1 and a finite Λ₀ ⊂ Λ*: ∫_I |Σ_{λ∈Λ₀} β(λ) e^{iλy}|² dy ≤ L K_δ² Σ_{λ,λ′∈Λ₀} w(λ) w(λ′) κ(λ − λ′).
(b) Q_T converges in L²(I) for every bounded I to a limit Q, and sup_{|I|=L≥1} (1/L) ∫_I |Q − Q_T|² ≤ K_δ² W_T, where W_T = Σ_{|λ|,|λ′|>T} w(λ) w(λ′) κ(λ − λ′) → 0.
(c) sup_{a≥0} |(1/L) ∫_a^{a+L} |Q|² dy − Σ_λ |β(λ)|²| → 0 as L → ∞.

*Proof.* (a) Expand the square. For real θ, |∫_I e^{iθy} dy| ≤ min(L, 2/|θ|) ≤ L κ(θ) when L ≥ 1. (b) Apply (a) to Λ₀ = {T < |λ| ≤ T′}; W_T → 0 is the tail of the convergent series of Lemma 1. (c) Fix I = [a, a + L] and put k_I(θ) = (1/L) ∫_I e^{iθy} dy, so |k_I(θ)| ≤ min(1, 2/(L|θ|)) ≤ κ(θ). Then (1/L) ∫_I |Q_T|² = Σ_{|λ|,|λ′|≤T} β(λ) β̄(λ′) k_I(λ − λ′). As T → ∞ the left side tends to (1/L) ∫_I |Q|² by (b). The right side tends to the full double series, which converges absolutely by Lemma 1. The diagonal λ = λ′ gives Σ |β|², because k_I(0) = 1. The off-diagonal part is at most ε(L) := Σ_{λ≠λ′} |β(λ)||β(λ′)| min(1, 2/(L|λ − λ′|)), which does not depend on a. Each term tends to 0 as L → ∞ and is dominated by K_δ² w(λ) w(λ′) κ(λ − λ′), which is summable. So ε(L) → 0 by dominated convergence. ∎

**Lemma 3 (identification with the primes; proved; uses R1).** Put F(x) = Δ_δ(x)/√x. For almost every y > 0, F(e^y) = −Q(y) + E(y), where

```
E(y) = e^(−y/2) · ½ [ln(1 − e^(−2y)) − ln(1 − (1 + δ)^(−2) e^(−2y))],        ∫₀^∞ |E(y)|² dy < ∞.
```

*Proof.* Apply (I1) at (1 + δ)x and at x and subtract. The main terms cancel, (1 + δ)x − x − δx = 0, and so does ln 2π. Since (1 + δ)^ρ x^ρ/ρ − x^ρ/ρ = c_ρ x^ρ,

```
Δ_δ(x) = −lim_T Σ_{|γ_ρ|≤T} c_ρ x^ρ + √x·E(ln x)
```

for x > 1 whenever neither x nor (1 + δ)x is a prime power. Under RH x^ρ = √x e^{iγ_ρ ln x}, and grouping equal ordinates gives Σ_{|γ_ρ|≤T} c_ρ x^ρ = √x Q_T(ln x). So Q_T(y) → −F(e^y) + E(y) for every y > 0 outside a countable set. By Lemma 2(b) a subsequence Q_{T_j} converges almost everywhere to Q, so Q = −F + E almost everywhere. Near y = 0, E(y) ≪ 1 + |ln y|; for y ≥ 1, E(y) ≪ e^(−5y/2). Both are square-integrable. ∎

**Lemma 4 (from the local mean square to R; proved).** Let V ≥ 0 be locally integrable on [0, ∞) with (1/L) ∫_a^{a+L} V → M_V uniformly in a ≥ 0. Put R̃(U) = (2/δ) ∫₀^U V(y) e^(−2(U−y)) dy. Then R̃ is bounded, and (1/L) ∫_a^{a+L} R̃ → M_V/δ uniformly in a ≥ 0.

*Proof.* Choose L₀ with sup_a (1/L₀) ∫_a^{a+L₀} V ≤ M_V + 1. Then ∫_v^{v+1} V ≤ A₁ := L₀(M_V + 1) for every v ≥ 0. Cutting [0, U] into unit pieces from the right, R̃(U) ≤ (2/δ) Σ_{k≥0} e^(−2k) A₁ ≤ 3A₁/δ. By Fubini, with t = U − y,

```
(1/L) ∫_a^{a+L} R̃(U) dU = (2/δ) ∫₀^∞ e^(−2t) A_L(a, t) dt,       A_L(a, t) = (1/L) ∫_{max(a−t, 0)}^{max(a+L−t, 0)} V.
```

For t ≤ t₀ the window has length between L − t₀ and L, so A_L = M_V + o(1) uniformly in a as L → ∞. For all t, 0 ≤ A_L ≤ 2A₁. Split at t₀, let L → ∞ and then t₀ → ∞: the limit is (2/δ) M_V ∫₀^∞ e^(−2t) dt = M_V/δ. ∎

**Proof of Theorem A (proved).** Put V(y) = F(e^y)² = Δ_δ(e^y)² e^(−y). The substitution x = e^y gives Δ_δ(x)² dx = V(y) e^{2y} dy, so

```
R(e^U, δ) = (2/(δ e^{2U})) ∫₀^U V(y) e^{2y} dy = R̃(U).
```

By Lemma 3, V = (Q − E)² almost everywhere. Q is real, because β(−λ) = β̄(λ) (c_δ(−λ) = c̄_δ(λ), m(−λ) = m(λ)). By Lemma 2(c), (1/L) ∫_a^{a+L} Q² → Σ_λ |β(λ)|² uniformly in a. Also (1/L) ∫_a^{a+L} E² ≤ ‖E‖²/L, and the cross term is O(L^(−½)) by Cauchy–Schwarz. So V has the uniform mean M_V = Σ_λ |β(λ)|² = Σ_λ m(λ)² |c_δ(λ)|², which is finite because m(λ) ≪ ln(|λ| + 2) and |c_δ(λ)| ≪ 1/(1 + |λ|). Lemma 4 gives sup_X R(X, δ) < ∞ and M(δ) = M_V/δ, uniformly in a. ∎

*Remarks.* (1) The mean is a Weyl-type mean, uniform in the position of the window. It implies the Cesàro limit (1/ln X) ∫₁^X R(x, δ) dx/x → M(δ) (take a = 0). (2) RH enters through (R1), (R2) and (R3) only. Without RH an off-line zero would add a term of size x^(2 Re ρ − 1) to V. That R is then unbounded is plausible, but it is not proved here (§9). (3) This is the increment version of Cramér's mean-value theorem (H. Cramér, Math. Z. 12 (1922) 147–153), which under RH gives the mean of (ψ(x) − x)²/x² over ln x as a sum over zeros; the exact form of Cramér's statement was not checked here.

---

## 4. Proof of Theorem B

Absolute convergence: for 0 < Re ρ < 1 and y ≥ 1, |y^ρ − 1| ≤ y + 1 and |y^(1−ρ) − 1| ≤ y + 1, so |H_δ(ρ)| ≤ (y + 1)²/|ρ(1 − ρ)|, and Σ 1/|ρ(1 − ρ)| < ∞ by (I2).

### 4.1 The Riesz-mean route (in full)

**Step 1 (proved, algebra).** Since (y^s − 1)(y^(1−s) − 1) = y + 1 − y^s − y^(1−s) and y + 1 = 2 + δ,

```
H_δ(s) = [(2 + δ) − y^s − y^(1−s)]/(s(1 − s)).
```

The multiset of zeros is invariant under ρ ↦ 1 − ρ (functional equation ξ(s) = ξ(1 − s), multiplicities included). So Σ_ρ y^(1−ρ)/(ρ(1 − ρ)) = Σ_ρ y^ρ/(ρ(1 − ρ)). With

```
J(y) := Σ_ρ y^ρ/(ρ(1 − ρ))        (absolutely convergent; real, because the zeros are also symmetric under ρ ↦ ρ̄)
```

this gives

```
Σ_ρ H_δ(ρ) = (2 + δ) J(1) − 2 J(1 + δ).                                                    (4.1)
```

**Step 2 (proved; the Riesz-mean explicit formula, as in fact 2048).** For y ≥ 1 let P(y) := Σ_{n≤y} Λ(n)(y/n − 1) and T(y) := Σ_{k≥1} y^(−2k)/(2k(2k + 1)). Then

```
P(y) = y ln y − (1 + γ) y + ln 2π + J(y) − T(y).                                            (4.2)
```

*The kernel.* For z > 0, (1/2πi) ∫_{(2)} z^s ds/(s(s − 1)) = max(z − 1, 0). For z < 1 close the contour to the right: there are no poles, so the integral is 0. For z ≥ 1 close it to the left: the residues are z at s = 1 and −1 at s = 0, total z − 1.

*The integral.* The Dirichlet series −ζ′/ζ(s) = Σ Λ(n) n^(−s) converges absolutely on Re s = 2, and ∫ |ds|/|s(s − 1)| < ∞ there. So, applying the kernel with z = y/n and summing,

```
P(y) = (1/2πi) ∫_{(2)} (−ζ′/ζ)(s) y^s ds/(s(s − 1)).
```

*The contour shift.* Fix an integer K ≥ 1 and take the rectangle with vertical sides Re s = 2 and Re s = −2K − 1 and horizontal sides Im s = ±T_j, with T_j from (I3). On the horizontal sides |y^s/(s(s − 1))| ≤ y²/T_j², and |ζ′/ζ| ≪ ln² T_j there by (I3); the two sides contribute ≪ (2K + 3) y² ln² T_j/T_j², which tends to 0 as j → ∞. On Re s = −2K − 1 the distance to the trivial zeros −2K and −2K − 2 is 1, so |ζ′/ζ| ≪ ln(|s| + 2), and that side contributes ≪ y^(−2K−1) ∫ ln(|t| + 2K + 3) dt/((2K + 1)² + t²) ≪ y^(−2K−1) ln K/K, which tends to 0 as K → ∞ because y ≥ 1.

*The residues* of (−ζ′/ζ)(s) y^s/(s(s − 1)):

- s = 1. Here −ζ′/ζ(s) = 1/(s − 1) − γ + O(s − 1) and y^s/s = y(1 + (s − 1)(ln y − 1) + …). The double pole gives y ln y − (1 + γ) y.
- s = 0. Here −ζ′/ζ(0) = −ln 2π, and y^s/(s(s − 1)) has residue −1. Together: ln 2π.
- s = ρ, a zero of multiplicity m. −ζ′/ζ has residue −m, so the residue is −m y^ρ/(ρ(ρ − 1)) = m y^ρ/(ρ(1 − ρ)). The zeros with |γ_ρ| < T_j give the partial sums of J(y), which converges absolutely since |y^ρ/(ρ(1 − ρ))| ≤ y/|ρ(1 − ρ)|.
- s = −2k (k ≥ 1). −ζ′/ζ has residue −1, and y^(−2k)/((−2k)(−2k − 1)) = y^(−2k)/(2k(2k + 1)). The residue is −y^(−2k)/(2k(2k + 1)); together they give −T(y).

Summing gives (4.2). The script checks the kernel residues and the residue at s = 1 symbolically, and ζ′/ζ(0) = ln 2π numerically.

**Step 3 (proved, elementary).** At y = 1, P(1) = 0 and T(1) = Σ_k (1/(2k) − 1/(2k + 1)) = 1 − ln 2, so (4.2) gives

```
J(1) = (1 + γ) − ln 2π + 1 − ln 2 = 2 + γ − ln 4π = C.
```

Without RH, J(1) = Σ_ρ (1/ρ + 1/(1 − ρ)) = 2 Σ_ρ Re(1/ρ), summed in conjugate pairs, using ρ ↦ 1 − ρ. It equals Σ_ρ 1/|ρ|² only under RH (then 1 − ρ = ρ̄).

For y > 1 put w = 1/y. Then Σ w^(2k)/(2k) = −½ ln(1 − w²) and Σ w^(2k)/(2k + 1) = (artanh w − w)/w, so

```
T(y) = Σ_k w^(2k) (1/(2k) − 1/(2k + 1)) = 1 − ½ ln(1 − y⁻²) − (y/2) ln((y + 1)/(y − 1)).
```

Insert (4.2) into (4.1) with y = 1 + δ, and sort the terms by origin. With 1 − y⁻² = δ(2 + δ)/(1 + δ)² and (y + 1)/(y − 1) = (2 + δ)/δ:

| origin | contribution to Σ_ρ H_δ(ρ) = (2 + δ)J(1) − 2J(1 + δ) |
|---|---|
| pole at s = 1 | (2 + δ)(1 + γ) − 2(1 + γ)(1 + δ) + 2(1 + δ) ln(1 + δ) = −(1 + γ)δ + 2(1 + δ) ln(1 + δ) |
| s = 0 (ζ′/ζ(0) = ln 2π) | −(2 + δ) ln 2π + 2 ln 2π = −δ ln 2π |
| trivial zeros | (2 + δ)T(1) − 2T(1 + δ) = δ − (2 + δ) ln 2 + (2 + δ) ln(2 + δ) − δ ln δ − 2 ln(1 + δ) |
| prime powers n ≤ 1 + δ | −2P(1 + δ) |

The sum is δ[ln(1/δ) − γ − ln 2π + 2 ln(1 + δ) + ((2 + δ)/δ) ln(1 + δ/2)] − 2P(1 + δ) = δ·𝔇(δ). ∎

Since J(y) = G(y) − ln 2π + T(y) with G as in fact 2048, (4.1) is the formula δ·𝔇(δ) = (2 + δ)C − 2[G(1 + δ) − ln 2π + T(1 + δ)] of §1.

*Checks.* The script verifies the identity (2 + δ)J(1) − 2J(1 + δ) = δ·𝔇(δ) symbolically with sympy, for all δ > 0: P(1 + δ) is a free symbol, and every logarithm is written in ln δ, ln(1 + δ), ln(2 + δ), ln 2 and ln π, so the expanded difference is identically 0. It also checks the closed form of T(y) against its series (coefficients of w² … w¹²), and the fact-2048 form against 𝔇 numerically at δ = 0.3, 1, 2, 3.7, 10, 57.2 (difference below 1e-26).

### 4.2 Weil's explicit formula (a check)

h_δ(r) = ∫ g(x) e^{irx} dx with g(x) = ∫ f(v) f(v − x) dv and f = e^(v/2)·1_[0,u](v), u = ln(1 + δ). Explicitly,

```
g(x) = (1 + δ) e^(−|x|/2) − e^(|x|/2)   for |x| ≤ u,        g(x) = 0 otherwise.
```

So g is even, continuous, piecewise smooth and supported in [−u, u], with g(0) = δ and g(±u) = 0. As an entire function of r, h_δ(r) = H_δ(½ + ir), so h_δ(γ_ρ) = H_δ(ρ) for every zero. Weil's formula reads

```
Σ_ρ h(γ_ρ) = h(i/2) + h(−i/2) − g(0) ln π + (1/2π) ∫ h(r) Re ψ(¼ + ir/2) dr − 2 Σ_n Λ(n) n^(−½) g(ln n).
```

Its terms:

- *Pole:* h(±i/2) = H_δ(1) = H_δ(0) = δ ln(1 + δ), from lim_{s→1} (y^(1−s) − 1)/(1 − s) = ln y. Together: **2δ ln(1 + δ)**.
- *π^(−s/2):* −g(0) ln π = **−δ ln π**.
- *Primes:* for n ≤ 1 + δ, g(ln n) = (1 + δ)n^(−½) − n^(½), so −2Λ(n) n^(−½) g(ln n) = −2Λ(n)((1 + δ)/n − 1). The sum is **−2P(1 + δ)**. It is empty for δ < 1; at δ = 1 it holds n = 2 with g(ln 2) = 0, so it contributes nothing.
- *Γ(s/2):* Gauss's integral ψ(z) = −γ + ∫₀^∞ (e^(−t) − e^(−zt))/(1 − e^(−t)) dt (Re z > 0), with t = 2x and (1/2π) ∫ h(r) cos(rx) dr = g(x), gives

  ```
  (1/2π) ∫ h(r) Re ψ(α + ir/2) dr = −γ g(0) + 2 ∫₀^∞ [g(0) e^(−2x) − e^(−2αx) g(x)]/(1 − e^(−2x)) dx.      (4.3)
  ```

  The interchange of integrals is justified by splitting the bracket as [e^(−2x) − e^(−2αx)] g(0) + e^(−2αx)(g(0) − g(x)). The first part is absolutely integrable. The second is ≥ 0 and is handled by Tonelli, using h ≥ 0 and ∫ h(r)(1 − cos rx) dr = 2π(g(0) − g(x)) = O(x).
  For α = ¼ on [0, u], substitute z = e^(−x). The integrand becomes (δz² − (1 + δ)z + 1)/(1 − z²) = (1 − δz)/(1 + z), and dx = −dz/z. Since (1 − δz)/(z(1 + z)) = 1/z − (1 + δ)/(1 + z),

  ```
  ∫₀^u = ∫_{1/(1+δ)}^1 (1/z − (1 + δ)/(1 + z)) dz = (1 + δ) ln(1 + δ/2) − δ ln(1 + δ).
  ```

  On [u, ∞) only g(0)e^(−2x) remains: ∫_u^∞ δ e^(−2x)/(1 − e^(−2x)) dx = −(δ/2) ln(1 − (1 + δ)⁻²). With ln(1 − (1 + δ)⁻²) = ln δ + ln(2 + δ) − 2 ln(1 + δ), (4.3) at α = ¼ equals **−γδ − δ ln 2δ + (2 + δ) ln(1 + δ/2)**.

The four bold terms add up to δ·𝔇(δ). The script checks (4.3) as an x-integral against these closed forms at α = ¼ and ¾ for seven δ in [0.01, 10] (to 1e-20), and the whole Weil sum, with the digamma integral done by quadrature in r, against 𝔇 at δ = 0.01, 0.1, 1, 3 (to 1e-17).

*Non-smooth g.* If one prefers the smooth-test-function version of Weil's formula, apply it to g_ε = g ∗ φ_ε, with φ_ε a smooth even mollifier supported in [−ε, ε], so that h_ε = h·φ̂_ε. As ε → 0: Σ_ρ h_ε(γ_ρ) → Σ_ρ h(γ_ρ), because |φ̂_ε(γ_ρ)| ≤ e^(ε/2) for |Im γ_ρ| ≤ ½ and Σ |h(γ_ρ)| < ∞; the pole and Γ terms converge by dominated convergence; and g_ε(ln n) → g(ln n) for the finitely many n involved. When 1 + δ is a prime power, the term n = 1 + δ tends to g(ln n) = 0.

### 4.3 What the two routes share, and where the constants come from

The two derivations are not independent in their inputs: both rest on the zero counting (I2), the bounds for ζ′/ζ (I3) and the functional equation, which also underlie Weil's formula. They differ in the test function and in the bookkeeping, so they check each other's algebra, not each other's foundations.

- *Where each constant sits depends on the route.* In the Weil split, 2 ln(1 + δ) is exactly the pole term. In the Riesz split, the pole gives (2(1 + δ)/δ) ln(1 + δ) − 1 − γ and the trivial zeros give back −(2/δ) ln(1 + δ). The route-free statement is the comparison with an even primitive character, which has the same gamma factor and no pole: for δ ≤ 1, 𝔇_χ − ln q = 𝔇 − 2 ln(1 + δ) (§5). In that sense "2 ln(1 + δ) is the pole".
- *Euler's γ* is the constant of Gauss's integral for ψ (Weil route), or the Laurent constant of ζ at 1 (Riesz route).
- *ln 2π* is ln π (from π^(−s/2)) plus ln 2 (from Re ψ(¼ + ir/2) ≈ ln(r/2)), i.e. the 2π of N(T) = (T/2π) ln(T/2πe).
- *ln(1/δ)* comes from the lags x > ln(1 + δ) in (4.3), i.e. from the density of zeros at height about 1/δ. This is the Goldston–Montgomery leading term.

---

## 5. Characters (round 102)

Let χ be a primitive character mod q > 1; then q is the conductor of χ. Let χ(−1) = (−1)^a with a ∈ {0, 1}. The completed L-function Λ(s, χ) = (q/π)^((s+a)/2) Γ((s + a)/2) L(s, χ) satisfies Λ(s, χ) = ε(χ) Λ(1 − s, χ̄) and has no pole. Define Δ_χ(x) = ψ((1 + δ)x, χ) − ψ(x, χ), with no main term, and R_χ(X, δ) = (2/(δX²)) ∫₁^X |Δ_χ(x)|² dx.

*Imprimitive characters.* If χ mod q is induced by the primitive χ* mod q*, then L(s, χ) = L(s, χ*) Π_{p|q} (1 − χ*(p) p^(−s)). The extra factors have zeros on Re s = 0, which are not zeros of L(s, χ*), and ψ(x, χ) differs from ψ(x, χ*) by O(ln x), which does not change the mean of R_χ. So for imprimitive χ every formula below holds with q* and χ* in place of q and χ.

The closed form:

```
𝔇_χ(δ) = ln(q/(2πδ)) − γ + Ω_a(δ) − (2/δ) Σ_{n ≤ 1+δ} Λ(n) Re χ(n) ((1 + δ)/n − 1),
Ω₀(δ) = ((2 + δ)/δ) ln(1 + δ/2),        Ω₁(δ) = (2(1 + δ)/δ) ln(1 + δ) − ((2 + δ)/δ) ln(1 + δ/2).
```

**Theorem B_χ (proved; no hypothesis).** Σ_ρ H_δ(ρ) over the nontrivial zeros of L(s, χ), with multiplicity, equals δ·𝔇_χ(δ). The sum is real.

**Theorem A_χ (proved; GRH for L(s, χ)).** sup_X R_χ(X, δ) < ∞, and the uniform mean of R_χ exists and equals (1/δ) Σ_λ m_χ(λ)² |c_δ(λ)|², summed over the distinct ordinates of the zeros of L(s, χ). Under GRH, 1 − ρ = ρ̄, so B_χ says (1/δ) Σ_λ m_χ(λ) |c_δ(λ)|² = 𝔇_χ(δ), and Corollary 1 holds word for word.

*Proof of A_χ.* That of Theorem A, with four changes. (I1) is the χ-version. The ordinates need not be symmetric for complex χ, and an ordinate λ = 0 (a zero at s = ½) is allowed; Lemma 2(c) does not need symmetry, and the bound |c_δ(λ)| ≤ K_δ/(1 + |λ|) holds at λ = 0. Q is now complex and V = |Q − E_χ|², with E_χ(y) = e^(−y/2)[A_χ((1 + δ)e^y) − A_χ(e^y)] in L²(0, ∞): for even χ the term −ln x of A_χ contributes the constant −ln(1 + δ), which is O(x^(−½)) after dividing by √x, and the logarithmic singularity at x = 1 is square-integrable. Lemma 1 uses N(t + 1, χ) − N(t, χ) ≪ ln(q(|t| + 2)). ∎

*Proof of B_χ: the Riesz-mean route (in full).* Since L(ρ, χ) = 0 implies L(1 − ρ̄, χ) = 0 (functional equation, then complex conjugation), the zeros are invariant under ρ ↦ 1 − ρ̄, and H_δ(1 − ρ̄) = conj H_δ(ρ). So with J_χ(y) := Σ_ρ y^ρ/(ρ(1 − ρ)),

```
Σ_ρ H_δ(ρ) = (2 + δ) J_χ(1) − 2 Re J_χ(1 + δ),        J_χ(1) real.
```

The contour shift of §4.1, with (−L′/L)(s, χ) y^s/(s(s − 1)) and P_χ(y) := Σ_{n≤y} Λ(n) χ(n)(y/n − 1), has these residues: s = 1 gives −y·L′/L(1, χ) (simple pole, L′/L is regular at 1); the zeros give J_χ(y); the trivial zeros give −T_o(y) (odd χ, zeros at −1, −3, …) or −T(y) (even χ, zeros at −2, −4, …), with

```
T_o(y) = Σ_{k≥1} y^(1−2k)/((2k − 1)2k) = ½ ln((y + 1)/(y − 1)) + (y/2) ln(1 − y⁻²),        T_o(1) = ln 2.
```

At s = 0: for odd χ, L(0, χ) ≠ 0 and the residue is L′/L(0, χ). For even χ, L(s, χ) has a simple zero at 0, −L′/L(s) = −1/s − B₀ + O(s) with B₀ = lim_{s→0} (L′/L(s, χ) − 1/s), and y^s/(s(s − 1)) = −(1/s)(1 + s(ln y + 1) + …); the residue is ln y + 1 + B₀. So

```
odd χ:   J_χ(y) = P_χ(y) + y L′/L(1, χ) − L′/L(0, χ) + T_o(y),
even χ:  J_χ(y) = P_χ(y) + y L′/L(1, χ) − ln y − 1 − B₀ + T(y).
```

The logarithmic derivative of the functional equation at s = 0, with ψ(½) = −γ − 2 ln 2 and ψ(1) = −γ, gives

```
odd χ:   L′/L(0, χ) + L′/L(1, χ̄) = −ln(q/π) + γ + ln 2,
even χ:  B₀ + L′/L(1, χ̄) = −ln(q/π) + γ + ln 2,
```

and L′/L(1, χ̄) is the complex conjugate of L′/L(1, χ). Inserting: J_χ(1) = 2 Re L′/L(1, χ) + ln(q/π) − γ for odd χ (minus a further 2 ln 2 for even χ), which is real, and in (2 + δ)J_χ(1) − 2 Re J_χ(1 + δ) the coefficient of Re L′/L(1, χ) is 2(2 + δ) − 2(1 + (1 + δ)) = 0. What remains is δ·𝔇_χ(δ). The script checks this symbolically for both parities (Re and Im of L′/L(1, χ) as free real symbols). ∎

*Check: the Weil route.* For L(s, χ) Weil's formula has no pole term, the archimedean factor ln(q/π) + Re ψ((½ + a)/2 + ir/2), and the prime term −Σ_n Λ(n) n^(−½)[χ(n) g(ln n) + χ̄(n) g(−ln n)] = −2 Σ_n Λ(n) Re χ(n) n^(−½) g(ln n). The ln(q/π) gives δ ln(q/π). For α = ¼ (even), (4.3) is as for ζ. For α = ¾ (odd), e^(−3x/2) g(x) = (1 + δ)e^(−2x) − e^(−x) on [0, u], so the integrand there is (e^(−x) − e^(−2x))/(1 − e^(−2x)) = e^(−x)/(1 + e^(−x)), whose integral is ln(1 + δ) − ln(1 + δ/2); the part on [u, ∞) is unchanged. Then (4.3) = −γδ − δ ln 2δ + 2(1 + δ) ln(1 + δ) − (2 + δ) ln(1 + δ/2) = δ(Ω₁(δ) − γ − ln 2δ). The script checks the α = ¾ digamma integral by quadrature in r at δ = 0.1 and 1 (to 1e-17).

*Consequences (proved).*

- For δ ≤ 1 and even χ: 𝔇_χ = 𝔇 − 2 ln(1 + δ) + ln q.
- Odd minus even, same q: (2/δ)[(1 + δ) ln(1 + δ) − (2 + δ) ln(1 + δ/2)] = δ/2 − δ²/4 + 7δ³/48 − ….
- Small δ: 𝔇_χ = ln(q/δ) + B + δ/4 − δ²/24 + … (even) or ln(q/δ) + B + 3δ/4 − 7δ²/24 + … (odd), against 𝔇 = ln(1/δ) + B + (9/4)δ + … for ζ. Here 9/4 = 2 (pole) + ¼ (even Γ factor).
- χ₄ and χ₃ are both odd, so 𝔇_χ₄ − 𝔇_χ₃ = ln(4/3) = 0.287682072452 exactly for δ ≤ 1. For 1 < δ ≤ 2, χ₃ gains +((δ − 1)/δ) ln 2 (χ₃(2) = −1) and χ₄ gains nothing (χ₄(2) = 0).
- The repo's values 2.346, 2.058 (δ = 0.1) and 4.584, 4.296 (δ = 0.01) for χ₄, χ₃ are reproduced: 2.3460172, 2.0583352, 4.5838428, 4.2961607.

---

## 6. Corollaries

**C1 (simple zeros; proved, RH).** This is Corollary 1 of §1. A zero of multiplicity m₀ ≥ 2 at ±λ₀ raises the mean by at least (2/δ) m₀(m₀ − 1)(√(1 + δ) − 1)²/(¼ + λ₀²), by the lower bound in (R3). The repo's sentence "a zero of multiplicity m would add (m² − m)|c_ρ|²" is right, provided 𝔇 counts zeros with multiplicity m.

**C2 (small δ; proved, exact rational coefficients by sympy).** For 0 < δ ≤ 1 the prime sum contributes nothing, and

```
𝔇(δ) − ln(1/δ) − B  = (9/4)δ − (25/24)δ² + (65/96)δ³ − (161/320)δ⁴ + …,
𝔇(δ) − HL_exact(δ)  = (7/4)δ − (7/8)δ² + (19/32)δ³ − (29/64)δ⁴ + (353/960)δ⁵ − …,
HL_exact(δ) := ln(1/δ) − γ − ln 2π + (1 + δ) ln(1 + δ)/δ          (the repo's Hardy–Littlewood value for finite δ).
```

The two series are different expansions and must not be mixed. The prime-free formula minus ln(1/δ) is analytic in |δ| < 1; its nearest singularity is that of ln(1 + δ) at δ = −1. By Corollary 1, 𝔇 is the mean of R as δ → 0 if the zeros are simple. This agrees with Goldston–Montgomery's leading term ln(1/δ) and with the constant B (§8). Note that Goldston–Montgomery average over x at fixed X with δ ≥ X^(−1+ε), while here δ is fixed and the average is over ln X.

**C3 (where 7/4 and 7/8 come from; proved under RH, 0 < δ ≤ 1).**

- 7/4 = 2 − ¼. The 2 is the pole: 2 ln(1 + δ) = 2δ − δ² + …. The −¼ is the Γ part beyond HL_exact: ((2 + δ)/δ) ln(1 + δ/2) − (1 + δ) ln(1 + δ)/δ = −δ/4 + δ²/8 − 7δ³/96 + …. The δ² coefficient −7/8 = −1 + ⅛ comes from the same two sources.

- *The link to N(T).* Write N(t) = θ(t)/π + 1 + S(t) for t > 0, with N(0+) = 0, so S(0+) = −1. Let N₀(t) = (t/2π) ln(t/2πe), D₀(δ) = (2/δ) ∫₀^∞ h_δ(t) N₀′(t) dt (the smooth-density diagonal of round 93) and ε(t) = θ(t)/π − N₀(t) + ⅛. ε is bounded, ε(0) = ⅛ and ε(t) = O(1/t). Under RH, 𝔇 = (2/δ) ∫_(0,∞) h_δ dN, and on (0, ∞), dN = N₀′ dt + ε′ dt + dS. Integrating ∫ h_δ ε′ dt by parts (h_δ → 0 at ∞):

  ```
  𝔇(δ) − D₀(δ) = (7/8)·(2/δ) h_δ(0)  +  [(2/δ) ∫_(0,∞) h_δ dS − (2/δ) h_δ(0)]  −  (2/δ) ∫₀^∞ h_δ′ ε dt,
  (2/δ) h_δ(0) = (8/δ)(√(1 + δ) − 1)² = 2δ − δ² + (5/8)δ³ − ….
  ```

  The S part: Weil's formula is exactly this split, its archimedean term being (2/π) ∫₀^∞ h_δ θ′ dt, so (2/δ) ∫_(0,∞) h_δ dS = 2 ln(1 + δ) − (2/δ)P(1 + δ). For δ ≤ 1 the bracket is 2 ln(1 + δ) − (8/δ)(√(1 + δ) − 1)² = δ³/24 − δ⁴/16 + …. The ε part: since (2/δ) ∫ h_δ θ′/π dt is the archimedean term divided by δ, −(2/δ) ∫ h_δ′ ε dt = [((2 + δ)/δ) ln(1 + δ/2) − (1 + δ) ln(1 + δ)/δ] + (1/δ)(√(1 + δ) − 1)² − (D₀ − HL_exact) = δ²/144 − δ³/576 + …, using the next item. Hence 𝔇 − D₀ = (7/8)(2/δ)h_δ(0) + O(δ²) = (7/4)δ + O(δ²): the 7/8 of N(T) − N₀(T) = 1 (pole) − ⅛ (Γ) is the coefficient, and S(t) enters only at δ³.

- *D₀ − HL_exact (proved; series computed exactly).* With ln|t| = ∫₀^∞ (e^(−x) − cos tx)/x dx (t ≠ 0), (1/2π) ∫ h_δ(t) dt = g(0) = δ and (1/2π) ∫ h_δ(t) cos(tx) dt = g(x), with g from §4.2:

  ```
  D₀(δ) = (1/(2πδ)) ∫_ℝ h_δ(t) ln(|t|/2π) dt = −ln 2π + (1/δ) ∫₀^u (δ e^(−x) − g(x))/x dx + E₁(u),       u = ln(1 + δ).
  ```

  The interchange of the t- and x-integrals is first done on [0, M] in x, where |e^(−x) − cos tx|/x ≤ 1 + min(2/x, t²x/2) and h_δ(t) ≪ 1/(1 + t²) make it absolutely integrable. Then M → ∞: on the t side, ∫₀^M (e^(−x) − cos tx)/x dx = ln|t| − Ci(|t|M) − E₁(M), and ∫ h_δ(t) Ci(|t|M) dt → 0 because Ci(z) = O(1/z) as z → ∞, Ci(z) = O(1 + |ln z|) near 0, and h_δ is bounded and integrable. On the x side g = 0 beyond u, which leaves δ E₁(u). Expanding the integrand in x, integrating term by term, and using E₁(u) = −γ − ln u + Σ_{k≥1} (−1)^(k+1) u^k/(k·k!), sympy gives

  ```
  D₀(δ) − HL_exact(δ) = −δ²/144 + δ³/144 − (703/115200) δ⁴ + O(δ⁵).
  ```

  Numerically (two independent quadratures of D₀, in t and in x, agree to 1e-15): D₀ − HL_exact = −6.306e-5, −1.653e-5, −2.723e-6 at δ = 0.1, 0.05, 0.02. The three-term series gives −6.311e-5, −1.653e-5, −2.723e-6; −δ²/144 alone gives −6.944e-5, −1.736e-5, −2.778e-6. This replaces the earlier "O(δ²) (proved)" claim, which had no computation behind it.

- *An integral of S(t) (proved under RH).* ∫₀^∞ h_δ′(t) S(t) dt = 4(√(1 + δ) − 1)² − δ ln(1 + δ) + P(1 + δ) for every δ > 0. This is the S part above, after one integration by parts ((2/δ)∫_(0,∞) h_δ dS = (2/δ)h_δ(0) − (2/δ)∫ h_δ′ S dt). For δ ≤ 1 the P term vanishes. It is the repo's "S(t) contributes exactly 2 ln(1 + δ) − (8/δ)(√(1 + δ) − 1)²".

**C4 (δ = 1; proved).** 𝔇(1) = 3 ln 3 − ln 4π − γ = ln(27/(4π)) − γ = 0.18759695413350542. The prime 2 sits at the edge of the support, where g(ln 2) = 0, so δ = 1 is still given by the prime-free formula.

**C5 (shape on (0, 1]; proved).** Differentiating,

```
𝔇′(δ) = 2/(1 + δ) − (2/δ²) ln(1 + δ/2) = 2φ(δ)/(δ²(1 + δ)),        φ(δ) = δ² − (1 + δ) ln(1 + δ/2).
```

φ(0) = 0, φ′(0) = −½ and φ″(δ) = 2 − 1/(2 + δ) − 1/(2 + δ)² ≥ 5/4 > 0. So φ has exactly one zero δ₀ in (0, ∞). Hence 𝔇 decreases strictly on (0, δ₀] and increases strictly on [δ₀, 1]. Numerically δ₀ = 0.7403218 and min 𝔇 = 𝔇(δ₀) = 0.1594331 (computed).

**C6 (δ > 1: the prime powers enter; proved).** 𝔇 is continuous on (0, ∞) and real-analytic between consecutive values δ = n − 1, n a prime power. At δ = n − 1 the slope jumps by −2Λ(n)/(n(n − 1)), the derivative of −(2/δ)Λ(n)((1 + δ)/n − 1) at its onset. Examples:

- 1 < δ ≤ 2: 𝔇(δ) = [prime-free formula] − ((δ − 1)/δ) ln 2.
- 2 < δ ≤ 3: subtract also (2/δ) ln 3 · (δ − 2)/3.
- At δ = 1: 𝔇′(1−) = 1 − 2 ln(3/2) = 0.18907 and 𝔇′(1+) = 0.18907 − ln 2 = −0.50408, so δ = 1 is a local maximum (a corner).

The prime term at δ = 1.5, 2, 3, 5, 10 is −0.231, −0.347, −0.706, −1.261, −2.162; the zero sums of §7 confirm the formula there.

**C7 (large δ and an RH criterion).** By (4.1), δ·𝔇(δ) = (2 + δ)C − 2J(1 + δ), with J(y) = G(y) − ln 2π + T(y) = Σ_ρ y^ρ/(ρ(1 − ρ)) the tube function of fact 2048 (G as in §1: G(y) = P(y) − y ln y + (1 + γ)y, P(y) = Σ_{n≤y} Λ(n)(y/n − 1)).

- *(a) The limit (proved, no hypothesis).* J(y)/y = Σ_ρ y^(ρ−1)/(ρ(1 − ρ)). Every term tends to 0 because Re ρ < 1 (ζ(1 + it) ≠ 0), and |y^(ρ−1)| ≤ 1 with Σ 1/|ρ(1 − ρ)| < ∞. By dominated convergence J(y) = o(y), so 𝔇(δ) → C = 0.0461914 as δ → ∞.
- *(b) The bounds (proved under RH).* Under RH |y^ρ| = √y and ρ(1 − ρ) = |ρ|², so |J(y)| ≤ √y Σ_ρ 1/|ρ|² = C√y. From δ𝔇 = (2 + δ)C − 2J(1 + δ) and (2 + δ) ∓ 2√(1 + δ) = (√(1 + δ) ∓ 1)²:

  ```
  C (√(1 + δ) − 1)² ≤ δ·𝔇(δ) ≤ C (√(1 + δ) + 1)²        for every δ > 0,
  ```

  so 𝔇(δ) = C + O(δ^(−½)).
- *(c) The converse (proved).* Fix any δ₁ > 0. Each of the two inequalities, required only for all δ ≥ δ₁, implies RH:

  ```
  (L)  δ·𝔇(δ) ≥ C (√(1 + δ) − 1)²  for all δ ≥ δ₁      ⟺   J(y) ≤ C√y   for all y ≥ y₁ = 1 + δ₁,
  (U)  δ·𝔇(δ) ≤ C (√(1 + δ) + 1)²  for all δ ≥ δ₁      ⟺   J(y) ≥ −C√y  for all y ≥ y₁.
  ```

  So RH ⟺ (L) ⟺ (U). *Proof for (L).* Put f(y) = C√y − J(y), which is ≥ 0 on [y₁, ∞).
  Step 1. For Re s > 1, ∫_{y₁}^∞ J(y) y^(−s−1) dy = Σ_ρ y₁^(ρ−s)/(ρ(1 − ρ)(s − ρ)) =: Φ_J(s). The interchange is allowed because |y^ρ| ≤ y, Σ 1/|ρ(1 − ρ)| < ∞ and ∫_{y₁}^∞ y^(−Re s) dy < ∞. On a compact set avoiding the zeros, |y₁^(ρ−s)| ≤ y₁^(1−Re s) and |s − ρ| is bounded below, so the series converges uniformly. Hence Φ_J is meromorphic on ℂ, with simple poles exactly at the distinct zeros ρ₀, residue m₀/(ρ₀(1 − ρ₀)) ≠ 0 (since 0 < Re ρ₀ < 1).
  Step 2. Since ∫_{y₁}^∞ C y^(½−s−1) dy = C y₁^(½−s)/(s − ½) for Re s > ½, for Re s > 1
  ```
  Φ_f(s) := ∫_{y₁}^∞ f(y) y^(−s−1) dy = C y₁^(½−s)/(s − ½) − Φ_J(s).
  ```
  The right side is meromorphic on ℂ with poles at s = ½ and at the zeros (ζ(½) ≠ 0, so they do not merge). ζ has no zero on the real segment (0, 1): there (1 − 2^(1−σ))ζ(σ) = Σ (−1)^(n−1) n^(−σ) > 0 and 1 − 2^(1−σ) < 0. So the right side is holomorphic on the real half-line (½, ∞).
  Step 3 (Landau's theorem). With y = e^v, Φ_f(s) = ∫_{ln y₁}^∞ f(e^v) e^(−sv) dv is the Laplace transform of a function ≥ 0. Its abscissa of convergence σ_c is at most 1. By Landau's theorem (for Laplace integrals of a function ≥ 0; e.g. D. V. Widder, *The Laplace Transform*, ch. II), if σ_c is finite, the function defined by the integral cannot be continued analytically to any neighbourhood of the real point σ_c. If σ_c > ½, the right side of Step 2 would be such a continuation. So σ_c ≤ ½: the integral converges absolutely and is holomorphic on Re s > ½.
  Step 4. The two sides agree on Re s > 1, hence on Re s > ½ minus the poles (a connected set), so the right side has no pole in Re s > ½. A zero ρ₀ with Re ρ₀ > ½ would be a pole there with residue −m₀/(ρ₀(1 − ρ₀)) ≠ 0, because C y₁^(½−s)/(s − ½) is holomorphic at ρ₀. So every zero has Re ρ ≤ ½, and by ρ ↦ 1 − ρ, Re ρ = ½.
  For (U) take f = J + C√y; the residue at ρ₀ is +m₀/(ρ₀(1 − ρ₀)), and the argument is the same. ∎
  This is the Mellin–Landau argument of RIEMANN.md §22, made one-sided. The value of C plays no role in the converse. Only large δ matter: dropping the hypothesis on any bounded range of δ does not change the conclusion. So a check of the bounds on a bounded range of δ says nothing about RH; it tests only the arithmetic.
- *(d) The approach to C (computed).* On a grid of step 0.05 on (0, 60], the smallest value of 𝔇 is 0.04320, at δ = 33.5, below C: 𝔇 oscillates around C for large δ. Both bounds of (b) hold on the grid of step 0.05 on (0, 100] (2000 points), as they must under RH.

**C8 (the Riesz-mean, "tube", versions; proved, checked against the zeros).**

- *Mean square of the tube function.* Under RH, J(e^v) e^(−v/2) = Σ_λ m(λ) e^{iλv}/(¼ + λ²) converges absolutely, so it is a Bohr almost periodic function, and the mean of its square over v is Σ_λ m(λ)²/(¼ + λ²)². With simple zeros this is

  ```
  K₂ = Σ_ρ 1/(ρ(1 − ρ))² = 6 + 2γ + 2γ² + 4γ₁ − π²/4 − 2 ln 4π = 7.42012728749297e-5.
  ```

  *Proof.* Repeat Step 2 of §4.1 with the kernel 1/(s(s − 1))², whose residues give (1/2πi) ∫_{(2)} z^s ds/(s(s − 1))² = (z + 1) ln z − 2(z − 1) for z ≥ 1 (and 0 for z < 1). With P₂(y) = Σ_{n≤y} Λ(n)[(y/n + 1) ln(y/n) − 2(y/n − 1)] and T₂(y) = Σ_{k≥1} y^(−2k)/(4k²(2k + 1)²) this gives J₂(y) := Σ_ρ y^ρ/(ρ(1 − ρ))² as

  ```
  J₂(y) = y[½ ln²y − (2 + γ) ln y + 3 + 2γ + γ² + 2γ₁] − ln 2π (ln y + 2) + 2γ₁ + γ² − π²/12 − T₂(y) − P₂(y).
  ```

  The first bracket is the residue at s = 1, from ζ′/ζ(s) = −1/(s − 1) + γ − (2γ₁ + γ²)(s − 1) + …. The middle terms are the residue at s = 0, from ζ′/ζ(0) = ln 2π and (ζ′/ζ)′(0) = π²/12 − 2γ₁ − γ² (from the functional equation). T₂(1) = Σ (1/(2k) − 1/(2k + 1))² = π²/6 − 3 + 2 ln 2. Then K₂ = J₂(1). Also K₂ = Σ (1/ρ + 1/(1 − ρ))² = 2σ₂ + 2C, with σ₂ = Σ_ρ ρ⁻² = 1 + γ² + 2γ₁ − π²/8 = −0.0461543172958. This is the classical value of σ₂; its attribution to Lehmer and Keiper is from memory and was not checked. The script checks the Laurent data of ζ′/ζ at 1 and 0, the residues and T₂(1). The RMS of the tube deviation is √K₂ = 0.00861 = 0.1865·C, against the sup bound C.
- *Increments.* Under RH the mean over ln x of |J((1 + δ)x) − J(x)|²/(δx) is (1/δ) Σ_λ m(λ)² |c_δ(λ)|²/(¼ + λ²). With simple zeros it equals Dt(δ) := (1/δ) Σ_ρ H_δ(ρ)/(ρ(1 − ρ)) = [(2 + δ)K₂ − 2J₂(1 + δ)]/δ; the last equality holds without hypothesis, by Step 1 of §4.1. Here too the primes enter only for δ > 1. Dt(δ)/δ → C as δ → 0, by dominated convergence (|c_ρ|²/δ² ≤ 1 + δ). The approach is slow: Dt(0.01)/0.01 = 0.0393.

**C9 (truncation of the zero sum; heuristic, checked numerically).** Truncate the zero sum at the height T of the last zero used, and add the smooth part and the boundary term as in §7. Then, exactly (RH),

```
closed − zeros = −(2/δ) ∫_T^∞ h_δ′(t) S(t) dt,        S(t) = N(t) − θ(t)/π − 1.
```

*Heuristic.* Insert the formal series S(t) = −(1/π) Σ_n Λ(n) sin(t ln n)/(√n ln n), and h_δ′(t) = 2√(1 + δ) ln(1 + δ) sin(ut)/t² + O(t⁻³), u = ln(1 + δ). The term of n gives

```
(2Λ(n)/(πδ)) · (√(1 + δ) ln(1 + δ)/(√n ln n)) · ∫_T^∞ [cos((u − ln n)t) − cos((u + ln n)t)]/t² dt.
```

Near a prime power, 1 + δ ≈ n, the middle factor is 1 and the cos((u + ln n)t) part is O(1/T²). So, to leading order,

```
closed − zeros ≈ (2Λ(n)/(πδ)) · I(T),        I(T) = ∫_T^∞ cos(wt)/t² dt = cos(wT)/T − |w|(π/2 − Si(|w|T)),        w = ln(1 + δ) − ln n.
```

At w = 0 this is 2Λ(n)/(πδT). For |w|T ≫ 1, I(T) ≈ −sin(|w|T)/(|w|T²): the peak divided by |w|T. So the resonance changes sign and has width about 1/T in ln(1 + δ). The truncated sums reach the level of the non-resonant rows of §7 (at most 2e-9) only when |ln(1 + δ) − ln n| ≫ 1/T for every prime power n. At δ = 2.999, |w|T = 8.3 and the error is still +4.98e-7. This is Landau's formula Σ_{0<γ≤T} x^ρ = −(T/2π)Λ(x) + O(ln T) seen in the variance. It explains why the repo's sums converge slowly at δ = 1.

For χ, the sum runs over both signs of the ordinates, truncated at T₊ and T₋. The Re χ(n) part of S_χ resonates; the Im χ(n) part gives cos(t ln n) terms, which do not resonate and cancel between the two signs. So

```
closed − zeros ≈ (Re χ(n) Λ(n)/(πδ)) · [I(T₊) + I(T₋)],        = 2 Re χ(n) Λ(n)/(πδT) at w = 0 when T₊ = T₋ = T.
```

The observed values are in §7: at every row with 1 + δ at or next to a prime power, the prediction matches the observed difference to 2e-9 for ζ, and to 1e-6 for characters at T ≈ 1000. χ₄ shows nothing at δ = 1 since χ₄(2) = 0, and neither does the complex character mod 5, since Re χ(2) = Re i = 0. (Heuristic: the series of S does not converge, and the remainder was not bounded.)

---

## 7. Numerical checks (computed)

**ζ: the closed form against 40,000 zeros.** Zero side: (2/δ)[Σ_{n≤N} h_δ(γ_n) + ∫_T^∞ h_δ θ′/π dt − h_δ(T) S(T+)], with T = γ_N (T = 33190.012 for N = 40,000) and S(T+) = N − 1 − θ(T)/π. The oscillatory part of the tail integral is done half-period by half-period with tanh-sinh quadrature, then an 8-term integration-by-parts series. The omitted term is −(2/δ) ∫_T^∞ h_δ′ S dt (C9). Two bounds for it:

- (S₁) |∫_{t₁}^{t₂} S(t) dt| ≤ 2.067 + 0.059 ln t₂ for t₂ > t₁ > 168π (T. Trudgian, Math. Comp. 80 (2011), Thm 2.2), after one more integration by parts, with s = √(1 + δ), u = ln(1 + δ) and |h_δ″| ≤ 2su²/t² + 8su/t³ + 6(2 + δ + 2s)/t⁴;
- (S) |N(t) − F(t)| < 0.137 ln t + 0.443 ln ln t + 1.588, F(t) = (t/2π) ln(t/2π) − t/2π + 7/8 (J. B. Rosser 1941, Thm 19), so |S(t)| is below that plus 1/(40t), with |h_δ′| ≤ 2su/t² + 2(2 + δ + 2s)/t³.

These bounds are rigorous only under three conditions: (a) all zeros above T lie on the line, which is verified up to 3·10¹² (Platt–Trudgian 2021) and formal beyond (fact 2047, not refereed); (b) the zero list below T is complete, for which there is no Turing-method proof here; (c) the constants of Trudgian and Rosser are as quoted, which was not rechecked against the papers.

The zero file: mpmath's `zetazero(n)`, which locates the n-th zero through Gram blocks, returns the file's ordinate at n = 1, 2, 1000, 20000, 40000 (agreement better than 1e-12), so the indexing is right at those points. S(γ_n+) at every 50th zero has mean 0.50 (0.51 and 0.49 on the two halves; a missing zero would shift it by −1 from there on).

| δ | 𝔇(δ) closed | zeros − closed | bound (S₁) | bound (S) | resonance term (C9) | residual |
|---|---|---|---|---|---|---|
| 0.001 | 5.49491150668117 | +1.2e-9 | 3.5e-7 | 5.4e-4 | – | – |
| 0.01 | 3.21247396010240 | −5.1e-10 | 3.3e-6 | 5.1e-4 | – | – |
| 0.1 | 1.10270616884989 | +2.3e-10 | 3.1e-5 | 5.1e-4 | – | – |
| 0.5 | 0.204702422036445 | −1.6e-9 | 1.3e-4 | 5.1e-4 | – | – |
| 1 | 0.187596954133505 | −1.330e-5 | 2.2e-4 | 5.0e-4 | −1.330e-5 (n = 2) | −1.6e-10 |
| 1.5 | 0.0867447359919386 | −8.1e-10 | 2.9e-4 | 4.9e-4 | – | – |
| 2 | 0.128705436305314 | −1.054e-5 | 3.5e-4 | 4.8e-4 | −1.054e-5 (n = 3) | −6.3e-10 |
| 2.999 | 0.0797441945246622 | +4.98e-7 | 4.2e-4 | 4.7e-4 | +4.97e-7 (n = 4) | +1.6e-9 |
| 3 | 0.0798007375292860 | −4.43e-6 | 4.2e-4 | 4.7e-4 | −4.43e-6 (n = 4) | −1.9e-9 |
| 3.001 | 0.0797419206506647 | +4.97e-7 | 4.2e-4 | 4.7e-4 | +4.97e-7 (n = 4) | +6.3e-11 |
| 5 | 0.0515093215827293 | +1.4e-9 | 5.2e-4 | 4.5e-4 | – | – |
| 6 | 0.0795872835034556 | −6.22e-6 | 5.5e-4 | 4.4e-4 | −6.22e-6 (n = 7) | +2.2e-10 |
| 10 | 0.0663954730327453 | −4.60e-6 | 6.3e-4 | 4.1e-4 | −4.60e-6 (n = 11) | +5.8e-10 |

Every difference is inside both bounds. Away from prime powers (δ = 0.001, 0.01, 0.1, 0.5, 1.5, 5) the difference is at most 2e-9. At δ = 1 the truncations at N = 10,000, 20,000 and 40,000 give −4.47e-5, −2.45e-5 and −1.330e-5 (predicted −4.47e-5, −2.45e-5, −1.330e-5): (zeros − closed)·T/(2Λ(n)/(πδ)) = −1.00 at all three, the 1/T law. The ratio is also −1.00 at δ = 2, 3, 6, 10 (n = 3, 4, 7, 11). At δ = 2.999 and 3.001 the difference has the opposite sign to δ = 3 and is much smaller, as C9 predicts.

**ζ: other checks.**

- C, σ₂ and K₂ from the 40,000 zeros plus the smooth tail and the boundary term agree with their closed forms to 1e-12, 1e-12 and 1e-18.
- Dt(δ) from the zeros against [(2 + δ)K₂ − 2J₂(1 + δ)]/δ at δ = 0.01, 0.1, 0.5, 1, 1.5, 3, 10: within 5e-15. The largest, −4.0e-15 at δ = 1, is the same resonance with an extra 1/t² in the weight: −2Λ(2)/(3πδT³) = −4.0e-15.
- The Weil route by quadrature and (4.3): §4.2. The Riesz route: §4.1, symbolic.
- 𝔇′ against a numerical derivative at δ = 0.3 (to 1e-20); the slope jumps at δ = 1, 2, 3, 4 against −2Λ(n)/(n(n − 1)) (to 1e-15).

**Characters, against zeros from PARI** (cypari2 `lfunzeros`, heights 0 to 1000, in blocks of 100). PARI can miss a close pair of zeros at its default subdivision, so the script counts the zeros against the Riemann–von Mangoldt formula block by block: S(γ_n+) = n − θ_χ(γ_n)/π + arg L(½, χ)/π must have mean near ½ in every block, and a block that fails is recomputed with a finer subdivision. For the complex character mod 5, two blocks needed this. Here θ_χ(t) = (t/2) ln(q/π) + Im ln Γ((½ + a + it)/2), and the term arg L(½, χ)/π comes from N(T, χ) = θ_χ(T)/π + S_χ(T) − arg L(½, χ)/π, where arg L(½, χ) is followed continuously along the real axis from +∞. It is 0 for real χ, +0.088 for the complex character mod 5 and −0.073 for the one mod 7 (and the negatives for χ̄, whose zeros are the negative ordinates). This fixes the boundary term S(T+) for complex χ, which the earlier script took as if it were real. The zero side is (1/δ) Σ_{±} [Σ h_δ(γ) + ∫_{T±}^∞ h_δ θ_χ′/π dt − h_δ(T±) S_±(T±+)] over both signs of the ordinates.

The characters: χ₄ (odd, mod 4), χ₃ (odd, mod 3), (·/5) (even, mod 5), the complex odd character mod 5 with χ(2) = i, and the even character mod 7 of order 3 with χ(3) = e^(2πi/3). The script checks the last two tables, their primitivity and parity against PARI.

| χ | δ = 0.01 | 0.1 | 0.5 | 1 | 1.5 | 2.5 | 3.5 |
|---|---|---|---|---|---|---|---|
| χ₄: 𝔇_χ | 4.5838428044 | 2.3460172309 | 0.9814217024 | 0.5273950277 | 0.3142687894 | 0.2494544590 | 0.3102969029 |
| zeros − closed | +1.7e-7 | −1.5e-7 | −1.8e-7 | −3.5e-7 | +4.8e-7 | −6.1e-7 | −2.9e-6 |
| χ₃: 𝔇_χ | 4.2961607320 | 2.0583351585 | 0.6937396300 | 0.2397129553 | 0.2576357772 | 0.2311790564 | 0.1543202212 |
| zeros − closed | −5.8e-7 | +4.9e-7 | +1.1e-6 | +4.42e-4 → +9.7e-7 | −1.1e-6 | +2.0e-6 | +1.7e-6 |
| (·/5): 𝔇_χ | 4.8020112108 | 2.5215237217 | 1.0032101183 | 0.4107405054 | 0.3256993051 | 0.3000987853 | 0.2907246295 |
| zeros − closed | −5.2e-7 | +4.6e-7 | +7.6e-7 | +4.42e-4 → +8.1e-7 | −8.2e-7 | +1.6e-6 | −9.7e-7 |
| mod 5, χ(2) = i: 𝔇_χ | 4.8069863557 | 2.5691607823 | 1.2045652538 | 0.7505385790 | 0.5374123407 | 0.3261163718 | 0.2690617417 |
| zeros − closed | −4.1e-7 | +4.0e-7 | −6.1e-8 | +3.6e-7 | −2.5e-7 | −6.2e-7 | +1.7e-6 |
| mod 7, order 3: 𝔇_χ | 5.1384834475 | 2.8579959583 | 1.3396823549 | 0.7472127421 | 0.5466470116 | 0.3553860485 | 0.2969654583 |
| zeros − closed | −3.8e-7 | +3.3e-7 | +5.4e-7 | +2.21e-4 → +3.3e-7 | −4.1e-7 | +2.0e-7 | +2.0e-6 |

The arrows give the residual after the resonance term of C9 at δ = 1 (Re χ(2) = −1, −1 and −½). Every other entry is within 3e-6, for the complex characters too. The alternatives are far off: adding the pole term 2 ln(1 + δ) would move 𝔇_χ by 0.020 (δ = 0.01) to 3.01 (δ = 3.5), and the wrong parity by 0.0050 to 0.69. The rows at δ = 1.5, 2.5, 3.5 test the χ-weighted prime term. Further checks, for all five characters: the sum over zeros of 1/|ρ|² (with tail and boundary term) agrees with the residue formula J_χ(1) = 2 Re L′/L(1, χ) + ln(q/π) − γ (minus 2 ln 2 for even χ) to 1e-9, a check of the residues of §5 for complex χ; and the two functional-equation relations of §5 hold to 1e-15 with L′/L computed by PARI.

**Primes against zeros, end to end.** R(X, δ) is computed from the primes by exact integration of the step function (sieve to 2·10⁶). At X = 10⁴, 3·10⁵ and 10⁶ it agrees with the zero double sum (1000 zeros and the diagonal tail) to between 1.3e-4 and 2.3e-3, the error of the truncated off-diagonal. Averages over ln X ∈ [ln 10², ln 10⁶] (20,000 points):

| δ | 0.1 | 0.5 | 1 |
|---|---|---|---|
| primes | 1.10865 | 0.20541 | 0.18819 |
| zeros (same window) | 1.11008 | 0.20561 | 0.18835 |
| 𝔇(δ) | 1.10271 | 0.20470 | 0.18760 |
| R range in the window | 0.69 – 1.71 | 0.13 – 0.29 | 0.12 – 0.27 |

The window averages differ from 𝔇 by the window's off-diagonal residue (+0.007, +0.001, +0.001), as Theorem A allows for a finite L.

---

## 8. Prior art

The δ → 0 constant 1 − γ − ln 2π is in Montgomery–Soundararajan (2004) and Chan (2003). The Guinand–Weil device for the variance diagonal is in de la Bretèche–Fiorilli (Math. Ann. 381, 2021, Lemma 2.2), evaluated there to O(δ). For a general L-function, Bui–Keating–Smith (J. London Math. Soc. 94, 2016, Theorem A) give the δ → 0 form with log q_F. The exact closed form for finite δ, its extension to δ > 1, and the equivalence with simple zeros were not found in print in those sources (a search of about 15 papers; textbooks such as Montgomery–Vaughan ch. 13 could not be opened). So these three items are an evaluation we did not find in print.

The leading term ln(1/δ) is Goldston–Montgomery's (1987). Theorem A is the increment version of Cramér's mean-value theorem (Math. Z. 12 (1922) 147–153; the exact form of his statement was not checked) and follows the classical proof pattern. The explicit formulas used (von Mangoldt's, the Riesz-mean form, Weil's) are classical. The tube function and its exact formula are fact 2048 of the repo.

---

## 9. Open points

1. *Without RH.* Theorem A assumes RH. That R(X, δ) is unbounded when a zero lies off the line is plausible (Landau-type), but not proved here.
2. *Simple zeros.* M(δ) = 𝔇(δ) is equivalent to the simplicity of all zeros (Corollary 1), which is open. Under RH alone only M(δ) ≥ 𝔇(δ) is proved.
3. *Finite windows.* Theorem A says nothing about the off-diagonal average over a window of fixed length L. The repo's three-decade averages (r96, r102) are such windows; in §7 the window [10², 10⁶] leaves +0.007 at δ = 0.1.
4. *The resonance (C9)* is heuristic: the series for S(t) does not converge, and the remainder of the leading-order formula was not bounded.
5. *The numerical tail* rests on the three conditions of §7 (RH above T, a complete zero list without a Turing-method proof, the quoted constants). For characters no rigorous remainder bound is given; completeness rests on the block-by-block count against the Riemann–von Mangoldt formula.
6. *Unchecked citations.* The exact form of Cramér's statement; the attribution of σ₂ to Lehmer and Keiper; Barner's conditions for Weil's formula with non-smooth test functions (not needed for the Riesz route); the chapter of Widder's book for Landau's theorem; the constants of Trudgian (Math. Comp. 80 (2011), Thm 2.2) and Rosser (1941, Thm 19) as quoted.
7. *Prior art.* The search covered about 15 papers; Montgomery–Vaughan ch. 13 and similar textbooks were not opened.
