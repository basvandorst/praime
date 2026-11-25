# The one-point sum Σ e^{2iθ_E(γ)} of an elliptic curve as a prime sum: statement, proof, consequences

Folder: `scratchpad/gl2proof/` (nothing in the repo was changed). Code: `gl2lib.py` (library), `exact_check.py` (exact
identities by quadrature), `num_tables.py` (tables, output `num_tables_out.txt`), `check_gl2.py` (independent check of
every identity; run as `check_gl2.py a` (A1–A5, A8–A10, 146 s) and `check_gl2.py b` (A6–A7, 65 s): 19/19 PASS, output
`check_gl2_out.txt`). Sources opened (LaTeX from arXiv) in `dl_*/`.

Labels: **proved** (complete proof here, given the cited lemmas), **proved (cited)** (the step is a precisely cited
lemma that was opened), **sketch**, **GAP**, **conjecture**, **heuristic**, **computed**.

## 0. Result in short

1. **Exact explicit formula (proved, unconditional).** For every analytic window h = (W ∗ Gaussian of width Δ), the sum
   over *all* nontrivial zeros of h(ρ)·X(ρ)^{−1} (X the functional-equation factor, X(½+iγ)^{−1} = e^{2iθ_E(γ)}) equals
   −(1/2π) Σ_n Λ_E(n) n^{−½} ∫ H(t) e^{i(2θ_E(t) − t ln n)} dt up to a left-edge term O(N^{−3/4}T^{−1/2}log NT). No root
   number appears. Checked to 10^{−9} on 5 curves × 3 windows (exact quadrature, no stationary phase).
2. **Smoothed theorem (proved, unconditional; uses Blomer–Khan–Young's stationary-phase lemma).** For dyadic windows
   [T₁,T₂] ⊂ [T,2T]:
   Σ_ρ h(ρ)X(ρ)^{−1} = −(e^{3πi/4}/√2) N^{−1/4} Σ_n Λ_E(n) n^{−1/4} e(−2√(n/N)) H(2π√(n/N)) + O_δ(√N T^{1/2}) for a
   Hann-type window smoothed at a scale Δ ∈ [T^{1/2+δ}, T^{1−δ}], and + O_δ(√N T^{1/2} + √N T^{3/2}Δ^{−1} log NT) for a box
   with Gaussian edges of width Δ. The constant −e^{3πi/4}/√2 of the stationary-phase sketch is **correct**.
3. **Sharp windows.** (i) Under GRH: Σ_{T₁<γ≤T₂} e^{2iθ_E(γ)} = (the same prime sum with Gaussian-edged weight, Δ =
   N^{1/4}T^{3/4}) + O(N^{1/4}T^{3/4} log NT) (proved). (ii) Sharp on both sides, unconditional, Arias de Reyna's method
   gives **O(√N T log²(NT)), not O(T^{1/2} log² T)** (sketch): degree 2 is ADR's κ = 2 case, whose error is T^{κ/2}. This
   is the size of the conjectured main term −w c_E T/2π, so the sharp-sharp form cannot decide the main term. The
   reviewer's error term is **wrong for degree 2**. Numerically the sharp-sharp discrepancy is 0.9–8.9 for T ≤ 400.
4. **Twist** r^{−iγ}: the prime sum runs with r^{1/4} e(−2√(nr/N)) H(2π√(nr/N)) (proved, same error).
5. **Voronoi over all integers (proved, exact identity):** Σ_n λ(n) n^{−½} I_n(r) = w Σ_m λ(m) m^{−½} Ĥ(ln(r/m)). For r = 1
   the m = 1 term gives the resonance (e^{3πi/4}/√2)N^{−1/4} Σ λ(n) n^{−1/4} e(−2√(n/N)) H(t_n) = (w/2π)∫H + O(√N T^{1/2}):
   "c = 1" for λ. Checked: exact identity to 10^{−9}; stationary-phase form to 0.003; sharp partial sums within 0.1–0.3 %.
6. **Consequences.** (a) zero-side law ⟺ prime-side law for smoothed sums (theorem); the sharp prime asymptotic implies
   both, the converse needs a Tauberian input (GAP). (b) Under GRH the smoothed prime sum is ≪ (x³/N)^{1/4} log x +
   (Nx)^{1/2} (theorem): the conjectured main term has the largest size GRH allows, up to log x. Unconditionally only Zhao's
   x^{5/6}log^{21}x (level 1). (d) The heuristic Λ_E = (λ log) ∗ μ_E + Voronoi gives c_E = F′(1) = B_E/L(1, Sym²E) and the
   recipe's line factors ρ(p), ρ(1/p) exactly (checked symbolically); the unjustified step is the type II range.
7. **Novelty:** the sharp theorem: KNOWN IN SUBSTANCE (ADR 2014 Thm 3.1 and Fujii 1992 for ζ with general κ; the degree-2
   transcription not found written). The smoothed theorem with explicit main-term constant: KNOWN IN SUBSTANCE (Weil
   explicit formula + stationary phase), not found written. The Voronoi resonance over integers: KNOWN (Ren–Ye; ILS).
   The constant c_E for primes: NOT FOUND (conjecture).

## 1. Setting and the phase

E/Q of conductor N, root number w, λ(n) = a_n/√n, L(s) := L(s,E) = Σ λ(n)n^{−s} (σ > 1); local factors
(1 − λ_p p^{−s} + p^{−2s})^{−1} (p ∤ N), (1 − λ_p p^{−s})^{−1} (p | N, λ_p ∈ {0, ±p^{−½}}). Λ(s) := N^{s/2}Γ_C(s + ½)L(s),
Γ_C(s) = 2(2π)^{−s}Γ(s), is entire of order 1 with Λ(s) = wΛ(1 − s) (modularity). −L′/L(s) = Σ Λ_E(n)n^{−s},
Λ_E(p^k) = (α_p^k + β_p^k) ln p (α_pβ_p = 1, α_p + β_p = λ_p) for p ∤ N, λ_p^k ln p for p | N; |Λ_E(n)| ≤ 2Λ(n).

**Lemma 1 (phase; proved).** Put X(s) := N^{½−s}Γ_C(3/2 − s)/Γ_C(s + ½) = N^{½−s}(2π)^{2s−1}Γ(3/2 − s)/Γ(s + ½). Then
L(s) = wX(s)L(1 − s), and for real t, X(½ + it)^{−1} = e^{2iθ(t)} with θ(t) := (t/2)ln N − t ln 2π + Im lnΓ(1 + it).
Hence L(½ + it) = w e^{−2iθ(t)} L(½ − it), and Z(t) := e^{iθ(t)}L(½ + it) satisfies Z = w·Z̄ (real for w = 1, imaginary
for w = −1). *Proof.* Divide the two forms of Λ; on s = ½ + it, X = N^{−it}(2π)^{2it}Γ(1 − it)/Γ(1 + it) and
Γ(1 − it) = conj Γ(1 + it). The last claim uses L(½ − it) = conj L(½ + it) (real coefficients). ∎ *Checked (A1):* 1.6·10^{−13}
(X^{−1} vs e^{2iθ}), 1.5·10^{−14} (FE with w via PARI lfun, Z = wZ̄).

**Lemma 2 (continuation; proved).** s ↦ X(s)^{−1} = N^{s−½}(2π)^{1−2s}Γ(s + ½)/Γ(3/2 − s) is meromorphic, with poles
only at s = −½ − k and zeros only at s = 3/2 + k (k ≥ 0); it is holomorphic and zero-free on −½ < σ < 3/2. For
s = ½ + iα it equals exp(2iθ(α)) with θ continued analytically, so **the summand attached to a zero ρ = ½ + iα_ρ,
α_ρ = γ − i(β − ½), is X(ρ)^{−1} = e^{2iθ(α_ρ)}**; under GRH it is e^{2iθ(γ)}. For σ ∈ [−¼, 3/2], t ≥ 1:
−X′/X(σ + it) = ln N − 2 ln 2π + ψ(σ + ½ + it) + ψ(3/2 − σ − it) = ln(Nt²/4π²) + O(1/t), hence
|X(σ + it)^{−1}| = (Nt²/4π²)^{σ−½} e^{O(1/t)} (Stirling). ∎

**Lemma 3 (Stirling for θ; classical, checked).** For t ≥ 1: 2θ(t) = t ln(Nt²/4π²) − 2t + π/2 − 1/(6t) + O(t^{−3});
2θ′(t) = ln N − 2 ln 2π + 2 Re ψ(1 + it) = ln(Nt²/4π²) + 1/(6t²) + O(t^{−4}); 2θ″(t) = −2 Im ψ′(1 + it) = 2/t − 1/(3t³) +
O(t^{−5}); θ^{(j)}(t) ≪_j t^{1−j} (j ≥ 2), with constants independent of N. (From lnΓ(z) = (z − ½)ln z − z + ½ln 2π +
1/(12z) + O(|z|^{−3}) at z = it, and Re ψ(1 + it) = ln t + 1/(12t²) + 1/(120t⁴) + …) *Checked (A2):* t³|2θ − Stirling| ≤ 0.006,
t⁴|2θ′ − …| and t⁵|2θ″ − …| bounded, t ∈ {20, 50, 100, 400}.

**Lemma 4 (zeros and L′/L; proved (cited)).** (i) N_E(t + 1) − N_E(t) ≪ log(N(|t| + 2)). (ii) For each Y ≥ 10 there is
Y′ ∈ [Y, Y + 1] with |L′/L(σ + iY′)| ≪ log²(N Y) for −1 ≤ σ ≤ 2. Source: Milinovich–Ng, *Simple zeros of modular
L-functions*, arXiv:1306.0854 (opened), Lemma 3.1 and Lemma 3.2, which derive them from Iwaniec–Kowalski Thm 5.38 and
Prop. 5.7 (IK not opened). MN state constants "depending on f"; the N-uniform form log(N…) is IK's (doubt 1, §9).

**Lemma 5 (stationary phase; proved (cited)).** Blomer–Khan–Young, *Distribution of mass of holomorphic cusp forms*,
Duke Math. J. 162 (2013), arXiv:1203.2573 (opened): **Lemma 8.1** (repeated integration by parts: if w is smooth, supported
on [α, β], w^{(j)} ≪_j XU^{−j}, |h′| ≥ R, h^{(j)} ≪_j YQ^{−j} (j ≥ 2), then ∫w e^{ih} ≪_A (β − α)X[(QR/√Y)^{−A} + (RU)^{−A}]);
**Proposition 8.2** (if 0 < δ < 1/10, Z := Q + X + Y + V₁ + 1, Y ≥ Z^{3δ}, V₁ ≥ V ≥ QZ^{δ/2}/Y^{½}, w supported on an interval
J of length V₁ with w^{(j)} ≪_j XV^{−j}, h with a unique t₀ ∈ J where h′(t₀) = 0, h″ ≫ YQ^{−2}, h^{(j)} ≪_j YQ^{−j}
(j ≥ 1), then ∫w e^{ih} = e^{ih(t₀)}h″(t₀)^{−½} Σ_{k≤3A/δ} p_k(t₀) + O_{A,δ}(Z^{−A}), p_k = √(2π)e^{πi/4}(k!)^{−1}
(i/(2h″(t₀)))^k G^{(2k)}(t₀), G = w e^{iΦ}, Φ(t) = h(t) − h(t₀) − ½h″(t₀)(t − t₀)²). The same lemmas restated in
Kıral–Petrow–Young, arXiv:1710.00916 (opened), Lemma 3.1.

## 2. The exact explicit formula

Let W be measurable, |W| ≤ 1, supported in [T₁, T₂] (0 < T₁ < T₂), Δ ≥ 1, g_Δ(z) := e^{−z²/Δ²}/(√π Δ), and
H(α) := ∫ W(u) g_Δ(α − u) du (entire in α), h(s) := H((s − ½)/i), so h(½ + it) = H(t). Since |g_Δ(t − iy)| =
e^{y²/Δ²}g_Δ(t): **|h(σ + it)| ≤ e^{(σ−½)²/Δ²}(|W| ∗ g_Δ)(t)** — bounded in vertical strips, Gaussian decay in |t|.

**Proposition 1 (proved, unconditional).** With I_n := ∫_ℝ H(t) e^{i(2θ(t) − t ln n)} dt,

    Σ_ρ m_ρ h(ρ) X(ρ)^{−1} = −(1/2π) Σ_{n≥2} Λ_E(n) n^{−½} I_n + E₀,
    E₀ := −(1/2πi) ∫_{(−¼)} h(s) X(s)^{−1} (L′/L)(s) ds,   |E₀| ≪ N^{−3/4} log(N T₂) (T₁^{−½} + e^{−T₁²/(4Δ²)}),

where ρ runs over **all** nontrivial zeros (multiplicity m_ρ; central zeros and γ < 0 included), and both sides converge
absolutely. With the twist r^{½−s} (r > 0) in the integrand: Σ_ρ m_ρ h(ρ)X(ρ)^{−1} r^{½−ρ} = −(1/2π) Σ_n Λ_E(n)n^{−½}
I_n(r) + E₀(r), I_n(r) := ∫H(t)e^{i(2θ(t) − t ln(nr))}dt, |E₀(r)| ≤ r^{3/4}·(bound above); on the line r^{½−ρ} = r^{−iγ}.

*Proof.* (1) F(s) := h(s)X(s)^{−1}(L′/L)(s) is meromorphic in −¼ ≤ σ ≤ 5/4 with poles only at the nontrivial zeros
(Lemma 2: X^{−1} is holomorphic there; L has no zeros for σ ≥ 1 (absolute convergence for σ > 1; non-vanishing on
σ = 1) nor for σ ≤ 0 by the functional equation; the trivial zeros −½ − k lie outside), residue m_ρ h(ρ)X(ρ)^{−1}.
(2) Integrate over the rectangle [−¼, 5/4] × [−Y′, Y′] with Y′ from Lemma 4(ii) (and −Y′ by conjugation symmetry). On the
horizontal sides |h| ≤ e^{1/Δ²}e^{−(Y′−T₂)²/Δ²}, |X^{−1}| ≪ (NY′²)^{3/4}, |L′/L| ≪ log²(NY′): the sides → 0. The number
of zeros below height Y is ≪ Y log(NY) and |h(ρ)X(ρ)^{−1}| ≪ e^{−(|γ|−T₂)²/Δ²}(N(|γ|+2)²)^{½}, so the zero sum
converges absolutely. (3) Right side σ₁ = 5/4: L′/L = −ΣΛ_E(n)n^{−s}; Σ|Λ_E(n)|n^{−5/4}∫|hX^{−1}|dt < ∞, so (Fubini)
(1/2πi)∫_{(σ₁)}F = −(1/2π)Σ_n Λ_E(n) J_n(σ₁), J_n(σ) := ∫h(σ+it)X^{−1}(σ+it)n^{−σ−it}dt. (4) For each n separately,
h X^{−1} n^{−s} is holomorphic on ½ ≤ σ ≤ σ₁ and tends to 0 uniformly as |t| → ∞ (Gaussian), so J_n(σ₁) = J_n(½) =
n^{−½} I_n (Lemma 1). Since J_n(σ₁) = n^{−½}I_n termwise and Σ|Λ_E(n)||J_n(σ₁)| ≤ Σ|Λ_E(n)|n^{−σ₁}∫|hX^{−1}| < ∞, the
series Σ Λ_E(n)n^{−½}I_n converges absolutely. (5) Left side = E₀. On σ = −¼: L′/L(s) = X′/X(s) − (L′/L)(1 − s),
|(L′/L)(5/4 − it)| ≤ 2|ζ′/ζ(5/4)|, |X′/X| ≤ log(N(|t|+2)²) + O(1), |X^{−1}(−¼ + it)| ≪ N^{−3/4}(1 + |t|)^{−3/2} (Lemma 2),
|h| ≤ e^{9/(16Δ²)}(1_{[T₁,T₂]} ∗ g_Δ)(t); integrate. The twisted case is identical (|r^{½−s}| = r^{½−σ}). ∎

*Computed:* (A4) and `exact_check.py`/`num_tables.py`: |zeros − primes| = 4·10^{−11} … 8·10^{−9} on 5 curves × 3 boxes
(Δ = 10, 20) and 2 more (Δ = 4, 5), I_n by quadrature over all n with |ln n − 2θ′(t)| ≤ 1.6 on the support. (This is
also a completeness check of the zero lists in those windows: one missing zero would give an error ≈ 1.)

## 3. Stationary phase: the smoothed theorem

Write φ_n(t) := 2θ(t) − t ln n, t_n the solution of 2θ′(t_n) = ln n (unique for t ≥ 1, since 2θ″ > 0 there; Lemma 3),
t_n* := 2π√(n/N) (= t_n + O(1/t_n)).

**Theorem 2 (proved given Lemmas 4–5; unconditional).** Fix 0 < δ < 1/10. Let T ≥ T₀(δ), [T₁, T₂] ⊂ [T, 2T], H = W ∗ g_Δ
as in §2 with one of
- (box) W = 1_{[T₁,T₂]}, T^{½+δ} ≤ Δ ≤ T^{1−δ};
- (Hann type) W(t) = w₀((t − T₁)/(T₂ − T₁)), T₂ − T₁ ≥ T/4, w₀ supported in [0,1], C¹ with ‖w₀″‖_∞ < ∞ (e.g. sin²(πu)),
  T^{½+δ} ≤ Δ ≤ T^{1−δ} (for w₀ ∈ C^∞_c(0,1) any 1 ≤ Δ ≤ T^{1−δ}).

Then

    Σ_ρ m_ρ h(ρ) X(ρ)^{−1} = −(e^{3πi/4}/√2) N^{−1/4} Σ_n Λ_E(n) n^{−1/4} e(−2√(n/N)) H(2π√(n/N)) + R,
    R ≪_{δ,w₀} √N T^{½}                                  (Hann type),
    R ≪_δ √N T^{½} + √N T^{3/2} Δ^{−1} log(NT)            (box),

uniformly in N. The same holds with the "exact" main term −(1/2π)Σ_n Λ_E(n)n^{−½}(2π/φ_n″(t_n))^{½}e^{i(φ_n(t_n)+π/4)}H(t_n)
(the difference is part of R). With the twist r^{½−ρ}: the main term becomes −(e^{3πi/4}/√2)N^{−1/4} r^{1/4}
Σ_n Λ_E(n)n^{−1/4}e(−2√(nr/N))H(2π√(nr/N)), same R (constants depend on r).

*Proof.* Start from Proposition 1 (E₀ is ≪ N^{−3/4}). Put D := 3Δ√(log T), x₊ := N(8T)²/4π², and ω ∈ C^∞_c with ω = 1 on
[T₁ − D, T₂ + D], supp ω ⊂ [T₁ − 2D, T₂ + 2D] ⊂ [3T/4, 5T/2], ω^{(j)} ≪ D^{−j}. Then |H(1 − ω)| ≤ T^{−9} e^{−(dist−D)²/Δ²}.

*(a) Tails of H.* Σ_{n≤x₊}|Λ_E(n)|n^{−½}∫|H(1 − ω)| ≪ √x₊·ΔT^{−9} ≪ √N T^{−7}.

*(b) n ≤ x₊ with t_n ∉ [T/2, 4T].* On supp(Hω), |φ_n′(t)| = |2θ′(t) − ln n| ≥ 2 ln(3/2) − O(T^{−2}); φ_n^{(j)} ≪ T^{1−j}.
Lemma 8.1 (X = 1, U = Δ, Q = Y = T, R ≍ 1): I_n(Hω) ≪_A T(T^{−A/2} + Δ^{−A}); summed with weights |Λ_E(n)|n^{−½} over
n ≤ x₊: ≪ √N T^{2−A/2}, negligible.

*(c) n with t_n ∈ [T/2, 4T].* Apply Prop. 8.2 on J = [T/4, 8T] (V₁ ≍ T) with h = φ_n (h″ = 2θ″ ≍ 1/T, h^{(j)} ≪ T^{1−j},
|h′| = |ln(t²/t_n²)| + O(T^{−2}) ≤ ln 256 on J, so Y = Q = T), w = Hω, X = 1, V = Δ (box, Hann type: H^{(j)} ≪_j Δ^{−j}), Z ≍ T; the
conditions hold for T ≥ T₀(δ) because Δ ≥ T^{½+δ}. Main term p₀: √(2π)e^{πi/4}w(t_n)·e^{iφ_n(t_n)}φ_n″(t_n)^{−½}. For k ≥ 1,
G^{(2k)}(t_n) = Σ_j C(2k, j) w^{(j)}(t_n)(e^{iΦ})^{(2k−j)}(t_n), and since Φ vanishes to third order at t_n, Faà di Bruno gives
(e^{iΦ})^{(m)}(t_n) = 0 for m = 1, 2 and ≪_m T^{⌊m/3⌋−m} for m ≥ 3. With |1/h″| ≍ T:

    |p_k(t_n)| ≪_k T^k|w^{(2k)}(t_n)| + Σ_{1≤j≤2k−3} T^{(2j−k)/3}|w^{(j)}(t_n)| + T^{⌊2k/3⌋−k}|w(t_n)|   (k ≥ 1),

and p₁ = √(2π)e^{πi/4}(i/(2φ_n″))w″(t_n) exactly. Multiply by (1/2π)|Λ_E(n)|n^{−½}φ_n″^{−½} ≍ |Λ_E(n)|n^{−1/4}N^{−1/4} and
sum over n. (i) Terms with w itself (only k ≥ 2): T^{⌊2k/3⌋−k} ≤ T^{−1}, and N^{−1/4}Σ_{n≤x(4T)}2Λ(n)n^{−1/4} ≪ N^{½}T^{3/2}
(Chebyshev ψ(x) ≪ x): total ≪ √N T^{½}. (ii) Terms with derivatives: they are supported where w^{(j)} ≠ 0. Summing over
unit t-intervals [k, k+1] (n-length ≍ NT, Σ_{n in it}|Λ_E(n)| ≤ 2(NT + 1)log x, n^{−1/4}N^{−1/4} ≍ N^{−½}T^{−½}),
Σ_n |Λ_E(n)|n^{−1/4}N^{−1/4}|f(t_n)| ≪ √N T^{½} log(NT)·(‖f‖_{L¹} + ‖f′‖_{L¹}) (and ≪ √N T^{½}·T‖f‖_∞ without log, by Chebyshev). Box: ‖H^{(j)}‖_{L¹} ≍
Δ^{1−j}; the largest term is k = 1: T‖H″‖₁ ≍ T/Δ; the others are smaller (T^kΔ^{1−2k} = (T/Δ)(T/Δ²)^{k−1}; T^{(2j−k)/3}
Δ^{1−j} ≤ T^{(2−k)/3} for k ≥ 2). Hann type: H^{(j)} = W″ ∗ g_Δ^{(j−2)} (j ≥ 2), ‖H^{(j)}‖_∞ ≪ T^{−2}Δ^{2−j}, ‖H′‖_∞ ≪ 1/T;
every term is ≪ T^{−1} pointwise on [T/2, 4T] (k = 1: T·T^{−2}; T^k T^{−2}Δ^{2−2k} ≤ T^{−1}; mixed T^{(2j−k)/3}T^{−2}Δ^{2−j}
≤ T^{−1}), so the total is ≪ √N T^{½} with Chebyshev only. Prop 8.2's remainder O(Z^{−A}) per n sums to ≪ √N T^{2−A}.
(iii) If t_n ∉ supp w, then p_k(t_n) = 0 for all k: only the O(Z^{−A}) remainder is left.

*(d) n > x₊.* Here use J_n(σ₁) = n^{−½}I_n (Proposition 1, step 4) with σ₁ := 1 + 1/log x₊. Split h = hχ₁ + hχ₂ with
χ₁ = 1 on [3T/4, 3T], supp χ₁ ⊂ [T/2, 4T], χ₁^{(j)} ≪ T^{−j}. The χ₂-part: |h(σ₁ + it)χ₂| ≤ 2e^{−T²/(16Δ²)}·(Gaussian in t),
so its contribution is ≪ √N T² e^{−T^{2δ}/16}Σ|Λ_E(n)|n^{−σ₁}, negligible. On supp χ₁ write X^{−1}(σ₁+it)n^{−σ₁−it} = e^{Ψ(t)},
Ψ′(t) = i(ℓ(σ₁ + it) − ln n), ℓ := −X′/X = ln(Nt²/4π²) + O(1/t) (Lemma 2), so |Ψ′| ≥ ln(n/x₊) + ln 4 − o(1) ≥ 1 and
Ψ^{(j)} ≪ T^{1−j}. Integrating by parts k times (the operator 𝒟f = −(f/Ψ′)′ of BKY's proof of Lemma 8.1; their bound uses
only |Ψ′| ≥ R and Ψ^{(j)} ≪ YQ^{−j}, also valid for complex Ψ) with weight derivatives ≪ Δ^{−j}:
|J_n(σ₁)| ≪_k T·(NT²)^{σ₁−½}n^{−σ₁}·R^{−k}(Δ^{−1} + T^{−½})^k, R = ln(4n/x₊). Since Σ_{n>x₊}Λ(n)n^{−σ₁}(ln(4n/x₊))^{−k} ≪ 1
(k ≥ 2), the total is ≪ √N T^{2−k/2}; k = 5 gives √N T^{−½}.

*(e) Main term.* From (c): −(1/2π)Λ_E(n)n^{−½}(2π/φ_n″(t_n))^{½}e^{i(φ_n(t_n)+π/4)}H(t_n). By Lemma 3: 2θ′(t_n) = ln n
⇒ t_n = t_n*(1 + O(t^{−2})); φ_n(t_n) = t_n[ln(Nt_n²/4π²) − ln n] − 2t_n + π/2 + O(1/t) = −4π√(n/N) + π/2 + O(1/t);
(2π/φ″)^{½} = (πt_n)^{½}(1 + O(t^{−2})) = π√2 (n/N)^{1/4}(1 + O(t^{−2})); H(t_n) = H(t_n*) + O(‖H′‖_∞/T). Hence the
term is −(1/2π)·π√2·e^{iπ/2}e^{iπ/4}·Λ_E(n)n^{−½}(n/N)^{1/4}e(−2√(n/N))H(t_n*)(1 + O(1/T)) =
**−(e^{3πi/4}/√2)N^{−1/4}Λ_E(n)n^{−1/4}e(−2√(n/N))H(t_n*)**·(1 + O(1/T)); the O(1/T) summed as in (c)(i) is ≪ √N T^{½}.
The twisted case: replace ln n by ln(nr) throughout; t_{n,r}* = 2π√(nr/N) and (2π/φ″)^{½} = π√2(nr/N)^{1/4}. ∎

*Remarks.* (1) The zeros with γ ≤ 0 contribute ≪ √N e^{−T²/Δ²}: under GRH the left side is Σ_{γ>0}H(γ)e^{2iθ(γ)} up to that.
(2) Second-order main term: adding p₁ and the k = 2, 3 terms with w itself, (i/(2φ″))(w″ − w′φ‴/φ″ − wφ⁗/(4φ″) +
5wφ‴²/(12φ″²)) (for constant w this is −i/(48t)·w), lowers the error (computed: §7). (3) Why Δ ≥ T^{½+δ}: the Fresnel
scale of each prime term is (πt)^{½}; a weight that changes faster cannot be evaluated term by term. For T ≤ 400 this
scale is ≈ 35, which is why short Hann windows (length 120) carry errors ≈ 1 at leading order (§7).
(4) **The bound is far from the data.** For N ≈ 300–900, T ≈ 400, √N T^{½} ≈ 350–600 exceeds the main term (6–100). The
measured stationary-phase error is 0.03–0.24 (Hann [20,400]): the trivial estimate of the correction terms is very lossy,
because each correction is itself a resonant prime sum with a weight whose ∫ vanishes (e.g. ∫(t/4)H″ = 0 for Hann).

## 4. Sharp windows

**Theorem 3(i) (proved, assumes GRH for L(s,E)).** For T₁, T₂ in [T, 2T], T ≥ N^{1/(1−4δ)}T₀(δ):

    Σ_{T₁<γ≤T₂} e^{2iθ(γ)} = −(e^{3πi/4}/√2)N^{−1/4}Σ_n Λ_E(n)n^{−1/4}e(−2√(n/N))H_Δ(2π√(n/N)) + O_δ(N^{1/4}T^{3/4}log(NT)),

with H_Δ = 1_{[T₁,T₂]} ∗ g_Δ, Δ = N^{1/4}T^{3/4}. *Proof.* Under GRH h(ρ) = H_Δ(γ) and |X(ρ)^{−1}| = 1, so the difference of the
sharp and the smoothed zero sums is ≤ Σ_γ |1_{(T₁,T₂]} − H_Δ|(γ) ≤ ½Σ_γ Σ_{e∈{1,2}} erfc(|γ − T_e|/Δ) ≪ (Δ + 1)log(NT)
(Lemma 4(i)). Theorem 2 (box) bounds the rest by √N T^{½} + √N T^{3/2}Δ^{−1}log(NT); Δ = N^{1/4}T^{3/4} balances. ∎
The prime side stays smoothed: the sharp prime sum differs from it by a sum over a t-window of width Δ around T₁, T₂ whose
trivial bound ≍ √N T^{½}Δ is ≫ T.

**Proposition 3(ii) (sketch; unconditional; transcription of Arias de Reyna 2014).**

    Σ_{0<γ≤T} X(ρ)^{−1} = −(e^{3πi/4}/√2)N^{−1/4}Σ_{n≤NT²/4π²}Λ_E(n)n^{−1/4}e(−2√(n/N)) + O(√N T log²(NT)).

ADR (arXiv:1403.7126, opened) prove, for ζ and every κ > 0 (Theorem 3.1), Σ_{0<Re α≤T}e^{2iκϑ(α)} = −(e^{πi(1−κ)/4}/√κ)
Σ_{n<(T/2π)^κ}Λ(n)n^{−½+1/(2κ)}e^{−2πiκn^{1/κ}} + O_κ(T^{(1−κ)/2}log T) + O_κ(T^{κ/2}log²T). Every error in their §4 comes from
|χ(σ + it)^{−κ}| ≍ (t/2π)^{κ(σ−½)} at σ₁ = 1 + 1/log T (top edge: Lemma 4.2; first-derivative-test parts V_{2,1}, V_{2,2};
Stirling correction; the endpoint remainders of their Lemma 2.4 summed in Lemma 4.5). For L(s,E), |X(σ + it)^{−1}| ≍
(√N t/2π)^{2(σ−½)}: the κ = 2 case with t/2π replaced by √N t/2π, the phase 2θ_E(t) − t ln n = 2π·2x ln(x/(e√(n/N))) + π/2
(x = t/2π, up to O(1/t)) is ADR's κ x ln(x/(er)) with κ = 2, r = √(n/N); their Lemma 2.4 main term κ^{−½}r^{a+½}e^{πi/4}
e^{−2πiκr} gives exactly the constant above. So every bound becomes √N T·log(NT) or √N T log²(NT). *Not written line by line*
(GAP 1); ADR's Lemma 2.4 is itself quoted there "inferred from Levinson and Gonek" without proof.

**The reviewer's statement** ("Theorem 3.1: … + O(T^{½}log²T)") is the κ = 1 corollary for ζ. For degree 2 the method gives
O(√N T log²(NT)), the same order as the conjectured main term −w c_E T/2π: the sharp-sharp version proves nothing about the
main term. A sharp-sharp error o(T) is **not proved** (GAP 2). Computed (§7): |sharp zeros − sharp primes| = 0.9–8.9 on 5
curves × 6 windows, T ≤ 400 (≤ 0.45 √T2), against main terms 4–200.

## 5. Voronoi over all integers

**Theorem 4 (proved, unconditional; exact identity).** For H as in §2 and r > 0, with Ĥ(ξ) := ∫H(t)e^{−iξt}dt,

    Σ_{n≥1} λ(n) n^{−½} ∫ H(t) e^{i(2θ(t) − t ln(nr))} dt = w Σ_{m≥1} λ(m) m^{−½} Ĥ(ln(r/m)),

both sides absolutely convergent. *Proof.* Consider (1/2πi)∫_{(5/4)} h(s)X(s)^{−1}L(s)r^{½−s}ds. Expanding L and shifting each
n-term to σ = ½ (as in Prop. 1, step 4) gives (1/2π)·(left side). On the other hand X^{−1}L = wL(1 − s) is entire, and
L(1 − s) is polynomially bounded in vertical strips, so the integral may be moved to σ = −¼, where L(1 − s) = Σλ(m)m^{s−1}
converges absolutely; shifting each m-term back to σ = ½ gives (w/2π)Σλ(m)m^{−½}∫H(t)(m/r)^{it}dt. ∎
This is the GL(2) Voronoi formula (equivalently the functional equation) for this test function (known in substance).

**Corollary 4.1 (resonance at dual index 1).** With the hypotheses of Theorem 2 and r = 1:
(e^{3πi/4}/√2)N^{−1/4}Σ_n λ(n)n^{−1/4}e(−2√(n/N))H(2π√(n/N)) = (w/2π)[∫H + Σ_{m≥2}λ(m)m^{−½}Ĥ(−ln m)] + O_E(√N T^{½}(…)),
where Ĥ(−ln m) = Ŵ(−ln m)e^{−Δ²ln²m/4} is negligible (box) or ≪ T^{−A} (Hann type). For integer r = q: main term
(w/2π)λ(q)q^{−½}∫H (dual index q); for r with no integer m near r: no main term. (Proof: Theorem 4 + the stationary-phase
analysis of §3 with λ in place of Λ_E, using Σ_{n≤x}|λ(n)| ≪_E x (Rankin–Selberg and Cauchy–Schwarz); the constant now
depends on E.) In the sharp form (heuristic constant, via partial summation of the density (w/√2)e^{−3πi/4}(nN)^{−1/4}):
Σ_{n≤x}λ(n)e(−2√(n/N)) ≈ (2√2/3)e^{−3πi/4}w(x³/N)^{1/4}. Literature: Ren–Ye, arXiv:1411.6022 (opened), Theorem 1.2 (full
level, GL(m), smooth weights: main term at α = m q^{1/m}, and "similar estimate for a sharp-cut sum"); level-N weight-2
version not seen there.

*Computed:* Theorem 4 exactly (quadrature) to 1.1·10^{−9} (A5: 8.3·10^{−10}); the stationary-phase left side against the
exact dual side: |difference| ≤ 0.003 (Hann windows, 5 curves; the m = 1 term (w/2π)·190 = ±30.239); sharp partial sums at
x ≈ 4100N: ratio 0.998–1.003, |S − main|/√x = 0.01–0.03 (§7); A10 (independent code, x = 4000N): |S/main − 1| = 0.0024, 0.0050.

## 6. Consequences

### (a) Zeros ⟺ primes for the main term (theorem, with the stated limits)

Put Z(H) := Σ_ρ h(ρ)X(ρ)^{−1}, P(H) := Σ_n Λ_E(n)n^{−1/4}e(−2√(n/N))H(2π√(n/N)). Theorem 2: Z(H) = −(e^{3πi/4}/√2)N^{−1/4}P(H)
+ O(E_H), E_H = √N T^{½} (Hann type). Hence, for any constant κ and any class of Hann-type weights at scale ≍ T,

    Z(H) = −(κ/2π)∫H + o(T)   ⟺   P(H) = (√2/2π)e^{−3πi/4} κ N^{1/4} ∫H + o(N^{1/4}T)        (T → ∞, N fixed),

because √N T^{½} = o(T). The project's law is κ = w c_E. **Sharp prime asymptotic ⇒ both.** If S(x) := Σ_{n≤x}Λ_E(n)e(−2√(n/N))
= A(x³/N)^{1/4} + o(x^{3/4}), then by partial summation P(H) = (3A/4)N^{−1/4}∫n^{−½}H(t_n*)dn + o(x^{½}) and
∫n^{−½}H(t_n*)dn = (√N/π)∫H dt, so P(H) = (3A/4π)N^{1/4}∫H + o(√N T): with A = (2√2/3)e^{−3πi/4}κ this is the P-law
((3/4π)(2√2/3) = √2/2π; checked symbolically and by quadrature, A9). Under GRH it also gives the sharp zero law with error
o(T) + O(N^{1/4}T^{3/4}log NT) (Theorem 3(i)). **Converse: GAP 3.** The smoothed laws (smoothing ≥ T^{½+δ} in t, i.e. ≥
x^{3/4+δ/2}N^{−1/4} in n) fix S(x) only up to the mass of a window of that width, trivially ≍ N^{1/4}x^{3/4+δ/2}, larger than
the main term; a Tauberian input (a bound for short sums of Λ_E(n)e(−2√(n/N))) is missing.

### (b) What is provable about the prime sum

- Over all integers: Theorem 4 / Corollary 4.1 (exact main term from the dual index 1).
- Over prime powers, unconditionally: trivially ≪ x. Zhao, *Oscillations of Hecke eigenvalues at shifted primes*, Rev. Mat.
  Iberoam. 22 (2006), arXiv:math/0508129 (opened), Theorem 2: Σ_{n≤x}λ(n)Λ(n)e(α√n) ≪_{α,f} x^{5/6}log^{21}(3x) for level 1,
  weight ≥ 12; his §6 says similar results hold for congruence subgroups "via analogous means" (not proved there), and that with
  Lindelöf or Montgomery's conjecture the bilinear forms give the majorant x^{3/4+ε}. (Λ_E(n) − λ(n)Λ(n) is supported on
  p^k, k ≥ 2, total ≪ x^{½}log x.) Zhao quotes the Iwaniec–Luo–Sarnak heuristic S(x) = Z x^{3/4} + O(x^{5/8+ε}), Z ≠ 0.
- **Under GRH (theorem, from Theorem 2):** |Z(H)| ≤ Σ_γ|H(γ)| ≪ (∫|H| + Δ)log(NT), hence for V smooth on [1, 2]
  (H(t) := (Nt²/4π²)^{1/4}V(Nt²/(4π²x)), sup H ≍ x^{1/4}):

      Σ_n Λ_E(n) e(−2√(n/N)) V(n/x) ≪_V x^{3/4} N^{−1/4} log x + (Nx)^{½}.

  (Apply Theorem 2, Hann type with w₀ ∈ C^∞_c, to W := H and Δ = 1; replacing (W ∗ g₁)(t_n) by W(t_n) costs ≪ x^{1/4}T^{−2}
  per term, ≪ N in total.) The conjectured main term (2√2/3)c_E(x³/N)^{1/4}·(…) is of exactly this order: GRH allows the resonance, and caps it at the
  conjectured size up to log x. (Known in substance from Zhao's Lindelöf remark, which gives x^{3/4+ε}.)
- **Type I/II (where it fails).** Λ_E = (λ log) ∗ μ_E (from −L′/L = (−L′)·(1/L)). A Vaughan/Heath-Brown decomposition gives
  type I sums Σ_{d≤D}μ_E(d)Σ_m λ(m)log m·f(dm), f(n) = e(−2√(n/N))V(n/x): each inner sum is evaluated exactly by Theorem 4
  with r = d (main term ∝ wλ(d)d^{−½}: dual index d). Type II sums Σ_{d∼D}Σ_{m∼M}a_d b_m e(−2√(dm/N)) with D, M in
  [x^ε, x^{1−ε}] have no Voronoi structure in either variable; Cauchy–Schwarz/large sieve give at best ≈ x^{5/6} (Zhao) and
  x^{3/4+ε} under Lindelöf. An asymptotic needs type II ≪ x^{3/4−ε} **and** no main term from type II: open, also under GRH.

### (c) A provable mean statement

The GRH bound above is the provable "main term cannot be larger than" statement. Unconditionally, Corollary 4.1 is the
provable average form: the main term exists for the full coefficient sequence λ(n) (constant 1), and at every integer
dual index q (constant λ(q)q^{−½}). A family average (root-number-weighted Petersson average over
weight-2 newforms of prime level N, where X is common to the family, the c = N Kloosterman term S(N, n; N) = −1 resonates with
e(−2√(n/N)) through J₁(4π√(n/N)), and Petersson's harmonic average of c_f is 1 + o(1)) is provable in principle with PNT and
Weil's bound; **not carried out (GAP 4).**

### (d) The heuristic for c_E = B_E/L(1, Sym²E), as a chain

1. Λ_E = (λ·log) ∗ μ_E, μ_E the coefficients of 1/L(s,E) (exact).
2. P(H) = Σ_d μ_E(d)d^{−1/4}Σ_m λ(m)ln m·m^{−1/4}e(−2√(dm/N))H(2π√(dm/N)) (exact rearrangement; absolutely convergent).
3. For each d, Theorem 4 with r = d: the inner sum without ln m equals (√2/2π)e^{−3πi/4}wN^{1/4}λ(d)d^{−3/4}∫H + (dual m ≠ d
   terms) + (stationary-phase error); with ln m, the stationary point gives ln m = 2θ′(t) − ln d, so the main term is
   (√2/2π)e^{−3πi/4}wN^{1/4}λ(d)d^{−3/4}∫H(t)(2θ′(t) − ln d)dt (theorem for each fixed d, as T → ∞).
4. **Sum over d of step 3 (the unjustified step).** The errors of step 3 are uniform only while the inner sum is long
   (d ≤ x^{½−ε}, say); for larger d the m-sums are short, the stationary-phase/Voronoi errors exceed their main terms, and the
   sum over d of the errors is not controlled. This is exactly the type II range of (b).
5. Granting step 4: Σ_d μ_E(d)λ(d)d^{−1}(2θ′(t) − ln d) = 2θ′(t)F(1) + F′(1), F(s) := Σ μ_E(d)λ(d)d^{−s} = B(s)/(ζ(s)L(s,Sym²E))
   (local factors F_p = (1 − X)(1 − λ(p²)X), p ∤ N; 1 − X/p, p | N; checked A8). F(1) = 0 and F′(1) = B_E/L(1,Sym²E) = c_E
   (the series at s = 1 converge by the prime number theorem for ζ·L(s,Sym²E) — theorem-level).
6. Hence P(H) "=" (√2/2π)e^{−3πi/4}w c_E N^{1/4}∫H and, by Theorem 2, Z(H) "=" −(w c_E/2π)∫H: the recipe's law.

**Lines.** The same chain with r = a/b ((a,b) = 1): only d = bd′ give an integer dual index ad′, so the line amplitude is
−(w/2π)(ab)^{−½}F′_{a,b}(1)∫H, F_{a,b}(s) := Σ_{d′}μ_E(bd′)λ(ad′)d′^{−s}. Its local factors are (p ∤ N, X = p^{−1}):
r = p: Σ_j μ(p^j)λ(p^{j+1})X^j = (1 − X)(λ_p − λ(p³)X); r = 1/p: Σ_j μ(p^{j+1})λ(p^j)X^j = −λ_p(1 − X). Divided by F_p(1/p)
these are **exactly the recipe's ρ(p) and ρ(1/p)** of `derivation.md` §2.1 (sympy, A8). So the Voronoi heuristic and the
ratios recipe give the same constants for c_E and for the lines; both rest on step 4.

## 7. Numerical table (5 fresh curves; zeros of round 144, a_p by PARI; `num_tables_out.txt`, `check_gl2_out.txt`)

Columns: Z = zero side, P = leading-order prime side (exact θ at the exact stationary point), P₂ = with the second-order term;
law = −w c_E∫H/2π (c_E frozen; the direct partial product to 2.5–7·10⁶ agrees to 2·10^{−4}).

Hann [20, 400] (the data's taper; ∫H = 190):

| curve | w | c_E | Z | P | \|Z−P\| | \|Z−P₂\| | law |
|---|---|---|---|---|---|---|---|
| N307 | +1 | 1.299 | −50.240 + 11.924i | −50.268 + 11.764i | 0.163 | 0.043 | −39.29 |
| N335 | −1 | 3.270 | 96.091 − 2.048i | 96.118 − 2.100i | 0.058 | 0.032 | 98.87 |
| N373 | −1 | 0.470 | 10.657 + 2.325i | 10.636 + 2.565i | 0.241 | 0.054 | 14.22 |
| N611 | +1 | 0.206 | −17.408 − 6.395i | −17.279 − 6.401i | 0.129 | 0.028 | −6.24 |
| N870 | +1 | 2.637 | −87.654 + 14.307i | −87.650 + 14.416i | 0.109 | 0.051 | −79.73 |

Six windows [20,400], [20,140], [140,260], [260,400], [60,230], [230,380], per curve (max over the windows):

| curve | Hann \|Z−P\| | Hann \|Z−P₂\| | **sharp** \|Z−P\| (range) | box Δ=10,10,20: exact \|Z−P\| | box: stationary \|Z−P\| / with 2nd order | twists r=2,½,3,3/2 \|Z−P\| |
|---|---|---|---|---|---|---|
| N307 | 1.22 | 0.33 | 1.2 – 7.7 | ≤ 4.0e−9 | 0.67, 0.72, 0.74 / 0.47, 3.68, 0.11 | ≤ 0.21 |
| N335 | 1.20 | 0.30 | 0.9 – 4.4 | ≤ 7.9e−9 | 0.07, 1.82, 0.93 / 0.38, 1.93, 0.17 | ≤ 0.19 |
| N373 | 1.12 | 0.52 | 4.0 – 8.9 | ≤ 1.5e−9 | 1.18, 0.77, 0.65 / 0.56, 0.46, 0.04 | ≤ 0.20 |
| N611 | 1.50 | 0.34 | 1.1 – 7.4 | ≤ 2.7e−9 | 1.03, 2.79, 0.89 / 1.10, 2.96, 0.17 | ≤ 0.20 |
| N870 | 0.96 | 0.37 | 1.9 – 7.7 | ≤ 6.8e−9 | 0.57, 1.25, 0.32 / 0.84, 1.52, 0.28 | ≤ 0.16 |

Boxes: [70,180] Δ = 10, [180,330] Δ = 10, [125,275] Δ = 20. "Exact" = Proposition 1 with I_n by quadrature (no
approximation): agreement to 10^{−9}. The stationary-phase errors behave as the theory says: small where the weight varies
slowly against (πt)^{½} ≈ 25–35 (Hann on [20,400]: 0.03–0.05 at second order), ≈ 1–3 for Δ = 10 at t ≈ 300, where the
second-order term even makes it worse (expansion parameter ≈ t/(4Δ²) ≈ 0.75, times O(1)), and 0.04–0.28 for Δ = 20 at second order.

Voronoi (all 5 curves): stationary-phase left side of Theorem 4 vs the exact dual side, Hann [20,400]: 0.003
(second order ≤ 0.004); [20,140], [140,260]: ≤ 0.003. Sharp Σ_{n≤x}λ(n)e(−2√(n/N)) at x = N·402²/4π² (1.26–3.56·10⁶): ratio
to (2√2/3)e^{−3πi/4}w(x³/N)^{1/4} = 1.002, 1.002, 1.003, 0.998, 0.999; |S − main|/√x = 0.01–0.03. The sharp prime sums
S(x) at the same x are far from c_E·main for small c_E (N373: −2913 + 616i vs 3098 + 3098i), as expected from §6(a).

## 8. Novelty (opened sources marked; mathematical queries only)

- Arias de Reyna, arXiv:1403.7126 (2014; **opened**, LaTeX): Theorem 3.1 for ζ and all κ > 0, unconditional, error
  O(T^{(1−κ)/2}log T) + O(T^{κ/2}log²T); corollary at κ = 1 (−ψ(T/2π) + O(T^{½}log²T)); "equivalent to one of the results of
  Schoißengeier" (1979, Monatsh. Math. 87; **not opened**, paywalled). Journal version not found.
- Fujii, *Some observations concerning the distribution of the zeros of the zeta functions III*, Proc. Japan Acad. 68 (1992)
  (**not opened**: Project Euclid and J-STAGE refused automated access): quoted by ADR as the RH-conditional formula
  Σ_{0<γ≤T}e^{iκγ log(κγ/2πea)} = −e^{πi/4}(√a/κ)Σ_{n≤(κT/2πa)^κ}Λ(n)n^{−½+1/(2κ)}e^{−2πian^{1/κ}} + O(T^{κ/2}log T loglog T) +
  O(T^{½−κ/2}log T). For κ = 2 this is the ζ-analogue of our kernel (n^{−1/4}e(−a√n), but with Λ of GL(1): no resonance).
- Hardy–Littlewood 1920 (via ADR): Σ e^{aiγ log(γθ)} = O(T^{½+a/2}) under RH.
- Blomer–Khan–Young (**opened**), Kıral–Petrow–Young (**opened**): the stationary-phase lemmas used.
- Milinovich–Ng, arXiv:1306.0854 (**opened**): Lemma 3.3 is the degree-2 Landau–Gonek formula Σ_{0<γ_f≤T}x^{ρ_f} =
  −(T/2π)Λ_f(x) + …, a linear twist; their approximate functional equation (Lemma 3.4) contains ε_f ψ_f(s)·(dual sum), the
  same X-factor, but no sum Σ X(ρ)^{−1} over zeros.
- Conrey–Ghosh–Gonek (search results only): the χ(1 − s) r^{−s} lemma (Gonek 1984 Lemma 5) gives Σ a(n)e(−n/r): degree 1,
  linear twists. Garunkštis–Steuding, arXiv:1011.5339 (**opened**, LaTeX): sums of ζ′(ρ_a) and ζ(ρ_a + iδ) over a-points, no
  Σ Δ(ρ)^{−1}. Sourmelidis–Steuding–Suriajaya (search snippet): Σ Δ(ρ_a) over a-points of ζ. Kalpokas–Steuding, Murty–Perelli,
  Murty–Zaharescu: **not opened** (no free text found in the time available).
- Zhao 2006 (**opened**) and Ren–Ye arXiv:1411.6022 (**opened**): resonance over integers (theorem, level 1) and over primes
  (bound x^{5/6}; ILS heuristic constant Z unspecified).
- Weishäupl arXiv:2201.08898 (opened in round 148): Gram-point means of Z for Hecke L-functions, not zeros.

**Verdicts.** Theorem (sharp form, ADR-type, Prop 3(ii)): **KNOWN IN SUBSTANCE** (ADR's method; the κ = 2 bookkeeping is
the degree-2 case; no written degree-2 version found). Smoothed Theorem 2 with the constant −e^{3πi/4}/√2 and the explicit
twist: **KNOWN IN SUBSTANCE** (Weil explicit formula + stationary phase), **written form NOT FOUND**. Exact Voronoi identity
and integer resonance: **KNOWN**. GRH bound x^{3/4}log x for smoothed prime sums: **KNOWN IN SUBSTANCE** (Zhao's Lindelöf
remark). c_E = B_E/L(1, Sym²E) and the line factors for primes: **NOT FOUND** (conjecture; heuristic chain §6(d)).

## 9. GAPs and doubts

GAPs: (1) Prop. 3(ii) is a sketch (ADR's proof with substitutions; ADR's own Lemma 2.4 is quoted without proof there).
(2) A sharp-sharp error o(T) is not proved; the method gives √N T log². (3) Smoothed laws ⇒ sharp prime asymptotic: needs a
Tauberian input. (4) The family average (§6(c)) is not carried out. (5) The prime-side asymptotic with c_E (type II).

Doubts: (1) Lemma 4 is quoted from Milinovich–Ng, whose constants "depend on f"; the N-uniform version is IK's Prop. 5.7 /
Thm 5.38 (not opened). Proposition 1 and Theorem 2 use Lemma 4(ii) only to let the horizontal sides tend to 0 (no
uniformity needed); the N-dependence in Theorem 3(i) and in the GRH bound of §6(b) relies on the uniform form of 4(i). (2) BKY Prop. 8.2 is used as stated; its
proof was not re-checked here; the bounds for p_k (Faà di Bruno) are mine and were checked numerically only through the
second-order term (A3: 9·10^{−3} → 1·10^{−4} relative, one n; tables). (3) The error bound √N T^{½} is far above the data's
agreement; the data's 0.03–0.24 (leading 0.06–0.24) is a computed fact, not a consequence of Theorem 2 at T ≤ 400.
(4) Theorem 2's hypotheses do not cover the data's literal configuration (pure Hann on [20, 400], 20-fold range, no
Gaussian smoothing): it covers dyadic pieces with Δ ≥ T^{½+δ}. The exact Proposition 1 covers any Δ ≥ 1 (the pure Hann
weight is the limit Δ → 0, legitimate under GRH; the n-series then converges only conditionally). (5) Theorem 3(i)
needs T ≥ N^{1+O(δ)}; at the data's T ≈ 400, N^{1/4}T^{3/4} ≈ 400 = T: not informative there. (6) Zero lists: complete up to
402 (round 144); Proposition 1 at 10^{−9} confirms completeness in the tested windows only. Central zeros and γ < 0 are absent
from the lists; their weight is e^{−T₁²/Δ²}. (7) The bad primes enter Λ_E exactly (λ_p^k ln p); additive reduction (λ_p = 0)
does not occur in the 5 curves (semistable), but the proofs cover it.

## 10. What a referee should check first

1. Lemma 1–2: the definition of X, the sign/branch of θ, and that the summand for an off-line zero is X(ρ)^{−1} (not
   e^{2iθ(γ)}); A1 checks the functional equation with w numerically.
2. Proposition 1, step 4: the per-n shift σ₁ → ½ and the termwise identity J_n(σ₁) = n^{−½}I_n that gives absolute convergence
   of the n-series on the critical line (this is what makes the formula unconditional and exact; A4 checks it to 10^{−9}).
3. Theorem 2 (c): the use of BKY Prop. 8.2 (parameters Y = Q = T, V = Δ, X = 1, Z ≍ T; the requirement V ≥ T^{½+δ}) and the
   p_k bounds via Faà di Bruno; (d): the tail n > x₊ on the σ₁ line.
4. Theorem 2 (e): the constant −(e^{3πi/4}/√2)N^{−1/4} (A3, one n at t = 250: the simplified term differs from the exact quadrature by 9.7·10^{−3} relative, of which ≈ 7.6·10^{−3} is the p₁ (H″) term of the Hann weight).
5. Prop. 3(ii) and the claim that the sharp-sharp error is √N T log², not T^{½}log² (the κ = 2 correspondence).
6. Theorem 4's two contour shifts (X^{−1}L = wL(1 − s) entire) and the absence of any other dual main term for r = 1.
7. §6(d) step 4 as the only non-theorem step, and A8 (the line factors equal the recipe's ρ(p), ρ(1/p)).
