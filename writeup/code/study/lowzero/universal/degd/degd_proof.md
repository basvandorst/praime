# The one-point sum Σ X_L(ρ)^{−1} in every degree: exact formula, stationary phase, dual side, singular series

Folder: `scratchpad/degd/` (nothing in the repository was changed). Code: `degd_lib.py` (phase θ_L, windows, oscillatory
sums), `check_degd.py` (sections A–F, every identity asserted; output `check_degd_out.txt`: **44/44 PASS, 35 s** on the busy
2-core machine; run as `(ulimit -v 1500000; timeout 900 nice .venv/bin/python -I check_degd.py)`), `explore.py`,
`explore2.py` (choice of the n-ranges). Sources opened in this round (LaTeX from arXiv): `dl_src/2203.09515` (Kaneko–Thorner),
`dl_src/1505.03122` (Lemke Oliver–Thorner), `dl_src/1911.10497` (Kaczorowski–Perelli); re-read from the degree-2 round:
`scratchpad/gl2proof/dl_bky` (Blomer–Khan–Young §8), `scratchpad/gl2proof/dl_1411_6022` (Ren–Ye).
Builds on the refereed notes `writeup/code/study/lowzero/ecross/proof/gl2_onepoint_proof.md` ("the degree-2 note":
Proposition 1, Theorem 2, Theorem 4) and `.../ecross/family/family_proof.md` §9a, and on `.../universal/general_law.md`.

Labels: **proved** (complete proof here, given the cited results), **proved (cited)** (a result of an opened source, quoted),
**sketch**, **GAP**, **heuristic**, **conjecture**, **computed**.

## 0. Result in short

Class: L = L(s, π), π a cuspidal automorphic representation of GL_d(𝔸_ℚ), d ≥ 1, unitary central character (normalised as in
Kaneko–Thorner), π ≠ 1; conductor Q, root number ε, archimedean parameters μ_1..μ_d; hypothesis **(H∞): Re μ_j ≥ 0 for all j**
(m := min_j Re μ_j). X(s) := Q^{½−s}∏_jΓ_R(1−s+μ̄_j)/Γ_R(s+μ_j), so L(s) = εX(s)L̃(1−s), L̃ = L(·, π̃), and X(½+it)^{−1} = e^{2iθ(t)}.

1. **G1 (exact explicit formula; proved, unconditional, all d; needs (H∞)).** For H = W ∗ g_Δ (W bounded, support [T₁, T₂]),
   h(s) = H(−i(s−½)), r > 0: Σ_ρ m_ρ h(ρ)X(ρ)^{−1}r^{½−ρ} = −(1/2π)Σ_n Λ_L(n)n^{−½}∫H(t)e^{i(2θ(t)−t ln(nr))}dt + E₀(r), all
   nontrivial zeros, both sides absolutely convergent. **Λ_L itself (not its conjugate) and no root number** on the prime side.
   If m > 0, E₀ is explicit: (1/π)∫Hθ′e^{2iθ}r^{−it} − (1/2π)Σ_n **conj Λ_L(n)** n^{−½}∫He^{i(2θ+t ln(n/r))} (Weil's explicit
   formula for this test function). Bound: |E₀(r)| ≪_π r^{3/4}Q^{−3/4}e^{9/(16Δ²)}log(Q(T₂+Δ+2)^d)[(T₂−T₁)(1+T₁/2)^{−3d/4}
   + (1+T₁)e^{−T₁²/(4Δ²)}] (m > ¼). Checked to 3·10⁻¹³–4·10⁻¹² on a complex character mod 5 (d = 1), Δ (d = 2), Sym² 11a (d = 3),
   with twists r = 3/2, 2; with conj Λ_L in place of Λ_L the identity fails by 8.0.
2. **G2 (stationary phase; proved given Blomer–Khan–Young, all d).** The n-th term is
   −A_d Q^{−1/(2d)}Λ_L(n)n^{−½+1/(2d)}e(−d(n/Q)^{1/d})(n/Q)^{ib/d}H(2π(n/Q)^{1/d})(1 + O(1/t)), **A_d = e^{i(φ₀+π/4)}/√d,
   φ₀ = (π/2)ΣRe μ_j − dπ/4**, b = ΣIm μ_j; first-order factor (real μ, or b = 0) **1 + iκ/t**, κ = (d²−1)/(24d) +
   ½ΣRe(μ_j − μ_j²), plus (i/2d)H′ + (it/2d)H″. Summed (smooth windows of length ≍ T): **error ≪ T^{(d−1)/2}** (≍ Q^{½}T^{(d−1)/2}
   under Ramanujan) with the closed-form leading kernel: **o(T) only for d ≤ 2** (d = 3: O(T), the size of the conjectured main
   term). With the K-term expansion as main term the error is ≪ T^{(d+1)/2−K}: **o(T) in every degree for K ≥ ⌈d/2⌉**, and
   (Corollary G2.2) the zero-side law and the prime-side law (closed-form kernel) are **equivalent in every degree** for smooth
   windows. Constants checked: κ = −0.8889 (Sym²), −0.34375 (genus 2), −30.1875 (Δ), 0 (d = 1), to 4 digits.
3. **G3 (exact dual side; proved, all d, no (H∞)).** Σ_n λ(n)n^{−½}∫H e^{i(2θ−t ln(nr))} = **ε Σ_m conj λ(m)** m^{−½}Ĥ(ln(r/m)).
   Log-weighted form G3′ (m > 0): Σ_b λ(b)ln b·b^{−½}I_b(a) = εΣ_m conj λ(m)m^{−½}∫H(2θ′ − ln m)e^{−it ln(a/m)}. Checked on
   Sym² 11a (r = 1, 2, 3/2; flipped ε fails by 3.4–15.8) to 1.5·10⁻¹³ and on χ mod 5 (complex ε; unconjugated λ or ε̄ fail).
   Over all integers the resonance has main term ε∫H at r = 1 and ε conj λ(q)q^{−½}∫H at r = q (known: Kaczorowski–Perelli's
   standard twist; Ren–Ye).
4. **Singular series (local algebra proved; the sum over moduli heuristic).** For every q the twisted sum
   Σ_n c_q(n)λ(n)n^{−½}I_n is given exactly by G3 (checked, q = 2, 3, 5, 6, d = 3, to 2.5·10⁻¹³); its resonant part is
   ε∫H·∏_{p|q}(p−1)(1 − G_p(1/p)), **G_p(X) := F_p(X)/(1−X)**, F_p(X) = Σ_k conj λ(p^k)μ_L(p^k)X^k, and F_p(1) = 0 at every
   unramified p (Ramanujan). Hardy–Littlewood: Σ_q μ(q)∏_{p|q}(1 − G_p(1/p)) = ∏_p G_p(1/p) = **c_L = B_L(1)/L(1, Ad π)**. The
   route Λ_L = (λ ln) ∗ μ_L + G3′ gives the same c_L **and** the same line factors ρ_L(r) as general_law.md, in every degree.
   Exact local factor: **d = 2: G_p = 1 − (|λ(p)|²−1)/p; d = 3: G_p = 1 − (|λ(p)|²−1)/p + ω(p)conj λ(p³)/p²**; in general
   1 − λ_{Ad}(p)/p + O(p^{−2}), λ_{Ad}(p) = |λ(p)|² − 1. Partial products at p ≤ 10⁶: 1.94469 (Sym² 11a; c = 1.94470), 3.27655
   (Sym² 15a; 3.27678), 2.36059 at 3·10⁵ (Δ; 2.36082), 1.25000 (χ mod 5; 5/4).
5. **Novelty.** G1: KNOWN IN SUBSTANCE (Weil's explicit formula with g = hX^{−1}); written form for d ≥ 2: NOT FOUND. G2: per-term
   KNOWN IN SUBSTANCE (asymptotics of the GL(d) Voronoi kernel, Ren–Ye Theorem 1.1, c₀ = −1/√m); the summed form with A_d, κ and
   the degree threshold: NOT FOUND. G3 and its integer resonance: KNOWN. The degree-d singular series of the prime resonance with
   G_p = F_p/(1−X): NOT FOUND. What stays open in every degree ≥ 2: the sum over moduli (minor arcs / type II).

## 1. Setting, cited facts

**Class.** 𝔉_d as in Kaneko–Thorner (arXiv:2203.09515, opened, §2): π cuspidal automorphic on GL_d(𝔸_ℚ), unitary central
character trivial on the diagonal ℝ_{>0}; we exclude π = 1 (ζ, which has poles). Analytic normalisation:
L(s) := L(s, π) = Σλ(n)n^{−s} = ∏_p∏_{j≤d}(1 − α_j(p)p^{−s})^{−1}, γ(s) := ∏_jΓ_R(s+μ_j), Γ_R(s) = π^{−s/2}Γ(s/2),
Λ(s) := Q^{s/2}γ(s)L(s). π̃: λ_π̃(n) = conj λ(n), α_j(p, π̃) = conj α_j(p), μ_j(π̃) = conj μ_j, Q(π̃) = Q.
L̃(s) := L(s, π̃) = conj L(s̄). −L′/L(s) = Σ_n Λ_L(n)n^{−s}, Λ_L(p^k) = Σ_jα_j(p)^k ln p; Λ_L̃(n) = conj Λ_L(n).

**Facts used (all from opened sources).**
- (F1) [Kaneko–Thorner (A), (B), Proposition "prop:class": π ∈ 𝔉_m ⇒ L(s, π) ∈ 𝔖(m)]. Dirichlet series and Euler product converge
  absolutely for σ > 1; Λ(s) is entire of order 1 (r_π = 0 for π ≠ 1) and Λ(s) = εΛ(1−s, π̃) with |ε| = 1.
- (F2) [KT (D)] Σ_n|Λ_L(n)|n^{−1−η} ≤ d/η + d log C(π) + O(d²) (η > 0), C(π) = Q∏(|μ_j| + 3); and
  Σ_{x<n≤xe^{1/T}}|Λ_L(n)| ≪ d x/T for x ≥ c d^{182d⁴}(C(π)T)^{144d³}. With T = 1: **Chebyshev** Σ_{n≤x}|Λ_L(n)| ≪_π x.
- (F3) [KT (E)] L(s) ≠ 0 for Re s ≥ 1.
- (F4) [KT Lemma "number-of-zeros"; via Soundararajan–Thorner and Iwaniec–Kowalski Prop. 5.7, not opened]
  #{ρ : 0 < β < 1, |γ − t| ≤ 1} ≪ d log(C(π)(2 + |t|)), every t ∈ ℝ; Hadamard product of Λ (KT eq. "Hadamard-product").
- (F5) [Lemke Oliver–Thorner, arXiv:1505.03122, opened, §2.1: Luo–Rudnick–Sarnak; Müller–Speh at ramified p]
  |α_j(p)| ≤ p^{½−1/(d²+1)} for **all** p; Re μ_j > −½ + 1/(d²+1); |Λ_L(n)| ≤ dΛ(n)n^{½−1/(d²+1)}.
- (F6) [KT §2] L(s, π × π̃) has a simple pole at s = 1; its local factor at p ∤ Q is ∏_{i,j}(1 − α_i(p)conj α_j(p)p^{−s})^{−1}.
- (F7) [Blomer–Khan–Young, arXiv:1203.2573, §8, LaTeX re-read]: **Lemma 8.1** (repeated integration by parts:
  w smooth on [α, β], w^{(j)} ≪ X U^{−j}, |h′| ≥ R, h^{(j)} ≪ Y Q_B^{−j} (j ≥ 2) ⇒ ∫w e^{ih} ≪_A (β−α)X[(Q_B R/√Y)^{−A} + (RU)^{−A}]);
  **Proposition 8.2** (0 < δ < 1/10, Z := Q_B + X + Y + V₁ + 1, Y ≥ Z^{3δ}, V₁ ≥ V ≥ Q_B Z^{δ/2}/Y^{½}, w supported on J of length V₁
  with w^{(j)} ≪ XV^{−j}, unique t₀ ∈ J with h′(t₀) = 0, h″ ≫ YQ_B^{−2}, h^{(j)} ≪ YQ_B^{−j} (j ≥ 1) ⇒
  ∫w e^{ih} = e^{ih(t₀)}h″(t₀)^{−½}Σ_{n≤3A/δ}p_n(t₀) + O(Z^{−A}), p_n = √(2π)e^{πi/4}(n!)^{−1}(i/(2h″(t₀)))^n G^{(2n)}(t₀),
  G = w e^{iΦ}, Φ = h − h(t₀) − ½h″(t₀)(t−t₀)², and **(d/dt₀)^j p_n ≪ X(V^{−j} + Q_B^{−j})((V²Y/Q_B²)^{−n} + Y^{−n/3})**).
  (BKY's "Q" is written Q_B here to avoid a clash with the conductor.)

**Hypothesis (H∞).** Re μ_j ≥ 0 for all j; m := min Re μ_j. It holds for every π tempered at ∞ in the usual sense used here:
holomorphic newforms (μ = (k−1)/2, (k+1)/2), Dirichlet characters (0 or 1), Artin representations, tempered Maass forms
(Re μ = 0), Sym² and genus-2 L-functions (μ ≥ ½), and it is the gamma-factor condition of the Selberg class. It fails only for
non-tempered archimedean components (exceptional eigenvalues), where some trivial zeros −μ_j lie inside the critical strip.
(GRC) |α_j(p)| ≤ 1 for all p (Ramanujan): used only for box windows in G2 and in the singular-series algebra.

**Windows.** W bounded measurable, |W| ≤ 1, supported in [T₁, T₂], 0 < T₁ < T₂; g_Δ(z) = e^{−z²/Δ²}/(√πΔ), Δ > 0; H := W ∗ g_Δ
(entire); h(s) := H(−i(s−½)), so h(½+it) = H(t). Since |g_Δ(t − iy)| = e^{y²/Δ²}g_Δ(t):
|h(σ+it)| ≤ e^{(σ−½)²/Δ²}(|W| ∗ g_Δ)(t) (bounded in vertical strips, Gaussian decay away from [T₁, T₂]).
Ĥ(ξ) := ∫H(t)e^{−iξt}dt = Ŵ(ξ)e^{−Δ²ξ²/4}. I_n(r) := ∫H(t)e^{i(2θ(t) − t ln(nr))}dt, I_n := I_n(1);
I_n^−(r) := ∫H(t)e^{i(2θ(t) + t ln(n/r))}dt.

## 2. Lemmas

**Lemma 1 (functional equation and phase; proved).** With X(s) := Q^{½−s}∏_jΓ_R(1−s+μ̄_j)/Γ_R(s+μ_j):
L(s) = εX(s)L̃(1−s), and X(½+it)^{−1} = e^{2iθ(t)} for real t with
θ(t) := (t/2)ln Q + Σ_j[Im lnΓ((½+μ_j+it)/2) − ((t + Im μ_j)/2)ln π]   (principal lnΓ, analytic on Re z > 0).
Hence L(½+it) = εe^{−2iθ(t)}conj L(½+it).
*Proof.* (F1): Q^{s/2}γ(s)L(s) = εQ^{(1−s)/2}∏Γ_R(1−s+μ̄_j)L̃(1−s); divide. On s = ½+it: Γ_R(½−it+μ̄_j) =
conj Γ_R(½+it+μ_j) (conj of π^{−w} is π^{−w̄}), so X(½+it)^{−1} = Q^{it}∏_j e^{2i arg Γ_R(½+it+μ_j)}, and
arg Γ_R(½+it+μ) = Im lnΓ((½+μ+it)/2) − ((t + Im μ)/2)ln π; Re((½+μ+it)/2) > 0 by (F5), so the branch is continuous. The
last claim uses L̃(½−it) = conj L(½+it). ∎ *Computed (A1, A1b):* 3.5·10⁻¹³ (real and complex μ, mpmath); FE with ε via PARI
`lfun` to 1.3·10⁻¹⁴ for χ mod 5 (non-self-dual, L̃ = L(χ̄)) and Sym² 11a; with −ε it fails by 2.0. For real μ_j this θ is
general_law.md's θ_L (the ln π terms combine to (t/2)ln(Q/π^d)).

**Lemma 2 (continuation, strip, growth; proved).** X^{−1} is meromorphic with poles only in 𝒫 := {−μ_j − 2k : j, k ≥ 0}
(poles of γ) and zeros only in {1 + μ̄_j + 2k}; it is **holomorphic and zero-free on −m < σ < 1 + m** (m = min Re μ_j; under
(H∞) this contains 0 < σ < 1; in general (F5) it contains |σ − ½| < 1/(d²+1)). For a zero ρ = ½ + iα_ρ (α_ρ = γ − i(β−½))
the summand is X(ρ)^{−1} = e^{2iθ(α_ρ)} with θ continued analytically; under GRH it is e^{2iθ(γ)}. For σ in a compact set and
|t| ≥ 1 + 2max_j|μ_j|:
|X(σ+it)^{−1}| = (Q∏_j|t + Im μ_j|/2π)^{σ−½}(1 + O(1/|t|)),   −X′/X(σ+it) = ln(Q(|t|/2π)^d) + O(1/|t|).
*Proof.* Γ has no zeros; the poles of Γ_R(s+μ_j) are at s + μ_j ∈ −2ℕ₀, those of Γ_R(1−s+μ̄_j) at 1 − s + μ̄_j ∈ −2ℕ₀.
Stirling |Γ(x+iy)| = √(2π)|y|^{x−½}e^{−π|y|/2}(1 + O(1/|y|)) for bounded x applied to (s+μ_j)/2 and (1−s+μ̄_j)/2 (same |y| =
|t + Im μ_j|/2) gives |Γ_R(s+μ_j)/Γ_R(1−s+μ̄_j)| = (|t+Im μ_j|/2π)^{σ−½}(1+O(1/t)); −X′/X = ln Q + Σ_j[ψ_R(s+μ_j) +
ψ_R(1−s+μ̄_j)] with ψ_R = Γ_R′/Γ_R = ½ψ(·/2) − ½ln π. ∎ *Computed (A2):* t·|ratio − 1| ≤ 0.02 at t = 50, 500, 5000, σ = −¼, 5/4.

**Lemma 3 (Stirling for θ; proved, checked).** With b := ΣIm μ_j, φ₀ := (π/2)ΣRe μ_j − dπ/4, c₁ := d/24 + ½ΣRe(μ_j − μ_j²):
2θ(t) = t ln(Q(t/2π)^d) − dt + b ln(t/2π) + φ₀ + c₁/t + O(t^{−2})   (O(t^{−3}) if all μ_j are real),
2θ′(t) = ln(Q(t/2π)^d) + b/t − c₁/t² + O(t^{−3}),  2θ″(t) = d/t − b/t² + O(t^{−3}),  θ^{(j)}(t) ≪_{j,μ} t^{1−j} (j ≥ 2).
*Proof.* lnΓ(z) = (z−½)ln z − z + ½ln 2π + 1/(12z) + O(|z|^{−3}) at z = a + iτ, a = (½+μ)/2, τ = t/2:
ln z = ln τ + iπ/2 − ia/τ + a²/(2τ²) + O(τ^{−3}), so
Im lnΓ(a+iτ) = τ ln τ − τ + Im(a)ln τ + (π/2)Re(a − ½) + Re[(a − a²)/2 − 1/12]/τ + O(τ^{−2}),
and (a − a²)/2 − 1/12 = 1/96 + μ/8 − μ²/8. Sum over j with the −((t+Im μ_j)/2)ln π terms. For real a the expansion of
lnΓ(a+iτ) has real coefficients times (iτ)^{−k}, whose imaginary parts vanish for even k. ∎ (d = 2, μ = {½, 3/2}: φ₀ = π/2,
c₁ = −1/6, the degree-2 note's Lemma 3. ζ: c₁ = 1/24, the Riemann–Siegel θ.) *Computed (A3):* the residual times t³ (real μ:
Sym², genus 2) and times t² (complex μ) converges at t = 100, 400, 1600. general_law.md's φ₀ = (π/2)Σμ_j − dπ/4 is
**confirmed** (with Re μ_j for complex μ_j); its "O(1/t)" is c₁/t, present **in every degree** (not only d ≥ 3), and for
b ≠ 0 there is an extra b ln(t/2π) (b = 0 for all L-functions tested here).

**Lemma 4 (zeros and good heights; proved given F1–F4).** Assume (H∞). (i) The zeros of Λ lie in 0 < σ < 1, and none lies
in 𝒫. (ii) #{ρ : |γ − Y| ≤ 1} ≪ d log(C(π)(|Y|+2)). (iii) For every real Y with |Y| ≥ Y₀(π) there is Y′ ∈ [Y, Y+1] with
|L′/L(σ + iY′)| ≪_π log²(|Y|+3) uniformly for −1 ≤ σ ≤ 2.
*Proof.* (i) L ≠ 0 for σ ≥ 1 (F1, F3), γ has no zeros, so Λ ≠ 0 there; Λ(s) = εΛ̃(1−s) gives Λ ≠ 0 for σ ≤ 0. 𝒫 ⊂ {σ ≤ 0}
by (H∞). (ii) is (F4). (iii) Hadamard: Λ′/Λ(s) = B + Σ_ρ[1/(s−ρ) + 1/ρ]. For s = σ+iY′, s₂ = 2+iY′: the zeros with
|γ − Y′| ≥ 1 contribute to Λ′/Λ(s) − Λ′/Λ(s₂) at most (2−σ)Σ_{k≥1}k^{−2}#{k ≤ |γ−Y′| < k+1} ≪ log(C(|Y|+2)) by (ii); those with
|γ − Y′| < 1 (n′ ≪ log(C(|Y|+3)) of them, counting |γ − Y − ½| < 3/2) contribute ≤ n′ + Σ1/|Y′−γ|. The union of the intervals
|Y′ − γ| < 1/(4n′+4) has measure ≤ ½, so some Y′ ∈ [Y, Y+1] has Σ1/|Y′−γ| ≤ n′(4n′+4) ≪ log²(C(|Y|+3)). Finally
|Λ′/Λ(s₂)| ≤ |L′/L(s₂)| + ½ln Q + Σ|ψ_R(s₂+μ_j)| ≪ log C(π)(|Y|+2) by (F2), and Λ′/Λ − L′/L = ½ln Q + Σψ_R(s+μ_j) ≪ log(|Y|+C)
for |Y′| > max|Im μ_j| + 1 (away from the poles of ψ_R). ∎ (The degree-2 note cited this from Milinovich–Ng/IK; here it is
derived from the opened Kaneko–Thorner lemma.)

## 3. Theorem G1: the exact explicit formula in degree d

**Theorem G1 (proved; unconditional; π ∈ 𝔉_d, π ≠ 1, (H∞)).** For W, Δ, H, h as in §1 and r > 0:

    Σ_ρ m_ρ h(ρ) X(ρ)^{−1} r^{½−ρ} = −(1/2π) Σ_{n≥2} Λ_L(n) n^{−½} I_n(r) + E₀(r),

where ρ runs over all zeros of Λ (multiplicity m_ρ; real zeros, central zeros and γ < 0 included) and both sides converge
absolutely. Here

(a) E₀(r) = −(1/2πi)∫_{(−¼)} h(s)X(s)^{−1}r^{½−s}(L′/L)(s)ds − Σ_{s₀∈𝒫, Re s₀ > −¼} Res_{s₀}[h X^{−1}r^{½−s}L′/L]
(the residue sum is empty when m > ¼; the residues are O_π(sup_{|s−s₀|≤ρ_π}|h(s)|), i.e. Gaussian-small in the distance from
−Im μ_j to [T₁, T₂]);

(b) if m > 0 (Weil form):
    E₀(r) = (1/π)∫H(t)θ′(t)e^{2iθ(t)}r^{−it}dt − (1/2π)Σ_{n≥2} conj(Λ_L(n)) n^{−½} I_n^−(r),   absolutely convergent;

(c) if m > ¼: |E₀(r)| ≤ C_π r^{3/4}Q^{−3/4}e^{9/(16Δ²)}log(Q(T₂+Δ+2)^d)·[(T₂−T₁)(1+T₁/2)^{−3d/4} + (1+T₁)e^{−T₁²/(4Δ²)}],
C_π depending on d, μ and Σ|Λ_L(n)|n^{−5/4} (≪ d log C(π) + d², F2).

Under GRH the left side is Σ_γ H(γ)e^{2iθ(γ)}r^{−iγ}. **Conjugations:** the stationary (resonant) prime sum carries Λ_L of L
itself; conj Λ_L = Λ_L̃ appears only in the non-stationary part E₀ (phase 2θ + t ln n, no stationary point at t > 0 large);
**ε does not appear** (it cancels in L′/L).

*Proof.* (1) F(s) := h(s)X(s)^{−1}r^{½−s}(L′/L)(s) is meromorphic on −¼ ≤ σ ≤ 5/4. Its poles there: the zeros of L, which are
the zeros of Λ (Lemma 4(i); residue m_ρh(ρ)X(ρ)^{−1}r^{½−ρ}, X^{−1} being holomorphic at ρ by Lemma 2) and the trivial zeros
(= 𝒫, since Λ = Q^{s/2}γL is entire and non-zero on σ ≤ 0); and the poles of X^{−1}, which also lie in 𝒫.
(2) Integrate over the rectangle [−¼, 5/4] × [−Y″, Y′] with good heights from Lemma 4(iii) (both signs). On the horizontal
sides |h| ≤ e^{1/Δ²}·½erfc((Y − T₂)/Δ), |X^{−1}| ≪ (QY^d)^{3/4} (Lemma 2), |L′/L| ≪ log²Y, |r^{½−s}| ≤ max(r, r^{−1})^{3/4}:
they tend to 0. The zero sum converges absolutely: |h(ρ)| ≤ e^{1/(4Δ²)}(|W| ∗ g_Δ)(γ), |X(ρ)^{−1}| ≪ (Q(|γ|+2)^d)^{½} for
|γ| ≥ γ₀ (finitely many others), and Lemma 4(ii).
(3) Right side σ₁ = 5/4: Σ|Λ_L(n)|n^{−5/4}∫|hX^{−1}|r^{−3/4}dt < ∞ (F2, Lemma 2), so by Fubini (1/2πi)∫_{(σ₁)}F =
−(1/2π)Σ_nΛ_L(n)J_n(σ₁), J_n(σ) := ∫h(σ+it)X^{−1}(σ+it)r^{½−σ−it}n^{−σ−it}dt.
(4) For each n, h X^{−1}r^{½−s}n^{−s} is holomorphic on ½ ≤ σ ≤ 5/4 (𝒫 ⊂ {σ ≤ 0}) with Gaussian decay uniform in σ, so
J_n(5/4) = J_n(½) = n^{−½}I_n(r) (Lemma 1). As Σ|Λ_L(n)||J_n(5/4)| < ∞ termwise, Σ_nΛ_L(n)n^{−½}I_n(r) converges absolutely.
(5) The left side is E₀ in form (a). For (b), take the left side at σ₀ ∈ (−min(m, ¼), 0) instead (𝒫 lies left of it). There
L′/L(s) = X′/X(s) − (L̃′/L̃)(1−s) = X′/X(s) + Σ_nΛ_L̃(n)n^{s−1}, absolutely convergent (1 − σ₀ > 1, F2 for π̃). Each piece is
holomorphic on σ₀ ≤ σ ≤ ½ (X′/X has poles only in 𝒫 and in Re s ≥ 1 + m) and is shifted to σ = ½ separately: on that line
X′/X(½+it) = −2θ′(t) (differentiate X(½+it)^{−1} = e^{2iθ(t)}), giving (1/π)∫Hθ′e^{2iθ}r^{−it}; and the n-th term gives
−(1/2π)Λ_L̃(n)n^{−½}I_n^−(r); absolute convergence of that series follows termwise as in (4).
(6) Bound (c), on σ₀ = −¼ (no point of 𝒫 on or right of it when m > ¼): |h| ≤ e^{9/(16Δ²)}(|W|∗g_Δ)(t),
|X^{−1}(−¼+it)| ≤ C_μQ^{−3/4}(1+|t|)^{−3d/4}, |L′/L(−¼+it)| ≤ ln(Q(|t|+2)^d) + C_μ + Σ|Λ_L(n)|n^{−5/4}, |r^{½−s}| = r^{3/4};
split |t| ≥ T₁/2 (factor (1+T₁/2)^{−3d/4}, ∫|W|∗g_Δ ≤ T₂−T₁) and |t| < T₁/2 ((|W|∗g_Δ)(t) ≤ e^{−T₁²/(4Δ²)}). ∎

*Remarks.* (1) Only the ℓ¹ bound (F2) is used for absolute convergence; the right and left lines can be 1 + η and −η for any
η > 0. With only the pointwise LRS bound |Λ_L(n)| ≤ dΛ(n)n^{θ_d}, θ_d = ½ − 1/(d²+1), one must take σ₁ > 1 + θ_d and
σ₀ < −θ_d; then form (b) needs m > θ_d and the general form (a) collects the trivial residues in (σ₀, 0]. (2) (b) is Weil's
explicit formula for the test function g(½+it) = H(t)e^{2iθ(t)}r^{−it}: its gamma term (1/2π)∫g·[ln Q + Σψ_R(½+it+μ_j) +
Σψ_R(½−it+μ̄_j)]dt equals (1/π)∫gθ′ because that bracket is 2θ′(t). (3) The bound (c) is crude: by (b) and
BKY Lemma 8.1 (the phases 2θ(t) − t ln r and 2θ(t) + t ln(n/r) have no stationary point at large t) E₀ is ≪ T^{−A} for
smooth windows (sketch); measured |E₀| ≤ 5·10⁻¹³ (quadrature floor) in all cases below.

*Computed (section C, `check_degd_out.txt`).* I_n by the trapezoidal rule on a grid of step 0.05 (the integrand is analytic in
a strip of half-width ½ + m and Gaussian-small at the grid ends), over all prime powers n ≤ N, N chosen from the decay of I_n
(`explore2.py`: |I_n| reaches the 10⁻¹³ floor at ln n ≈ 8.5, 9.8, 14 for the three windows); E₀ from (b) with n ≤ 2000;
zeros from `universal/zeros/` (negative-height zeros: −γ for self-dual L, −γ(L̃) for χ):

| L | d | window | zeros | prime side | \|E₀\| | \|zeros − primes − E₀\| | #zeros / #prime powers |
|---|---|---|---|---|---|---|---|
| χ mod 5 (X5) | 1 | [400,440], Δ = 4 | −7.8585015897 − 3.4117778311i | same to 10 digits | 1.9e−13 | **1.8e−12** | 82 / 1280 |
| X5, r = 3/2 | 1 | same | −8.3772659570 − 2.9945268042i | same | 5.5e−13 | **4.3e−12** | 82 / 2328 |
| Δ (D12) | 2 | [200,260], Δ = 5 | −20.1209720495 − 2.3972294850i | same | 3.7e−13 | **6.4e−13** | 137 / 3320 |
| Sym² 11a (S11a) | 3 | [40,50], Δ = 4 | −4.0456530563 + 1.2091698906i | same | 4.3e−14 | **2.9e−13** | 97 / 100283 |
| S11a, r = 2 | 3 | same | −4.8931589263 − 0.9693024473i | same | 6.8e−14 | **3.0e−13** | 97 / 100283 |

With conj Λ_L (= Λχ̄) in place of Λ_L for X5: |difference| = 8.025 (the conjugation is detected). The S11a coefficients
(built from a_p(11a) and PARI's local factor at 11) equal PARI `lfunan` of the object whose zeros were computed to 2·10⁻¹⁵ for
n ≤ 3000; Δ likewise. Genus 2 (degree 4) was **not run**: at every height the needed range is n ≲ 2·10⁶ (Fresnel scale
√(t/d) small), and PARI produced the genus-2 coefficients to 2·10⁵ in 62 s, so the run would exceed the time budget (GAP 7).

## 4. Theorem G2: stationary phase in degree d

Write φ_n(t) := 2θ(t) − t ln n, t_n the solution of 2θ′(t_n) = ln n (unique for t ≥ t₀(μ), since θ″ > 0 there),
t_n* := 2π(n/Q)^{1/d} (t_n = t_n* − b/d + O(1/t)), A_d := e^{i(φ₀+π/4)}/√d, κ := c₁ − 1/(24d) = (d²−1)/(24d) + ½ΣRe(μ_j − μ_j²).

**Lemma 5 (one term; proved from F7 and Lemma 3).** Let H be one of the windows of Theorem G2 at height T (so
[T₁, T₂] ⊂ [T, 2T]) and t_n ∈ [T/2, 4T]. Then
I_n = 2πA_d (n/Q)^{1/(2d)} e(−d(n/Q)^{1/d}) (n/Q)^{ib/d} H(t_n*) (1 + O(1/T)) + (derivative terms),
and, if b = 0 (e.g. all μ_j real),
I_n = 2πA_d (n/Q)^{1/(2d)} e(−d(n/Q)^{1/d}) [H(1 + iκ/t) + (i/2d)H′ + (it/2d)H″](t_n*) + O(T^{½}·T^{−2}) (smooth windows).
*Proof.* BKY Prop. 8.2 at t₀ = t_n: p₀ gives √(2π)e^{πi/4}H(t_n)e^{iφ_n(t_n)}φ_n″(t_n)^{−½}. Lemma 3: φ_n(t_n) =
−dt_n* + (b/d)ln(n/Q) + φ₀ + O(1/t) and, for b = 0, = −dt_n* + φ₀ + c₁/t_n* + O(t^{−3}) (from 2θ′(t_n) = ln n,
t_n = t_n*(1 + c₁/(dt²) + …)); (2π/φ_n″)^{½} = (2πt_n*/d)^{½}(1 + O(t^{−1})) (O(t^{−2}) if b = 0), and (2πt*/d)^{½} =
2π(n/Q)^{1/(2d)}/√d. The order-1/t terms of p₁, p₂, p₃ are (i/2φ″)H″, −iφ‴H′/(2φ″²), −iφ⁗H/(8φ″²) and 5iφ‴²H/(24φ″³); with
φ″ = d/t, φ‴ = −d/t², φ⁗ = 2d/t³ they are (it/2d)H″, (i/2d)H′ and (−1/4 + 5/24)(i/dt)H = −iH/(24dt); together with e^{ic₁/t}
this is the bracket. ∎
*Computed (B, one n, smooth window, t = 100…800):* |I/lead − 1|·t → |κ|: 0.8889 (Sym², κ = −8/9), 0.34375 (genus 2,
κ = −11/32), 30.18 (Δ, κ = −30.1875 = −k(k−2)/4 − 3/16 at k = 12), 0 for d = 1 (κ = 0 for μ ∈ {0, 1}: the degree-1 kernel has no
power corrections, only edge terms); after the first-order factor |I/first − 1|·t² is bounded (0.24, 0.016, 448). Complex μ
with b = 0 (μ = 0.4 ± 3i, 1): κ = 9.3511 confirmed. Complex μ with b = 1.3: the leading constant with (n/Q)^{ib/d} holds
(error O(1/t)); the first-order formula above is **not** claimed for b ≠ 0.

**Theorem G2 (proved given F2, F7; unconditional except where (GRC) is stated).** Fix δ ∈ (0, 1/10), T ≥ T₀(δ, π),
[T₁, T₂] ⊂ [T, 2T]. Windows:
- (smooth) W(t) = w₀((t−T₁)/(T₂−T₁)), T₂ − T₁ ≥ T/4, w₀ ∈ C_c^∞((0,1)), 1 ≤ Δ ≤ T^{1−δ};
- (Hann type) as above with w₀ ∈ C¹, ‖w₀″‖_∞ < ∞, T^{½+δ} ≤ Δ ≤ T^{1−δ};
- (box, assuming GRC) W = 1_{[T₁,T₂]}, T^{½+δ} ≤ Δ ≤ T^{1−δ}.

(a) **Closed-form leading kernel.**

    Σ_ρ m_ρ h(ρ)X(ρ)^{−1} = −A_d Q^{−1/(2d)} Σ_n Λ_L(n) n^{−½+1/(2d)} e(−d(n/Q)^{1/d}) (n/Q)^{ib/d} H(2π(n/Q)^{1/d}) + R,
    R ≪_{π,w₀,δ} T^{(d−1)/2}                                        (smooth, Hann type),
    R ≪_{d,μ,δ} Q^{½}T^{(d−1)/2} + Q^{½}T^{(d+1)/2}Δ^{−1}log(QT)        (box, GRC).

Under GRC the smooth/Hann bound is ≪_{d,μ,w₀,δ} Q^{½}T^{(d−1)/2} uniformly in Q. With the twist r^{½−ρ}: n → nr in the
kernel, factor r^{1/(2d)}, same R (constants depend on r).

(b) **K-term expansion (smooth windows).** For every K ≥ 1,
Σ_ρ m_ρh(ρ)X(ρ)^{−1} = −(1/2π)Σ_n Λ_L(n)n^{−½}e^{iφ_n(t_n)}φ_n″(t_n)^{−½}Σ_{k<3K}p_k(t_n) + O_{π,w₀,K}(T^{(d+1)/2−K}).

(c) **When is R = o(T)** (the conjectured main term −εc_L∫H/2π is ≍ T): (a) gives o(T) **iff d ≤ 2** (d = 1: O(1); d = 2:
O(T^{½}), the degree-2 note's √N T^{½}); for **d = 3 the leading-kernel error is O(T), the size of the main term**, and for
d ≥ 4 it is larger. (b) gives o(T) for every d once **K ≥ ⌈d/2⌉** (d = 3, 4: the first-order bracket of Lemma 5 suffices
when b = 0).

*Proof.* As the degree-2 note's Theorem 2, with these changes. Start from G1 (E₀ ≪ Q^{−3/4}); cut-off ω around the window,
x₊ := Q(8T/2π)^d.
(i) Tails of H and n ≤ x₊ with t_n ∉ [T/2, 4T]: on supp(Hω) ⊂ [3T/4, 5T/2], |φ_n′| ≥ d ln(3/2) − O(1/T), φ_n^{(j)} ≪ T^{1−j};
BKY Lemma 8.1 (X = 1, U = Δ or T, Q_B = Y = T, R ≍ 1) gives ≪ T(T^{−A/2} + Δ^{−A}) per n, ≪ (QT^d)^{½}T^{1−A/2} in total.
(ii) t_n ∈ [T/2, 4T]: BKY Prop. 8.2 on J = [T/4, 8T] with h = φ_n (h″ = 2θ″ ≍ d/T, h^{(j)} ≪ T^{1−j}, |h′| ≤ d ln 16 + o(1)
on J, so Y = Q_B = T), w = Hω, X = 1, V = Δ (box, Hann) or T (smooth), Z ≍ T; the conditions are those of degree 2 (V ≥ T^{½+δ}
suffices; the degree enters only the constants). Each term after p₀ is ≪ T^{−1} relative to p₀, pointwise, for smooth and
Hann-type windows: the degree-2 note's Faà di Bruno bound |p_k| ≪ T^k|w^{(2k)}| + Σ_{1≤j≤2k−3}T^{(2j−k)/3}|w^{(j)}| +
T^{⌊2k/3⌋−k}|w| uses only Φ^{(m)}(t_n) ≪ T^{1−m} and 1/φ_n″ ≍ T/d, so it holds in every degree (BKY's own uniform bound
(V²Y/Q_B²)^{−n} + Y^{−n/3} is weaker for n = 1, 2 and is used only in (b)); so their sum over n is ≪ Σ_{n≤x₊}|Λ_L(n)|n^{−½}T^{½}T^{−1} ≪ (QT^d)^{½}T^{−½} = Q^{½}T^{(d−1)/2} by Chebyshev (F2). For
the box the k = 1 term (t/2d)H″ is concentrated on the two edges; summing over unit t-intervals, each containing ≍ dn/t ≥
n^{1−1/d} integers, GRC and Brun–Titchmarsh (d ≥ 2; trivially with a log for d = 1) give Σ_{t_n∈[k,k+1]}|Λ_L(n)| ≪ dn/t, hence
Σ_n|Λ_L(n)|n^{−½}t_n^{½}|f(t_n)| ≪ Q^{½}T^{(d−1)/2}log(QT)(‖f‖₁ + ‖f′‖₁), and f = tH″ gives Q^{½}T^{(d+1)/2}Δ^{−1}log. (Without
GRC the short-interval bound of F2 is available only for x ≥ (C(π)T)^{144d³}, far beyond x ≍ QT^d: GAP 5.) The O(Z^{−A})
remainders sum to ≪ (QT^d)^{½}T^{−A}.
(iii) n > x₊: on σ₁ = 1 + 1/log x₊, |X^{−1}(σ₁+it)n^{−σ₁−it}| has phase derivative ≥ ln(n/x₊) + d ln 2 − o(1); k integrations by
parts give |J_n(σ₁)| ≪ T(QT^d)^{σ₁−½}n^{−σ₁}(ln(4n/x₊))^{−k}(Δ^{−1} + T^{−½})^k, summed (F2 on e-adic intervals) ≪ Q^{½}T^{1+d/2−k/2};
k = d + 3.
(iv) Leading term: Lemma 5; the replacements t_n → t_n*, φ″(t_n) → d/t_n*, H(t_n) → H(t_n*) cost O(1/T) relatively, summed
≪ Q^{½}T^{(d−1)/2} as in (ii).
(b): keep the BKY terms k < 3K exactly; the rest is ≪ T^{½}·T^{−K} per n (BKY's derivative bound with j = 0), summed ≪
(QT^d)^{½}T^{½−K}. (c): compare exponents: (d−1)/2 < 1 ⟺ d < 3; (d+1)/2 − K < 1 ⟺ K > (d−1)/2. ∎

*Why the threshold is real.* The first correction in (a) is itself a resonant prime sum, −A_dQ^{−1/(2d)}Σ_nΛ_L(n)n^{−½+1/(2d)}
e(…)·[iκH/t + (i/2d)H′ + (it/2d)H″](t_n*): if the law holds it is O(1) (its weight integrates to O(1)), but the only
unconditional bound is the triangle inequality over the ≍ QT^d terms, Q^{½}T^{(d−1)/2}. So in degree ≥ 3 the leading-kernel
prime sum and the zero sum are not provably within o(T) of each other; either the correction sums are kept (b), or the law is
assumed for them (Corollary G2.2).

*Computed (F; illustration, not a test of the theorem's range: heights 20–400).* Hann [a, b] ∗ g₁, zeros against the
leading-kernel prime sum and against the first-order bracket:

| L | d | window | zeros | leading kernel | \|diff\| | first order | \|diff\| | ∫H/2π |
|---|---|---|---|---|---|---|---|---|
| X5 | 1 | [200,400] | −15.917 − 10.384i | −16.137 − 10.646i | 0.343 | −15.843 − 10.392i | 0.075 | 15.9 |
| Δ | 2 | [100,300] | −36.284 − 1.512i | −35.491 − 7.110i | 5.654 | −36.684 − 1.579i | 0.405 | 15.9 |
| Sym² 11a | 3 | [20,60] | −6.636 + 2.933i | −6.453 + 2.320i | 0.639 | −6.629 + 2.783i | 0.150 | 3.18 |

**Corollary G2.2 (zero law ⟺ prime law in every degree, smooth windows; proved, with the expansion bookkeeping in outline).**
Fix [a, b] ⊂ (0, ∞), 𝒱 := C_c^∞((a, b)); for v ∈ 𝒱 put v_T(t) := v(t/T), Z_T(v) := Σ_ρ h(ρ)X(ρ)^{−1} for H = v_T ∗ g₁, and
P_T(v) := Σ_nΛ_L(n)n^{−½+1/(2d)}e(−d(n/Q)^{1/d})(n/Q)^{ib/d}v_T(t_n*). For κ ∈ ℂ, as T → ∞:

    [∀v ∈ 𝒱: Z_T(v) = −(κ/2π)T∫v + o(T)]   ⟺   [∀v ∈ 𝒱: P_T(v) = κ′T∫v + o(T)],   κ′ = (√d/2π)e^{−i(φ₀+π/4)}Q^{1/(2d)}κ.

*Proof.* Expansion lemma (outline): collecting in Lemma 5's BKY expansion the asymptotic series of θ^{(i)}(t) in powers of 1/t
(Lemma 3), of t_n − t_n*, φ_n(t_n) and φ_n″(t_n), and the heat expansion of v_T ∗ g₁ (moments of g₁), one gets for t_n* ∈ [aT/2, 2bT]
I_n = 2πA_d(n/Q)^{1/(2d)}e(−d(n/Q)^{1/d})(n/Q)^{ib/d}Σ_{j<K}T^{−j}(𝒟_jv)(t_n*/T) + O(T^{½−K}), with 𝒟₀ = id and 𝒟_j linear maps
𝒱 → 𝒱 independent of T (𝒟_jv = Σ_{a≤2j}ψ_{j,a}v^{(a)}, ψ_{j,a} Laurent monomials in τ). With G1 this gives
Z_T(v) = −A_dQ^{−1/(2d)}[P_T(v) + Σ_{1≤j<K}T^{−j}P_T(𝒟_jv)] + O(T^{(d+1)/2−K}).  (∗)
(⇐) P_T(𝒟_jv) = O(T), so the j ≥ 1 terms are O(T^{1−j}) = o(T); take K > (d−1)/2. (⇒) Put E_T(v) := P_T(v) − κ′T∫v; (∗) and
the zero law give E_T(v) = −Σ_{j≥1}T^{−j}E_T(𝒟_jv) + O(1) + o(T) (the main terms of the corrections are κ′T^{1−j}∫𝒟_jv = O(1)).
Iterate L times (finitely many functions, the zero law applied to each) and use the trivial bound |E_T(w)| ≪_w Q^{1/(2d)}T^{(d+1)/2}
on the last generation: E_T(v) = o(T) + O(T^{(d+1)/2−L}) = o(T) for L > (d−1)/2. ∎ The sharp prime law
Σ_{n≤x}Λ_L(n)e(−d(n/Q)^{1/d}) = κ e^{−i(φ₀+π/4)}(2√d/(d+1))(x^{d+1}/Q)^{1/(2d)} + o(x^{(d+1)/(2d)}) implies the smoothed one by
partial summation (x^{(1−d)/d}dx = (d/2π)Q^{1/d}dt); the converse needs a Tauberian input (GAP 3, as in degree 2).

## 5. Theorem G3: the exact dual side in degree d

**Theorem G3 (proved; π ∈ 𝔉_d, π ≠ 1; no (H∞) needed).** For W, Δ, H as in §1 and r > 0,

    Σ_{n≥1} λ(n) n^{−½} I_n(r) = ε Σ_{m≥1} conj(λ(m)) m^{−½} Ĥ(ln(r/m)),     both sides absolutely convergent.

*Proof.* Consider 𝒥 := (1/2πi)∫_{(5/4)}h(s)X(s)^{−1}L(s)r^{½−s}ds. Expanding L and shifting each n-term to σ = ½ (holomorphic
on ½ ≤ σ ≤ 5/4: the poles of X^{−1} have real part ≤ −Re μ_j < ½ by F5) gives 𝒥 = (1/2π)·(left side), absolutely as in G1(4).
On the other hand X^{−1}L = εL̃(1−s) is entire; |L̃(w)| ≤ exp(C|w|^{1+η}) in vertical strips (Λ̃ of order 1, 1/Γ of order 1),
which the Gaussian decay of h dominates, so the line moves to σ = −¼, where L̃(1−s) = Σconj λ(m)m^{s−1} converges absolutely
(F1); shifting each m-term (entire integrand h(s)m^{s−1}r^{½−s}) back to σ = ½ gives (ε/2π)Σconj λ(m)m^{−½}∫H(t)(m/r)^{it}dt. ∎

**Theorem G3′ (log-weighted dual side; proved; needs m > 0).** For a > 0:
Σ_b λ(b) ln b · b^{−½} I_b(a) = ε Σ_m conj(λ(m)) m^{−½} ∫H(t)(2θ′(t) − ln m) e^{−it ln(a/m)} dt.
*Proof.* As G3 with −L′(s) = Σλ(b)ln b·b^{−s} and X^{−1}(−L′) = ε[−(X′/X)(s)L̃(1−s) + L̃′(1−s)] (differentiate L = εXL̃(1−s));
X′/X is holomorphic on −m < σ < 1 + m, so the lines 1 + η and −η (η < m) are admissible; on σ = ½, −X′/X = 2θ′. ∎
The resonant (m = a) term is ε conj λ(a)a^{−½}∫H(t)(2θ′(t) − ln a)dt.

**Corollary G3.1 (integer resonance; proved for smooth windows, GRC for the coefficient sums).** For r = q ∈ ℕ the term m = q
of G3 is ε conj λ(q)q^{−½}∫H; for smooth windows Ĥ(ln(q/m)) ≪_A T^{1−A} for m ≠ q. With Lemma 5 for λ in place of Λ_L:
2πA_dQ^{−1/(2d)}q^{1/(2d)}Σ_nλ(n)n^{−½+1/(2d)}e(−d(nq/Q)^{1/d})(nq/Q)^{ib/d}H(t_{nq}*) = ε conj λ(q)q^{−½}∫H + R_λ, with
R_λ ≪ Q^{½}T^{(d−1)/2}(log T)^{d−1} under GRC (Σ_{n≤x}|λ(n)| ≤ Σd_d(n)), o(T) for d ≤ 2 and, with the K-term expansion, in
every degree. In the normalisation of the law: "constant 1" at r = 1 and ε conj λ(q)q^{−½} at r = q.
This is Kaczorowski–Perelli's standard twist F(s, α) = Σa(n)n^{−s}e(−αn^{1/d}) (arXiv:1911.10497, opened, §1.2: spectrum
{d(m/q_F)^{1/d}}, n_α = q_Fd^{−d}α^d, first pole at s₀* = (d+1)/(2d) − iθ_F, residue ∝ ω_F·conj a(n_α); our α = d(r/Q)^{1/d}
gives n_α = r, and θ_F = b/d is our (n/Q)^{ib/d}), and Ren–Ye's Corollary 1.1 (level 1, opened).

*Computed (D).* Sym² 11a (ε = 1), window [36,44], Δ = 4, all n ≤ 8·10⁵:

| case | left side | ε·dual side | \|diff\| | with −ε |
|---|---|---|---|---|
| r = 1 | 7.9134259406 + 0.0559431377i | same | **1.45e−13** | 15.83 |
| r = 2 | 6.3949419681 + 0.3887276212i | same | **1.32e−13** | 12.81 |
| r = 3/2 | −0.5989979912 − 1.5954843061i | same | **1.38e−13** | 3.41 |
| G3′, a = 2 | 60.9101152496 + 4.0621581050i | same | **1.33e−12** | (resonant term 54.55) |
| G3′, a = 3/2 | −5.9331348696 − 12.9867872061i | same | **1.35e−12** | |
| χ mod 5, r = 2 (ε = 0.851 + 0.526i) | 13.7287685621 − 24.0490653811i | same | **4.1e−12** | no conj: 54.3; ε̄: 29.1 |
| χ mod 5, r = 1 | 34.1620099560 + 21.2824922846i | same | **4.4e−12** | no conj: 0.58; ε̄: 42.3 |

## 6. The singular series in degree d

### 6.1 Local algebra (proved)

At a prime p let P_p(X) = ∏_j(1 − α_jX) = Σ_kμ_L(p^k)X^k (so μ_L(p^k) = (−1)^ke_k(α)), λ(p^k) = h_k(α), and
F_p(X) := Σ_k conj λ(p^k)μ_L(p^k)X^k (general_law.md's F_p; c_L = ∏_pF_p(1/p)/(1−1/p) as a conditionally convergent product,
= B_L(1)/L(1, Ad π)).

**Lemma 6 (proved).** (i) If conj α_j = α_j^{−1} for all j (p unramified, GRC at p), then F_p(1) = 0, so
G_p(X) := F_p(X)/(1 − X) is a polynomial of degree d − 1 and the local factor of c_L is exactly c_{L,p} = G_p(1/p).
(ii) G_p(X) = 1 − (|λ(p)|² − 1)X + O(X²), and |λ(p)|² − 1 = tr Ad(diag α) =: λ_{Ad}(p).
(iii) d = 2: G_p(X) = 1 − (|λ(p)|² − 1)X (= 1 − λ(p²)X for trivial central character);
d = 3: **G_p(X) = 1 − (|λ(p)|² − 1)X + ω(p)conj λ(p³)X²**, ω(p) = α₁α₂α₃ (= 1 − 2|a|² + ωā³ with a = λ(p));
in general the top coefficient is (−1)^{d+1}ω(p)conj λ(p^d).
(iv) R_{p,1}(X) := P_p(X)Σ_{j≥0}λ(p^{1+j})X^j = (1 − P_p(X))/X, i.e. its coefficients are r_{p,1,k} = −μ_L(p^{k+1}); and
R_{p,a}(X) := P_p(X)Σ_jλ(p^{a+j})X^j is a polynomial of degree ≤ d − 1 for every a ≥ 0 (no unitarity needed).
*Proof.* (i) Σ_k(−1)^ke_k(α)h_k(α^{−1}) is the constant term of ∏_j(1 − α_jz)·∏_j(1 − α_j^{−1}z^{−1})^{−1} = ∏_j(−α_jz) = (−1)^de_dz^d,
whose constant term is 0 for d ≥ 1. (ii) coefficient of X in F_p is −|λ(p)|², and F_p = (1−X)G_p. (iii) For d = 3 write
e₁ = a, e₂ = ωā, e₃ = ω (unitarity: conj e_k(α) = e_{d−k}(α)/e_d(α)); F_p = 1 − |a|²X + (ωā³ − |a|²)X² + (2|a|² − ωā³ − 1)X³
and h₃(ᾱ) = ω̄(1 − 2|a|² + ωā³). (iv) P_pΣ_jλ(p^{a+j})X^j = (1 − P_p·Σ_{j<a}λ(p^j)X^j)/X^a. ∎
*Computed (E1):* (i) |F_p(1)| ≤ 5·10⁻¹⁴, (iv) 7·10⁻¹⁵, (iii) 2·10⁻¹⁵, for 200 random unitary Satake sets in each d = 1..5;
S11a: F₃(1/3) = 566/729 = 0.776406 (constants file 0.776406035665).

So, under GRC at the unramified primes, **c_L = ∏_p G_p(1/p) = ∏_p(1 − λ_{Ad}(p)/p + O_d(p^{−2}))**: the first-order shape
asked for is right, and the exact factor is the degree-(d−1) polynomial G_p. For d = 3:
c_L = ∏_p[1 − (|λ(p)|² − 1)/p + ω(p)conj λ(p³)/p²] (unramified p; exact F_p(1/p)/(1−1/p) at ramified p).

### 6.2 Route A: Hardy–Littlewood major arcs (one step proved, one step heuristic)

**Step A1 (heuristic).** On primes Λ_L(p) = λ(p)ln p; prime powers contribute O(T^{½}) to the smoothed prime sum. Model
Λ(n) by its major-arc expansion Σ_{q≤R}(μ(q)/φ(q))c_q(n) (Hardy's Ramanujan expansion of Λ(n)φ(n)/n; classical, not opened),
so Σ_nΛ_L(n)n^{−½}I_n "≈" Σ_q(μ(q)/φ(q))T_q, T_q := Σ_nc_q(n)λ(n)n^{−½}I_n.

**Step A2 (theorem; proved, all d).** For every q ≥ 1, with c_q(n) = Σ_{e|(q,n)}eμ(q/e) and λ(ek) = Σ_{fj=k}r_e(f)λ(j),
r_e the coefficients of 𝓡_e(s) := ∏_{p^a‖e}R_{p,a}(p^{−s}) (Lemma 6(iv)):

    T_q = ε Σ_{e|q} e μ(q/e) Σ_f r_e(f)(ef)^{−½} Σ_m conj λ(m) m^{−½} Ĥ(ln(ef/m))      (G3 with r = ef, finitely many terms).

The resonant terms (m = ef) sum to ε∫H·Σ_{e|q}eμ(q/e)𝔊(e), 𝔊 multiplicative, 𝔊(p) = Σ_k r_{p,1,k}conj λ(p^{1+k})p^{−1−k} =
1 − F_p(1/p) (Lemma 6(iv)); for squarefree q:

    resonant part of T_q = ε∫H · ∏_{p|q}(p𝔊(p) − 1) = ε∫H · φ(q) ∏_{p|q}(1 − c_{L,p}),   c_{L,p} := F_p(1/p)/(1 − 1/p).

So the modulus q resonates at the dual indices ef (e | q, f ∈ supp r_e: for q = p, e = p, f = p^k, k < d: m ∈ {p, …, p^d}; e = 1:
m = 1) with total weight φ(q)∏_{p|q}(1 − c_{L,p}). For d = 2, level 1: weight λ(p²)(p−1)/p at q = p, the family note's §9a
("dual index q², weight (φ(q)/q)λ(q²)"). *Computed (D):* the exact identity for Sym² 11a, q = 2, 3, 5, 6, to 3.5·10⁻¹⁴–2.5·10⁻¹³
(at Δ = 4 the non-resonant dual terms are not small: the resonance needs Δ·ln(1 + 1/m) ≫ 1).

**Step A3 (heuristic: the sum over q, i.e. the minor arcs).** Σ_q(μ(q)/φ(q))·φ(q)∏_{p|q}(1 − c_{L,p}) = ∏_p c_{L,p} = c_L
(Euler product in increasing p, conditionally convergent). Hence the prime sum "=" εc_L∫H and, by G1, Σ_ρ "=" −εc_L∫H/2π:
**the singular series of the GL(d) resonance over primes is c_L = B_L(1)/L(1, Ad π)**, in every degree. Unjustified: that the
major arcs give the whole main term (Steps A1, A3).

*Remark (Voronoi picture; sketch, sources not opened).* For (q, Q) = 1 and a primitive character χ mod q,
X_{π⊗χ}(s) = q^{d(½−s)}X_π(s)·(a factor of bounded phase variation), so the twisted sum resonates at dual index q^dr, where
conj(λχ)(q^dr) = 0: primitive characters never resonate, and the whole resonance of c_q-twists comes from the divisibility
pieces used in A2. In the additive-twist (Voronoi) form, GL(2): dual index q² (family §9a); GL(3) (Miller–Schmid; not
opened): the resonance condition n₁²n₂ = q³ selects (n₁, n₂) = (q, q) and (1, q³), whose coefficients A(q, q) = |λ(q)|² − 1
and A(q³, 1) ∝ conj λ(q³) are exactly the two non-trivial coefficients of G_p in Lemma 6(iii). Not used anywhere.

### 6.3 Route B: Λ_L = (λ·ln) ∗ μ_L and G3′ per divisor (exact for each divisor; the sum over divisors heuristic)

1. (exact) −L′/L = (−L′)·(1/L): Λ_L(n) = Σ_{ab=n}μ_L(a)λ(b)ln b.
2. (exact) Σ_nΛ_L(n)n^{−½}I_n = Σ_aμ_L(a)a^{−½}Σ_bλ(b)ln b·b^{−½}I_b(a) — absolutely convergent (via the σ₁-line representation
   n^{−½}I_n = J_n(σ₁): Σ_a|μ_L(a)|a^{−σ₁}·Σ_b|λ(b)|ln b·b^{−σ₁} < ∞ for σ₁ > 1, using F5 for |μ_L(p^k)| ≪ p^{kθ_d}).
3. (theorem for each fixed a; G3′) the inner sum is εΣ_m conj λ(m)m^{−½}∫H(2θ′ − ln m)e^{−it ln(a/m)}; resonant term m = a.
4. **(GAP 1: the unjustified step)** Sum the resonant terms over all a and drop the off-diagonal m ≠ a: for a ≍ x the m-sums near
   m ≈ a(1 ± 1/T) are not small; this is the type II range, open in every degree ≥ 2.
5. Granting 4: Σ_aμ_L(a)conj λ(a)a^{−1}(2θ′(t) − ln a) = 2θ′(t)F_L(1) + F_L′(1) with F_L(s) = Σconj λ(a)μ_L(a)a^{−s} =
   B_L(s)/L(s, π×π̃) (values of the analytic continuation; the series converge conditionally at s = 1, not proved here: GAP 9);
   F_L(1) = 0 (F6) and F_L′(1) = c_L. Hence Σ_ρ "=" −(εc_L/2π)∫H.
6. **Lines** r = u/v (lowest terms): only a = va′ give an integer dual index ua′; the amplitude is ε(uv)^{−½}F′_{u,v}(1),
   F_{u,v}(s) = Σ_{a′}μ_L(va′)conj λ(ua′)a′^{−s}, whose p-factor (r = p^v) is p^{|v|/2}S_p(v) with
   S_p(v) = Σ_{m−h=v}conj λ(p^m)μ_L(p^h)p^{−(m+h)/2}: divided by F_p(1/p) this is **general_law.md's ρ_L(p^v) = p^{|v|/2}S_p(v)/F_p(1/p)**.

**Agreement.** Route A's line version (the same computation with I_n(r)) has p-factor (p/(p−1))p^{−v/2}Σ_hμ_L(p^h)conj λ(p^{h+v})p^{−h},
which equals route B's/the recipe's S_p(v)/(1 − 1/p) identically (E1c: 5·10⁻¹⁶ for v = −3..3, d = 1..5). So both heuristics
give the constant c_L **and** the line factors ρ_L(r) of general_law.md in every degree, with conj λ in the numerator slot and
μ_L of L itself; the root number enters as ε (G3), the conjugate of λ through the dual side.

### 6.4 Numerical comparison (E2)

| L | c_L (constants file) | ∏_{p≤10³} G_p(1/p) | 10⁴ | 10⁵ | 3·10⁵ | 10⁶ | Σ_{q≤2·10⁵} μ(q)∏_{p\|q}(1−G_p) |
|---|---|---|---|---|---|---|---|
| Sym² 11a (d = 3) | 1.94470 | 1.93614 | 1.94395 | 1.94499 | 1.94338 | **1.94469** | 1.94170 |
| Sym² 15a (d = 3) | 3.27678 | 3.24125 | 3.27073 | 3.27510 | 3.27650 | **3.27655** | 3.27157 |
| Δ (d = 2) | 2.36082 | 2.34798 | 2.35773 | 2.36017 | **2.36059** | – | 2.35988 |
| χ mod 5 (d = 1) | 1.25 = 5/4 | 1.25000 | 1.25000 | 1.25000 | 1.25000 | 1.25000 | 1.25000 |

The products converge (slowly, like the prime number theorem for Ad π) to the constants computed as B(1)/L(1, Ad) with PARI
L-values; the q-sums converge more slowly. For d = 1, G_p = 1 at p ∤ q and G_q(1/q) = q/φ(q) (q prime): c = q/φ(q).

## 7. What is provable about the main term in degree d (summary)

- Exact: G1 (all zeros ↔ all prime powers, kernel X^{−1}(½+it)n^{−it}), G3, G3′, the T_q identity.
- Asymptotic: G2 — the closed-form resonance e(−d(n/Q)^{1/d}) represents the zeros within o(T) only for d ≤ 2; in every degree
  with ⌈d/2⌉ terms of the expansion; zero law ⟺ prime law (Corollary G2.2) for smooth windows, every degree.
- Over all integers the resonance has constant 1 (G3.1; known). Over primes the constant is c_L only heuristically (GAP 1).
- Sharp windows: not treated (Arias de Reyna's method gives errors T^{d/2}, already ≍ T at d = 2: degree-2 note §4).

## 8. Novelty (opened sources marked; mathematical queries only)

- **G1.** Weil's explicit formula for GL(d) with the test function g = hX^{−1}r^{½−s} (form (b)): **KNOWN IN SUBSTANCE**
  (Weil; for GL(m) e.g. Rudnick–Sarnak 1996, Iwaniec–Kowalski Thm 5.12: **not opened**). The written statement — all zeros, the
  phase X(ρ)^{−1}, absolute convergence on the critical line, Λ_L without conjugation and no root number — for d ≥ 2: **NOT FOUND**
  (degree 2 in fact 2052 of this project).
- **G2.** Per term: the stationary phase of ∫H X^{−1}(½+it)n^{−it}dt is the asymptotics of the GL(d) Voronoi (Meijer-G) kernel:
  **KNOWN IN SUBSTANCE** (Ren–Ye, arXiv:1411.6022, **opened**, Theorem 1.1: full asymptotic expansion, c₀ = −1/√m, level 1;
  Kaczorowski–Perelli's standard twist, poles s_ℓ = (d+1)/(2d) − ℓ/d, arXiv:1911.10497, **opened**). Degree 1: Gonek's lemma,
  Arias de Reyna 2014, Durkan–Hughes–Pearce-Crump arXiv:2601.18025 (abstract **opened**: ζ only). The summed statement over zeros
  with A_d = e^{i(φ₀+π/4)}/√d, the first-order constant κ = (d²−1)/(24d) + ½ΣRe(μ_j−μ_j²), the threshold d ≤ 2 for the closed-form
  kernel and the every-degree equivalence: **NOT FOUND**.
- **G3, G3′, Corollary G3.1.** The functional equation for this test function (a smoothed GL(d) Voronoi formula at a
  multiplicative twist): **KNOWN**; the integer resonance at α = d(q/Q)^{1/d} with weight ω·conj a(q): **KNOWN** (Kaczorowski–Perelli
  standard twist, Selberg class; Ren–Ye Corollary 1.1 for GL(m), level 1).
- **Singular series.** The identification of c_L = B_L(1)/L(1, Ad π) with the Hardy–Littlewood singular series of
  Σ_nΛ_L(n)e(−d(n/Q)^{1/d}), the exact local factor G_p = F_p/(1 − X) of degree d − 1 (d = 3: 1 − (|λ(p)|²−1)X + ω(p)conj λ(p³)X²),
  and the agreement of the major-arc and type-I routes with the ratios recipe's lines: **NOT FOUND** (web searches with
  mathematical terms only; MathSciNet not available). The degree-2 analogue (Iwaniec–Luo–Sarnak/Zhao heuristic) has no explicit
  constant in the sources read in the degree-2 round.

## 9. GAPs and doubts

GAPs:
1. The prime-side asymptotic with constant c_L (sum over moduli / type II) — open in every degree ≥ 2; both heuristic routes
   rest on it.
2. Corollary G2.2: the expansion lemma (structure 𝒟_jv = Σψ_{j,a}v^{(a)} with T-independent coefficients) is written in outline,
   not term by term.
3. Smoothed laws ⇒ sharp prime asymptotic: Tauberian input missing (as in degree 2).
4. (H∞) excludes archimedean parameters with Re μ_j < 0; then G1 needs Λ(−μ_j − 2k) ≠ 0 and the trivial zeros inside the strip
   as extra residues (not written).
5. G2 for box windows uses GRC (Brun–Titchmarsh on intervals of length ≍ n/t); the unconditional short-interval ℓ¹ bound (F2)
   holds only for x ≥ (C(π)T)^{144d³}.
6. Uniformity in Q: constants depend on π except under GRC (smooth windows); the E₀ constant depends on μ.
7. No degree-4 numerical check (genus-2 coefficients to ≈ 2·10⁶ would take ≫ 10 min with PARI on this machine).
8. The Voronoi remark of §6.2 (GL(3) dual indices) is from memory of Miller–Schmid, not opened; nothing depends on it.
9. Conditional convergence of Σμ_L(a)conj λ(a)a^{−1}(ln a)^k at s = 1 is not proved; Step 5 of route B uses the analytic
   continuation of F_L.

Doubts:
1. Facts F1–F4 are quoted from Kaneko–Thorner's summary of 𝔉_m (their Proposition "prop:class" cites Soundararajan–Thorner and
   IK, not opened); F5 from Lemke Oliver–Thorner's summary of LRS and Müller–Speh (not opened).
2. BKY Prop. 8.2 is used as stated (its proof not re-checked); the first-order bracket of Lemma 5 was derived here and checked
   numerically (B, F) only for constant and slowly varying windows.
3. The numerical G1/G3 checks use Gaussian-smoothed boxes at low heights (t ≤ 460) with Δ = 4–5: they test the exact identities
   (G1, G3, G3′, T_q), not the asymptotic regime of G2; G2 itself is checked per term (B) and illustrated (F).
4. Zero lists: completeness is confirmed only in the tested windows (a missing zero would show as an error of size ≈ |H(γ)| ≈ 1);
   zeros are float64 (≈10⁻¹³ absolute), which sets the 10⁻¹³–10⁻¹² agreement.
5. b = ΣIm μ_j vanishes under Kaneko–Thorner's normalisation of the central character for principal-series components (argued,
   not proved here); the G2 leading term is stated with b, the first-order bracket only for b = 0.
6. The genus-2 constants in `universal/` are direct products (not B/L(1, Ad)); they were not re-derived here.

## Referee report (25 November 2026)

An independent referee (own code in [referee/](referee/)) gave these verdicts. **Theorems as stated:** Lemmas 1–4, **G1** (exact
explicit formula, forms (a)–(c); unconditional, needs Re μ_j ≥ 0), Lemma 5 (one term: A_d, φ₀, κ; the constants confirmed for
d = 1…5, κ to 10⁻⁸–10⁻⁵), **G3** (no condition on μ), G3′, Corollary G3.1, **Lemma 6**, the per-q identity of Route A and its
resonant part. **Heuristic, labelled as such:** the sum over q and Route B step 4. Independent checks: the functional equation
with conjugations on complex characters mod 7 and on Sym² of 11a (10⁻¹⁶; the wrong conjugation fails by ≥ 0.7); G1 on new
windows and twists (10⁻¹²–10⁻¹³; conj Λ on the resonant side fails by 9.9–18), including a low window where E₀ = 0.07 and 0.35
so that its Weil form is really tested (10⁻¹¹); G3 and G3′ on Sym² of 11a (10⁻¹²–10⁻¹³; flipped ε fails by 8–10). Corrections,
which override the text above where they differ:

1. **G2, error E₀.** The proof of G2 quotes "E₀ ≪ Q^{−3/4}", which G1(c) gives only with the factor T^{1−3d/4}·log and only when
   min Re μ_j > ¼. Fix: move the left line to σ₀ = −K and collect the trivial-zero residues (each Gaussian-small in the
   distance of −Im μ_j from the window): E₀ ≪ T^{1−d(K+½)}·log T + Gaussian terms. No o(T) conclusion changes.
2. **G2(c), wording.** "o(T) only for d ≤ 2; d = 3: O(T)" is a statement about the provable (triangle-inequality) bound on R,
   not about R. Measured for Sym² of 11a on windows [20,60], [40,120], [80,240]: the sum of |leading terms| is 386, 1540, 6164,
   the triangle bound of the correction 59, 121, 242, the actual |zeros − leading kernel| 0.73, 0.46, 0.16.
3. **Lemma 6 needs less.** F_p(1) = 0 needs only unitarity at an unramified p, not Ramanujan; G_p has degree ≤ d − 1; the
   parenthetical for d = 3 should read ω(p)·conj λ(p³) = 1 − 2|a|² + ω·ā³. b = Σ Im μ_j = 0 for every π in the class (central
   character trivial on ℝ_{>0}), so the factor (n/Q)^{ib/d} is 1 and the first-order bracket holds throughout.
4. **Corollary G2.2** is a theorem once the expansion lemma is written out (routine; not done). Both directions are correct.
5. Smaller: in G1(a), when Re μ_j = ¼ take a line in (−½, 0) that avoids the poles; Route B step 2 should take σ₁ = 2;
   Corollary G3.1 does not need Ramanujan (Rankin–Selberg mean square). In check_degd.py the thresholds in sections C–E are
   looser than the agreement reached, and no window has E₀ large enough to test its explicit form (the referee's r3 does).

## 10. What a referee should check first

1. Lemma 1–2: the definition of X with μ̄_j in the numerator, θ with the −((t+Im μ_j)/2)ln π terms, the strip (−m, 1+m), and that
   the summand at an off-line zero is X(ρ)^{−1}.
2. G1 step (1)/(5): that under (H∞) the only poles in the rectangle are the zeros of Λ and points of 𝒫 (none right of σ₀ in
   form (b)); the conjugation in E₀ (Λ_L̃ = conj Λ_L) and the absence of ε.
3. Lemma 4(iii): the good-height argument from the opened zero count (replaces the IK citation of the degree-2 note).
4. Lemma 5 and Theorem G2(ii): BKY parameters (Y = Q_B = T, V = Δ or T) are degree-independent; the sum over n of the per-term
   errors gives Q^{½}T^{(d−1)/2} — the source of the d ≤ 2 threshold; the constant A_d and κ (check B: 0.8889, 0.34375, 30.18).
5. Corollary G2.2: the induction (⇒ direction) and the expansion lemma.
6. G3/G3′: the entire-function shift with order-1 growth (no good heights needed), the conj λ and ε on the dual side.
7. Lemma 6(i) (F_p(1) = 0 under unitarity) and the T_q identity's resonant part ∏(p𝔊(p) − 1) = φ(q)∏(1 − c_{L,p}).

## 11. Files

| file | content |
|---|---|
| `degd_proof.md` | this note |
| `degd_lib.py` | θ_L, θ_L′, windows H and Ĥ, the oscillatory sums I_n (trapezoid, recurrence), prime powers, local algebra |
| `check_degd.py` | sections A (phase, FE, growth, Stirling), B (one-term stationary phase, A_d, κ), C (G1 exact on zeros), D (G3, G3′, T_q), E (Lemma 6, Euler products vs constants), F (G2 illustration) |
| `check_degd_out.txt` | output: 44/44 PASS, 35 s |
| `explore.py`, `explore2.py` | decay of I_n in ln n (choice of N) |
| `dl_src/2203.09515`, `dl_src/1505.03122`, `dl_src/1911.10497` | LaTeX of the opened sources |
