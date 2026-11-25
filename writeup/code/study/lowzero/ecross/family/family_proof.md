# The one-point sum Σ e^{2iθ_f(γ)} on average over a Petersson family: theorems, the constant, and one form

Folder: `scratchpad/family/` (nothing in the repository was changed). Code: `famlib.py` (Kloosterman sums, Petersson series,
PARI family data, stationary-phase weights), `check_family.py` (sections a–d; outputs `check_family_out_{a,b,c,d}.txt`).
Builds on the refereed note `writeup/code/study/lowzero/ecross/proof/gl2_onepoint_proof.md` ("the note": Proposition 1,
Theorem 2, Theorem 4 there). Source opened: Iwaniec–Luo–Sarnak, arXiv:math/9901141 (DVI source, converted to text by
`dvitext.py`; `dl_ils/ils.txt`). Labels: **proved**, **proved (cited)** (a precisely cited result is used; "not opened" is said
when the source was not read), **sketch**, **GAP**, **conjecture**, **computed**.

## 0. Result in short

Family: N prime, k ∈ K := {2, 4, 6, 8, 10, 14} fixed (S_k(SL₂(ℤ)) = 0, so the Hecke basis H_k(N) of S_k(Γ₀(N)) consists of
newforms), harmonic weights ω_f, Σ^h X_f := Σ_{f∈H_k(N)} ω_f X_f. The phase θ(t) = (t/2)ln N − t ln 2π + Im lnΓ(k/2 + it) is
the same for every f. Z_f(H) := Σ_ρ m_ρ h(ρ)X(ρ)^{−1} (all zeros; = Σ_γ H(γ)e^{2iθ(γ)} under GRH for L(s,f)), H = W(·/T) ∗ g₁,
W ∈ C_c^∞((1,2)) fixed, ∫H = T∫W.

1. **Theorem A (exact decomposition; proved, unconditional).** Σ^h ε_f Z_f(H) = −(1/2π)∫H + 𝓔₁ + Σ_{m≥2}𝓚_m + O(T^{1/2}),
   where 𝓔₁ is the prime-number-theorem remainder of one explicit smooth function and 𝓚_m is a sum over primes of the
   Kloosterman sum S(Np, 1; Nm). The main term comes from the modulus c = N of the root-number-twisted Petersson formula
   and is exact: it rests on the identity ∫₀^∞ u^{−1/2−it}J_{k−1}(4π√(u/N))du = (√N/2π)e^{−2iθ(t)} (Bessel–Mellin; the J-Bessel
   kernel of Petersson is the inverse of the functional-equation factor), so no stationary phase is needed for the main term.
2. **Theorem B (GRH for ζ and Dirichlet L-functions; proved).** Error O(T^{1/2} + N^{−1/2}T^{5/2}log³(NT)): the average
   equals −(1/2π)∫H(1 + o(1)) for T → ∞, T ≤ N^{1/3}(log N)^{−3}. **Theorem C (unconditional; proved, with Page's theorem
   cited, not opened):** the same for T → ∞, T ≤ exp(c₀√log N).
3. **The constant.** The family constant is **M = 1**. The conjecture predicts Σ^h ε_f·(−ε_f c_f∫H/2π), i.e. M_pred = Σ^h c_f.
   **Proposition 5 (GRH for L(s, Sym²f), as assumed in ILS):** Σ^h c_f = 1 + O(N^{−1/4+ε}). **They match.** The match is
   structural: the harmonic measure makes λ_f(p) Sato–Tate distributed, and E_ST[1 − λ(p²)/p] = 1 for every p. Natural
   weights: the mean of c_f tends to 6/π², and the natural-weight family average of ε_f Z_f is −(6/π²)∫H/2π (sketch).
   Computed at small prime levels: Σ^h c_f scatters around 1 (k = 4, 10 levels 13–47: mean 1.06; k = 2, 9 levels 43–97: mean
   1.04), not a test of the limit.
4. **Theorem D (twisted moments and lines; GRH, proved as B):** Σ^h ε_f λ_f(ℓ) Z_f(H; r) = −(∫H/2π)·A(ℓ, r) + o(T) with
   A(ℓ, r) = μ(m)c_m(ℓ)/(φ(m)√ℓ) when r = ℓ/m², (m, N) = 1, and A = 0 for every other r (c_m = Ramanujan sum).
   **Proposition D′ (proved; sympy check a7):** A(ℓ, r) equals the Sato–Tate average of the conjecture's amplitude
   λ_f(ℓ)c_f ρ_f(r)/√(ab), prime by prime. So every Hecke-twisted first moment of the family, at every rational line,
   agrees with the conjecture: the conjecture's c_f and ρ_f(r) have exactly the right Hecke–Fourier coefficients.
5. **Second moment (C):** the diagonal is (1/2π)∫H²ln(Nt²/4π²); the off-diagonal is a bilinear form in S(p, p′; Nm) over pairs
   of primes with moduli up to ≈ NT², out of reach of GRH. "Almost all f" follows from one stated second-moment inequality
   (S2) (trivial proof); (S2) itself is a **GAP**.
6. **One form (D).** For Δ (S₁₂(1) = ℂΔ) Petersson gives λ_Δ(n) exactly as a Kloosterman series. Its c = 1 term carries the
   amplitude 1/ω_Δ = 0.3521; the conjecture says c_Δ = B_Δ/L(1, Sym²Δ) = 2.3608. **The zeros of L(Δ) decide: computed slope
   1.016 (sharp) and 0.989 (Hann) on [20, 400] against c_Δ; 1/ω_Δ is excluded** (−21 against −145). Numerically 85 % of the
   amplitude comes from moduli 5 ≤ c ≤ 60 at heights ≤ 400 (c growing with T). Within the Petersson decomposition (not as a statement about primes: the same
   split occurs for the sum over all integers, a theorem; see §9a) a bound beating the square-root scale uniformly in c ≤ T^{1+η}
   would contradict the conjecture; the constant is the Hardy–Littlewood singular series of the resonance (§9a). **Theorem E (proved):** the ratios conjecture for
   L(s, f) in the t-aspect (two statements R1, R2 below) implies the one-point law for one form, error O(T^{1/2+ε}).
7. **Δ data point (new, outside elliptic curves):** zeros of L(Δ) to 420 (431 zeros; complete: the exact explicit formula
   of Proposition 1, weight 12, holds to 2·10^{−12}); Σ_{20<γ≤400} e^{2iθ(γ)} = −145.03 − 1.20i against −c_Δ·380/2π = −142.78;
   Hann −70.61 − 0.03i against −71.39. Over primes alone, Σ_{n≤x}Λ_Δ(n)e(−2√n)/((2√2/3)e^{−23πi/4}x^{3/4}) = 2.30 + 0.08i at
   x = 2·10⁵ (c_Δ = 2.36; over all integers the same ratio is 0.998 + 0.03i, theorem value 1).

## 1. Setting

k ∈ K, N prime, H_k(N) the newforms (a_f(1) = 1), λ_f(n) = a_f(n)n^{−(k−1)/2}, L(s,f) = Σλ_f(n)n^{−s},
Λ(s,f) = N^{s/2}Γ_C(s + (k−1)/2)L(s,f) = ε_f Λ(1−s,f). Hecke: λ_f(m)λ_f(n) = Σ_{d|(m,n),(d,N)=1}λ_f(mn/d²) (ILS (113),
opened); λ_f(N)² = 1/N. **Root number:** ε_f = i^k μ(N) λ_f(N)√N = −i^k λ_f(N)√N (ILS (115), (118) for k = 2, opened;
general k: checked against PARI's root number for 9 forms with k = 2, 4, check b5).
X(s) := N^{½−s}Γ_C(1 − s + (k−1)/2)/Γ_C(s + (k−1)/2), L(s,f) = ε_f X(s)L(1−s,f), X(½+it)^{−1} = e^{2iθ(t)} with θ as above
(note, Lemma 1, with Γ(1+it) replaced by Γ(k/2+it); X^{−1} is holomorphic and zero-free in −(k−1)/2 < σ < (k+1)/2).
Λ_f(p^j) = (α_p^j + β_p^j)ln p (p ≠ N), λ_f(N)^j ln N; |Λ_f(n)| ≤ 2Λ(n).
**Harmonic weights:** ω_f := Γ(k−1)/((4π)^{k−1}⟨f,f⟩_{Γ₀(N)\H}); then (Lemma 1)
ω_f = 2π²/((k−1)·N·L(1, Sym²f)) with L_N(s, Sym²f) = (1 − N^{−1−s})^{−1} (for k = 2 this is ILS (124), (126), (127):
12/((N − N^{−1})Σλ_f(t²)/t), and Σλ_f(t²)t^{−s} = L(s,Sym²f)/ζ^{(N)}(2s); for general k computed: three independent
values agree to 10^{−11} for Δ (check b1), and the Petersson sum rule holds with these ω_f at (N,k) = (7,4), (13,4), (11,2),
(37,2), check b3). Σ^h 1 = Δ_{k,N}(1,1) = 1 + O(N^{−3/2+ε}).
**The conjecture (round 148, per form):** Z_f(H) = −ε_f c_f(1/2π)∫H + o(T) as T → ∞, c_f = B_f/L(1,Sym²f),
B_f = ∏_{p≠N}B_p·(1 − 1/N)^{−1}, B_p = (1 − e₁/p)/(1 − e₁/p + e₁/p² − 1/p³), e₁ = λ_f(p)² − 1; equivalently
c_f = "∏_{p≠N}(1 − λ_f(p²)/p)·(1 + 1/N)" (regularised) = G_f(1), G_f(s) := Σ_{(d,N)=1}μ(d)λ_f(d²)d^{−s}·(1 − N^{−1−s})/(1 − N^{−s}).
Lines: Z_f(H; r) := Σ_ρ h(ρ)X(ρ)^{−1}r^{½−ρ} = −ε_f c_f ρ_f(r)(ab)^{−½}(1/2π)∫H + o(T), r = a/b (derivation.md §2.1).
**Test functions:** H := W_T ∗ g₁, W_T(t) = W(t/T), W ∈ C_c^∞((1,2)) real, g₁(z) = e^{−z²}/√π; h(s) := H((s−½)/i). H is entire,
|h(σ+it)| ≤ e^{(σ−½)²}(|W_T| ∗ g₁)(t), H^{(j)} ≪_j T^{−j}, H(t) ≪ e^{−d(t)²/2} at distance d(t) from [T, 2T];
Ĥ(ξ) = ∫H e^{−iξt}dt = TŴ(Tξ)e^{−ξ²/4} ≪_A T^{1−A}|ξ|^{−A}. Put I_u := ∫H(t)e^{i(2θ(t) − t ln u)}dt (u > 0): **independent of f**.
Put t_u := 2π√(u/N), x₋ := N(T/4)²/4π², x₊ := N(4T)²/4π² (so t_u ∈ [T/4, 4T] ⟺ u ∈ [x₋, x₊]).

## 2. Lemmas

**Lemma 1 (Petersson; proved (cited)).** For m, n ≥ 1: Σ_{f∈H_k(N)} ω_f λ_f(m)λ_f(n) = δ(m,n) + 2πi^{−k}Σ_{c≡0 (N)} S(m,n;c)c^{−1}
J_{k−1}(4π√(mn)/c). Source: ILS (125)–(126) for k = 2 (opened); general k: Iwaniec–Kowalski Prop. 14.5 (not opened).
*Computed:* Δ (k = 12, N = 1): ω_Δλ_Δ(n) = Δ(1,n) for 14 values of n to 2.9·10^{−15} (b2); (7,4), (13,4): to 6·10^{−9};
(11,2), (37,2): to 3·10^{−4} (the weight-2 series converges like c^{−3/2}; truncation at c ≤ 1500N) (b3).

**Lemma 2 (root-number Petersson; proved).** For n ≥ 1:
    Σ^h ε_f λ_f(n) = −2πN^{−½} Σ_{m≥1} S(Nn, 1; Nm) m^{−1} J_{k−1}(4π√(n/N)/m),
absolutely convergent. *Proof.* ε_fλ_f(n) = −i^k√N λ_f(N)λ_f(n) = −i^k√N λ_f(Nn) (Hecke, the d = 1 term only since (d,N) = 1);
apply Lemma 1 with (Nn, 1): δ = 0, c = Nm, J_{k−1}(4π√(Nn)/(Nm)) = J_{k−1}(4π√(n/N)/m). Convergence: |S(a,b;c)| ≤
τ(c)(a,b,c)^{½}c^{½} (Weil; cited, not opened) and |J_ν(y)| ≤ (y/2)^ν/ν!: Σ_m τ(Nm)(Nm)^{½}m^{−1−(k−1)} < ∞. ∎
*Computed (b4):* (7,4), (13,4) to 7·10^{−8} (m ≤ 600; to 8·10^{−10} with m ≤ 3000), (11,2), (37,2) to 6·10^{−4} (m ≤ 1500).

**Lemma 3 (Kloosterman sums; proved).** (i) For (m,N) = 1: S(Nn, 1; Nm) = μ(N)S(n, N̄; m), NN̄ ≡ 1 (m) (ILS (157), opened;
a1). In particular S(Nn,1;N) = μ(N) = −1. (ii) For (b,m) = 1: Σ_{a∈(ℤ/m)^×} S(ℓa, b; m) = μ(m)c_m(ℓ). *Proof:*
Σ_a S(ℓa,b;m) = Σ_{x∈(ℤ/m)^×}e(bx̄/m)Σ_{a∈(ℤ/m)^×}e(ℓax/m) = c_m(ℓ)Σ_x e(bx̄/m) = c_m(ℓ)c_m(b) = μ(m)c_m(ℓ), since
c_m(ℓx) = c_m(ℓ) for (x,m) = 1. (a2.) (iii) For K(a) := S(Na,1;Nm), a ∈ (ℤ/m)^×, write K = Σ_χ ĉ_χ χ (χ mod m). Then
Σ_χ|ĉ_χ| ≤ m. *Proof:* Σ_χ|ĉ_χ|² = φ(m)^{−1}Σ_a|K(a)|² ≤ φ(m)^{−1}Σ_{a mod m}|S(a,N̄;m)|² = m (orthogonality:
Σ_a|S(a,b;m)|² = mφ(m)); Cauchy–Schwarz. (a3.) The principal coefficient is ĉ_{χ₀} = −μ(m)²/φ(m) (by (i), (ii)).

**Lemma 4 (Bessel–Mellin; proved (cited) + computed).** For ν ≥ 1, −ν < Re s < 3/2: ∫₀^∞ y^{s−1}J_ν(y)dy =
2^{s−1}Γ((ν+s)/2)/Γ((ν−s)/2 + 1) (Gradshteyn–Ryzhik 6.561.14, not opened; a4: to 2·10^{−21} at Re s = 0, 0.6, 1; absolutely
convergent only for Re s < ½). With y = a√u: ∫₀^∞u^{w−1}J_ν(a√u)du = (2/a)^{2w}Γ(ν/2 + w)/Γ(ν/2 + 1 − w). For ν = k−1,
a = 4π/(m√N), w = ½ − it:

    ∫₀^∞ u^{−½−it} J_{k−1}(4π√(u/N)/m) du = (√N m/2π)·(Nm²/4π²)^{−it}Γ(k/2 − it)/Γ(k/2 + it) = (√N m/2π) e^{−2iθ(t)} m^{−2it}   (★)

(a5: to 4·10^{−20}). On Re w = 0 the kernel has modulus exactly |a′ + iy|^{−1}, a′ = (k−1)/2 (a6).

**Lemma 5 (the f-free integral I_u; proved, given BKY Lemma 8.1 as in the note).** (i) |I_u| ≤ C T^{½}, |∂_u I_u| ≤ C T^{3/2}/u
for all u (second-derivative test, |∫_a^b G e^{iF}| ≤ 8(sup|G| + ∫|G′|)/√r when |F″| ≥ r, with F″ = 2θ″ ≥ c/T on [T/2, 3T],
plus the Gaussian tails of H; ∂_uI_u = −(i/u)∫tH(t)e^{iF}). (ii) For u ∉ [x₋, x₊]: on supp W_T ± 3√(log T) the phase derivative
2θ′(t) − ln u has modulus ≥ 1, so BKY Lemma 8.1 (X = 1, U = T, Q = Y = T, R = 1) gives |I_u| + u|∂_uI_u|/T ≪_A T^{−A};
for u ≥ e^{10}x₊ or u ≤ e^{−10}x₋ the contour shift I_u = −i∫_{(½±σ)}h(s)X(s)^{−1}u^{½−s}ds (0 < σ < k/2 on the left, any σ on
the right; |X^{−1}(½+σ+it)| ≍ (N(t² + k²))^{σ}) adds the factor (x₊/u)^{σ}, resp. (u/x₋)^{σ}. ∎

**Lemma 6 (the smooth functions Φ_m; proved).** Φ_m(u) := u^{−½}J_{k−1}(4π√(u/N)/m)I_u. With μ_m := min(1, (m/T)^{½}, (T/m)^{k−1}):
|Φ_m(u)| ≪ u^{−½}T^{½}μ_m, |Φ_m′(u)| ≪ u^{−3/2}T^{3/2}μ_m on [x₋, x₊] (|J_ν(y)| ≪ min(1, y^{−½}, (y/2)^ν/ν!), the same for
J_ν′, and the J-factor's relative derivative ≪ (t_u/m + 1)/u ≪ T/u), rapid decay outside (Lemma 5(ii)). Hence
∫₀^∞√u|Φ_m′|du ≪ T^{3/2}μ_m, ∫₀^∞u|Φ_m′|du ≪ √N T^{5/2}μ_m, Σ_p ln p|Φ_m(p)| ≪ √N T^{3/2}μ_m (Chebyshev).

## 3. Theorem A (exact decomposition; proved, unconditional)

For N prime, k ∈ K, T ≥ 3:

    Σ^h ε_f Z_f(H) = −(1/2π)∫H + 𝓔₁ + Σ_{m≥2} 𝓚_m + 𝓟,      |𝓟| ≪ T^{½},
    𝓔₁ := −N^{−½}[Σ_n Λ(n)Φ₁(n) − ∫₀^∞Φ₁(u)du],     𝓚_m := N^{−½}m^{−1}Σ_p ln p·S(Np, 1; Nm)·Φ_m(p).

*Proof.* (1) Note, Proposition 1 (refereed: unconditional theorem; its proof only uses the functional equation, which holds
for weight k with the X above): Z_f(H) = −(1/2π)Σ_n Λ_f(n)n^{−½}I_n + E₀(f), |E₀(f)| ≪ N^{−3/4}log(NT)T^{−½}, both uniform in
f. Multiply by ε_fω_f and sum over the finitely many f: Σ^hε_fZ_f = −(1/2π)Σ_n c̄(n)n^{−½}I_n + Ē₀, c̄(n) := Σ^hε_fΛ_f(n).
(2) Prime powers n = p^j, j ≥ 2: |c̄(n)| ≤ 2 ln p·Σ^h1 ≤ 4 ln p; by Lemma 5 and Mertens Σ_{p²∈[x₋,x₊]}ln p/p ≪ 1, so this part
is ≪ T^{½} (outside [x₋, x₊] negligible): no Petersson needed here. (Through Petersson the squares would resonate too, at
m = 1, but with total weight ≪ N^{−½}log T.) Ē₀ ≪ N^{−3/4}: both go into 𝓟.
(3) Primes: c̄(p) = ln p·Σ^hε_fλ_f(p) = −2πN^{−½}ln p Σ_m S(Np,1;Nm)m^{−1}J_{k−1}(4π√(p/N)/m) (Lemma 2), so the prime part is
Σ_m N^{−½}m^{−1}Σ_p ln p S(Np,1;Nm)Φ_m(p) (the double series converges absolutely by Weil, |J_ν(y)| ≤ (y/2)^ν/ν! and
Lemma 5). (4) m = 1: S(Np,1;N) = −1 for every p (also p = N), so the term is −N^{−½}Σ_p ln pΦ₁(p) = −N^{−½}∫₀^∞Φ₁ + 𝓔₁ +
N^{−½}Σ_{p^j,j≥2}ln pΦ₁(p^j), the last ≪ N^{−½}log(NT) (into 𝓟). (5) **∫₀^∞Φ₁ = (√N/2π)∫H.** For U < ∞, Fubini on
[0, U] × ℝ gives ∫₀^UΦ₁ = ∫H(t)e^{2iθ(t)}[∫₀^U u^{−½−it}J_{k−1}(4π√(u/N))du]dt. The inner integral tends to (★) as U → ∞
and is bounded by C(1+|t|) uniformly in U (split at y = 1 + 4|t| in the y-variable; beyond it the phase y ∓ 2t ln y of
J's asymptotic has derivative ≥ ½: first-derivative test). Dominated convergence (H has Gaussian tails) and
∫₀^∞|Φ₁| < ∞ (Lemma 5) give ∫₀^∞Φ₁ = ∫H(t)e^{2iθ(t)}(√N/2π)e^{−2iθ(t)}dt = (√N/2π)∫H. Hence the m = 1 term is
−(1/2π)∫H + 𝓔₁ + O(N^{−½}log NT). ∎

*Remark.* For each m the same argument gives ∫₀^∞Φ_m = (√N m/2π)Ĥ(2 ln m) (★): for m ≥ 2 this is ≪_A √N m T^{1−A}(ln m)^{−A}:
**only m = 1 is resonant**; every m ≥ 2 contributes only through the fluctuation of primes in progressions mod m.

## 4. Theorem B (GRH; proved) and Theorem C (unconditional; proved with a cited theorem)

**Theorem B.** Assume RH for ζ and GRH for L(s,χ), χ mod m, m ≤ M₀ (M₀ = T⁵ for k = 2, T³ for k ≥ 4). For 3 ≤ T ≤ N^{1/3}:

    Σ^h ε_f Z_f(H) = −(T/2π)∫W + O_{W,k}( T^{½} + N^{−½}T^{5/2}log³(NT) ).

In particular Σ^hε_fZ_f(H) = −(1 + o(1))(1/2π)∫H when T → ∞ and T ≤ N^{1/3}(log N)^{−3}.
*Proof.* (i) 𝓔₁: ψ(u) − u ≪ √u log²u (RH, von Koch; cited) and partial summation, Σ_nΛ(n)Φ₁(n) − ∫Φ₁ =
−∫₁^∞(ψ(u) − u)Φ₁′(u)du + O(|Φ₁(1)| + ∫₀^1|Φ₁|): |𝓔₁| ≪ N^{−½}log²(NT)∫√u|Φ₁′| ≪ N^{−½}T log²(NT) (Lemma 6, μ₁ ≍ T^{−½}).
(ii) 2 ≤ m ≤ M₀, (m,N) = 1: K(p) := S(Np,1;Nm) depends on p mod m (Lemma 3(i)); for p ∤ m, K(p) = Σ_χĉ_χχ(p). So
Σ_p ln p K(p)Φ_m(p) = Σ_χĉ_χΣ_nΛ(n)χ(n)Φ_m(n) − Σ_χĉ_χΣ_{p^j,j≥2}χ(p^j)ln pΦ_m(p^j) + Σ_{p|m}ln pK(p)Φ_m(p). The last sum:
p ≤ m and |J_{k−1}(4π√(p/N)/m)| ≤ (2π/√(Nm))^{k−1}: total over m ≪ N^{−½}. The prime powers: ≤ m·Σ_{p^j}ln p|Φ_m(p^j)| ≪
m T^{½}μ_m log(NT); with the factor N^{−½}m^{−1} and summed over m: ≪ N^{−½}T^{3/2}log²(NT). Main part: under GRH,
|ψ(u,χ) − δ_{χ₀}u| ≤ C√u log²(mu) for every χ mod m (imprimitive χ differ from the primitive one by ≤ log m·log u; cited),
and partial summation: |Σ_nΛ(n)χ(n)Φ_m(n) − δ_{χ₀}∫Φ_m| ≪ log²(mNT)T^{3/2}μ_m. The χ₀-term δ_{χ₀}ĉ_{χ₀}∫Φ_m is ≪_A T^{1−A}
(Remark after Theorem A; |ĉ_{χ₀}| ≤ 1). With Σ_χ|ĉ_χ| ≤ m (Lemma 3(iii)):
|𝓚_m| ≪ N^{−½}T^{3/2}μ_m log²(mNT). Σ_{m≤T}μ_m ≪ T, Σ_{T<m≤M₀}μ_m ≪ T log M₀ (k = 2), ≪ T (k ≥ 3): Σ_{m≤M₀}|𝓚_m| ≪
N^{−½}T^{5/2}log³(NT). (iii) m ≡ 0 (N): c = N²m′, |S| ≤ τ(c)c^{½}, μ_m = (T/m)^{k−1} (T ≤ N): Σ ≪ N^{ε}T^{3/2}(T/N)^{k−1} ≤
N^{−½}T^{5/2}. (iv) m > M₀, (m,N) = 1: |𝓚_m| ≤ N^{−½}m^{−1}τ(m)m^{½}Σ_p ln p|Φ_m(p)| ≪ τ(m)T^{3/2}m^{−½}(T/m)^{k−1};
Σ_{m>M₀} ≪ T^{k+½}M₀^{3/2−k}log M₀ ≪ log T (k = 2, M₀ = T⁵), ≪ T^{5−2k} (k ≥ 3, M₀ = T³). (v) Theorem A collects. ∎

**Theorem C (unconditional).** There is an absolute c₀ > 0 with Σ^hε_fZ_f(H) = −(1/2π)∫H·(1 + o(1)) for T → ∞,
T ≤ exp(c₀√log N). *Proof.* As B, with Page's theorem (Davenport, *Multiplicative Number Theory*, ch. 20; **cited, not
opened**): for m ≤ exp(c₂√log u), ψ(u,χ) = δ_{χ₀}u − δ_{χ=χ₁}u^{β₁}/β₁ + O(u e^{−c₃√log u}). The error now enters through
∫u|Φ_m′|e^{−c₃√log u}du ≪ √N T^{5/2}μ_m e^{−c₃√log x₋}: |𝓔₁| ≪ T²e^{−c₃√log N}, Σ_{m≤M₀}|𝓚_m| ≪ T^{7/2}e^{−c₃√log N}log M₀ with
M₀ := exp(c₂√log x₋); the tail m > M₀ by Weil as (iv): ≪ T^{5/2}M₀^{−½}log M₀. The exceptional term contributes
−ĉ_{χ₁}Φ̃_m(β₁), Φ̃_m(s) = ∫u^{s−1}Φ_m: by Lemma 4 (convergent for Re w = β₁ − ½ < 3/4, Fubini as in Theorem A(5))
Φ̃_m(β₁) = ∫H(t)e^{2iθ(t)}(Nm²/4π²)^{w}Γ((k−1)/2 + w)/Γ((k+1)/2 − w)dt, w = β₁ − ½ − it, whose phase has t-derivative
−2 ln m + O(1/t) (no stationary point for m ≥ 2): A integrations by parts give ≪_A √N m T^{2−A}. Choose c₀ ≤ min(c₂, c₃)/8. ∎
(Siegel zeros are real; the resonance lives at height ≈ T, where only the classical zero-free region is needed.)

**Remark (sharper range; sketch, GAP 1).** Replacing partial summation by the explicit formula for L(s,χ) with the exact
transform Φ̃_m(½+iγ) = ∫H(t)e^{2iθ(t)}(Nm²/4π²)^{i(γ−t)}Γ((k−1)/2 + i(γ−t))/Γ((k+1)/2 − i(γ−t))dt (kernel modulus
|(k−1)/2 + i(γ−t)|^{−1}, stationary point γ = t(1 ± 1/m)) gives |Φ̃_m(ρ)| ≪ √(m/T) per zero, hence N^{−½}T²log²(NT) and the
range T ≤ N^{½−ε}: the analogue of ILS's support (−2, 2) (x = NT² ≤ N^{2−ε}). The uniform bound for the zero sum (degenerate
stationary point at γ ≈ 0, zeros far from [T, 2T]) is not written out.

## 5. The constant M = 1 and the harmonic average of c_f

**Proposition 5 (assumes GRH for L(s, Sym²f) for all f ∈ H_k(N), as ILS (129); proved (cited)).**
Σ^h c_f = 1 + O(N^{−1/4+ε}).
*Proof.* c_f = G_f(1) (§1; F(1+x) = x c_f + O(x²) and ζ(1+x) = 1/x + O(1)). G_f(s) = B_f(s)/L(s,Sym²f)·(N-factor), B_f(s)
absolutely convergent for Re s > ½ (B_p = 1 + O(p^{−2s})). Under GRH(Sym²f), 1/L(s,Sym²f) ≪ (N(|t|+2))^{ε} on Re s ≥ ½ + ε
(standard consequence of GRH, e.g. Iwaniec–Kowalski §5.9; not opened). Mellin: Σ_n g_f(n)n^{−1}e^{−n/y} =
(1/2πi)∫_{(1)}G_f(1+w)Γ(w)y^w dw = c_f + O(y^{−½+ε}N^{ε}), g_f the coefficients of G_f (n = dN^j, g_f = μ(d)λ_f(d²)·(1 − 1/N)^{[j≥1]}).
Average with Lemma 1: Σ^hλ_f(d²) = δ_{d,1} + O(dN^{−3/2+ε}) (Weil and J_{k−1}(y) ≪ y). Hence Σ^hc_f = (1 + 1/N + O(e^{−N/y}))
+ O(yN^{−3/2+ε}) + O(y^{−½+ε}N^ε); y = N^{½}. ∎
**Comparison.** Theorem B/C give M = 1 (the c = N Kloosterman term, a Ramanujan sum, times the Bessel–Mellin identity).
The conjecture gives M_pred = Σ^h ε_f²c_f = Σ^h c_f → 1. **The conjecture's constant is right on average.** Mechanism:
under harmonic weights λ_f(p) becomes Sato–Tate distributed (Σ^hλ_f(n) → δ_{n,1}), and c_f's local factor 1 − λ_f(p²)/p has
Sato–Tate mean 1 (a8). (The note's §6(c) anticipated this; it is now carried out.)
**Natural weights (sketch, needs GRH(Sym²) for the ILS conversion (136)–(137)).** ω_f^{−1} ∝ L(1,Sym²f) ∝ Σ_tλ_f(t²)/t, so
the natural mean of c_f is Σ^h B_f/Σ^h L(1,Sym²f) → ∏_p E_ST[B_p]/E_ST[L_p(Sym²,1)] = ∏_{p≠N}(1 − p^{−2})·(1 + 1/N) =
(6/π²)/(1 − 1/N) (a8: E_ST[(1 − λ(p²)/p)Σ_jλ(p^{2j})p^{−j}] = 1 − p^{−2}); with Theorem D at ℓ = t² (t ≤ N^{δ}):
|H_k(N)|^{−1}Σ_f ε_fZ_f(H) = −(6/π²)(1/2π)∫H(1 + o(1)). The natural-weight family constant is 6/π², not 1.
**Computed (c, P = 5000 for B_f, all embeddings via PARI):** Σ^h c_f at k = 4, N = 13…47 (10 levels): 1.03, 1.36, 1.44, 0.96,
0.78, 1.20, 0.81, 0.97, 0.87, 1.16 (mean 1.06); k = 2, N = 43…97 (9 levels): 1.09, 1.10, 0.71, 1.34, 0.97, 1.09, 0.99, 1.07,
1.02 (mean 1.04). Natural means 0.54–0.87 (k = 4) and 0.67–0.88 (k = 2) against the limit 0.61. With 3–11 forms per
level these scatter widely; they are consistent with the limits, not a test of them (the proved rate is N^{−1/4}).

## 6. Twisted first moments and lines (B)

**Theorem D (GRH as in Theorem B; proved as Theorem B).** Fix ℓ ≥ 1 with (ℓ,N) = 1 and r > 0. Then
Σ^hε_fλ_f(ℓ)Z_f(H; r) = −(1/2π)∫H·A(ℓ, r) + O_{ℓ,r}(T^{½} + N^{−½}T^{5/2}log³(NT)), where
A(ℓ, r) = μ(m)c_m(ℓ)/(φ(m)√ℓ) if r = ℓ/m² with m ≥ 1 an integer (and (m, N) = 1, automatic for N large), else A = 0.
*Proof.* Note, Proposition 1 with the twist (I_n → I_{nr}). For p ∤ ℓ, ε_fλ_f(ℓ)λ_f(p) = −i^k√Nλ_f(Nℓp) (p | ℓ: finitely
many primes, O(T^{½})). Lemma 1: the m-th term is N^{−½}m^{−1}Σ_p ln p S(Nℓp,1;Nm)Φ_{m}^{(ℓ,r)}(p) with
Φ_m^{(ℓ,r)}(u) = u^{−½}J_{k−1}(4π√(ℓu/N)/m)I_{ur}. Its χ₀-part has ĉ_{χ₀} = μ(N)μ(m)c_m(ℓ)/φ(m) (Lemma 3(i),(ii)) and, by
(★) with N → Nm²/ℓ, ∫Φ_m^{(ℓ,r)} = (√N m/(2π√ℓ))Ĥ(ln(rm²/ℓ)): = (√Nm/2π√ℓ)∫H if rm² = ℓ, ≪_A T^{1−A} otherwise
(|ln(rm²/ℓ)| ≥ c(r,ℓ) > 0). Product: −(1/2π)∫H·μ(m)c_m(ℓ)/(φ(m)√ℓ). The rest is bounded exactly as in Theorem B. ∎
Special cases: ℓ = 1: lines exactly at r = 1/m², amplitude μ(m)²/φ(m) (r = 1: 1); r = 1, ℓ = m²: μ(m)/m; r = ℓ = p: 1/√p.

**Proposition D′ (proved; a7).** For p ∤ N let e = v_p(ℓ), v = v_p(r). With the conjecture's local factor
(c_fρ_f)_p = p^{|v|/2}S_p(v)/(1 − 1/p), S_p(v) = Σ_{m−h=v}λ(p^m)μ_f(p^h)p^{−(m+h)/2} (derivation.md §2.1), and Sato–Tate
expectation E_ST (λ(p^j) = U_j(cos θ), orthonormal):
E_ST[λ(p^e)S_p(v)]/(1 − 1/p) = p^{−e/2} (v = e); −p^{−e/2} (v = e − 2, e ≥ 1); 1/(p−1) (v = −2, e = 0); 0 otherwise,
which is the p-factor of μ(m)c_m(ℓ)/(φ(m)√ℓ)·(ab)^{½}·(ab)^{−½}. *Proof.* μ_f(p) = −λ_p, μ_f(p²) = 1; U₁U_j = U_{j−1} + U_{j+1};
E[U_iU_j] = δ_ij. For v = e: E[U_e(U_e − U₁U_{e+1}/p + U_{e+2}/p²)]p^{−e/2} = (1 − 1/p)p^{−e/2}. For v = e − 2: the h = 1, 2 terms
give (−1 + 1/p)p^{−e/2}. Parity or v > e or v ≤ e − 4: every product has zero mean. ∎ (sympy: all e ≤ 3, |v| ≤ 4.)
**Consequence.** With Proposition 5's argument (GRH(Sym²)) Σ^hλ_f(ℓ)c_fρ_f(r)/√(ab) → A(ℓ, r): **every Hecke-twisted first
moment, at every rational line, is the one the conjecture predicts.** Equivalently the residual R_f := Z_f + ε_fc_f∫H/2π
satisfies Σ^h ε_fλ_f(ℓ)R_f(H; r) = o(T) for all fixed ℓ, r.
**Slope (sketch).** Approximating c_f by Σ_{d≤N^δ}μ(d)λ_f(d²)/d (GRH(Sym²)) and using Theorem D with ℓ = d² (uniformity in
ℓ ≤ N^{2δ} not written; needs |Z_f| ≪ T log NT, i.e. GRH for L(s,f)): Σ^hε_fc_fZ_f(H) = −(∫H/2π)Σ^hc_f² + o(T), with
Σ^hc_f² → ζ(2)/ζ(4)-type constant (a8: E_ST[(1 − λ(p²)/p)²] = 1 + p^{−2}). So the harmonic regression slope of ε_fZ_f on
the law −c_f∫H/2π is 1 + o(1): the family form of the project's slope test (0.989 ± 0.011 on 158 curves).

## 7. Second moment and "almost all f" (C)

Expanding |Z_f|² with Proposition 1 and Lemma 1 (λ_f(p)λ_f(p′), no root number):
Σ^h|Z_f(H)|² = 𝓓 + 𝓞 + (prime powers, E₀), 𝓓 = (1/4π²)Σ_p ln²p·p^{−1}|I_p|²·Σ^h1, and by stationary phase
(|I_p|² = πt_pH(t_p)²(1 + O(1/T))) and the prime number theorem **𝓓 = (1/2π)∫H(t)²ln(Nt²/4π²)dt·(1 + o(1))** (proved, as the
note's Theorem 2(c)): the random-phase variance, of order T log(NT). 𝓞 = (1/2π)i^{−k}Σ_{m≥1}(Nm)^{−1}Σ_{p≠p′}ln p ln p′(pp′)^{−½}
I_pĪ_{p′}S(p,p′;Nm)J_{k−1}(4π√(pp′)/(Nm)): here the Bessel argument is ≈ t_pt_{p′}/(πm) ≍ T²/m, oscillatory for all m ≲ T²,
and |S(p,p′;Nm)| ≍ (Nm)^{½} has no Ramanujan reduction. Trivially 𝓞 ≪ N^{½}T⁴. For "almost all" one needs
(S2) Σ^h|Z_f(H)|² ≤ (Σ^hc_f²)(∫H/2π)² + o(T²),
i.e. 𝓞 must carry a main term ≍ T² and cancel most of 𝓓 (the data's rigidity, scatter 13.9 against √(#zeros) ≈ 29, says
𝓞 indeed lowers the variance). **Proposition (proved, trivial):** (S2) and the slope statement of §6 give
Σ^h|Z_f + ε_fc_f∫H/2π|² = o(T²), hence for every η > 0 the harmonic measure of {f : |Z_f(H) + ε_fc_f∫H/2π| > ηT} → 0.
(S2) is **GAP 2**: it is a statement about Σ_{p,p′}S(p,p′;c) with c up to ≈ x = NT² (character expansions mod c have ≈ x
terms; GRH for Dirichlet L-functions does not reach it). It would follow from the ratios conjecture for the family (two
L-functions, root-number twisted); not attempted. Under GRH for every L(s,f) only |Z_f| ≪ T log(NT), so Σ^h|Z_f|² ≪ T²log².

## 8. One form (D)

**8.1 Δ and the Kloosterman expansion (proved; d6).** For N = 1, k = 12, S₁₂(1) = ℂΔ, Lemma 1 gives the exact identity
λ_Δ(n) = (2π/ω_Δ)Σ_{c≥1}S(n,1;c)c^{−1}J₁₁(4π√n/c), ω_Δ = 2.840287375167 = 2π²/(11L(1,Sym²Δ)), L(1,Sym²Δ) = 0.631792945721,
ε_Δ = i^{12} = 1. Inserting it into Proposition 1, the c = 1 term is −(1/ω_Δ)Σ_p ln pΦ₁(p) with ∫Φ₁ = (1/2π)∫H by (★):
**the c = 1 term carries the amplitude 1/ω_Δ = 0.3521** (its PNT remainder is o(T) under RH by the zero-sum argument of the
Remark in §4: sketch; computed: −10.61 against −(1/ω_Δ)∫H/2π = −10.65 on Hann [20, 400]). The conjecture: c_Δ = B_Δ/L(1,Sym²Δ)
= 1.4915509730/0.6317929457 = **2.360822** (direct partial product to 2·10⁵: 2.36043). The zeros (§9) give 2.40 (sharp) and 2.33
(Hann). So the moduli c ≥ 2 carry (c_Δ − 1/ω_Δ)/c_Δ = **85 %** of the amplitude. Each fixed c ≥ 2 is non-resonant (its
χ₀-integral is (m/2π)Ĥ(2 ln c)·…, negligible), so the main term is built by moduli c → ∞ with T; with J₁₁(y) negligible for
y ≲ 6, only c ≤ 4π√p/6 ≈ 2T/6 matter. Computed (Hann [20, 400], cumulative Σ_{c≤C} of the prime side, fraction of the law):
C = 1: 0.149; 2: 0.147; 5: 0.164; 10: 0.240; 20: 0.465; 30: 0.666; 40: 0.853; 60: 1.011; ≥ 80: 1.022 (stable): **the main term
is built across 5 ≲ c ≲ 60 at heights t ≤ 400 (c between T/80 and T/7), not by any fixed modulus**.

**8.2 What a Hypothesis-S-type bound can and cannot do (proved, given the conjecture for Δ).** ILS (arXiv version,
Hypothesis 4, = Hypothesis S of the journal version; opened): Σ_{p≤X, p≡a (c)}e(2√p/c) ≪_ε X^{½+ε}; "all that we need … is a
bound X^δc^A with some δ < 3/4" (δ = 7/8 Vinogradov, δ = 3/4 from the density hypothesis). In our problem the c-th term is
c^{−1}Σ_p ln p S(p,1;c)Φ_c(p), Φ_c ≈ p^{−½}(c/T)^{½}T^{½}e(2√p(±1/c − 1)): the analytic phase has frequency 1 ∓ 1/c, not 1/c,
and the arithmetic weight is a Kloosterman sum. (a) If for some δ > 0, uniformly in c ≤ T/3,
|Σ_p ln p S(p,1;c)Φ_c(p)| ≪ c^{½}x^{−δ}sup_u √u|Φ_c(u)| (more than square-root cancellation; the exponent x^{½−δ} of the first version was the trivial bound, referee's correction; the range should be c ≤ T^{1+η}), then Σ_{2≤c≤T/3}|c-th term| ≪
T^{1−2δ} and Z_Δ(H) = −(1/ω_Δ)∫H/2π + o(T) (under RH for the c = 1 remainder): **this contradicts the conjecture, since
c_Δ = 2.3608 ≠ 0.3521 = 1/ω_Δ**. So the conjecture for Δ implies that no such bound holds. (b) Square-root cancellation
(δ = 0, an "S-type" hypothesis) gives each c-term O(T^{ε}) and the sum O(T^{1+ε}): it neither proves nor refutes the law.
This holds **within the Petersson decomposition only** (referee's correction, §9a): the split over c is a property of the
Kloosterman expansion of λ_Δ, not of primes, and in the Hardy–Littlewood decomposition an upper bound on the minor arcs
would suffice. GRH alone gives only the trivial bound (note §6(b), refereed).

**8.3 Theorem E (conditional; proved).** Let f be one form (any level, weight k), H as in §1. Assume, uniformly for complex
α, γ with |α|, |γ| ≤ 2/log T, Re γ ≥ 1/log T, and smooth w supported in [T/2, 3T] with w^{(j)} ≪_j T^{−j}:
(R1) ∫w(t)L(½ − α − it, f)/L(½ + γ + it, f)dt = F_f(1 − α + γ)∫w + O(T^{½+ε}),
(R2) ∫w(t)e^{2iθ(t)}L(½ + α − it, f)/L(½ + γ − it, f)dt = O(T^{½+ε}),
with F_f(s) = Σλ_f(n)μ_f(n)n^{−s} (μ_f the coefficients of 1/L(s,f)). (R1) is the CFKRS ratios conjecture for L(s,f) (one
L at −it over one at +it; the recipe keeps the diagonal n = m and drops the swapped term, which carries X(½ − α − it)(nm)^{−it});
(R2) is the recipe's statement that the "e^{4iθ}" term of derivation.md §2.1 has no main term. Then
Z_f(H) = −ε_f c_f(1/2π)∫H + O(T^{½+ε}).
*Proof.* Let δ ∈ [1/log T, 2/log T] with no zero on Re s = ½ ± δ. By the residue theorem (horizontal edges → 0 as in the
note's Proposition 1), Z_f(H) = (1/2πi)[∫_{(½+δ)} − ∫_{(½−δ)}]h(s)X(s)^{−1}(L′/L)(s,f)ds. *Right line:* (L′/L)(½+δ+it) =
∂_α[L(½+α+it)/L(½+δ+it)]_{α=δ} and X^{−1}(½+δ+it)L(½+α+it) = ε_fκ_α(t)L(½−α−it), κ_α(t) := X(½+α+it)/X(½+δ+it)
(smooth, κ_δ = 1, ∂_jκ ≪ (log NT)^j). By (R1) with w = h(½+δ+i·)κ_α (smooth after a partition of unity; Gaussian tails
negligible) and Cauchy's formula for ∂_α on |α − δ| = 1/log T (error T^{½+ε}log T):
right = (ε_f/2π)∂_α[∫hκ_αF_f(1 − α + δ)]_{α=δ} = (ε_f/2π)[F_f(1)∫h∂_ακ − F_f′(1)∫h(½+δ+it)dt] = −(ε_fc_f/2π)∫H, using
F_f(1) = 0, F_f′(1) = c_f (derivation.md §2.1) and ∫h(½+δ+it)dt = ∫H(t − iδ)dt = ∫H (H entire, Gaussian decay).
*Left line:* (L′/L)(s) = (X′/X)(s) − (L′/L)(1−s); the X′/X part is ∫H(t+iδ)X^{−1}(½−δ+it)(X′/X)dt = O(1) (phase 2θ, 2θ′ ≍ log(NT²),
first-derivative test); the other part is (1/2π)∫H(t+iδ)X^{−1}(½−δ+it)(L′/L)(½+δ−it)dt with X^{−1}(½−δ+it) = e^{2iθ(t)}·
(smooth) — O(T^{½+ε}) by (R2) differentiated as above. ∎
Status: (R1)+(R2) is the cleanest hypothesis found under which the one-form law is a theorem; it is strong (a ratios
conjecture). By 8.2 no Hypothesis-S-type upper bound can replace it, and GRH does not (note §6(b)). Weaker natural
candidates (a type II hypothesis for Λ_f = (λ_f log) ∗ μ_f with prescribed main terms; note §6(d) step 4) are not clean.

## 9. Numerical table (check_family_out_*.txt)

**Δ (k = 12, N = 1, ε = 1), zeros of L(Δ) to 420 (PARI lfunzeros; 431 zeros; N(T) vs θ(T)/π: 59/58.96, 159/159.36,
277/276.41, 404/404.27).** c_Δ = 2.360822, 1/ω_Δ = 0.352077.

| window | Σ e^{2iθ(γ)} sharp | law −c_ΔT/2π | ratio | Hann | law/2 | ratio | c = 1 term only (1/ω_Δ) |
|---|---|---|---|---|---|---|---|
| [20,400] | −145.03 − 1.20i | −142.78 | 1.016 | −70.61 − 0.03i | −71.39 | 0.989 | −21.29 |
| [20,140] | −40.56 − 0.25i | −45.09 | 0.900 | −20.50 − 0.15i | −22.54 | 0.909 | −6.72 |
| [140,260] | −45.43 − 6.53i | −45.09 | 1.008 | −21.10 + 1.15i | −22.54 | 0.936 | −6.72 |
| [260,400] | −59.04 + 5.58i | −52.60 | 1.122 | −28.32 + 3.56i | −26.30 | 1.077 | −7.84 |
| [60,230] | −63.43 + 0.65i | −63.88 | 0.993 | −33.28 − 2.89i | −31.94 | 1.042 | −9.53 |
| [230,380] | −57.38 + 2.75i | −56.36 | 1.018 | −30.04 + 6.06i | −28.18 | 1.066 | −8.41 |

Noise scale: |Im| ≤ 6.5 sharp, ≤ 6.1 Hann on the sub-windows; the full-window deviation from the law is 2.3 (sharp) and 0.8
(Hann); the alternative 1/ω_Δ is off by 124 (sharp). Exact explicit formula (Proposition 1, weight 12, Gaussian boxes
[100,200] Δ = 5 and [250,330] Δ = 6, I_n by quadrature over all prime powers): |zeros − primes| = 1.2·10^{−12}, 2.4·10^{−12}
(completeness of the zero list in those windows). Prime side at the exact stationary point, Hann [20,400]: −70.461 − 0.132i
against the zeros −70.613 − 0.034i (with the simplified kernel t_n* = 2π√n the phase error O(k²/t) is large for k = 12:
−72.12 − 8.51i). Resonance sums S(x)/((2√2/3)e^{−23πi/4}x^{3/4}): over all integers 0.992 + 0.056i, 0.996 + 0.041i, 0.998 + 0.029i
at x = 5·10⁴, 10⁵, 2·10⁵ (theorem: 1; this also checks the weight-k phase e^{−i(2k−1)π/4}); over prime powers 2.271 + 0.250i,
2.292 + 0.186i, 2.301 + 0.076i (law 2.361; c = 1 term 0.352).
**Petersson and root numbers:** b1–b5 (§2). **Family averages:** §5 (c1, c2). **Identities:** a1–a8 (§2, §6).
Totals: a 8/8, b 11/11, c 2/2, d 9/9 PASS (see the output files).

## 9a. Referee report (25 November 2026) and what it adds

An independent referee (own code in [referee/](referee/)) confirmed: Lemmas 1–6, the weight ω_f = 2π²/((k−1)·N·L(1, Sym²f)) (12
digits on four curves), the root number ε_f = −i^k·λ_f(N)·√N (42 embeddings, both signs), **Theorem A (unconditional, exact; main
term exactly −(1/2π)∫H)**, **Theorem B (GRH; T ≤ N^{1/3}(log N)^{−3})**, Proposition 5 (GRH for Sym²f), Theorem D for fixed ℓ, r,
Proposition D′. Independent tests: N = 61, k = 4 from eigenform data: family sum/(−∫H/2π) = 0.86, 0.92, 1.15, 1.03 at T = 4, 6, 8,
10; N = 10007, k = 4: family sum minus the predicted line sum ≤ 0.026 in ten (ℓ, r) cases at T = 6. Corrections:

1. **Theorem C** needs |ĉ_χ| = |τ(χ)|²/φ(m) ≤ q₁/φ(m) and the gain (T ln m)^{−A} from integration by parts for the exceptional
   character; with these the total is ≪ T^{1−A}. An exceptional zero is harmless (no stationary point for m ≥ 2).
2. **Theorem E:** hypothesis (R2) is redundant. Under GRH the left line moves to −¼ with no zeros crossed and is the term E₀ of
   Proposition 1 of the first note. So: **the ratios conjecture (R1, the one-over-one case of Conrey–Farmer–Zirnbauer, with
   |α| ≤ 3/log T) and GRH for L(s, f) imply the one-point law.** This is the recipe written as a hypothesis.
3. **§8.2** as first written used the trivial bound in place of square-root cancellation and drew a conclusion about primes
   that holds only inside the Petersson decomposition (corrected above). For the sum over all integers (constant 1, a
   theorem) the c = 1 term also carries exactly 1/ω_Δ = 0.352, and moduli 10 < c ≤ 60 carry the other 65 %.
4. **Theorem D** should be stated in its uniform form −(1/2π)·Σ_m μ(m)c_m(ℓ)/(φ(m)√ℓ)·Ĥ(ln(rm²/ℓ)): at N = 10007 and r = ½ the
   family value is as large as the main term at T = 6. Proposition 5 alone does not test c_f (c_f ≡ 1 also has mean 1); the
   discriminating statements are Theorem D at ℓ = m², r = 1 and Proposition D′.
5. Smaller points: for N | m the Kloosterman sum S(Np, 1; Nm) is exactly 0; Theorem B's range N^{1/3} comes from partial
   summation, the explicit formula per character should give N^{½−ε}, and T ≈ N^{½} is the support-2 barrier of
   Iwaniec–Luo–Sarnak; no assertion of check_family.py tests Theorems A, B or D directly (referee/r3_family.py and
   referee/r6_largeN.py do).

**The constant is a Hardy–Littlewood singular series (referee, level 1).** Expand Λ(n) ≈ Σ_q μ(q)c_q(n)/φ(q). By the twisted
functional equation each sum Σ_n c_q(n)λ(n)·(kernel) resonates only at dual index q², with weight (φ(q)/q)·λ(q²) times the
q = 1 value (computed for Δ: q = 2, 3, 5, 7: 10.866, 12.928, 12.689, 22.08 against 10.867, 12.933, 12.633, 22.24). So the
constant is Σ_q μ(q)λ_f(q²)/q = Π_p (1 − λ_f(p²)/p) = c_f (for Δ: partial sums 1, 1.915, 2.250, 2.335, 2.356, 2.3595, 2.3605 at
Q = 1, 10, …, 10⁶; c_Δ = 2.3608). The one-point law for one form is therefore the singular-series conjecture for the Voronoi
resonance over primes: the major arcs give the constant, and what is open is a minor-arc bound o(x^{3/4}) (Zhao 2006 has
x^{5/6}). It does not follow from GRH for Dirichlet L-functions: the classical main terms of primes in progressions, summed
over the Petersson moduli, give 1/ω_Δ = 0.352, not 2.3608. Σ^h c_f at six further levels (k = 4, N = 61–151): mean 1.015.

## 10. Prior art (opened sources marked)

- **Iwaniec–Luo–Sarnak** (arXiv:math/9901141; **opened**, DVI text): §3 (large level, N prime, k = 2): Petersson with harmonic
  weights (125)–(127), ε_f = −μ(N)λ_f(N)√N (118), and in §3.2 the term IV_B: A*[λ_f(Nn²)λ_f(p)] via S(Nn²,p;lN) =
  μ(N)S(n², N̄p; l), characters mod l, GRH, the mean μ(l)Σ_{d|(n²,l²)}…, PNT; their main term (−½∫_{|ξ|>1}φ̂) comes from
  ∫J₁(y)dy = 1 at the transition p ≈ N, with a non-oscillating weight. Hypothesis 4 (§1; Hypothesis S in the journal
  version) as quoted in 8.2. **Theorem A/B: the method is KNOWN (ILS IV_B); the statement — the root-number-weighted
  one-point sum at height T, resonance through (★), constant 1 — NOT FOUND.** Lemma 2, Lemma 3(i): KNOWN (ILS).
- Proposition 5 (harmonic mean of B_f/L(1,Sym²f)): **KNOWN IN SUBSTANCE** (ILS (129) is the same GRH(Sym²) approximation;
  harmonic moments of L(1,Sym²f)^{±1}: Royer 2001, Cogdell–Michel 2004, Lau–Wu: **not opened** (the arXiv API returned
  nothing in the time available)).
- Theorem D/D′ (twisted lines = Sato–Tate average of the conjecture's ρ_f(r)): **NOT FOUND**.
- Hughes–Miller (n-level densities with root-number splitting beyond support 1), Ricotta–Royer (lower-order terms of
  one-level densities, Sym^r families), "phase sums / Gram points for automorphic L-functions on average": **not opened**;
  no source seen that averages Σ_γ e^{2iθ_f(γ)} over a family.
- Theorem E (ratios ⇒ law): **KNOWN IN SUBSTANCE** (Conrey–Snaith-type deductions from the ratios conjecture; the specific
  twisted ratio (R2) not seen). 8.2 (Hypothesis S cannot give the one-form law): **NOT FOUND**.

## 11. GAPs and doubts

GAPs: (1) Remark §4: range N^{½−ε} under GRH (zero-sum bounds not written); Theorem B as proved stops at N^{1/3}.
(2) (S2), hence "almost all f" (§7). (3) Unconditional Proposition 5 (needs a zero-density input for the Sym² family; cited
works not opened). (4) The natural-weight corollary and the slope statement of §6 are sketches (uniformity in ℓ ≤ N^{δ}).
(5) 8.1's statement that the c = 1 PNT remainder is o(T) for N = 1 uses the zero-sum argument of GAP 1 (with N = 1 partial
summation gives only O(T log²T)); the computed value −10.61 vs −10.65 supports it. (6) The one-form law itself (Theorem E
is conditional).
Doubts: (1) Lemma 1 for general k and Page's theorem are cited without opening; GRH bounds for 1/L(s,Sym²f) likewise.
(2) Theorems B/C need T ≤ N^{1/3}: the family statement is at heights T small compared with the level, the opposite of the
per-curve data (T ≈ 400, N ≈ 300–3000); it is the average law, not the law for fixed f as T → ∞. (3) ε_f = i^kμ(N)λ_f(N)√N
for k ≥ 4 was checked numerically (b5), not derived here (Atkin–Lehner theory). (4) Theorem A–C need k ∈ K (no oldforms);
for k ≥ 12 or composite N the oldforms and the newform-Petersson formula (ILS Prop. 2.8) would enter. (5) Section c uses
L(1,Sym²f) from PARI with numerically generated Dirichlet coefficients (lfuncheckfeq ≈ −46…−50) and B_f truncated at
p ≤ 5000 (relative error ≈ 10^{−4}); c_f can be negative (λ_f(2)² > 3). (6) The Δ result is one form at T ≤ 400; the
sub-window ratios range 0.90–1.12. (7) The c-decomposition of 8.1 uses the leading stationary phase at the exact point (its
total matches the zeros to 0.18).

## 12. What a referee should check first

1. (★) and Theorem A(5): the conditional convergence of the Bessel–Mellin integral and the Fubini/dominated-convergence step.
2. Lemma 2 and Lemma 3(i): S(Np,1;N) = μ(N) for every p, and the sign chain ε_f = −i^kλ_f(N)√N → −2πN^{−½}Σ_m.
3. Theorem B (ii): Σ_χ|ĉ_χ| ≤ m, the bound ∫√u|Φ_m′| ≪ T^{3/2}μ_m, and the sum over m (where the exponent 5/2 comes from).
4. Proposition 5 and the comparison M = 1 = lim Σ^hc_f (harmonic, not natural weights).
5. Theorem D's χ₀ computation and Proposition D′ (the Sato–Tate identity, a7).
6. 8.2(a): the implication "power-saving bound for Kloosterman–prime sums ⇒ amplitude 1/ω_Δ", and the numbers 1/ω_Δ = 0.3521,
   c_Δ = 2.3608 (three computations of ω_Δ agree; c_Δ against the direct product).
7. Theorem E: the left-line term and hypothesis (R2).
