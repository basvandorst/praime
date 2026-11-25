# Major arcs of the radial resonance of Gaussian primes: the axis zone (theorem), the Pythagorean rays (theorem given a cited zero-free region with growth bounds for Hecke L-functions with modulus; no GRH), the singular series as a sum over rays, the open minor arcs

Folder: `scratchpad/arcs/`. Nothing in the repository was changed. Labels: **proved** (complete proof here from the
cited statements), **literature (opened)** (statement quoted from a source opened in this session),
**literature (secondary)** (statement as quoted by another opened source), **not opened**, **computed**,
**heuristic**, **GAP** (a step that is not proved; what is missing is said exactly).

## 0. Result in short

| item | status |
|---|---|
| Lemma 1: exact factorisation e(−\|z\|)ψ_m(z) = u^m·e(−Re(z w̄)/h)·[smooth zone factor] | **proved** (algebra) |
| Lemma 2: ζ ↦ Re(ζ w̄) mod h maps (ℤ[i]/h)^* onto (ℤ/h)^* with equal fibres; mean of e(−·/h) = μ(h)/φ(h); A_h = (4/π)μ(h)/φ(h) for every hypotenuse h (0 if h is not squarefree) | **proved**; checked h ≤ 145 |
| Lemma 3: ∫G dA = W_V Φ_M X^{3/4} + O(X^{1/4}); Φ_M = e^{−iπ/4} + O_N(M^{−N}) for a cutoff ≡ 1 near 0 | **proved** (stationary phase); computed |
| **Theorem A (axis zone, h = 1)**: Z_axis(X; M) = (4/π)·I(X; M) + o(X^{3/4}) | **proved, unconditional**, from Stucky's Theorem 1.1 (= Ricci's theorem at θ = 1) |
| **Theorem B (Pythagorean ray, hypotenuse h)**: Z_u(X; M) = A_h u^m I(X; M) + o(X^{3/4}) | **proved assuming G_h** (Coleman 1990, Theorems 1 and 2 with modulus h: zero-free region and growth bounds; known here only from secondary quotations), through R_h, which is proved from G_h in [rh/rh_proof.md](rh/rh_proof.md) (refereed, [rh/referee/referee_rh.md](rh/referee/referee_rh.md)); no GRH |
| Theorem B′: same with error O(X^{1/2} log² X) under GRH for the Hecke L-functions of ℚ(i) of modulus h | **proved under GRH** (explicit formula; constants in M not optimised) |
| R_h for h > 1 | **proved up to G_h** ([rh/rh_proof.md](rh/rh_proof.md)): sectors of width x^(−θ) in residue classes for θ < 1/(2 + 8κ), κ the exponent of a pointwise bound for the Hecke L-functions in the full conductor: θ < 64/231 with Wu's bound (read), θ < 3/10 with Söhne's (secondary), θ < 1/(4 − 8δ) with Michel–Venkatesh (read). In print before: 5/24 (Panidapu–Thorner via Coleman's zero density). Not read: Coleman 1990 (the input G_h), Koval'chik 1974, Ricci's thesis: whether R_h is already in print is unknown |
| Corollary 1: Σ_{rays, h ≤ H} A_h ξ^m = (16/π)·S_H, S_H → (π/4)c (the identity c = (4/π)Π(1 − 2cos 2mθ_p/(p−1)) is the sum over rays) | **proved** (algebra exact; convergence by the standard Perron + zero-free-region argument, sketched); computed S_{10⁶} − (π/4)c = 2.8·10⁻⁵ |
| Corollary 2: the conjecture (smooth form) ⟺ the minor-arc statement (MA), given R_h for every h | **proved** (equivalence) |
| (MA) itself | **open** (trivial bound O(X); GRH gives only O(X^{3/4} log X) for the full sum) |
| numerics to N z ≤ 10⁷ (`check_arcs.py`, 85 PASS; `check_arcs_indep.py`, independent) | **computed**: axis zone vs Theorem A: 1.0 sd; Pythagorean zones: own term + neighbouring rays within 0.1–0.8 sd; own term alone 1.5–3.4 sd off (finite-X cross-talk, quantified); class statement of R_h: mean of e(−Re(zw̄)/h) within 0.4–1.2 sd of μ(h)/φ(h), never 0 mod p \| h |

## 1. Setting

- z = a + bi runs over Gaussian prime **elements** (all four associates). Λ(z) := ln N z (and ln N π for z = unit·π^k
  when prime powers are allowed; their total weight up to x is O(√x ln x) and never matters below).
- m ≡ 0 (mod 4), ψ_m(z) = (z/|z|)^m (a function of the ideal (z)). e(y) = e^{2πiy}.
- **Normalisation of the conjecture.** Summed over prime *ideals* (one generator each), P_m(x) = Σ_{N𝔭 ≤ x}
  ψ_m(𝔭) ln N𝔭 e(−√N𝔭) ~ K x^{3/4}, K = (2/3)c e^{−iπ/4} (primeside.md §1–§2; c = 2.026446 for m = 4). Summed over
  **elements** the sum is 4P_m. In smooth form (V ∈ C_c^∞((0, ∞))):

      Σ_full(X) := Σ_z Λ(z) ψ_m(z) e(−|z|) V(N z/X)  =  4c·W_V·e^{−iπ/4}·X^{3/4} + o(X^{3/4})      (Conjecture C)

  with W_V := (1/2)∫_0^∞ v^{−1/4} V(v) dv. (Check: 4·(3/4)·K·∫ y^{−1/4}V(y/X) dy = 4cW_V e^{−iπ/4}X^{3/4}.)
- **Directions.** D = the rational points of the unit circle: ξ = w/h with w = (x, y) ∈ ℤ², gcd(x, y) = 1,
  x² + y² = h² (h = 1: the four axes ±1, ±i). For h > 1, h is odd and every prime factor of h is ≡ 1 (mod 4); the
  primitive w with |w| = h are exactly w = ε·π₀² (ε a unit, π₀ = Π_{p|h} π_p^{(±)} a primitive Gaussian integer of norm
  h): 4·2^{ω(h)} vectors (checked for h ≤ 145 and, through the identity of §6, for every h ≤ 10⁴).
- **Zone sum.** Fix V ∈ C_c^∞((0, ∞)) with supp V ⊂ [a₀, b₀], 0 < a₀, and an even F ∈ C_c^∞(ℝ), F = 1 on [−1, 1],
  supp F ⊂ [−2, 2]; F_M(y) := F(y/M). For u = w/h ∈ D put t = Re(z ū) = Re(z w̄)/h, s = Im(z ū) (the coordinates of
  z along and across the ray) and

      Z_u(X; M) := Σ_{z prime, t > 0} Λ(z) ψ_m(z) e(−|z|) V(N z/X) F_M(s/√t).

  The zone is |s| ≤ 2M√t, a parabola around the ray; its angular half-width at radius r is ≈ 2M/√r ≍ M X^{−1/4}.
- **The common main-term integral.** With α = arg(z ū),

      G(z) := e^{imα} · e(−(|z| − t)) · V(N z/X) · F_M(s/√t)   (t > 0; G := 0 for t ≤ 0),      I(X; M) := ∫_{ℝ²} G dA.

  I does not depend on u (rotation invariance). The theorem's main term "of the axis zone with weight 1" is I.

## 2. Statements

**Lemma 1 (exact factorisation; proved).** For u = w/h ∈ D and every z with t > 0:
ψ_m(z) e(−|z|) = u^m · e(−Re(z w̄)/h) · G(z) / (V F_M). Proof: ψ_m(z) = (z/|z|)^m = u^m (zū/|z|)^m = u^m e^{imα};
e(−|z|) = e(−t) e(−(|z| − t)) and t = Re(z w̄)/h. ∎ No Taylor expansion of |z| = t + s²/2t − s⁴/8t³ + … is used:
it sits inside G, which is smooth.

**Lemma 2 (the weight of a ray; proved).** Let h > 1 and w primitive with |w| = h. For each p^k ∥ h let π_p be the
prime above p with π_p^{2k} | w (exactly one of π_p, π̄_p divides π₀). Then
(i) for ζ ∈ ℤ[i] the residue Re(ζ w̄) mod p^k depends only on ζ mod π_p^k, and ζ mod π_p^k ↦ Re(ζ w̄) mod p^k is the
map ζ ↦ (w̄/2)·ζ under ℤ[i]/π_p^k ≅ ℤ/p^k, a bijection of the unit groups;
(ii) hence ζ ↦ Re(ζ w̄) mod h maps (ℤ[i]/h)^* onto (ℤ/h)^* with all fibres of size Φ(h)/φ(h) = φ(h), where
Φ(h) = |(ℤ[i]/h)^*| = φ(h)²;
(iii) (1/Φ(h)) Σ_{ζ ∈ (ℤ[i]/h)^*} e(−Re(ζ w̄)/h) = c_h(1)/φ(h) = μ(h)/φ(h) (c_h the Ramanujan sum);
(iv) for every Gaussian prime z ∤ h, Re(z w̄) is not ≡ 0 modulo any prime factor of h.
Proof. 2Re(ζw̄) = ζw̄ + ζ̄w and π_p^{2k} | w, so 2Re(ζw̄) ≡ ζ w̄ (mod π_p^k); w̄ is prime to π_p (π_p ∤ π̄₀) and 2 is
a unit mod p. ℤ/p^k → ℤ[i]/π_p^k is an isomorphism (p splits). This gives (i); the CRT over p^k ∥ h and the
projection (ℤ[i]/p^k)^* → (ℤ[i]/π_p^k)^* (fibres |(ℤ[i]/π̄_p^k)^*| = φ(p^k)) give (ii); (iii) is the mean of
e(−r/h) over r ∈ (ℤ/h)^*; (iv) is (i) with k ≥ 1 for z coprime to p. ∎
Hence the weight of the ray of u is

    A_h := (4/π)·μ(h)/φ(h) = (4/π)·Π_{p|h} (−1/(p−1))  (h squarefree),   A_h = 0  (p² | h, e.g. 7 + 24i, h = 25),   A_1 = 4/π.

Example h = 5, w = 3 + 4i = (2+i)²: Re(z w̄) = 3a + 4b ≡ 3·(a − 2b) ≡ 3·(z mod (2+i)) (mod 5), so 3a + 4b ≢ 0 (mod 5)
for every Gaussian prime other than the associates of 2 ± i, and A_5 = (4/π)(−1/4) = −1/π.
**Which (w, class) belongs to which ray:** the zone of u = w/h sees exactly one additive character of ℤ[i]/h,
ζ ↦ e(−Re(ζ w̄)/h), i.e. the Fourier coefficient of the prime density at the frequency u mod ℤ². Distinct rational
unit vectors are distinct mod ℤ² except the four axes, which all reduce to 0 (if ξ − ξ′ ∈ ℤ², |ξ| = |ξ′| = 1, then
ξ − ξ′ ∈ {0, ±1, ±i, ±1 ± i, ±2, ±2i} and the only rational solutions are axis pairs). For composite h, e.g.
h = 65, the 16 vectors (±33, ±56), (±56, ±33), (±63, ±16), (±16, ±63) are 16 rays, each with its own character and
the same weight A_65 = +(4/π)/48, and its own phase (w/65)^m. In the Hardy–Littlewood bookkeeping by moduli q of
primeside.md §4, the ray of 3 + 4i collects the dual vectors k·(3 + 4i) of every modulus q = 5k; the 2-D bookkeeping
by frequencies gives the single coefficient μ(h)/φ(h) at once.

**Lemma 3 (the main-term integral; proved).** With W_V = (1/2)∫ v^{−1/4}V(v) dv and Φ_M := ∫ F(y/M) e^{−iπy²} dy,

    I(X; M) = W_V·Φ_M·X^{3/4} + O_{V,F}(M³(m² + M²)·X^{1/4}),        Φ_M = e^{−iπ/4} + O_{N,F}(M^{−N})  for every N.

Proof. In the coordinates (t, y = s/√t), dA = √t dy dt; α = arctan(y/√t) is odd in y, and
|z| − t = √(t² + ty²) − t, V((t² + ty²)/X), F(y/M) are even, so ∫ e^{imα}(even) dy = ∫ cos(mα)(even) dy, and
cos(mα) = 1 + O(m²y²/t), √(t² + ty²) − t = y²/2 + O(y⁴/t), V((t² + ty²)/X) = V(t²/X) + O(y²/t) on |y| ≤ 2M,
t ≥ √(a₀X)/2. The leading term is ∫√t V(t²/X) dt · ∫F(y/M) e(−y²/2) dy = X^{3/4}W_V Φ_M; the error is
∫√t·M·(m²M² + M⁴)/t dt ≪ M³(m² + M²)X^{1/4}. For Φ_M: Φ_M = M∫F(v)e^{−iπM²v²}dv; stationary phase at v = 0 with
large parameter M² gives the full expansion e^{−iπ/4}Σ_k c_k F^{(2k)}(0) M^{−2k}; F ≡ 1 near 0 kills every k ≥ 1. ∎
(The Fresnel tail of a sharp cutoff 1_{|y|≤M} is ≍ 1/(πM), oscillating; of a Gaussian cutoff e^{−πy²/M²} it is
e^{−iπ/4}[(1 − i/M²)^{−1/2} − 1] = O(M^{−2}).) Computed (`check_arcs.py` C): |Φ_M − e^{−iπ/4}| = 5.0·10⁻², 1.0·10⁻³,
2.1·10⁻⁵ at M = 1, 2, 3; at X = 10⁷, M = 3, m = 4: I = 53296.86 − 53230.17i by two independent quadratures (agree to
10⁻⁹ relative), closed form W_VΦ_M X^{3/4} = 53264.68 − 53262.47i (relative difference 6·10⁻⁴).

**Hypothesis R_h (Ricci–Stucky in the residue classes mod h).** For every ζ ∈ (ℤ[i]/h)^*, every fixed φ ∈ ℝ and
c > 0, and I_x = (φ, φ + c x^{−1/4}] or (φ − c x^{−1/4}, φ]:

    Σ_{z ≡ ζ (mod h), N z ≤ x, arg z ∈ I_x} Λ(z) = (2/π)·|I_x|·x/Φ(h)·(1 + o(1))      (x → ∞).

It follows from any statement of the form "uniformly for sectors of width ≥ x^{−θ}, some θ > 1/4". R_1 for φ = 0 is
a case of Stucky's Theorem 1.1 (below), with conjugation z ↦ z̄ for the orientation (φ − δ, φ].

**Theorem A (axis zone; proved unconditionally).** For fixed m ≡ 0 (4), V, F, M as above,

    Z_1(X; M) = (4/π)·I(X; M) + o(X^{3/4}) = (4/π)·W_V Φ_M X^{3/4} + o(X^{3/4}),

and the same for the zones around i, −1, −i (they are the images under multiplication by a unit). The inert primes
z = q ≡ 3 (mod 4) lie on the axis; their contribution, Σ_q 2 ln q V(q²/X) ≍ √X, is part of the o(X^{3/4}).

**Theorem B (Pythagorean ray; proved assuming R_h).** Let u = w/h ∈ D. If R_h holds, then for fixed m, V, F, M

    Z_u(X; M) = A_h·u^m·I(X; M) + o(X^{3/4}),      A_h = (4/π) μ(h)/φ(h).

In particular Z_{(7+24i)/25}(X; M) = o(X^{3/4}).

**Theorem B′ (GRH; proved under GRH, sketch-level constants).** Assume GRH for L(s, ξ) for every Hecke
Grössencharacter ξ of ℚ(i) of modulus dividing h. Then Z_u(X; M) = A_h u^m I(X; M) + O_{m,M,V,F}(h·X^{1/2} log²(hX)).

**Elementary form (Theorem B for h = 5; conditional on R_5, proved under GRH).** Let r = |z|. Among Gaussian primes
z with N z ∈ [X, 2X] lying within M√r of the ray through 3 + 4i, the projection t = (3a + 4b)/5 of z on the ray
lies in (1/5)ℤ and is **never** ≡ 0 (mod 1) (Lemma 2 (iv); unconditional, algebra), and the four values
t ≡ 1/5, 2/5, 3/5, 4/5 (mod 1) each carry 1/4 + o(1) of the ln N z weight (R_5). Since |z| = t + s²/(2t) + O(M⁴/r),
it is |z| − s²/(2t), not |z| itself, that sits at k/5 mod 1 (|z| mod 1 is shifted by the Fresnel term s²/2t ∈
[0, M²/2]; the wording "the distance to the origin mod 1 avoids 0" is right for t, not for |z|). The
Fresnel-weighted consequence: Σ ln N z·e(−|z|)·V F_M over the zone = (mean of e(−k/5) over k = 1..4 = −1/4)·(4/π)·
u^m·I + o(X^{3/4}) = −(1/π) u^m W_V Φ_M X^{3/4} + o(X^{3/4}); for all lattice points (weight 1) t mod 1 is uniform
on all five values and the same zone sum is O_N(X^{−N}) (Poisson summation: the phase −|z| − ⟨z, n⟩ is stationary
only on the axes, and the weight is smooth on scales ≥ M^{−1}X^{1/4}; sketch).
For general h: t ∈ (1/h)ℤ, h·t = Re(z w̄) is equidistributed over (ℤ/h)^* (R_h), and the mean of e(−t) is μ(h)/φ(h).

## 3. Cited theorems (what is used, in what form)

| source | statement (as opened) | used? |
|---|---|---|
| **Stucky**, *Gaussian primes in narrow sectors*, arXiv:2008.11325 (final version 2021), **Theorem 1.1** (opened, LaTeX source) | ψ(x, y; φ, δ) = Σ_{x−y < N𝔞 ≤ x, φ < arg 𝔞 ≤ φ+δ} Λ(𝔞) (𝔞 ideals, arg 𝔞 defined mod π/2). "For any ε > 0, φ ∈ ℝ, x sufficiently large, θ > 7/10, and δx^θ ≥ x^{7/10+ε}, we have ψ(x, x^θ, φ, δ) ∼ 2δx^θ/π." Proof (§3): φ enters only through the Fourier coefficients c_m of the angular cutoff, and only \|c_m\| is used, so the bound is uniform in φ and in δ in the range. | **yes**: θ = 1, φ = 0, δ(x) = c·x^{−1/4} (c fixed; δx = c x^{3/4} ≥ x^{7/10+ε} for ε < 1/20). Only finitely many instances are needed, so no uniformity is used. |
| **Ricci**, *Local distribution of primes*, PhD thesis, Univ. of Michigan, **1976** (not 1933) (not opened; secondary: Stucky intro, Huang–Liu–Rudnick intro) | "for all ε > 0 and δ ≥ x^{−3/10+ε}, ψ(x, x; φ, δ) ∼ 2δx/π" (Stucky); HLR: (1.1) #{𝔭: X < N𝔭 ≤ 2X, θ_𝔭 ∈ I} ∼ (\|I\|/(π/2))·X/log X holds for \|I\| > X^{−3/10+ε}; zero density N(σ; T, K) ≪ T K^{(10/3)(1−σ)}(log K)^B for T = o(K) for the family L(s, Ξ_k), 0 < k ≤ K; Weyl bound L(1/2 + it, λ^m) ≪ (m² + t²)^{1/6} log³(m² + t²) (Stucky Lemma 6.1, "Chapter 2, Theorem 4"). | the statement is the θ = 1 case of Stucky's theorem; used through Stucky. Whether the thesis treats non-trivial moduli: **unknown** (not opened). |
| **Huang–Liu–Rudnick**, *Gaussian primes in almost all narrow sectors*, arXiv:1903.04005 (opened, LaTeX) | almost all sectors with \|I\| > X^{−ρ}, ρ < 3/5; quotes Ricci (3/10), GRH: \|I\| > X^{−1/2+ε}; Kubilius' zero-free region 1 − β ≫ (log V log log V)^{−3/4}; Harman–Lewis 2001: existence in narrower sectors by a sieve, "without an asymptotic expression"; Coleman (Proc. LMS 1990): small balls, individual results on GRH and unconditionally. | background; almost-all results do not give the specific boxes we need; Harman–Lewis gives no main term (useless for an asymptotic zone sum). |
| **Kubilius** 1950/1955 (not opened; secondary: HLR) | zero-free region above; existence/asymptotics in sectors of width x^{−c} with c < 1/4 (exponent not opened) | no |
| **Maknys**, Acta Math. Hungar. 42 (1983) (not opened; secondary: Stucky) | claimed polar boxes of area x^{11/16}; Heath-Brown's review found an error (corrected exponent 0.7349…, itself relying on a flawed zero-density estimate) | no |
| **Coleman**, *The distribution of points at which binary quadratic forms are prime*, Proc. LMS (3) 61 (1990) 433–456 (not opened; secondary: Stucky, Panidapu–Thorner) | hybrid large sieve Thm 6.2: Σ_{\|m\|≤M}∫_{−T}^{T}\|Σ_{N𝔞≍N} c(𝔞)λ^m(𝔞)N𝔞^{−it}\|² dt ≪ (M² + T² + N)Σ\|c\|² (arbitrary c); large values Thm 7.3; zero density Cor 8.2: N(σ, M, T) ≪ (M² + T²)^{(12/5+ε)(1−σ)}(log MT)^B, "proved in this level of generality" for Grössencharacters with non-trivial modulus (Panidapu–Thorner) | would be the input for R_h; gives only exponent 5/24 (below) |
| **Coleman**, *A zero-free region for the Hecke L-functions*, Mathematika 37 (1990) 287–304 (not opened; secondary) | σ ≥ 1 − C(log V)^{−2/3}(log log V)^{−1/3}, V = (4m² + t²)^{1/2} (Theorems 1, 2); general Grössencharacters per Panidapu–Thorner | as above |
| **Coleman**, *The distribution of points at which norm-forms are prime*, J. Number Theory 41 (1992) 359–378 | **not opened** (no abstract found) | possibly relevant (number fields, small regions); exponent unknown |
| Coleman, *Relative norms of prime ideals in small regions*, Mathematika 43 (1996) 40–62 (extract only) | set-up with Grössencharaktere mod a fixed ideal 𝔣 | not checked further |
| **Panidapu–Thorner**, *Short-interval sector problems for CM elliptic curves*, arXiv:2105.11093 (opened, LaTeX) | for the Grössencharacter of a CM curve (non-trivial conductor): if δ + θ < 5/24, x^{1−δ} ≤ h ≤ x, \|I\| ≥ x^{−θ}, then (1/h)Σ_{x<p≤x+h, θ_p∈I} log p ∼ μ(I), via Coleman's ZFR (Thm 2) and ZDE (Cor 8.2) | shows the **published exponent for non-trivial modulus is 5/24** (angles of primary primes on the full circle) |
| Kaufman 1978/79; Söhne, J. Number Theory 66 (1997) 225–250 (not opened; secondary: Watt arXiv:1307.5333 Lemma 10) | ζ(1/2 + it, λ^d) ≪ T^{1/6} log^C(T + 2) (Kaufman, Gaussian field, conductor 1); "similar sub-convexity bounds … for general number fields" (Söhne; whether non-trivial modulus and the hybrid (d, t) aspect are covered: not verified) | would close one step of the GAP (§5) |
| Iwaniec–Kowalski, *Analytic Number Theory*, Thm 5.15 (quoted from memory, not re-opened) | under GRH, ψ(f, x) = r x + O(√x log x log(x^d q(f))), q(f) the analytic conductor | Theorem B′ |
| Hecke 1918/1920 | Grössencharacters, functional equation, L(s, ξ) entire for ξ ≠ 1 | throughout |

**Area comparison.** Our boxes (§4) are polar boxes of norm length Y ≍ X/K and angle β = X^{−1/4}/K: area
≍ X^{3/4}/K². Short radial length is **not** needed (Y is a fixed fraction of X), so for h = 1 Stucky's θ = 1 case
(= Ricci, sectors of width ≥ x^{−3/10+ε}) suffices with room 1/20 in the exponent; the area condition
x^{7/10+ε} ≤ X^{3/4}/K² holds. For h > 1 the angle X^{−1/4} is *narrower* than the published threshold
x^{−5/24} for non-trivial modulus (5/24 < 1/4): that is the GAP.

## 4. Proofs

**Proposition (box lemma; proved).** Let h ≥ 1, ζ ∈ (ℤ[i]/h)^*, u = e^{iφ_u}, and suppose R_h holds at φ = φ_u
(for h = 1, φ_u = 0: Stucky). Then for G as in §1

    Σ_{z ≡ ζ (h)} Λ(z) G(z) = (4/π)/Φ(h)·∫ G dA + o(X^{3/4}).

Proof. Rescale z = r e^{i(φ_u + α)}, r = √X ρ, α = X^{−1/4} ā. Then G = γ(ρ, ā) + O(X^{−1/4}) uniformly with all
derivatives, where γ(ρ, ā) = e(−ρā²/2) V(ρ²) F(√ρ ā/M) is a fixed function in C_c^∞((0, ∞) × ℝ) (checked: r(1 − cos α)
= ρā²/2 + O(X^{−1/2}); √r sin α/√cos α = √ρ ā + O(X^{−1/2}); e^{imα} = 1 + O(mX^{−1/4})). Fix an integer K and cut the
(ρ, ā) plane into boxes [j/K, (j+1)/K) × [l/K, (l+1)/K); in z they are polar boxes
B_{jl} = {N z ∈ [X j²/K², X(j+1)²/K²), arg z − φ_u ∈ [lβ, (l+1)β)}, β = X^{−1/4}/K. Only C(M, V)·K² boxes meet
supp γ. On each, |G(z) − γ(j/K, l/K)| ≤ η_K + O(X^{−1/4}), η_K := 2‖∇γ‖_∞/K. Write P_ζ(B) := Σ_{z∈B, z≡ζ} Λ(z).
Each B is a difference of at most four sectors (φ_u, φ_u + c x^{−1/4}] or (φ_u − c x^{−1/4}, φ_u] at the two norms
x = Xj²/K², X(j+1)²/K², with finitely many constants c (c = ρ^{1/2}·l/K); R_h gives
P_ζ(B) = (4/π)Φ(h)^{−1}·area(B) + o_K(X^{3/4}) (area(B) = (β/2)·(norm length) ≍ X^{3/4}/K², and the errors o(|I_x|·x)
are o(X^{3/4}) because x ≍ X). Boundary rays carry Λ-mass ≤ √X log X. Hence

    Σ_{z≡ζ} Λ G = Σ_B γ(j/K, l/K) P_ζ(B) + E,   |E| ≤ (η_K + o(1))·Σ_B P_ζ(B) ≤ C η_K X^{3/4} + o(X^{3/4}),

and Σ_B γ(j/K, l/K)(4/π)Φ^{−1}area(B) = (4/π)Φ^{−1}∫G dA + O(η_K X^{3/4}) (Riemann sum, dA = X^{3/4}ρ dρ dā).
So limsup_X X^{−3/4}|Σ_{z≡ζ}ΛG − (4/π)Φ^{−1}∫G| ≤ C′η_K for every K; let K → ∞. ∎

**Proof of Theorem A.** h = 1, w = 1, e(−Re(z w̄)) = e(−a) = 1 exactly. By Lemma 1, Z_1 = Σ_z Λ(z)G(z) (all z with
a > 0 near the axis). Boxes with l ≥ 0 are differences of sectors (0, δ], δ = c x^{−1/4}: Stucky's Theorem 1.1 with
θ = 1, φ = 0 (for a sector of width ≤ π/2 around the axis, elements z with arg z ∈ (0, δ] are in bijection with
ideals with arg 𝔞 ∈ (0, δ] mod π/2). Boxes with l < 0 are mapped onto boxes with l ≥ 0 by z ↦ z̄, which preserves
primality, norm and Λ. The Proposition gives Z_1 = (4/π)∫G + o(X^{3/4}); Lemma 3 evaluates ∫G. ∎

**Proof of Theorem B.** By Lemma 1, Z_u = u^m Σ_{ζ mod h} e(−Re(ζ w̄)/h) Σ_{z ≡ ζ} Λ(z)G(z), exactly. Classes not
prime to h contain only associates of divisors of h (norm ≤ h² < a₀X for large X): no contribution. The Proposition
(under R_h) gives each class (4/π)Φ(h)^{−1}∫G + o(X^{3/4}), and Lemma 2 (iii) sums the phases:
Z_u = u^m (4/π)(μ(h)/φ(h))∫G + o(X^{3/4}). ∎

**Proof of Theorem B′ (under GRH; constants not optimised).** Expand G in angular Fourier series,
G(re^{iφ}) = Σ_k g_k(r) e^{ikφ}, and the class weight e(−Re(ζw̄)/h) = Σ_χ ĉ(χ)χ(ζ) over the characters of
(ℤ[i]/h)^*; Parseval gives Σ_χ|ĉ(χ)|² = 1, so Σ_χ|ĉ(χ)| ≤ Φ(h)^{1/2} ≤ h. Summing over the four associates,
Σ_z Λ(z)χ(z)e^{ik arg z}g_k(|z|) vanishes unless χ(ε)ε^k = 1 for the units, and then equals
4Σ_𝔭 Λ(𝔭)ξ(𝔭)g_k(√N𝔭) with ξ = χ·(·/|·|)^k a Grössencharacter of modulus h, analytic conductor ≪ h²(|k| + 3)². Under
GRH, ψ(y, ξ) = δ_ξ y + O(y^{1/2} log y log(y h(|k|+3))) (IK Thm 5.15, degree 2). Partial summation:
|Σ_𝔭 Λξ g_k(√N𝔭) − δ_ξ∫g_k(√y)dy| ≪ X^{1/2} log² (hX)·(‖g_k‖_∞ + TV(g_k)). The φ-support of G at radius r has
width ≤ 4M/√r and its φ-derivatives are ≪ (M + 1/M + m/√r)√r per order, so Σ_k ‖g_k‖_∞ ≪ M² + 1 + m; the
r-derivative of G is ≪ (M² + 1)/r, so Σ_k TV(g_k) ≪ (M² + 1)(M² + 1 + m). The only pole is ξ = 1 (χ = χ₀, k = 0),
giving 4∫g_0(√y)dy = (4/π)∫G dA. Total error ≪ h(M² + 1)(M² + m + 1)X^{1/2} log²(hX). ∎ (The same argument proves
R_h for all sectors of width ≥ x^{−1/2+ε}.)

## 5. The GAP for h > 1 (exactly what is missing)

**Missing input:** R_h for h > 1: an asymptotic count of Gaussian primes in a fixed residue class mod h in sectors of
angular width c·x^{−1/4} and norm ≤ x. Equivalently: Ricci's theorem (any exponent θ > 1/4) for the Hecke
Grössencharacters ξ = χ·(·/|·|)^k with χ a character mod h, all k. For h = 5 this is very concrete: since the units
map onto (ℤ[i]/(2+i))^* (i ≡ 3, of order 4 mod 5), each prime ideal 𝔭 ∤ 5 has exactly one generator z_𝔭 ≡ 1
(mod 2 + i), and by Lemma 2 the class of the zone generator is the quarter-turn position of z_𝔭; so R_5 says: **the
angles mod 2π of the "primary" generators z_𝔭 ≡ 1 (mod 2+i) are equidistributed in sectors of width c x^{−1/4}**
(Ricci/Stucky give this for angles mod π/2 only, i.e. for the frequencies k ≡ 0 mod 4).

**By how much it is missing.** Published asymptotics for Grössencharacters with non-trivial modulus reach angular
width x^{−θ} for θ < 5/24 = 0.2083 (Panidapu–Thorner via Coleman's zero density with exponent 12/5 in (M² + T²)); we
need θ > 1/4 = 0.25. With a zero-density estimate N(σ; K) ≪ K^{A(1−σ)+ε} for the family {ξ k}_{|k|≤K} at bounded
height, the Hoheisel argument gives θ < 1/A: Coleman's A = 24/5 gives 5/24; we need A < 4; Ricci's A = 10/3 (modulus 1)
gives 3/10.

**How to close it (not done; expected routine, not verified).** Stucky's proof (θ = 1 case) uses: (a) Coleman's hybrid
large sieve and large-values theorem — for *arbitrary* coefficients c(𝔞), so a fixed character ξ₀ = χ·(·/|·|)^{r} mod h
is absorbed into the coefficients (in the Halász–Huxley duality the kernel ⟨λ^{m}N^{−it}, λ^{m′}N^{−it′}⟩ never sees
ξ₀, because |ξ₀(𝔞)| = 1 on ideals prime to h) — no change; both signs ±m must now be kept (L(s, ξ₀λ^m) ≠ L(s, ξ₀λ^{−m})); (b) the Heath-Brown
identity — valid for Λξ, μξ since ξ is completely multiplicative — no change; (c) Coleman's zero-free region — stated
for general Grössencharacters (per Panidapu–Thorner) — no change; (d) the main term — prime ideal theorem in ray
classes; (e) **Ricci's Weyl-type bound L(1/2 + it, λ^m) ≪ (m² + t²)^{1/6} log³** (Stucky Lemma 6.1) used for the long
smooth factors, whose coefficients cannot absorb ξ₀: needed as L(1/2 + it, ξ₀λ^m) ≪_h (m² + t²)^{1/6+ε}. Its proof
(2-D exponential sums over lattice points via the approximate functional equation) should go through on the cosets
ζ + hℤ[i] with h-dependent constants; Söhne 1997 may contain it (not opened). Also: Ricci's thesis and Coleman's 1992
norm-form paper were not opened and might already contain R_h. **Until (e) is written or found, Theorem B is
conditional (on R_h, or on GRH via Theorem B′).**

## 6. Corollary 1: the singular series is the sum over the rays (proved; computed)

Put g(h) := μ(h)Π_{p|h}2cos(2mθ_p)/φ(h) = Π_{p|h}(−2cos(2mθ_p)/(p−1)) for squarefree h with all p | h ≡ 1 (mod 4),
g(h) = 0 otherwise, θ_p = arg π_p (π_p = a + bi, a² + b² = p).

(a) **Term by term (proved, algebra).** Σ_{ξ ∈ D, h(ξ) = h} ξ^m = Σ_ε Σ_{π₀} ε^m (π₀/|π₀|)^{2m} = 4·Π_{p|h} 2cos(2mθ_p)
(m ≡ 0 mod 4; π₀ over the 2^{ω(h)} primitive elements of norm h up to units; each prime contributes e^{±2imθ_p}).
Hence

    Σ_{ξ ∈ D, h(ξ) ≤ H} A_{h(ξ)} ξ^m = (4/π)·4·S_H,      S_H := Σ_{h ≤ H} g(h),

and the h = 1 term is the four axes. Checked: ray enumeration = multiplicative formula for every h ≤ 10⁴ (max
difference 6·10⁻¹⁷), and independently from 316 Pythagorean triples (`check_arcs_indep.py`, S_1000 = 1.587324).

(b) **The limit (proved modulo a standard step).** G(s) := Σ g(h)h^{−s} = Π_{p≡1(4)}(1 − λ_p/((p−1)p^s)),
λ_p = 2cos 2mθ_p, and G(s)·L(s + 1, ψ_{2m}) = E(s) with E an Euler product absolutely convergent in Re s > −1/2
(each factor is 1 + O(p^{−σ−2} + p^{−2σ−2}); ψ_{2m}((1+i)) = i^m = 1, ψ_{2m}((q)) = 1). L(s, ψ_{2m}) is a Hecke
L-function, entire, with the classical zero-free region and 1/L ≪ log(|t| + 3) there; Perron's formula and a contour
in that region give S_H = G(0) + O(exp(−c√log H)) — the argument that proves Σμ(n)/n = 0, written for L(s, ψ_{2m}) (not
written out here: **standard step**). Abel's theorem and the prime ideal theorem for ψ_{2m} (Σ_p g(p) converges) give
G(0) = Π_{p≡1(4)}(1 + g(p)) in increasing p. So

    lim_{H→∞} S_H = Π_{p ≡ 1 (4)} (1 − 2cos(2mθ_p)/(p − 1)) = (π/4)·c,

and the sum over all rays of the theorem's main terms is 4c·I(X; M) → 4c W_V e^{−iπ/4} X^{3/4} (M → ∞): exactly
Conjecture C. The convergence is **conditional**: Σ|g(h)| = Π(1 + 2|cos|/(p−1)) diverges (Σ 1/p over p ≡ 1 mod 4).

(c) **Numbers (computed, m = 4; `series_conv.py`).** (π/4)c_exact = 1.591567 (c_exact = 2.026446 from
constants.json, B/L(1, Ad)).

| H | S_H | S_H − (π/4)c | Π_{p≤H}(1 + g(p)) |
|---|---|---|---|
| 10² | 1.567841 | −2.4·10⁻² | 1.603482 |
| 10³ | 1.587324 | −4.2·10⁻³ | 1.587433 |
| 10⁴ | 1.587387 | −4.2·10⁻³ | 1.588835 |
| 10⁵ | 1.591130 | −4.4·10⁻⁴ | 1.591115 |
| 10⁶ | 1.591595 | +2.8·10⁻⁵ | 1.591484 |
| 2·10⁶ | 1.591732 | +1.7·10⁻⁴ | 1.591541 |

The product over p ≤ 10⁷ (from the Gaussian primes of `check_arcs.py`) is 1.591505 (6·10⁻⁵ from (π/4)c_exact).
The fluctuations of S_H are of size ≈ H^{−1/2} (as from the zeros of L(s, ψ_8) on Re s = 1/2), so the plateau
S_{10³} ≈ S_{10⁴} is a coincidence, not a different limit.

## 7. Corollary 2 and the open remainder (the minor arcs)

For H, M fixed and X large the zones of the directions with h ≤ H are pairwise disjoint (two such directions differ in
angle by ≥ 1/(hh′) ≥ 1/H², the zones have half-width ≤ 2.1M(a₀X)^{−1/4}). Define the minor-arc sum

    R_{H,M}(X) := Σ_z Λ(z)ψ_m(z)e(−|z|)·[V(N z/X) − Σ_{ξ ∈ D, h(ξ) ≤ H} V(N z/X)F_M(s_ξ/√t_ξ)]     (weight ≥ 0).

**Open statement (MA):** lim_{H, M → ∞} limsup_{X→∞} X^{−3/4}|R_{H,M}(X)| = 0.

**Corollary 2 (proved).** Assume R_h for every h (e.g. GRH). Then Conjecture C ⟺ (MA). Proof: Σ_full − 4cW_V
e^{−iπ/4}X^{3/4} = R_{H,M} + [Σ_{h(ξ)≤H} Z_ξ − (16/π)S_H I] + (16/π)S_H(I − W_VΦ_M X^{3/4}) +
(16/π)W_V X^{3/4}(S_HΦ_M − (π/4)c e^{−iπ/4}); the second and third terms are o(X^{3/4}) (Theorems A, B; Lemma 3), the
last is ε(H, M)X^{3/4} with ε → 0 (Corollary 1, Lemma 3). So |limsup|Σ_full − main|/X^{3/4} − limsup|R_{H,M}|/X^{3/4}|
≤ ε(H, M) for all H, M. ∎ Unconditionally, Theorem A alone reduces nothing: the rays h ≥ 5 are part of the
remainder until R_h is proved.

**What is known about the remainder.**
- Trivial: |R_{H,M}| ≤ Σ Λ(z)V ≍ X — a factor X^{1/4} too large.
- Under GRH for L(s, ψ_m) alone, the explicit formula bounds the **full** sum Σ_full by the number of zeros in the
  stationary window (|γ| ≈ π√(vX), ≈ √X log X zeros, each contributing ≈ X^{1/2}·X^{−1/4}): O(X^{3/4} log X) — misses
  by a logarithm (sketch; the same statement for elliptic curves is in STUDY.md round 148). The asymptotic needs the
  zeros' correlation with the phase (fact 2055: GRH + the one-over-one ratios conjecture).
- Unconditional power savings for the classical analogue Σ_p λ_f(p) log p e(α√p): x^{5/6} (Zhao 2006, level 1, as
  quoted in STUDY.md; not re-opened) — x^{1/12} short of x^{3/4}.
- Why Vinogradov-type bounds fall short (heuristic reasoning about the method): on a box of the minor arcs the phase
  e(−|z|) is a linear character e(−Re(z ξ̄)) plus a quadratic term, with ξ the local direction; approximating ξ by a
  rational direction of height q ≤ Q, a Vinogradov/Mitsui-type bound saves at most q^{−1/2} on that box. Even with
  perfect savings and q up to X^{1/4} (beyond that the zones of width X^{−1/4} no longer resolve directions) the
  saving is ≥ X^{−1/8} per box, while (MA) needs X^{−1/4} on average: square-root cancellation across boxes is
  required. Note also that R_{H,M} still contains the Hardy–Littlewood main terms of all rays with h > H (whose sum
  converges only conditionally, Corollary 1), so (MA) mixes "tails of the singular series" with genuine minor arcs.
- Under GRH (Theorem B′) the rays are handled uniformly up to h ≤ X^{1/8−ε} (error Σ_{h≤H} h·2^{ω(h)}X^{1/2+ε}); the
  zones stay disjoint up to about the same height (1/H² ≥ 4M X^{−1/4}). This is the major-arc range a circle-method
  proof under GRH would start from; the complement is (MA) with H = X^{1/8−ε}.

## 8. Numerical check (computed; `check_arcs.py`, output `check_arcs_out.txt`, 85 PASS, 140 s, < 300 MB)

Setting: all Gaussian prime elements with Re z > 0 and N z ≤ X = 10⁷ (1,328,947 elements); m = 4 (and m = 0);
M = 3; F(y) = 1 on |y| ≤ 1, 0 on |y| ≥ 2 (C^∞ step built from e^{−1/x}); V(v) = smooth bump with support [0.05, 1],
flat on [0.30, 0.75] (the theorem allows any V ∈ C_c^∞((0, ∞)); a bump on [1, 2] would use only N z ∈ [X, 2X] and give
a weaker test at this height). sd = random-phase noise per component, √(Σ(ΛW)²/2); 88 control sums with e(−c|z|),
c ∈ {0.6713, 0.8219, 1.1597, 1.3711}, have rms/sd 0.89 (real), 0.74 (imag).

**Zone sums, m = 4** (Z without the inert primes; "own" = Theorem's main term A_h u^m I; "nb" = heuristic
Hardy–Littlewood contribution of the other rational rays with h′ ≤ 1000 passing through the zone,
Σ A_{h′}∫ψ_m e(Re(zξ̄′) − |z|)W_u dA; it is not part of the theorem):

| zone | h | A_h | Z | own (theorem) | (Z−own)/sd | nb | (Z−own−nb)/sd | sd |
|---|---|---|---|---|---|---|---|---|
| axis | 1 | +1.27324 | 66816 − 65624i | 67860 − 67775i | 1.00 | 1004 − 936i | 1.55 | 2395 |
| 3+4i | 5 | −0.31831 | 29464 − 6410i | 23414 − 5167i | 2.57 | 4729 − 951i | 0.56 | 2399 |
| 4+3i | 5 | −0.31831 | 5094 − 29568i | 5196 − 23407i | 2.57 | 1104 − 4774i | 0.77 | 2399 |
| 5+12i | 13 | −0.10610 | 7846 + 8592i | 5695 + 5608i | 1.53 | 3809 + 3149i | 0.69 | 2400 |
| 12+5i | 13 | −0.10610 | −8135 − 9499i | −5600 − 5702i | 1.90 | −3067 − 3783i | 0.22 | 2400 |
| 8+15i | 17 | −0.07958 | 10730 + 4063i | 5528 + 2318i | 2.29 | 4528 + 2620i | 0.46 | 2398 |
| 15+8i | 17 | −0.07958 | −5153 − 10436i | −2311 − 5531i | 2.36 | −2679 − 4638i | 0.13 | 2398 |
| 7+24i | 25 | 0 | 2027 + 6398i | 0 | 2.80 | 2834 + 5695i | 0.45 | 2399 |
| 24+7i | 25 | 0 | −5130 − 4023i | 0 | 2.72 | −5946 − 2574i | 0.69 | 2399 |
| 20+21i | 29 | −0.04547 | 8358 − 7253i | 2648 − 2173i | 3.18 | 6847 − 5120i | 0.47 | 2401 |
| 21+20i | 29 | −0.04547 | 5884 − 9357i | 2176 − 2645i | 3.19 | 4147 − 7223i | 0.28 | 2401 |

m = 0 (3+4i and 4+3i coincide by z ↦ i z̄, etc.): axis 66773 − 65516i vs own 67818 − 67816i (1.05 sd); 3+4i
−20930 + 21935i vs own −16955 + 16954i (2.66 sd; with nb 0.61 sd); 5+12i 1.87 → 0.48; 8+15i 2.48 → 0.27;
7+24i 2.84 → 0.54; 20+21i 3.37 → 0.37 sd. Inert primes in the axis zone: +1645 (real); with them the axis reads
68461 − 65624i, 0.93 sd from the theorem. Pooled amplitude of the 8 Pythagorean zones with A_h ≠ 0 (m = 4, Z regressed
on own + nb): 1.014 ± 0.048. The nb model is not converged in h′: |nb| at cuts h′ ≤ 100 / 300 / 1000 is 1.29 / 1.61 /
2.01 sd for 3+4i, 4.24 / 4.05 / 3.56 sd for 20+21i (so the model has a systematic uncertainty of ~0.5 sd).

**Reading.** The axis zone agrees with Theorem A at 1 sd. For h ≥ 5 the own term alone is 1.5–3.4 sd off at
X = 10⁷: with M = 3 the zones (half-width 0.11–0.23 rad) still contain other rational rays (for 3+4i: 20+21i,
8+15i, 33+56i, …) whose resonance the theorem only makes o(X^{3/4}) as X → ∞. Adding them heuristically brings every
zone within 0.13–0.77 sd. How fast the neighbours leave (model only, h′ ≤ 1000, no primes; section F of the output):

| zone | \|nb\|/\|own\| at X = 10⁷ | 10⁹ | 10¹¹ | 10¹³ |
|---|---|---|---|---|
| axis | 0.014 | 0.000 | 0 rays | 0 rays |
| 3+4i | 0.201 | 0.012 | 0.000 | 0 rays |
| 5+12i | 0.618 | 0.046 | 0.001 | 0.000 |
| 20+21i | 2.496 | 0.048 | 0.002 | 0.000 |
| 7+24i (own = 0): \|nb\|/X^{3/4} | 0.036 | 0.0037 | 0.0001 | 0.0000 |

So the theorem's regime (own term dominant) begins around X ≈ 10⁹ for these zones at M = 3; at 10⁷ the data test
"theorem + neighbouring rays", not the theorem alone.

**The class statement (the actual input R_h), flat part |s| ≤ M√t, N z ∈ [0.05X, X]:**

| ray | h | primes | Re(z w̄) ≡ 0 mod p \| h | χ² on (ℤ/h)^* (dof) | mean e(−Re(zw̄)/h) (Λ-weighted) | μ(h)/φ(h) | sd |
|---|---|---|---|---|---|---|---|
| 3+4i | 5 | 53496 | 0 | 0.6 (3) | −0.2510 − 0.0019i | −0.2500 | 0.0031 |
| 5+12i | 13 | 53586 | 0 | 5.6 (11) | −0.0829 − 0.0011i | −0.0833 | 0.0031 |
| 8+15i | 17 | 53532 | 0 | 4.2 (15) | −0.0600 − 0.0026i | −0.0625 | 0.0031 |
| 7+24i | 25 | 53528 | 0 | 6.2 (19) | +0.0009 − 0.0007i | 0 | 0.0031 |
| 20+21i | 29 | 53776 | 0 | 11.6 (27) | −0.0378 − 0.0022i | −0.0357 | 0.0031 |

Lattice points in the same 3+4i zone: 3a + 4b mod 5 has shares 0.20003, 0.19999, 0.19999, 0.19999, 0.20000 on
0, 1, 2, 3, 4. Every χ² is *below* its degrees of freedom (z-scores −0.98 … −2.10): prime counts in residue classes
are under-dispersed relative to Poisson (as in the Hooley–Montgomery variance, ≈ (log h/log X) of Poisson; heuristic
explanation, not tested further).

**Independent check** (`check_arcs_indep.py`, plain Python and sympy, different enumeration): Lemma 2 for h = 5, 13, 25,
65; the axis and 3+4i zone sums (68460.81 − 65623.82i and 29463.81 − 6409.83i: the same as `check_arcs.py` to the
printed digits); S_1000 = 1.587324 from Pythagorean triples. All pass.

## 9. Novelty verdict (searches with mathematical terms only; two web searches in this round plus primeside.md §6)

| statement | verdict | basis |
|---|---|---|
| Gaussian primes in narrow sectors (Theorem A's input) | **KNOWN** | Ricci 1976; Stucky (arXiv:2008.11325) |
| Theorem A (axis zone with the radial phase) | **NOT FOUND** as stated; a direct consequence of Ricci/Stucky (KNOWN IN SUBSTANCE once Lemma 1 is seen) | no source with e(\|z\|) found |
| Gaussian primes and Pythagorean directions: the weight A_h = (4/π)μ(h)/φ(h) of the ray of a primitive triple, zero for non-squarefree hypotenuse | **NOT FOUND** | searches: no hit on rays/resonance |
| equidistribution of \|z\| mod 1 for Gaussian primes | **KNOWN** in the trivial sense (\|z\| = √p; √p mod 1 is equidistributed over primes: classical, Vinogradov; refinements Balog 1983 — from memory, not re-opened) | the zone statement (t mod 1 on {k/h}, never 0, equidistributed) is NOT FOUND, and its equidistribution part is R_h (open unconditionally) |
| radial exponential sums over Gaussian primes Σ Λ ψ_m e(−\|z\|) | **KNOWN IN SUBSTANCE** as the prime resonance of a CM form (Iwaniec–Luo–Sarnak heuristic; integer version Kaczorowski–Perelli) | the major-arc decomposition by rays and the identity c = Σ_rays (Corollary 1) NOT FOUND |
| angles of Grössencharacters with modulus in short sectors | **KNOWN** up to 5/24 (Panidapu–Thorner via Coleman) | R_h at exponent > 1/4: NOT FOUND |

## 10. What a referee should check first

1. The application of Stucky's Theorem 1.1 at θ = 1, φ = 0, δ = c x^{−1/4} (hypotheses: θ > 7/10, δx^θ ≥ x^{7/10+ε});
   the element/ideal bijection for sectors inside (−π/4, π/4]; the conjugation trick for the negative side.
2. Lemma 2 (i): 2Re(ζw̄) ≡ ζw̄ (mod π_p^k) because π_p^{2k} | w; and that the zone of u sees exactly the character
   ζ ↦ e(−Re(ζw̄)/h) (Lemma 1 is exact).
3. Whether R_h (h > 1) is in fact in Ricci's thesis, Coleman 1992 (J. Number Theory 41) or Söhne 1997; if not,
   step (e) of §5 (Weyl bound for L(1/2 + it, ξ₀λ^m), fixed modulus) is the one real piece of work left.
4. Theorem B′: the bounds Σ_k‖g_k‖ ≪ M² + 1 + m and Σ_k TV(g_k) ≪ (M² + 1)(M² + m + 1), and the use of IK Thm 5.15
   uniformly in the infinity type k (q(ξ) ≍ h²k²).
5. Corollary 1 (b): the Perron/zero-free-region step (standard, not written out) and Abel's theorem for the order of
   the product.
6. That the numerics at 10⁷ test "theorem + neighbouring rays", not the theorem alone (§8, table of |nb|/|own|).

## 11. Doubts

1. Theorem B is unconditional only for h = 1. For h > 1 it rests on R_h (GRH suffices); the gap is the exponent 1/4
   versus the published 5/24 for non-trivial modulus. I did not open Ricci's thesis, Coleman 1992, or Söhne 1997.
2. The uniformity of Stucky's Theorem 1.1 in φ is read from the proof; it is not needed (φ = 0 and finitely many δ(x)).
3. Theorem B′'s constants (powers of M, the factor h) are crude; IK Thm 5.15 is quoted from memory.
4. Corollary 1 (b) uses a standard argument that is sketched, not written in full.
5. The neighbour model of §8 is heuristic (Hardy–Littlewood frequencies + quadrature) and not converged in h′ at the
   0.5 sd level; the agreement 0.13–0.77 sd should not be read as a precision test.
6. The GRH remark on the full sum (O(X^{3/4} log X)) is a sketch; the x^{5/6} bound (Zhao 2006) is for level-1
   holomorphic forms and was not re-opened; its CM analogue is assumed, not checked.
7. The caller's wording "the distance to the origin mod 1 avoids 0" holds for the projection t, not for |z|.
8. The caller's S_m "over Gaussian primes" is 4 × the ideal sum P_m ~ K x^{3/4}; all constants here are per element.
9. Ricci's result is from 1976 (thesis), not 1933.

## 12. Files (scratchpad/arcs/)

- `arcs_theorem.md` — this note.
- `check_arcs.py`, `check_arcs_out.txt` — sections A (Lemma 2), B (identity), C (Fresnel, ∫G), D (zone sums),
  F (cross-talk vs X), E (class statement); 85 PASS.
- `check_arcs_indep.py`, `check_arcs_indep_out.txt` — independent recomputation.
- `series_conv.py`, `series_conv_out.txt` — S_H to 2·10⁶.
- Sources opened (LaTeX, in their own folders): `../src_2008.11325/` (Stucky), `../src_1903.04005/` (Huang–Liu–Rudnick),
  `../src_2105.11093/` (Panidapu–Thorner), `../src_1307.5333/` (Watt).

## Referee report (25 November 2026)

An independent referee (own code in [referee/](referee/)) checked the note. **Confirmed as stated:** Lemma 1 (exact factorisation; in box
coordinates G is a fixed smooth function up to O(m·X^(−1/4))), Lemma 2 (proved again for h = p, p², pq; brute force for h = 5, 13, 25,
65, 85; the zone weight depends only on the class mod the ideal of norm h), Lemma 3, Theorem B (conditional on R_h), Theorem B′
(under GRH; constants at sketch level), Corollary 1 (the Perron step needs the classical zero-free region for L(s, ψ_2m), not only
L(1, ψ_2m) ≠ 0) and Corollary 2. check_arcs.py has 86 PASS lines, not 85.

**Theorem A** is a correct deduction from the published sector theorem (Stucky, Q. J. Math. 72 (2021) 1357–1377, Theorem 1.1 with
θ = 1, which is Ricci's theorem; hypotheses checked against the source: θ = 1 allowed, sector width c·x^(−1/4) admitted with room
1/20, any anchor angle, sharp cutoffs, prime-power ideals; the theorem gives no rate). Corrections: the theorem is applied with a
width depending on x, which its statement allows; each ideal has exactly one generator in the zone.
**A caveat on the cited proof, reported by the referee and not verified further here:** in the cited paper the mean-value
lemma for the case T ≤ M (attributed there to Ricci's thesis) is stated for arbitrary coefficients, and the referee gives a
numerical counterexample to that generality (coefficients supported on rational integers); the proof uses it only for
divisor-bounded coefficients, for which the referee expects the bound to hold but did not find it written. Theorem A is
therefore "a theorem from the literature (Ricci; Stucky)", with this remark attached.

**The gap R_h** (primes in a residue class in sectors of width x^(−1/4)): not found in print (in print: exponent 5/24 through a
zero-density estimate). The referee's analysis of the cited method: a Weyl-strength bound is not needed; with a pointwise bound
(m² + t²)^κ for the twisted Hecke L-functions the method reaches the exponent 1/(2 + 8κ), so any subconvex bound κ < 1/4 suffices,
and one is in print for all Hecke characters (Michel–Venkatesh 2010); the other steps carry over with the character absorbed in
the coefficients, with the same dependence on the mean-value lemma as Theorem A. So R_h holds at the level of a proof sketch; it
is not a published result and is not claimed here as a theorem.


## Second referee report: the hypothesis R_h (25 November 2026)

R_h is no longer a gap of method. [rh/rh_proof.md](rh/rh_proof.md) proves it from two inputs: a pointwise bound for
L(1/2 + it, ξ) in the full analytic conductor (Michel–Venkatesh, read; Wu, read, with an explicit exponent) and the
zero-free region with growth bounds for Hecke L-functions with modulus (Coleman, Mathematika 37 (1990), Theorems 1 and 2;
not opened: Theorem 2 is quoted word for word with the modulus in Sun–Zhang, arXiv 2008.09677, and Theorem 1 is described
with the modulus in its zbMATH review). An independent referee ([rh/referee/referee_rh.md](rh/referee/referee_rh.md), own
code) re-derived the exponent 1/(2 + 8κ), checked every lemma and found no failing step; R_h comes out in the exact form
used in section 2 (norm ≤ x, sharp sector, fixed anchor, uniform in the anchor and in the width). Consequences for this note:
- **Theorem B** is conditional on that citation, not on GRH. Theorem B′ remains the only version with a power saving.
- **Section 5 (e)** (Ricci's Weyl bound for twisted L-functions) is not needed: any subconvex bound in the full conductor suffices.
- **Section 3, the mean-value lemma of the cited proof** (Stucky, Lemma 4.2): false as printed, in both of its ranges (for
  T ≤ M: c = 1 on the ideals generated by rational integers gives a factor of order √N (2.06·√N at T = 1, M = 10N, fact 2058); for T > M: c = 1 on the ideals of one
  norm). For the coefficients that occur in the proof it holds up to powers of the logarithm (rh_proof.md, Lemmas 5 and 6),
  so the published theorem and Theorem A stand.
