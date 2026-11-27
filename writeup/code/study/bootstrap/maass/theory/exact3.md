# The odd real character problem: q* ≥ 3 − ε from Ξ_χ²·B, unconditionally (27 November 2026)

**In short.** Problem of [../maass.md](../maass.md), sections 2-3: Gamma factor Γ_R(s+1), no pole, conductor q, one probability measure μ_p on [−1, 1] per prime, coefficient ⟨x^m⟩_p log p at p^m. Take ĝ = Ξ_χ(t)²·B(t), χ = χ₋₃, B ≥ 0 an even polynomial.

1. **Proved.** For this g the explicit formula of L(s, χ₋₃) holds with zero side identically 0. No hypothesis on the zeros is needed: Λ·Λ′·B is entire and odd under s ↦ 1 − s (Lemma 1).
2. **Proved.** A family positive at conductor q for all f ∈ C_c^∞(ℝ) is positive on g (Lemma 3). Hence log q ≥ log 3 − 2·Defect(g)/g(0) (Theorem 4).
3. **Computer-assisted (one code, arb ball arithmetic).** For the best B found at each degree, q*_∞ ≥ 2.7400 (B = 1), 2.99731 (degree 4), 3 − 3.1·10⁻⁶ (degree 8), 3 − 4.5·10⁻¹⁰ (degree 12), 3 − 9.1·10⁻¹⁷ (degree 16) and 3 − 3.81·10⁻²² (degree 22). Defect/g(0) is already below 10⁻⁶ at degree 8 (5.14·10⁻⁷). Table in section 5.
4. **Not proved.** Exactness, q* = 3. It needs Defect → 0 along some sequence of B. Computed: Defect falls, with stalls (section 6). The upper bound q* ≤ 3 needs GRH for L(s, χ₋₃).

## 1. Normalisation (that of [../maass.md](../maass.md), section 2, with κ = 1)

Test functions f ∈ C_c^∞(ℝ), complex valued; g(w) = ∫ f(v+w) f̄(v) dv; ĝ(t) = ∫ g(u) e^{itu} du = |f̂(t)|² ≥ 0. For a family μ = (μ_p) of Borel probability measures on [−1, 1], one per prime, and q > 0:

`Q_{μ,q}(f) = (log q) g(0) + A₁(g) − 2 Re Σ_p Σ_{m≥1} ⟨x^m⟩_{μ_p} (log p) p^{−m/2} g(m log p)`

`A₁(g) = −(γ + log π) g(0) + ∫₀^∞ [2e^{−2u} g(0) − e^{−3u/2}(g(u) + g(−u))] / (1 − e^{−2u}) du`

This is the repo's form with Gamma factor Γ_R(s+1), c(p^m) = ⟨x^m⟩ log p; for complex f the two prime sums of maass.md combine to 2 Re. Q depends only on Re g, and Re(f*f̃) = U*Ũ + V*Ṽ for f = U + iV, so positivity for all real f and for all complex f are the same condition.

**Definitions.** μ is *feasible at q* if Q_{μ,q}(f) ≥ 0 for every f ∈ C_c^∞(ℝ), and *feasible at q and support X* if this holds for every f supported in [−a, a], X = e^{2a}. q*(X) = inf of the feasible q at support X; q*_∞ = inf of the feasible q. Feasibility at q implies feasibility at every q′ > q (the term (log q) g(0) grows, g(0) = ‖f‖² ≥ 0), and feasibility at support X′ implies it at X < X′. So q*(X) increases with X and q*(X) ≤ q*_∞. The repo's LP computes q*(X) in a finite trigonometric basis; this note is about q*_∞.

## 2. The objects (derived here; checked numerically, section 7)

- χ = χ₋₃: χ(n) = 0, 1, −1 for n ≡ 0, 1, 2 mod 3. Λ(s) = 3^{s/2} Γ_R(s+1) L(s, χ), Γ_R(s) = π^{−s/2} Γ(s/2).
- θ(x) = Σ_{n≥1} n χ(n) e^{−πn²x/3}. **Classical (C1)**, theta transformation with root number +1: θ(1/x) = x^{3/2} θ(x) for x > 0 (Poisson summation; Davenport, *Multiplicative Number Theory*, ch. 9, recalled, not reopened here; checked at x = 0.6 and 1.7 to 25 digits).
- **Mellin (proved).** Λ(s) = 3^{−1/2} ∫₀^∞ θ(x) x^{(s+1)/2} dx/x for every s ∈ ℂ. For Re s > 1 this is the termwise integral of Γ(w) a^{−w} = ∫ e^{−ax} x^{w−1} dx (absolutely convergent). By (C1) θ(x) = O(x^{−3/2} e^{−π/(3x)}) as x → 0 and θ(x) = O(e^{−πx/3}) as x → ∞, so the integral is entire in s. Substituting x ↦ 1/x and using (C1) gives Λ(1 − s) = Λ(s).
- With x = e^{2u}: Ξ(t) := Λ(1/2 + it) = ∫ φ(u) e^{itu} du, φ(u) = (2/√3) e^{3u/2} θ(e^{2u}). By (C1), φ(−u) = φ(u). So Ξ is real on ℝ, even and entire. φ is smooth, and every derivative is O(exp(−e^{2|u|})). Integrating by parts k times gives (1 + |t|)^k |Ξ(t)| ≤ C_{k,Y} on |Im t| ≤ Y, for every k and Y; the same holds for Ξ′.
- **The family.** B(t) = Σ_{k=0}^{D} b_k t^{2k}, real, B ≥ 0 on ℝ, B ≢ 0. g_B := (1/2π) ∫ Ξ(t)² B(t) e^{−itu} dt = Σ_k b_k (−1)^k (φ*φ)^{(2k)}. g_B is real, even and smooth, every derivative is O(exp(−e^{|u|})), and g_B(0) = (1/2π) ∫ Ξ² B dt > 0.
- **Closed form (Lemma 2, proved).** For every real w, (φ*φ)(w) = (4/3) e^{3w/2} Σ_{N≥1} χ(N) d(N) N K₀(2πN e^w/3). Proof: write φ(v+w) with θ(e^{2(v+w)}) and φ(v) = φ(−v) with θ(e^{−2v}). Each pair (n, m) gives ∫ exp(−αe^{2v} − βe^{−2v}) dv = K₀(2√(αβ)) with 2√(αβ) = 2πnm e^w/3. Tonelli applies (all terms positive after taking absolute values, Σ nm K₀(c nm e^w) < ∞), and χ(n)χ(m) = χ(nm), Σ_{nm=N} 1 = d(N). Derivatives: with z = cNe^w, c = 2π/3, and D = z d/dz, D^j[z^{3/2}K₀(z)] = z^{3/2}[p_j(z)K₀(z) + q_j(z)K₁(z)], where p_{j+1} = (3/2)p_j + z p_j′ − z q_j, q_{j+1} = (1/2)q_j + z q_j′ − z p_j, p₀ = 1, q₀ = 0 (from K₀′ = −K₁, zK₁′ = −zK₀ − K₁). The differentiated series converge locally uniformly, by the K bound below.
- Constants: 4/3 = (2/√3)² and c = 2π/3 come from the conductor 3 in θ. d(N) comes from χ * χ = χ·d, since χ is completely multiplicative. The exponent 3/2 comes from Γ_R(s+1).
- **Local polynomials.** P_p(x) = Σ_{m≥1} x^m (log p) p^{−m/2} g(m log p) for x ∈ [−1, 1]; Defect_p = P_p(χ(p)) − min_{[−1,1]} P_p ≥ 0; Defect(g) = Σ_p Defect_p. For p = 3, χ(3) = 0 and P_3(0) = 0, so Defect_3 = −min P_3. The series converge absolutely, by the decay of g.

## 3. Lemma 1 (the identity; proved, unconditional)

For every B as above:

`(log 3) g_B(0) + A₁(g_B) − 2 Σ_{n≥2} χ(n) Λ(n) n^{−1/2} g_B(log n) = 0`, that is, `(log 3) g_B(0) + A₁(g_B) = 2 Σ_p P_p(χ(p))`.

**Proof (route 1: contour, no zeros).** Put B̃(s) = B(−i(s − 1/2)), a polynomial with B̃(1 − s) = B̃(s), because B is even. Let H(s) = Λ(s)² B̃(s), so that H(1/2 + it) = ĝ_B(t), and let F(s) = H(s) Λ′(s)/Λ(s) = Λ(s) Λ′(s) B̃(s).

(i) F is entire. A zero of Λ, on or off the critical line and of any multiplicity, is a zero of H of twice its order. It is therefore never a pole of F. *This is the only place where the zeros enter, and they do no harm.*

(ii) F(1 − s) = −F(s): differentiate Λ(1 − s) = Λ(s).

(iii) On every vertical strip a ≤ Re s ≤ b, |F(s)| ≤ C(1 + |Im s|)^{−2}, by the bounds for Ξ and Ξ′ in horizontal strips.

By Cauchy on the rectangle [−1, 2] × [−T, T] with T → ∞, ∫_{(2)} F ds = ∫_{(−1)} F ds. The substitution s ↦ 1 − s and (ii) turn the right side into −∫_{(2)} F ds. Hence ∫_{(2)} F ds = 0.

On Re s = 2, Λ′/Λ(s) = (1/2) log 3 + ψ_R(s+1) + L′/L(s, χ), with ψ_R = Γ_R′/Γ_R, and −L′/L(s, χ) = Σ χ(n)Λ(n) n^{−s} converges absolutely and uniformly there (Euler product, classical (C2)). Summation and integration may be exchanged, since Σ Λ(n) n^{−2} ∫ |H(2 + iy)| dy < ∞. Each of the three pieces H·n^{−s}, H·ψ_R(s+1) and H is holomorphic on 1/2 ≤ Re s ≤ 2: the poles of ψ_R(s+1) are at s = −1, −3, … (the trivial zeros, which Λ does not have). Each piece also decays in that strip, so it can be moved to Re s = 1/2, where s = 1/2 + it and ds = i dt:

- (1/2πi) ∫ H(s) n^{−s} ds = (1/2π) ∫ ĝ(t) n^{−1/2} e^{−it log n} dt = n^{−1/2} g(log n);
- (1/2πi) ∫ H(s) ds = g(0);
- (1/2πi) ∫ H(s) ψ_R(s+1) ds = (1/2π) ∫ ĝ(t) Re ψ_R(3/2 + it) dt, since ĝ is even and real and Im ψ_R(3/2 + it) is odd in t.

So 0 = (1/2)(log 3) g(0) + (1/2π) ∫ ĝ Re ψ_R(3/2+it) dt − Σ χ(n)Λ(n) n^{−1/2} g(log n). It remains to show (1/π) ∫ ĝ(t) Re ψ_R(3/2+it) dt = A₁(g). Since ψ_R(s) = −(1/2) log π + (1/2) ψ(s/2), Gauss's integral ψ(z) = −γ + ∫₀^∞ (e^{−x} − e^{−zx})/(1 − e^{−x}) dx (Re z > 0, classical (C3)) at z = 3/4 + it/2, with x = 2u, gives

`2 Re ψ_R(3/2+it) = −log π − γ + ∫₀^∞ 2[e^{−2u} − e^{−3u/2} cos(tu)]/(1 − e^{−2u}) du`.

The integrand is at most C(1 + min(1/u, t²u)) for u ≤ 1 and at most Ce^{−u} for u ≥ 1, so its u-integral is O(log(2 + |t|)). Since ∫ |ĝ| log(2 + |t|) dt < ∞, Fubini and g(u) = (1/2π) ∫ ĝ(t) cos(tu) dt give the formula for A₁. □

**Route 2 (zeros; classical, cited).** Weil's explicit formula for L(s, χ₋₃) with test function h = ĝ_B reads Σ_ρ h(γ_ρ) = (log 3) g(0) + A₁(g) − 2 Σ χΛ n^{−1/2} g(log n), with ρ = 1/2 + iγ_ρ over all nontrivial zeros (γ_ρ complex when ρ is off the line). It is valid for h holomorphic on |Im t| ≤ 1/2 + δ with h ≪ (1 + |t|)^{−1−δ} there, which our h satisfies (Iwaniec–Kowalski, *Analytic Number Theory*, Thm. 5.12; recalled, **not reopened** in this session). Every term on the left is h(γ_ρ) = Λ(ρ)² B(γ_ρ) = 0. Route 1 needs neither the zeros nor the theorem.

**Checks (computed, section 7).** The identity holds to 4·10⁻³¹ relative (the working precision; 10⁻¹⁸ asserted) for B = 1 and for B = 1 + t²/7, with g(0) and A₁ computed on the Ξ side from mpmath's Dirichlet L-function and the prime sum from the Bessel form. These are two separate codes; the check tests the code, not the identity.

## 4. Feasible families are positive on g_B; the bound

**Lemma 3 (proved).** If μ is feasible at q, then Q_{μ,q}(g_B) := (log q) g_B(0) + A₁(g_B) − 2 Σ_p ∫ P_p dμ_p ≥ 0.

*Proof.* (a) Fejér–Riesz: B ≥ 0 on ℝ gives B = |b|² on ℝ for a polynomial b with complex coefficients (pair the complex roots; real roots have even multiplicity). Put F := (1/2π) ∫ Ξ(t) b(t) e^{−itu} dt = b(i d/du) φ. F is smooth, |F| + |F′| ≤ C e^{−4|u|}, F̂ = Ξ b, so F*F̃ has transform Ξ² B and F*F̃ = g_B. (b) Take η ∈ C_c^∞, 0 ≤ η ≤ 1, η = 1 on [−1, 1], and F_M(u) = F(u) η(u/M) ∈ C_c^∞. By feasibility, Q_{μ,q}(F_M) ≥ 0. (c) Q_{μ,q}(F_M) → Q_{μ,q}(g_B) as M → ∞. For g_M = F_M*F̃_M, |g_M(u) − g_B(u)| ≤ C e^{−2M} e^{−2|u|}, using |F − F_M| ≤ |F| 1_{|u|>M} and e^{−4|v+u|−4|v|} ≤ e^{−2|u|} e^{−2|v+u|−2|v|}. So g_M(0) → g_B(0), and the prime sums converge, since |⟨x^m⟩| ≤ 1 and Σ Λ(n) n^{−5/2} < ∞. For A₁ write A₁(h) = (1/2π) ∫ ĥ Ω dt, Ω(t) = 2 Re ψ_R(3/2 + it) = O(log(2 + |t|)) (section 3). Cauchy–Schwarz gives ∫ |ĝ_M − ĝ_B| log(2 + |t|) dt ≤ ‖(1 + |t|)(F̂_M − F̂)‖₂ ‖(|F̂_M| + |F̂|) log(2 + |t|)/(1 + |t|)‖₂ ≤ C ‖F_M − F‖_{H¹} → 0. □

**Theorem 4 (proved, given Lemmas 1-3).** If μ is feasible at q, then for every B as above

`log q ≥ log 3 − 2 Defect(g_B)/g_B(0)`,  that is,  `q*_∞ ≥ 3 exp(−2 Defect(g_B)/g_B(0))`.

*Proof.* Subtracting Lemma 1 from Lemma 3 gives 0 ≤ Q_{μ,q}(g_B) = (log q − log 3) g_B(0) − 2 Σ_p [∫ P_p dμ_p − P_p(χ(p))]. Each bracket is ≥ min P_p − P_p(χ(p)) = −Defect_p. The sums converge absolutely (∫ Σ_m = Σ_m ∫ by absolute convergence), and g_B(0) > 0. □

**Corollary 5 (exactness, conditional).** (a) If inf_B Defect(g_B)/g_B(0) = 0, then q*_∞ ≥ 3, unconditionally. (b) Under GRH for L(s, χ₋₃) the point masses δ_{χ(p)} are feasible at q = 3: then Q(f) = Σ_ρ |f̂(γ_ρ)|² by the explicit formula for f ∈ C_c^∞. So q*_∞ ≤ 3. (a) and (b) together give q*_∞ = 3. Neither hypothesis of (a) or (b) is proved here.

**Hypotheses, and where each is used.** (C1) theta transformation: evenness of φ (g real and even), Λ entire with strip decay (Lemma 1 (iii), Lemma 3 (a)), Λ(1−s) = Λ(s) (Lemma 1 (ii)). (C2) Euler product: Lemma 1, prime side. (C3) Gauss's integral: the form of A₁ (Lemma 1, last step; Lemma 3 (c)). Fejér–Riesz: Lemma 3 (a). B ≥ 0: Lemma 3 (so that the competitor's positivity applies) and g_B(0) > 0 in Theorem 4; Lemma 1 holds for every even B. B even: Lemma 1 (ii). Not used anywhere: GRH, any computation of zeros, the Boas–Kac–Krein theorem, the truncation lemma of [../../theory/convex/convex_limit.md](../../theory/convex/convex_limit.md) (not needed, since Lemma 3 cuts F, not g).

## 5. The computation and the certificate

**Search (computed, not part of the proof).** g_B is linear in b, and Defect(g_B) is convex in b (P_p(χ(p)) is linear, min_x P_p is concave). It was minimised with g_B(0) = 1 as an LP (HiGHS) with cutting planes in x.
- Degrees 0-10: plain LP in the orthonormal basis of Ξ²dt ([c3opt.py](c3opt.py)).
- Degrees 12-22: iterative refinement ([c3refine2.py](c3refine2.py)). The step is centred at the current point and scaled by the current defect, in the coordinates of the SVD of the map b ↦ (g(0), (log p) p^{−m/2} g(log p^m)). All data are in mpmath at 250 digits; B ≥ 0 is imposed on a t-grid together with b_D ≥ 0.
- Basis data: [c3basis.py](c3basis.py) at 900 bits ([c3core.py](c3core.py)).

The plain float LP stops at about 10⁻¹⁰. The polynomials found are stored in [cert_d*.json](.) as exact decimals: 40 digits up to degree 10, 70 digits from degree 12.

**Certificate (computer-assisted, [validate_exact3.py](validate_exact3.py), arb at 600 bits, sympy).** For each stored B:
1. B > 0 on ℝ: exact Sturm count of Σ b_k s^k on [0, ∞) (no root) and b_0 > 0.
2. g_B(0) and g_B(log n) for every prime power n ≤ 400, as balls from the Bessel form. The N-tail is bounded by |p_j|, |q_j| ≤ (Σ|coeff|) z^j for z ≥ 1, K_ν(z) ≤ √(π/2z) e^{−z+1/(2z)} (from K_ν(z) = ∫₀^∞ e^{−z cosh s} cosh(νs) ds and cosh s ≥ 1 + s²/2) and d(N) ≤ N. K₀ and K₁ come from arb's bessel_k or, for large z, from the asymptotic series with the remainder bound of DLMF 10.40(ii) (read: for x > 0 and ℓ ≥ max(|ν| − 1/2, 1) terms, the remainder is at most the first neglected term).
3. For each prime p ≤ 400, a rigorous lower bound of min_{[−1,1]} P_p (powers p^m ≤ 400), by branch and bound. On an interval where the ball of P_p′ has a sign, the minimum is at an end. Otherwise a convex Taylor model P(x) ≥ P(c) + P′(c)(x−c) + m(x−c)²/2 is used, with m ≤ P″ on the interval. Upper bound of Defect_p = upper(P_p(χ(p))) − lower(min).
4. Everything above 400 (the primes > 400 and the powers p^m > 400 of smaller primes) is bounded together by 2 Σ_{n>400} log n · n^{−1/2} |g(log n)|. A geometric series of the same envelope bounds it: below 10⁻²⁹⁷ at every degree.
5. The bound q*_∞ ≥ 3 exp(−2 upper(Defect)/lower(g(0))) as a ball's lower end. The sum of the per-prime lower bounds P_p(χ(p)) − (an attained value of P_p) shows that the upper bound is tight to the digits given.

| deg B | Defect(g)/g(0), certified interval | q*_∞ ≥ 3 − … (certified) | largest share | signs of b_k |
|---|---|---|---|---|
| 0 | [4.531, 4.532]·10⁻² | 0.2600 | p = 3 | + |
| 2 | [6.130, 6.131]·10⁻³ | 3.656·10⁻² | 5 | all + |
| 4 | [4.470, 4.471]·10⁻⁴ | 2.681·10⁻³ | 7 | all + |
| 6 | [1.380, 1.381]·10⁻⁵ | 8.284·10⁻⁵ | 3 (then 11) | all + |
| 8 | [5.136, 5.137]·10⁻⁷ | 3.082·10⁻⁶ | 13 | all + |
| 10 | [2.751, 2.752]·10⁻⁹ | 1.651·10⁻⁸ | 17 | all + |
| 12 | [7.479, 7.480]·10⁻¹¹ | 4.488·10⁻¹⁰ | 19 | all + |
| 14 | [5.252, 5.253]·10⁻¹³ | 3.152·10⁻¹² | 23 | all + |
| 16 | [1.510, 1.511]·10⁻¹⁷ | 9.062·10⁻¹⁷ | 31 | all + |
| 18 | [1.509, 1.510]·10⁻¹⁷ | 9.060·10⁻¹⁷ | 31 | all + |
| 20 | [1.509, 1.510]·10⁻¹⁷ | 9.060·10⁻¹⁷ | 31 | all + |
| 22 | [6.349, 6.350]·10⁻²³ | 3.810·10⁻²² | 41 | b_9 < 0, rest + |

(B = 1 gives the bound 2.740045; degree 4 gives 2.997319.) What is certified: every entry, given arb's ball arithmetic, sympy's root count, the DLMF bound and the correctness of one code. What is only computed: that these B are optimal at their degree. They are the outputs of a numerical search, and the bound needs no optimality.

## 6. Structure, and what enters as the degree grows (computed; the readings are heuristic)

- **Where the cost sits.** At every degree one prime carries almost all of Defect: 3, 5, 7, 3/11, 13, 17, 19, 23, 31, 31, 31, 41 at degrees 0 to 22. Before it, g_B has (near-)zeros or the right sign at the prime powers: a zero at log 3 (first order, so that P_3 has its minimum at x = 0), and g(log p) of sign −χ(p) or ≈ 0 beyond. For example, at degree 16, g/g(0) = 0.219, 9·10⁻¹⁶, −0.0036, 9·10⁻¹⁵, −6·10⁻²⁴ at n = 2, 3, 4, 5, 7. The carrier is a prime with χ(p)·g(log p) > 0, where P_p has its minimum at −χ(p), so that Defect_p ≈ 2 log p · p^{−1/2} |g(log p)|.
- **p = 2 never carries cost.** g(log 2) > 0 and g(log 4) < 0 make P_2 increasing on [−1, 1], with minimum at −1 = χ(2).
- **The ramified prime is cheap.** Defect_3 ≈ c₁²/(4c₂) is second order in g(log 3), with c₁ = log 3·g(log 3)/√3 and c₂ = log 3·g(log 9)/3 > 0.
- **Decay.** log₁₀(Defect/g(0)) is −1.34, −2.21, −3.35, −4.86, −6.29, −8.56, −10.13, −12.28, −16.82, −16.82, −16.82, −22.20. That is roughly one prime per two degrees up to degree 16, with no stall below degree 16.
- **Something new enters at degree 18: positivity of B and the sign of g at infinity.** From the Bessel form, g_B(u) ~ B(i·(2π/3)e^u)·(4/3)e^{3u/2}K₀((2π/3)e^u) as u → ∞, so the eventual sign of g_B is (−1)^D sign(b_D). B ≥ 0 forces b_D ≥ 0, so the eventual sign is (−1)^D. Beyond the last controlled prime, the cost falls on the class χ(p) = (−1)^D. This holds for the largest share outside p = 3 at every degree except 18 and 20, where the degree-16 pattern (class +1) persists: 7, 5, 7, 11, 13, 17, 19, 23, 31, 41 at degrees 0, 2, …, 16, 22. At degrees 18 and 20 the search found no gain over degree 16: the new top coefficients stay at b_9/b_0 < 10⁻²⁴ and b_10/b_0 < 10⁻⁴⁰. In a run without the constraint B ≥ 0 (log not kept), the search at degree 20 lowered Defect only with b_9, b_10 < 0, that is with B < 0 for large t, which is inadmissible. At degree 22 (b_11 > 0, one negative inner coefficient, B > 0 certified by Sturm) the defect falls again, to 6.3·10⁻²³, with cost at 41. So up to degree 16 the constraint ĝ ≥ 0 is slack (all b_k > 0), and from degree 18 it binds through the leading coefficient. This is the analogue of the parity rule of the ζ case (bootstrap.md 5f, where degree 4k+2 equals degree 4k), but here it acts on a sign pattern, not on positivity of g. Degree 24: the search was cut by its 10-minute limit while still improving ([c3refine2_d24.log](c3refine2_d24.log)); nothing stored, nothing certified.
- **The repo's LP dual (support 24, [dual_odd_out.txt](../../../../../data/study/maass/dual_odd_out.txt)).** It has near-zeros at n = 3, 5, 7, 13, like g_B here. It is positive at n = 4 and 16, which dual.py flags as a failed sign condition. That per-power test is not the right criterion for m ≥ 2. The right one is that P_2 has its minimum at x = −1. From the printed values, P_2′(−1) = −4.05·10⁻⁵ and P_2″(−1) = 0.0271 (units of G(1), factor log 2 omitted), so the minimum sits just inside −1 and the loss is second order. The construction does not fail there. In the B found here g(log 4) < 0 instead: the extremal function is not unique.

## 7. Checks of the set-up (computed; [validate_exact3.py](validate_exact3.py) section A)

- (C1) at x = 0.6 and 1.7 to 10⁻²⁵; Ξ(t) = ∫ φ cos at t = 0 and 6.3 against mpmath's Dirichlet L, to 10⁻²⁰.
- Closed form of φ*φ against direct quadrature at log 1, log 2, log 3, to 20 digits.
- g_k(0) = (1/2π) ∫ Ξ² t^{2k} dt for k ≤ 2 (L-function side) against the closed form, to 20 digits.
- K₀ and K₁ (arb, or the asymptotic series with the DLMF bound) against mpmath at 130 digits, z from 2.09 to 838: 5·10⁻¹²⁵ relative.
- Lemma 1 for B = 1 and B = 1 + t²/7, with g(0) and A₁ from the L-function side: residual 3.7·10⁻³² and 2.2·10⁻³² (g(0) = 0.0885 and 0.1395).

These test the code, not the theorems.

## 8. The theorem as it stands, and the gap

**Theorem (computer-assisted; Lemmas 1-3 and Theorem 4 proved, the numbers certified by one code).** Every family of probability measures (μ_p) on [−1, 1], one per prime, that keeps the Weil form with Gamma factor Γ_R(s+1) non-negative at conductor q for all f ∈ C_c^∞(ℝ), has q > 3 − 3.81·10⁻²². This is unconditional: no Riemann hypothesis, no zeros, no unproved input beyond the certified arithmetic (arb, sympy's root count, DLMF 10.40(ii), one code). With the degree-8 polynomial alone (five coefficients), q > 3 − 3.1·10⁻⁶.

**Gap 1: exactness.** q*_∞ = 3 needs (a) inf_B Defect(g_B)/g_B(0) = 0, which is unproved; computed through degree 22, falling but with a stall at degrees 18-20; and (b) GRH for L(s, χ₋₃) for q*_∞ ≤ 3. Two routes to (a), neither attempted: an explicit limit ĝ = Ξ²H, H ≥ 0 entire, with Defect = 0, that is g(log 3) = 0 and the sign pattern sign g(log p) = −χ(p) at every prime together with the polynomial conditions at small p (the analogue of Johns' interpolation function for ζ); or a quantitative approximation argument showing that each extra pair of degrees controls one more prime.

**Gap 2: finite support.** The theorem is about q*_∞. The repo's numbers 2.9999977254 (X = 8) and 2.9999999996 (X = 14) are q*(X) at finite support, and nothing here bounds q*(X) from below. A finite-support version needs Lemma 3 with explicit constants: a bound for F = b(i d/du)φ outside [−a + δ, a − δ] (F decays like exp(−(π/3)e^{2|u|}), but b is the Fejér–Riesz factor of B, so its roots must be enclosed rigorously) and the resulting change of Q. Not done.

**Gap 3: review.** One code, no referee; the second-code checks of section 7 cover the closed form and Lemma 1, not the branch and bound or the tails.

## 9. Classical and derived; prior art

Classical: (C1)-(C3), Fejér–Riesz, the Bessel integral for K₀, the explicit formula (route 2). Derived here: Lemma 1 by route 1 for this family, the closed form for χ₋₃ (the ζ analogue with d(N) and K₀, K₁ is in [../../theory/convex/certs/README.md](../../theory/convex/certs/README.md)), Lemma 3, Theorem 4 and the computation. The idea that Ξ² kills the zero side is Johns', Part I, for ζ (as read by the repo, [../../bootstrap.md](../../bootstrap.md) section 5f). The structured version with one measure per prime and the local defect was not searched for in print in this run: no literature search was made. The searches of [../maass.md](../maass.md), section 1 did not find "exactness of the structured bound".

## 10. Files

- [validate_exact3.py](validate_exact3.py): the one script that recomputes and asserts every number above and ends with "N checks, 0 failed". `quick` skips section A. Run as `(ulimit -v 1200000; timeout 900 nice .venv/bin/python -I validate_exact3.py)`. Output: [validate_exact3_out.txt](validate_exact3_out.txt).
- Logs: [c3refine2_d24.log](c3refine2_d24.log) (the interrupted degree-24 search).
- Search: [c3core.py](c3core.py), [c3basis.py](c3basis.py) (→ c3basis_D24.json, c3basis_D8.json), [c3opt.py](c3opt.py) (→ c3opt_d*.json), [c3refine2.py](c3refine2.py) (→ c3ref2_d*.json), [c3mkcert.py](c3mkcert.py) (→ cert_d*.json); [c3refine.py](c3refine.py) is the first refinement, superseded.
