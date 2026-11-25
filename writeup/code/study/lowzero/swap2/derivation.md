# swap2: a registered re-test of the swapped (Gram-lattice) term between two L-functions

Folder: this scratch folder (nothing in the repo was changed). No phase statistic was computed on a reserved zero set; of
the reserved files only the keys were read (`freeze.py` reads `z.files` only). Zeros used: r138_Fn (used set, validation),
lzeros/q7_k1.txt and q13_k1.txt (non-reserved complex characters, pilot of the conventions), synthetic sets (dry runs).

Labels: **proved** (algebra from stated formulas), **conjecture** (ratios recipe / Hardy–Littlewood), **heuristic**,
**computed**, **classical** (in print), **retrospective** (on data where the answer was known).

| file | output | what it does |
|---|---|---|
| `hlib.py` | | H₁₂(u) Euler factors, J₁₂(σ) = ∫H K for all pairs (real fast path; general complex path), character tables (PARI), κ, ε |
| `check_exact.py` | `check_exact_out.txt`, `.json` | J from the Euler product cut at 3·10³, 3·10⁴, 3·10⁵ against PARI `lfun` for L(1 − iu, ψ*) (26 pairs: real, complex, (χ, χ̄), complex × real) |
| `grid_check.py` | `grid_check_out.txt` | u-quadrature (step 0.05 to 6) against step 0.025 to 8 and 0.1 to 6 |
| `fn_validate.py` | `fn_validate_P3000.txt`, `fn_validate_P30000.txt` | own X₁₂ and predictor on Fn against the review's numbers |
| `fn_onepoint.py` | `fn_onepoint_out.txt` | one-point sums on Fn against the mean law and the per-function resonance prime sum |
| `pilot_complex.py` | `pilot_complex_out.txt` | conventions for complex χ on non-reserved zeros (mod 7, mod 13, both half-lines) |
| `complete_calib.py` | `complete_calib_out.txt` | completeness rule for sets stored from height 100 only (complex characters) |
| `freeze.py` | `pred_OO.npz`, `pred_C.npz`, `pred_FN.npz`, `onepoint.npz`, `freeze_manifest.json`, `freeze_*_out.txt` | the frozen predictions |
| `noise_model.py`, `power.py` | `noise_model.json`, `power_{OO,HI,C}.json`, `power_*_out.txt` | expected value of every reading under each hypothesis; expected se |
| `onepoint_expected.py` | `onepoint_expected_out.txt`, `.json` | the one-point readings computed from the primes alone (frozen expectations) |
| `swap2_test.py` | `run_validation_FN.txt`, `dryrun_poisson.txt`, `dryrun_gram.txt` (+ `run_*/`) | the registered test; validation on Fn; dry runs |
| `round_draft.txt` | | registration draft |

## Summary

1. **One-point law (derived; theorem in substance).** For primitive χ mod q, parity a, ε_χ = τ(χ)/(i^a √q),
   κ = q/φ(q), θ_χ(t) = (t/2) ln(q/π) + Im lnΓ((1/2 + a + it)/2) (odd in t), and the zeros 1/2 + iγ of L(s, χ):

       Σ_{T₁<γ≤T₂} e^{2iθ_χ(γ)} = −ε_χ κ (T₂ − T₁)/2π + o(T),     c_χ = −ε_χ κ,

   on **both** half-lines (for γ < 0 the same c_χ with θ_χ extended as an odd function; equivalently the mirrored
   negative zeros are the zeros of χ̄ with c_χ̄ = −ε̄_χ κ). For real χ, ε = 1 and c = −κ (the review's −k/φ(k)).
   Per function at finite T it equals, up to a remainder that is 0.2 % of the main term, the prime sum
   R_w = −(i^a/√q) Σ_n χ(n)Λ(n) e(−n/q) w(2πn/q) (Hann taper w: per-function rms |G − R| = 0.23 on Fn, against
   |law| 40–90; slope 1.0000 ± 0.0001). So the one-point readings are **data checks** of the explicit formula plus the
   prime number theorem in progressions, not tests of a conjecture (the ζ case is Arias de Reyna's Theorem 3.1).
2. **Two-function law (conjecture: ratios recipe).** With ψ = χ₁χ̄₂ (mod lcm, ψ(p) = 0 for p | q₁q₂) and K(u) = e^{−u²/2σ²},

       X₁₂ = (1/T) Σ_{γ₁,γ₂} e^{2iθ₁(γ₁) − 2iθ₂(γ₂)} K(γ₂ − γ₁),   E X₁₂ = (1/4π²) κ₁κ₂ ε₁ ε̄₂ J₁₂,
       J₁₂ = ∫ H₁₂(u) K(u) du,   H₁₂(u) = |L(1 − iu, ψ)|² A₁₂(−iu) = ∏_{p∤q₁q₂} (1 − 2/p + y)/((1 − ȳ)(1 − 1/p)²),  y = ψ(p) p^{iu}/p.

   Real pairs: ε₁ε̄₂ = 1 and J real. Real × complex: the factor is ε̄_χ (or ε_χ). Complex pairs: ε₁ε̄₂, J complex
   (median |Im J|/|J| = 0.44 for complex × complex, 0.24 for complex × real). (χ, χ̄): factor ε_χ², ψ = χ²: the
   correlation of the zeros of one L(s, χ) at +t with those at −t (not the self term, which needs ψ principal).
   For distinct primitive characters ψ is never principal; H has no pole.
3. **J is evaluated exactly enough.** Because K smooths H at the scale σ, the Euler product converges fast for J
   (not for H pointwise): cut at 3·10⁴ it agrees with PARI's L(1 − iu, ψ*) to ≤ 7.7·10⁻⁵ absolute (4·10⁻⁵ relative)
   on 26 pairs; the frozen values agree with the exact ones to ≤ 4.3·10⁻⁵ on all 42 compared (pair, σ) values.
4. **Validation on the used set (retrospective, as the review).** Own code reproduces every printed digit of the
   review (σ = 0.5: 1.0037 ± 0.0067, R² 0.911; H = 1: 0.9645 ± 0.0239; σ = 1: 1.0090 ± 0.0083; one-point 1.0084);
   the registered pipeline with its margin window gives 1.0049 ± 0.0069 (σ = 0.5) and 1.0112 ± 0.0085 (σ = 1).
   Pilot on non-reserved complex characters: one-point law with ε holds on both half-lines; X for the 6 pairs follows
   ε₁ε̄₂J, not J (no root numbers) and not conj J.
5. **Frozen**: all 148,785 OO pairs, 16,290 HI pairs, 25,960 complex-design pairs (15,400 CC incl. 88 (χ, χ̄), 10,560
   CR), σ = 0.5 and 1.0, complex values, with the alternatives; one-point laws and prime-side resonance sums for every
   function. sha256 in `freeze_manifest.json` and in `swap2_test.py`.
6. **Test**: `swap2_test.py all --designs OO,HI,C` (about 50 s wall, 0.52 GB); dry runs: Poisson → every reading 0
   (no term); jittered Gram points → exactly the H ≡ 1 values (slope 0.156 on the law, H − 1 part 0.001).
7. **Power** (noise model from Fn): OO slope se ≈ 0.004 (H_const at 0.155: z ≈ 210; H_half z ≈ 125); complex slope
   se ≈ 0.006–0.007 (H_noroot at 0.115: z ≈ 130); root-number coefficient se ≈ 0.004; direction reading se ≈ 0.01
   (law 1 vs conj J −1).
8. **Prior art**: the ζ one-point theorem (Arias de Reyna 2014, Fujii 1992, Schoißengeier 1979); Landau–Gonek for
   L(s, χ) (Banks 2023) without the phase-weighted term; no one-point χ-version and no cross Gram-lattice statement or
   test found. Verdicts in §12.

## 1. Conventions (computed checks)

- Functional equation: Λ(s, χ) = (q/π)^{(s+a)/2} Γ((s+a)/2) L(s, χ) = ε_χ Λ(1 − s, χ̄), ε_χ = τ(χ)/(i^a√q),
  τ(χ) = Σ_n χ(n) e^{2πin/q}. PARI's `lfunrootres(lfuncreate(Mod(m, N)))[3]` equals τ/(i^a√q) to 2.5·10⁻¹⁴ for all 88
  r121 characters (`freeze_complex_out.txt`); `lfuncreate(Mod(m,N))` is the Conrey character whose values `chareval`
  gives (L(2) to 10⁻⁸ against the Dirichlet series).
- On the line, X_χ(1/2 + it) := G(1/2 − it)/G(1/2 + it) = e^{−2iθ_χ(t)}, so Z_χ(t) = ε_χ^{−1/2} e^{iθ_χ(t)} L(1/2 + it, χ)
  is real for t of both signs: checked with PARI for four r121 characters at t = ±123…±480 (|Im Z| < 10⁻¹² |Z|).
  Hence L(1/2 + it) = ε^{1/2} e^{−iθ} Z(t): zeros lie between the "Gram points" of χ (L real) where
  e^{2iθ} ε̄ = −1, which is the sign of the one-point law.
- Negative heights: L(s̄, χ̄) = conj L(s, χ), so the zeros of L(s, χ) at −γ are the conjugates of the zeros of L(s, χ̄)
  at +γ; θ_χ̄ = θ_χ (same q, a) and θ is odd. ε_χ̄ = conj ε_χ (τ(χ̄) = χ(−1) conj τ(χ)). In the designs each complex
  character gives two one-sided functions: χ (γ > 0) and χ̄ (the mirrored γ < 0), both with phases e^{2iθ_χ(|γ|)}.

## 2. One-point law

**Recipe (conjecture-level derivation, as derivation.md of round 141 §3.5).** The zero density is
d(t) = (1/π)[θ′(t) + Re (L′/L)(1/2 + δ + it, χ)] (δ → 0⁺): each zero contributes Re 1/(δ + i(t − γ)) → πδ(t − γ),
and Re of (1/2)ln(q/π) + (1/2)ψ_Γ((1/2 + a + it)/2) is θ′. So

    Σ_γ e^{2iθ(γ)} = (1/π)∫ e^{2iθ}θ′ dt + (1/2π)∫ e^{2iθ}(L′/L)(1/2 + δ + it, χ) dt + (1/2π)∫ e^{2iθ}(L′/L)(1/2 + δ − it, χ̄) dt.

The first term is a boundary term (e^{2iθ}/2πi, O(1)). In the second the weight e^{2iθ} cancels the phase of the
swapped L(s + α, χ) = ε X_χ(s + α) L(1 − s − α, χ̄); the diagonal of L(1 − s − α, χ̄)/L(s + γ, χ) is 1/L(1 + γ − α, χ₀)
= (γ − α)κ + O((γ − α)²), and ∂_α[ε (qt/2π)^{−α}(γ − α)κ] at γ = α = δ → 0 is −εκ. The third term carries e^{4iθ}
after its swap and averages out. Hence **c_χ = −ε_χ κ_χ** per unit (T₂ − T₁)/2π.

**Resonance form (heuristic stationary phase; the method of Arias de Reyna's Theorem 3.1).** Term by term in
−Σ χ(n)Λ(n) n^{−1/2−it}, the phase 2θ(t) − t ln n is stationary at t_n = 2πn/q (2θ′ = ln(qt/2π)), with
2θ(t) = t ln(qt/2π) − t − π/4 + aπ/2 + O(1/t) and 2θ″ = 1/t, so

    (1/2π)∫ e^{i(2θ(t) − t ln n)} n^{−1/2} w(t) dt ≈ w(t_n) (1/2π)√(2π t_n) n^{−1/2} e^{−it_n − iπ/4 + iaπ/2 + iπ/4} = w(t_n) i^a e(−n/q)/√q,

    G_w(χ) := Σ_γ e^{2iθ_χ(γ)} w(γ) ≈ R_w(χ) := −(i^a/√q) Σ_n χ(n)Λ(n) e(−n/q) w(2πn/q).

Its mean (proved from PNT in progressions, fixed q): Σ_{b mod q} χ(b)e(−b/q) = χ(−1)τ(χ), so
R_w ≈ −(i^a/√q)(q∫w/2π)χ(−1)τ(χ)/φ(q) = −κ ε (∫w)/2π since i^a χ(−1) = i^a(−1)^a = 1/i^a. For ζ (q = 1) this is
Arias de Reyna's −ψ(T/2π). **Computed:** on Fn with a Hann taper on [100, 600], G_w − R_w has rms 0.23 per function
(max 0.57, against |law| ≈ 40), and the regression of G − law on R − law has slope 0.999 and correlation 1.000; with the
sharp window the edge (Fresnel) terms leave an rms of 3.4–3.8 (`fn_onepoint_out.txt`). The 1.008 of the review is the
primes' own value (prime-side slope 1.0067, sharp; 1.0111 Hann, equal to the zeros' 1.0111). Pilot, complex: χ₇ = χ₇(3, ·)
G/law = 0.979 − 0.045i (sharp), |G − R| = 0.18 (Hann); χ₁₃(2, ·) and its mirror: |G − R| = 0.17 and 0.08; the mirrored
χ₇ side gives 1.76 and is the one side flagged incomplete (count range −2.78, step 2.07: a missing pair).
Consequence: the one-point "law" is a theorem in substance; its readings are registered as data checks with frozen
prime-side values (`onepoint_expected_out.txt`): OO 0.9983 (sharp)/0.9950 (Hann), HI 1.0038/1.0040,
C 1.0157 − 0.0023i; even/odd OO 0.9993/0.9974; the κ regression c₁ = −0.9927.

## 3. Two-function law

**From the recipe (derivation of round 141, §3.3–3.4; re-derived by the review by two routes).**

    R₁₂_swap(t₁, t₂) = (1/4π²){κ₁κ₂ ε₁ε̄₂ e^{2i(θ₂(t₂) − θ₁(t₁))} H₁₂(t₂ − t₁) + c.c.},  H₁₂(u) = L(1 − iu, ψ) L(1 + iu, ψ̄) A₁₂(−iu),

with u = γ₂ − γ₁ (zero of χ₂ minus zero of χ₁: α + β = −iu from s + α = 1/2 + iγ₁, 1 − s + β = 1/2 − iγ₂).
Euler factor (proved): for p ∤ q₁q₂, with y = ψ(p)p^{−1+iu}, L-factors (1 − y)^{−1}(1 − ȳ)^{−1} and A-factor
(1 − y)(1 − 2/p + y)/(1 − 1/p)² give h_p = (1 − 2/p + y)/((1 − ȳ)(1 − 1/p)²); for p | q₁q₂ the factor is 1. Since
L(1 + iu, ψ̄) = conj L(1 − iu, ψ), H = |L(1 − iu, ψ)|²·A₁₂(−iu) with A absolutely convergent
(per factor (1 − 2/p + 2y/p − y²)/(1 − 1/p)² = 1 + O(p^{−2})), L(s, ψ) = L(s, ψ*)∏_{p|q₁q₂, p∤f}(1 − ψ*(p)p^{−s}).

**The statistic.** X₁₂ := (1/(√π s T))∫_W G₁ conj G₂ dt with G_j(t) = Σ e^{2iθ_j(γ)}e^{−(t−γ)²/2s²}, s = σ/√2, equals
(1/T)ΣΣ e^{2iθ₁(γ₁) − 2iθ₂(γ₂)} K(γ₂ − γ₁) w(m) with w the window indicator smoothed at width s at the midpoint m
(proved: ∫e^{−(t−γ₁)²/2s² − (t−γ₂)²/2s²}dt = √π s K(u) × [Gaussian in m]). Inserting R₁₂_swap: the first term's phase
cancels exactly, its integral is T·(1/4π²)κ₁κ₂ε₁ε̄₂∫H K du for every pair whatever the conductors; the c.c. term
oscillates like e^{4iθ} and the smooth d̄₁d̄₂ term is killed by K: (1/T)∫∫e^{2iθ₁ − 2iθ₂}K d̄₁d̄₂ ~ d̄² σ√2π e^{−σ²L²/2}
(estimate: ≤ 2·10⁻³ at σ = 0.5, L = 7.4, relevant only for equal conductors, where the t-average does not remove it;
nil at σ = 1). The diagonal term (1/2π²)Re D_ψ(−iu) contributes only through orbits n ≈ e^L (the Heisenberg range),
mean zero over pairs, of order 5·10⁻⁴ per pair (heuristic estimate). Hence **E X₁₂ = (1/4π²)κ₁κ₂ε₁ε̄₂J₁₂** (conjecture). X₂₁ = conj X₁₂ (K even) and
J₂₁ = ∫H_ψ̄K = conj J₁₂ (proved: H_ψ̄(u) = conj H_ψ(−u)), so the law is consistent under exchange.
With H ≡ 1 it is the product of the two one-point terms: (−ε₁κ₁/2π)(−ε̄₂κ₂/2π)·σ√2π.

**Cases.**
- (i) real × real (any parities): ε = 1, ψ real, H(−u) = conj H(u), J real; parities enter only through θ. Im X: 0.
- (ii) real × complex: factor ε̄_χ (complex second) or ε_χ (complex first); J complex.
- (iii) complex × complex: factor ε₁ε̄₂; J complex; |L(1 − iu, ψ)|² ≠ |L(1 + iu, ψ)|², the even K symmetrises, and the
  convention of direction/role enters only through Im J: the alternative pairing (ψ ↔ ψ̄) gives conj J (proved above);
  reversing u alone gives the same J.
- χ₂ = χ̄₁: ψ = χ₁², factor ε_χ², J with L(1 − iu, χ²). Since the zeros of χ̄ at +γ are the mirrors of those of χ at
  −γ, this is the correlation of one L-function's zeros across the real axis, e^{2iθ_χ(γ₁) + 2iθ_χ(γ₂)}K(γ₁ + γ₂) over
  γ₁ > 0 > γ₂. It is **not** the one-function swapped term under two labels: that term correlates nearby zeros on one
  side and has ψ = χ₀ (pole); here ψ = χ² is non-principal for complex χ (it may be imprimitive: χ₁₂₄(3, ·)² has
  conductor 31, χ₃₈₈(3, ·)² conductor 97) and H has no pole. For real χ the two labels coincide and the case is the self case.
- Equal modulus, distinct characters: ψ non-principal mod q, possibly imprimitive with small conductor (for real
  ±D with 8 | D: ψ* = χ₋₄ or χ±₈); handled by ψ(p) = 0 at p | q₁q₂ in the product (checked against PARI, (−104, 104)).
- ψ principal happens only for the same primitive character (χ₁ = χ₂): then |L(1 − iu, χ₀)|² ~ (φ/q)²/u², the formula
  diverges at u = 0, and the pole cancels against the double pole of the diagonal term (1/2π²)Re(L′/L)′ ~ −1/(2π²u²)
  and the δ(u) (same-zero) term: (1 − e^{−iLu})²/(4π²u²) near 0; X becomes d̄ plus the K-smoothed form factor at the
  Heisenberg time (the BK/GUE kink). No such pair is in any design.

**Heuristic meaning (BK correspondence; not computed here).** By the explicit formula in frequency form, G_j at the
Heisenberg frequency is carried by Λ(n)χ_j(n)e(−n/q_j); E X₁₂ is then a correlation of the two progressions at
n/q₁ ≈ m/q₂. H ≡ 1 is the product of their means (PNT in progressions, as in §2); H − 1 is the Hardy–Littlewood singular
series of the linear form q₂n − q₁m, as Bogomolny–Keating obtained ζ's off-diagonal term from Hardy–Littlewood. So the
conjectural content tested is a prime-pair (HL-type) statement seen through the zeros; unlike round 141's diagonal law,
it is not a consequence of the explicit formula plus PNT. At σ = 0.5 the Gaussian is far narrower than the
stationary-phase width √(2πt) ≈ 45, so X₁₂ is not a local prime sum and no per-pair prime-side predictor was built.

## 4. Evaluating J₁₂

J = Σ_u w_u [H(u) + H(−u)] on u = 0, 0.05, …, 6 (trapezoid, folded; K(6σ) = e^{−18}); H from the Euler product to
P = 3·10⁴ (p = 2, 3 multiplied directly: h₃ vanishes where ψ(3)3^{iu} = −1, e.g. u = 0 for real ψ(3) = −1; p ≥ 5 via
logarithms). All 546² real pairs at once: ln H_ij = Σ_p χ_iχ_j(a_p − b_p)/2 + χ_i²χ_j²(a_p + b_p)/2 (a, b = ln h for
ψ(p) = ±1), two matrix products per u (36 s). Complex: p ≤ 100 directly per pair, p > 100 by the series of ln h_p in y
to order 3 as matrix products (125 s for 236 functions); the two code paths agree to 1.6·10⁻⁸ on the 60 × 60 real block.
**Why the product converges for J:** the tail ∏_{p>P} h_p = exp(Σ 2Re ψ(p)p^{iu}/p + O(p⁻²)) oscillates at
frequencies ln p > ln P where K̂ = e^{−σ²ln²p/2} is negligible; it acts only through H_P's content near ln p, whose
coefficients are O(1/(pm²)) — random-sign estimate 5·10⁻⁵, systematic O(1/(P ln P)) (heuristic).
**Computed** (`check_exact_out.txt`): against PARI `lfun` (lfuninit on |Im s| ≤ 6.5; ψ* from `znchartoprimitive`; removed
factors; A to 2·10⁶): max |ΔJ| 1.35·10⁻³ (P = 3·10³), 7.7·10⁻⁵ (3·10⁴), 1.1·10⁻⁵ (3·10⁵) for σ = 0.5; 6.3·10⁻⁴,
4.1·10⁻⁵, 8.8·10⁻⁶ for σ = 1; |J| from 0.12 to 14. Quadrature: step 0.025 to 8 or 0.1 to 6 change J by ≤ 2·10⁻⁸.
Frozen J vs exact on 42 values: ≤ 4.3·10⁻⁵. On Fn the slope moves by 1·10⁻⁴ between P = 3·10³ and 3·10⁴.
Ranges (frozen): OO J ∈ [−9.3, 25.0] at σ = 0.5 ([−10.4, 29.5] at σ = 1); predicted X has mean 0.094, sd 0.137 (σ = 0.5),
22 % of the pairs negative (min −0.254).

## 5. Validation on the used set Fn (retrospective)

`fn_validate.py` (own θ with exact loggamma, own Kronecker table, review's window [100, 600] and grid): completeness
flags D = 1153 only (1.51, as in round 138); 185 sets, 17,020 pairs.

| | review s8 | here, P = 3·10³ | here, P = 3·10⁴ | swap2_test.py (window [100, 594]) |
|---|---|---|---|---|
| σ = 0.5 slope (icpt) | 1.0037 ± 0.0067 | 1.0037 ± 0.0067 | 1.0038 ± 0.0067 | 1.0049 ± 0.0069 |
| σ = 0.5 H = 1 slope | 0.9645 ± 0.0239 | 0.9645 ± 0.0239 | 0.9645 ± 0.0239 | |
| σ = 0.5 H − 1 part | 1.0038 ± 0.0072 | 1.0038 ± 0.0072 | 1.0039 ± 0.0072 | 1.0055 ± 0.0074 |
| σ = 1.0 slope | 1.0090 ± 0.0083 | 1.0090 ± 0.0083 | 1.0090 ± 0.0083 | 1.0112 ± 0.0085 |
| R² (σ = 0.5) | 0.911 | 0.911 | 0.911 | |
| one-point slope | 1.008 ± 0.007 | 1.0084 (jackknife 0.0057) | | 1.0084 ± 0.0057 |

Im X on Fn: slope −0.0049 ± 0.0046 (σ = 0.5), −0.0111 ± 0.0054 (σ = 1; 2.1 se), mean −0.0008. Bins of |ln(q₂/q₁)|
(Fn spans 0–0.26 only): 1.0058 ± 0.0079, 1.0034 ± 0.0078. Parities: 1.0136 ± 0.0106 (even-even), 0.9984 ± 0.0101,
1.0036 ± 0.0075 (`run_validation_FN.txt`). Residual per pair after the law 0.0425 against sd(Im X) 0.0456: the law
leaves Re X at the noise level of Im X.

## 6. Pilot on non-reserved complex zeros (computed; disclosed)

q7_k1.txt is χ₇(3, ·) (order 6, odd; |L(1/2 + iγ)| < 10⁻¹² at its zeros), q13_k1.txt is χ₁₃(2, ·) (order 12, odd).
Window [106, 784], T = 678 (`pilot_complex_out.txt`):

| pair | X (σ = 0.5) | law ε₁ε̄₂J | no root numbers | conj J |
|---|---|---|---|---|
| (χ₇, χ₁₃) | 0.070 + 0.102i | 0.054 + 0.099i | 0.068 + 0.090i | 0.081 − 0.079i |
| (χ̄₇, χ̄₁₃) | 0.061 − 0.097i | 0.054 − 0.099i | 0.068 − 0.090i | 0.081 + 0.079i |
| (χ₁₃, χ̄₁₃) | −0.027 + 0.023i | −0.016 + 0.026i | 0.031 + 0.003i | −0.012 + 0.029i |
| (χ₇, χ̄₇) | 0.019 + 0.039i | 0.013 + 0.019i | 0.004 − 0.022i | −0.019 − 0.012i |

At σ = 1 (χ₁₃, χ̄₁₃): −0.089 + 0.010i against law −0.069 + 0.021i, no root numbers +0.050 + 0.052i. Noise per pair ≈ 0.03.
Six pairs (two characters) cannot measure the constant; they fix the conventions (ε₁ε̄₂ and J, not conj J). The χ̄₇
side is incomplete (flagged), so its pairs are indicative only.

## 7. Frozen predictions

`freeze.py` (steps real, complex, fn, onepoint, summary; each < 3 min). Characters: real Kronecker symbols (PARI);
complex: the r121 keys N%d_m%d (Conrey label, PARI `chareval`), one-sided functions χ and χ̄ (the mirror).
- `pred_OO.npz` (sha256 884866ca…c4a8): the 546 real functions in the order R120 (r117 + r120, 181) + B (184) + Dn (181);
  all 148,785 pairs I < J; J_σ, pred_σ = κκJ/4π², pred1_σ (H ≡ 1), H(0). HI = the R120 block (16,290 pairs).
- `pred_C.npz` (f1fa89ed…249a): 176 complex one-sided functions + 60 real (seeded draw, `default_rng(20261125)`,
  from R120); CC pairs (15,400: 3,828 χχ, 3,828 χ̄χ̄, 7,656 χχ̄′, 88 (χ, χ̄)), CR pairs (10,560); for each σ: J,
  pred (law), pred1 (H ≡ 1), prednoroot (κκJ/4π²), predconjJ (ε₁ε̄₂ conj J); q, a, ε, κ per function.
- `onepoint.npz` (2ee4adf4…905d): per function of OO (100–600), HI (600–1000), C (100–500), FN: the law
  −εκ(T₂ − T₁)/2π (sharp) and half of it (Hann), and the prime-side sums R (sharp, Hann).
- `pred_FN.npz` (75098b6a…ac21): the same for the used set (validation design only).
- `freeze_manifest.json`: hashes (also of hlib.py, freeze.py), the lists of discriminants and characters, the seed.

## 8. The test script

`swap2_test.py STEP --designs OO,HI,C` with STEP = complete | onepoint | x | readings | all; cached in `--out`
(default `run_data/`), resumable, refuses frozen files whose sha256 differs. Windows: OO [100, 594], HI [600, 994],
C [106, 494] (6 units inside every band edge without stored zeros, so every pair with weight above e⁻¹⁸ is present);
one-point bands OO (100, 600], HI (600, 1000], C (100, 500]. θ by scipy `loggamma`. Completeness: real sets the
round-138 rule on the band (count from 0 against θ/π ≤ 1.5 at 200 heights; largest unfolded gap ≤ 3.0); complex sides,
stored from 100 only: centred count ≤ 1.5, step statistic max_k |mean(s[:k]) − mean(s[k:])| ≤ 0.7 (k = 20…180), gap
≤ 3.0. Calibration (`complete_calib_out.txt`): on the 185 complete Fn sets cut at 100, no false flag (step ≤ 0.33);
a removed zero is found in 93.5 %, a removed pair in 100 % (the round-138 rule from 0: 97.8 %, 100 %). A flagged set is
excluded with all its pairs and reported. Jackknife: delete one unit (an L-function; for complex characters both
sides together) with all its pairs, exact via normal equations.
Time: validation (FN) 5 s; dry runs on all three designs 46–48 s wall, 0.52 GB peak; the real run reads the zero files
in addition (seconds).
Dry runs (`dryrun_poisson.txt`, `dryrun_gram.txt`):
- Poisson points with each function's density: every reading within ~2 se of 0 (OO slope 0.0010 ± 0.0013; complex
  slope −0.0002 − 0.0017i ± 0.004; one-point 0.002 ± 0.006): **no term**. All sets are flagged by completeness
  (counts off by ~30; synthetic flags are reported, not applied).
- Jittered Gram points (one point per Gram interval of χ at θ − arg ε/2 = π(k + ½) + πjξ with e^{−2π²j²} = κ/L, which
  reproduces the one-point law exactly): OO slope 0.156 ± 0.017 on the law (expected under H ≡ 1: 0.155), split
  [H = 1 part, H − 1 part] = [1.000, 0.001]; complex slope 0.335 (H ≡ 1 expected 0.340), H − 1 part −0.0006; one-point
  1.0006 and 0.9996: the test separates the recipe's H from independently modulated Gram lattices.

## 9. Power (noise model from Fn; `power_*_out.txt`)

Fn residuals decomposed into per-function additive and multiplicative parts and a pair part (sd 0.009, 0.044, 0.040
for Re at σ = 0.5); simulated on Fn the model gives slope sd 0.0077 against the data jackknife 0.0067 (conservative);
scaled by √(L_iL_j/L_Fn² · 500/T). Expected (law) se and the alternatives' values (σ = 0.5; σ = 1 similar):

| reading | se | law | half | none | H ≡ 1 | no root | conj J |
|---|---|---|---|---|---|---|---|
| OO R1 slope | 0.0039 | 1 | 0.5 | 0 | 0.155 | | |
| OO R2 H − 1 coefficient | 0.0035 | 1 | 0.5 | 0 | 0 | | |
| OO R3 bins (5) | 0.004–0.006 | 1 | 0.5 | 0 | 0.147–0.164 | | |
| OO R6 law < −0.01, no icpt | 0.010 | 1 | 0.5 | 0 | −0.537 | | |
| HI R1 slope | 0.008 | 1 | 0.5 | 0 | 0.153 | | |
| C1 complex slope | 0.0065 | 1 | 0.5 | 0 | 0.340 | 0.115 | 0.728 |
| C2 no-root coefficient | 0.004 | 0 | 0 | 0 | −0.031 | 1 | 0.038 |
| C3 i Im J coefficient | 0.008–0.011 | 1 | 0.5 | 0 | 0.001 | 0.016 | −1 |
| C4 (χ, χ̄) slope | 0.026–0.035 | 1 | 0.5 | 0 | 0.351 | −0.027 | 0.768 |
| C5 complex × real | 0.008 | 1 | 0.5 | 0 | 0.376 | 0.172 | 0.777 |

## 10. Prior art (mathematical queries; opened sources marked)

- Arias de Reyna, arXiv:1403.7126 (**opened**, ar5iv): Theorem 3.1 for κ = 1: Σ_{0<γ≤T} e^{2iϑ(γ)} = −ψ(T/2π) +
  O(T^{1/2}log²T), by Cauchy's theorem on ζ′/ζ·χ(s)^{−κ} and stationary phase (Levinson, Gonek); cites Fujii 1992
  (an RH-conditional formula of this type), Schoißengeier 1979 (similar results). No Dirichlet L-functions.
- Banks, "Comparing zeros of distinct Dirichlet L-functions", arXiv:2302.07073 (**opened**): Gonek's uniform Landau
  formula for L(s, χ) (main term −(T₂ − T₁)/2π·Λ(x)χ(x)); the functional-equation (Gauss-sum) part is an error term;
  no phase-weighted sum, no correlations.
- Akbary–Murty (uniform distribution of (αγ) mod 1 for a class of L-series), Weishäupl (weak Gram law for Hecke
  L-functions, arXiv:2201.08898), Saraeb (Gram points of real-character L-functions, arXiv:2412.13438): search snippets
  only; none states the χ-version Σ e^{2iθ_χ(γ)} ~ −ε_χ κ T/2π.
- Cross term: Conrey–Simons–Smolyarenko, "Pair correlation of different L-functions" (aimath.org/~kaur/publications/u8.pdf:
  HTTP 404 again; the indexed snippet says it derives the correlation of the zeros of two L-functions from a ratios
  conjecture, same level for simplicity); Bogomolny–Leboeuf 1994 set the cross term to zero (round 141); ratios tests
  (Goes–Jackson–Miller et al., the quadratic-character experiment, Fiorilli–Miller) concern one-level densities of
  families, not cross terms between two given L-functions.

Verdicts. One-point χ-law with ε_χ and κ: **KNOWN IN SUBSTANCE** (the ζ theorem's proof transfers, the mean coming
from Σ_b χ(b)e(−b/q) = χ(−1)τ(χ); the χ statement and the per-function prime-sum form with i^a e(−n/q)/√q were not
found in print). Swapped cross term with ε₁ε̄₂ and H₁₂: **KNOWN IN SUBSTANCE** for equal levels (very likely in the
unpublished CSS preprint, unopened); for unequal conductors with the phase-cancelling statistic, the (χ, χ̄) mirror case
and the HL reading: **NOT FOUND**. Any numerical test of a cross term on zeros: **NOT FOUND**.

## 11. Doubts

1. The constant 1/4π² and the "+ c.c." structure come from the recipe; only the H ≡ 1 limit is anchored in a theorem
   (the one-point law, §2). The retrospective 1.004 ± 0.007 on Fn is the only data support; the review chose X₁₂ after
   seeing the window statistic fail, so Fn is not a test.
2. The recipe's finite-height error is unknown. Known small terms: d̄d̄ leakage ≤ 2·10⁻³ per pair (σ = 0.5, equal
   conductors only), Heisenberg-range diagonal orbits ~5·10⁻⁴ per pair with mean 0, J's numerics ≤ 5·10⁻⁵. Precedent:
   the self two-point amplitude 0.998 ± 0.014 (round 136). A 1 % allowance is stated outside the bands.
3. The noise model is fitted on real characters with conductors 1000–1300; complex characters (smaller conductors,
   two sides from one L-function, possibly correlated) may be noisier. Decisions use the data's own jackknife.
4. The jackknife with dyadic (pair) data is standard but not proved unbiased here; deleting both sides of a complex
   character is the conservative unit.
5. The completeness rule for complex sides is weaker than round 138's (93.5 % of single missing zeros); the per-side
   one-point data check |G − R| (Hann) is a second, very sensitive diagnostic (a missing zero moves G by up to 1, the
   complete-set rms is 0.23); it is reported, not used to exclude.
6. Height dependence of the law is not modelled (leading order has none); HI tests it.
7. The 546 OO sets were used in earlier rounds for other statistics (counts, Σ², round 141's noise computation on
   window covariances). X₁₂ and the one-point sums were never computed on them, nor on r134_C or r121.
8. The convention pilot (6 pairs) used non-reserved data and supports ε₁ε̄₂J over the alternatives; it slightly
   un-blinds the convention, not the constant.
9. GRH is assumed (a zero off the line would not be in the lists). θ is computed exactly (scipy `loggamma`); the
   leakage estimates in doubt 2 are order-of-magnitude estimates, not computations.
10. The HL interpretation (§3) is heuristic; a direct test of the prime-pair statement on primes was not built.
