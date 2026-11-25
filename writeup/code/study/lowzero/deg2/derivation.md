# Degree 2: the two-point function of the zeros of an elliptic-curve L-function, and its number variance

Files (all in this folder): `deg2lib.py` (the law), `deg2_curves.py` → `curves_all.json` (candidates), `deg2_select.py` →
`dev_set.json`, `fresh_set.json`; `deg2_zeros.py` (zeros, completeness, repair); `deg2_predict.py` (frozen per-curve
predictions); `deg2_test.py` (readings); `validate.py` → `validate_out.txt` (independent recomputation); `contrasts.py` →
`contrasts_out.txt` (per-prime contrasts by traversal); `cue_noise.py` → `cue_noise_840.json`; `make_noise_model.py` →
`noise_model.json`; `power.py` → `power/power_out.txt`; `registry.py` → `registry_out.txt`; `deg2cm.py`, `cm_pilot.py` →
`cm_pilot_out.txt` (CM pilot); development outputs `dev_zeros/`, `dev_predictions.json`,
`dev_readings.json`, `dev_test_out.txt`; the draft `round_draft.txt`; the frozen `fresh_predictions.json`.

Labels: **proved** (algebra from the stated recipe), **conjecture** (the recipe itself), **computed**, **retrospective**
(compared with data that were looked at: the development set).

## 0. Sources, read

- **D. J. Smith, thesis** (Bristol 2016, arXiv:1607.03297; TeX source read, chapters 1–3 and the numerics). Section
  "The 2-point correlation statistic for the Selberg class": the Conrey–Snaith contour method applied to the ratios
  conjecture of a function F in the Selberg class gives (Heuristic 1)
  R₂ = R₂ᵈ + R₂ᵒ, R₂ᵈ(x) = (1/2π²)[((F⊗F̄)′/(F⊗F̄))′(1 + ix) − B(ix; F)],
  R₂ᵒ(x) = (1/(2π² r²))·(F⊗F̄)(1 + ix)(F⊗F̄)(1 − ix)·e^(−2πix·d̄)·A(ix; F) (real part understood),
  with F⊗F̄ = ∏_p exp(Σ_l l|b(p^l)|² p^(−ls)), r its residue at 1, A an Euler product over the four-fold local sums
  Σ_(h+m=k+n) a(p^m)a*(p^n)μ(p^h)μ*(p^k)p^(…) times exp(Σ_l l|b(p^l)|²(2p^(−l) − p^(−l(1−s)) − p^(−l(1+s)))), and B the
  matching local sum. For E it lists α_p = e^(iθ_p), β_p = e^(−iθ_p) (good p) and α_p = λ(p), β_p = 0 (p | N).
  Two typographical slips in the thesis: the bracket around (… − B) is missing in Heuristic 1, and the functional
  equation of L(s, E) has L(s, E) on both sides. Neither affects the formula.
- **H. M. Bui, J. P. Keating, D. J. Smith**, "On the variance of sums of arithmetic functions over primes in short
  intervals and pair correlation for L-functions in the Selberg class", J. London Math. Soc. 94 (2016), arXiv:1506.03741
  (TeX source read, sections 4–5): Conjecture RCF (ratios conjecture for F ∈ S) and Conjecture RCd (the product of two
  logarithmic derivatives) are exactly the formula above; the paper uses it for the prime side only.
- Both works check the formula **on the prime side only** (variance of Σ Λ_F(n) in short intervals for ζ, Ramanujan τ,
  and the curves 27a, 37a, 37b; arithmetic correlations for 11a). Neither computes zeros. Confirmed by reading the
  numerics sections and figure captions.

## 1. Normalisation, density and unfolding

L(s, E) = Σ λ(n)n^(−s), λ(n) = a_n/√n; Λ(s) = N^(s/2) Γ_C(s + ½) L(s, E) = w Λ(1 − s), Γ_C(s) = 2(2π)^(−s)Γ(s). Zeros
½ + iγ. Local roots: good p: α_p, β_p = e^(±iθ_p), α_p + β_p = λ_p = a_p/√p, α_pβ_p = 1; multiplicative p (p ∥ N):
α_p = a_p/√p = ±p^(−1/2), β_p = 0; additive p: none.

- θ_E(t) := arg[N^(it/2) Γ_C(1 + it)] = (t/2)ln N − t ln 2π + Im lnΓ(1 + it). The argument principle on the
  rectangle (Λ entire, zeros symmetric under γ → −γ, r = order of the central zero) gives
  **N⁺(t) = #{0 < γ ≤ t} = θ_E(t)/π − r/2 + S(t)** (**proved** up to the standard S(t) bookkeeping). Checked: the mean
  of N⁺(t) − θ_E(t)/π + r/2 over 202 heights is between −0.02 and +0.05 on the 20 development sets (r = 0, 1, 2);
  `lfunzeros` returns the central zero r times (checked on 37a, 389a, 5077a, r = 1, 2, 3).
- Density: d̄(t) = θ_E′(t)/π = (1/2π)·L(t), **L(t) = ln N − 2 ln 2π + 2 Re ψ(1 + it) = ln(N t²/(2π)²) + O(t⁻²)**.
  So d̄ = 2 × (1/2π) ln(√N·t/2π): twice the density of a degree-1 L-function of conductor √N (**proved**). At t = 20–400
  and N = 300–3000, L = 8.6–15.9; the band averages are 12.9–14.7 (development set).
- Unfolding u = θ_E(γ)/π. Then e^(−iεL) = e^(−2πir) exactly in unfolded units, r = d̄ε (as for degree 1).

## 2. The ratios recipe for one curve (conjecture), and the diagonal term (proved from it)

R(α, β; γ, δ) := ⟨L(½+α+it)L(½+β−it)/(L(½+γ+it)L(½+δ−it))⟩_t = Y(α, β; γ, δ) + q^(−α−β)·Y(−β, −α; γ, δ) + error,
q = N t²/(2π)² (from X(½+α+it)X(½+β−it) = q^(−α−β)(1 + O(1/t)); w² = 1 drops out), and
Y = ∏_p Y_p, Y_p = Σ_(m+h=n+k) a(p^m)a(p^n)μ(p^h)μ(p^k) p^(−m(½+α)−n(½+β)−h(½+γ)−k(½+δ)), a, μ the coefficients of
L_p and L_p^(−1) (all real: E is self-dual).

**Lemma 1 (proved).** Y_p(α′, β; α, β) = 1 for all α′. *Proof:* summing over n + k = j at δ = β gives
Σ_(n+k=j) a(p^n)μ(p^k) = [X^j] L_p(X)L_p(X)^(−1) = δ_(j0), which forces j = 0, i.e. m = h = n = k = 0.
Hence ∂_α ln Y = ∂_β ln Y = 0 on the diagonal γ = α, δ = β, and

  J(α + β) := ∂_α∂_β Y|_(γ=α,δ=β) = Σ_p Σ_j (ln p)² (Σ_(m+h=j) m a(p^m)μ(p^h))² p^(−j(1+α+β)).

Since Σ_j X^j Σ_(m+h=j) m a(p^m)μ(p^h) = X L_p′/L_p = Σ_j (α_p^j + β_p^j)X^j:

  **J(x) = Σ_p Σ_(k≥1) (ln p)² W_p(k) p^(−k(1+x)),  W_p(k) = (α_p^k + β_p^k)²** — the orbit p^k with weight
  |α_p^k + β_p^k|²: good p: W_p(k) = (2 cos kθ_p)² (k = 1: a_p²/p; k = 2: (a_p²/p − 2)²); multiplicative p:
  W_p(k) = p^(−k); additive p: 0.

This is Smith's ((F⊗F̄)′/(F⊗F̄))′ − B: the first is Σ k(ln p)²W_p(k)p^(−ks), and B is the Σ(k − 1)(…) part (checked
term by term; `validate.py` (b) to 1e-15). The statement "a(p) = 0 is missing at k = 1 and doubled at k = 2" is
W_p(1) = 0, W_p(2) = (0 − 2)² = 4 (amplitude |α² + β²| = 2) (**proved**); |a_p| → 2√p gives W_p(1) → 4. For p = 2:
|a_2| = 0, 1, 2 give (W(1), W(2)) = (0, 4), (½, 9/4), (2, 0); for p = 3: |a_3| = 0, 1, 2, 3 give (0, 4), (⅓, 25/9),
(4/3, 4/9), (3, 1). The weights depend on |a_p| only (θ → π − θ). Sato–Tate means: E W(1) = 1, E W(k ≥ 2) = 2.

## 3. The off-diagonal term (proved from the recipe)

Differentiating q^(−α−β)Y(−β, −α; γ, δ) at γ = α, δ = β: only the term that differentiates the two vanishing
factors 1/ζ(1 − α + γ), 1/ζ(1 − β + δ) survives (as for ζ). With x = α + β:

  K(x) = ∏_p O_p(x)/(1 − 1/p)²,  **O_p(x) = Y_p(−β, −α; α, β) = Σ_(m+h=n+k) a(p^m)a(p^n)μ(p^h)μ(p^k) p^(−(m+h)) p^(x(m−k))**

(the exponent depends on α, β only through x, because m − k = n − h). With c_j(y) = Σ_(m+h=j) a(p^m)μ(p^h)y^m:
O_p(x) = Σ_j p^(−j(1+x)) c_j(p^x)², and by partial fractions of (1 − αX)(1 − βX)/((1 − αyX)(1 − βyX)) a closed form:

  O_p = 1 + A₁²u/(1 − u) + 2A₁A₂v/(1 − v) + A₂²w/(1 − w), u, v, w = (α², 1, β²)·p^(x−1),
  A₁ = (1 − 1/y)(αy − β)/(y(α − β)), A₂ = (1 − 1/y)(βy − α)/(y(β − α)), y = p^x;

  multiplicative p: O_p = 1 + (1 − p^(−x))² p^(x−2)/(1 − p^(x−2)); additive p: O_p = 1.

Checks: for ζ (a(p^m) = 1, μ(p) = −1) this gives O_p = (1 − 2/p + p^(−1−x))/(1 − p^(x−1)) and
O_p·(1 − p^(x−1))(1 − p^(−1−x))/(1 − 1/p)² = Conrey–Snaith's A_p(x) (**proved**, by hand); for a Dirichlet prime p | k,
O_p = 1 and the factor is (1 − 1/p)^(−2), Bogomolny–Keating's frozen bad-prime factor. Series, closed form and the
literal four-fold sum agree to 2e-15 (`validate.py` (a)).

Polar part. O_p = 1 + λ_p²(p^(x−1) + p^(−x−1) − 2/p) + O(p⁻²), so

  **K(x) = ζ(1 + x)ζ(1 − x) · L(1 + x, Sym²E) L(1 − x, Sym²E)/L(1, Sym²E)² · A⊗(x)**,

A⊗ = ∏_p O_p(x)·L⊗_p(1)²/(L⊗_p(1 + x)L⊗_p(1 − x)), L⊗_p = ζ_p·L_p(Sym²), an absolutely convergent product (terms
O(p⁻²)). L(1, Sym²E) is the residue of ζ(s)L(s, Sym²E) at 1. This is Smith's (1/r²)(F⊗F̄)(1 + x)(F⊗F̄)(1 − x)A
(his F⊗F̄ omits ζ_p at bad primes; K is the same). **The off-diagonal amplitude of an elliptic curve carries
|L(1 + iε, Sym²E)|²/L(1, Sym²E)²**: a global, curve-specific factor with no degree-1 analogue. Sym² local factors:
good p: [(1 − α²X)(1 − X)(1 − β²X)]^(−1); p ∥ N: (1 − X/p)^(−1), X = p^(−s) (PARI `lfunsympow`: a_S(p) = a_S(p²) = 1 in
the arithmetic normalisation; L(2, Sym²E) agrees with the Euler product to 1e-9, `validate.py` (c);
`lfuncheckfeq` −49 to −51 bits).

Similarly J(x) = (ζ′/ζ)′(1 + x) + (L′/L)′(1 + x, Sym²E) + D(x), D = Σ_p D_p absolutely convergent:
good p: D_p = −(ln p)²[ρ(α²z) + 2ρ(z) + ρ(β²z)], ρ(u) = u²/(1 − u)², z = p^(−1−x);
p ∥ N: D_p = −(ln p)²[ρ(z/p) + z/(1 − z)²].

## 4. The law

  **R(r) − 1 = (2/L²)·Re[ J(iε) + e^(−iεL) K(iε) ],  ε = 2πr/L,  L = ln(N t²/(2π)²)**  (conjecture: the ratios
  recipe; Bui–Keating–Smith Conjecture RCd, Smith Heuristic 1, specialised to E with all local factors).

- GUE limit (**proved**): J(iε) = −1/ε² + O(1), K(iε) = 1/ε² + O(1) (A⊗(0) = 1, O_p(0) = 1), and with ε = 2πr/L,
  (2/L²)(−1 + cos εL)/ε² = −(sin πr/πr)².
- No 1/ε term: O_p′(0) = ln p·Σ_j p^(−j)[2(Σ maμ)(Σ aμ) − j(Σ aμ)²] = 0 by Lemma 1's identity (**proved**), so
  K(iε) has no iκ/ε term that would give an L-linear R(0).
- R(0) = 0: the O(1) terms cancel; computed |C(r)| ≤ 5e-7 at r = 10⁻⁴ (0.12-size terms cancelling; N = 1906, 5077;
  `validate.py` (d), `test_lib.py`). **Computed**, not proved here.
- Root number and rank do not enter (w² = 1); the rank enters only the count offset −r/2 and the central zeros,
  outside the band.
- **Symmetry type:** one L-function in the t-aspect has unitary (GUE) local statistics whatever the symmetry type of a
  family containing it (Rudnick–Sarnak 1996 for the n-level correlations of any primitive automorphic L-function,
  restricted support). The orthogonal symmetry of elliptic-curve families concerns low-lying zeros averaged over the
  family (Katz–Sarnak); it does not apply here. The arithmetic terms above are the lower-order corrections to GUE.
- Number variance: Σ²(ℓ) = ℓ + 2∫₀^ℓ(ℓ − r)(R(r) − 1)dr = Σ²_GUE(ℓ) + 2∫₀^ℓ(ℓ − r)C(r)dr, C = R − 1 + sinc²(r).
  Averaged over the band with the zero-density weight (10 Gauss nodes in t). **Proved** for a stationary window;
  local stationarity is heuristic (as in degree 1).

Numerics (`deg2lib.py`): ζ parts on an ε grid of step 0.0025 (mpmath, 30 digits); L(s, Sym²E), its first and second
derivative from PARI (`lfuninit(…, [2, 0.05, 8.5], 2)`, ε step 0.01, 801 points, analytic s = PARI s − 1); D and ln A⊗
over p < 10⁵ (vectorised); cubic splines; Gauss–Legendre 8 nodes on pieces of 0.025 up to ℓ = 10. Checks
(`validate.py`): no-spline trapezoid recomputation with the closed-form O_p agrees to 7e-7 (ℓ ≤ 4) and 5.5e-5 (ℓ = 10);
p < 3·10⁵ against p < 10⁵ moves Σ² by ≤ 9e-6. About 18 s per curve (PARI 5 s, Euler products 12 s).

## 5. What a prime does (per-prime decomposition)

Changing a_p at one prime with everything else held changes Σ² by
- diagonal, traversal k: **ΔΣ²_k(ℓ) = (2/π²)·ΔW_p(k)·sin²(πkℓ ln p/L)/(k²p^k)** (exact; **proved** as in degree 1);
- off-diagonal: 2∫₀^ℓ(ℓ − r)(2/L²)Re[e^(−2πir)K(iε)(T′_p/T_p − 1)]dr, T_p = O_p/(1 − 1/p)² (host-dependent, small).

Single-prime contrasts at L = 14 (host 2413a-type curve [1,−1,1,−10,14]; `contrasts_out.txt` has L = 12 and a second
host; off-diagonal parts differ between hosts by ≤ 0.004):

```
l                       1.0     2.0     4.0     5.0     6.0     8.0    10.0
p=2 |a|=0 vs 2   k=1  -0.0049 -0.0190 -0.0688 -0.0998 -0.1308 -0.1818 -0.2026
                 k=2  +0.0047 +0.0172 +0.0455 +0.0506 +0.0464 +0.0187 +0.0000
                 k>=3 -0.0012 -0.0039 -0.0051 -0.0029 -0.0010 -0.0015 -0.0062
                 off  -0.0000 -0.0000 -0.0029 -0.0048 -0.0029 -0.0044 -0.0022
                 total-0.0013 -0.0057 -0.0313 -0.0567 -0.0884 -0.1690 -0.2109
p=2 |a|=0 vs 1   total-0.0002 -0.0009 -0.0045 -0.0081 -0.0136 -0.0386 -0.0634  (k=2: +0.0021 .. +0.0222)
p=2 |a|=1 vs 2   total-0.0011 -0.0048 -0.0269 -0.0487 -0.0747 -0.1304 -0.1476
p=2 2|N vs |a|=1 total-0.0047 -0.0152 -0.0337 -0.0334 -0.0265 -0.0119 -0.0098  (k=1: 0, since W(1) = 1/2 both)
p=3 |a|=0 vs 3   total-0.0081 -0.0334 -0.1292 -0.1752 -0.2017 -0.1666 -0.0645  (k=2: +0.0038 .. +0.0167)
```

- At L = 12 the p = 2 contrast |a_2| = 0 vs 2 is −0.0018, −0.048, −0.204 at ℓ = 1, 4, 8 (k = 1 and k = 2 alone:
  −0.0002, −0.039, −0.197: the scout's numbers, which omit k ≥ 3 and the off-diagonal term).
- **The second traversal is not small**: for p = 2 it cancels the first at ℓ ≤ 1 and is +0.046 to +0.051 at ℓ = 4–6
  (L = 14), half the size of the first-traversal contrast there. "Repetitions flat" (k ≥ 2 at the Sato–Tate weight 2)
  differs from the law by up to ±0.05 per curve.
- **Less rigid than GUE**: with |a_2| = 2, |a_3| = 3 and large |a_5|, |a_7| the curve's Σ²(ℓ) exceeds GUE for
  ℓ ≈ 2–9 (e.g. 5077a: 0.72 against 0.56 at ℓ = 8, band 20–400; development curve 1441: law 0.610, measured 0.610,
  GUE 0.557). In degree 1 the arithmetic terms only lower Σ² below GUE at these ℓ.
- The between-curve spread of the predicted Σ²(ℓ) is 0.017 (ℓ = 1), 0.072 (ℓ = 4), 0.082 (ℓ = 8) (development set):
  15–20 times the per-curve noise.

## 6. Data: zeros, completeness, statistic

- PARI `lfunzeros(lfuncreate(E), 402)` (default divz = 8): 10–18 s per curve at N = 400–2700 on this busy machine
  (lfunzeros cost grows ≈ t³: to 200: 4 s, 300: 10 s, 400: 23 s for N = 1906; zeros per curve in [20, 400]: 767–879).
- Completeness (registered): s(t) = N⁺(t) − θ_E(t)/π + r/2 at 202 heights; development: max|s| = 0.95–1.49, mean
  −0.02 … +0.05, largest unfolded gap 2.08–2.73; nothing flagged, nothing added. Degree-2 S(t) fluctuates more than
  degree 1, so the raw-|s| limit 1.5 of round 138 would fire falsely; a missed close pair shifts s by −2 for the rest of
  the band, so the rule uses the running mean of s over 21 heights (limit 1.0) and max|s| ≤ 2.0, plus gap > 3.0 →
  re-search with divz = 64, then 256.
- Σ²_i(ℓ): u = θ_E(γ)/π, windows every 0.01 unfolded units, np.var of counts, ℓ = 0.5 … 10.

## 7. Development set (20 curves, used; retrospective)

One curve per cell (|a_2| ∈ {0, 1, 2, 2 | N}) × (|a_3| ∈ {0, 1, 2, 3, 3 | N}); band [20, 400]; 16,808 zeros.
`dev_test_out.txt`:

| reading | development (20 curves) | law | alternatives |
|---|---|---|---|
| level, data − law, 20 lengths | −0.0015 … +0.0010 (se 0.0005–0.0012) | 0 | GUE: −0.003 … −0.17 away |
| γ_A (between-curve amplitude) | 1.002 ± 0.005 | 1 | GUE 0 |
| γ_deg1 (weights vs "all weights 1") | 1.002 ± 0.005 | 1 | 0 |
| γ_k1 (all repetitions vs Sato–Tate-flat) | 1.048 ± 0.024 | 1 | 0 |
| γ_k1p2 (doubled orbit of p = 2) | 1.051 ± 0.033 | 1 | 0 |
| γ_w1p2 (first traversal of 2, a_2²/2 vs 1) | 0.990 ± 0.008 | 1 | 0 |
| γ_offzeta (Sym² factor of the off-diagonal term) | 0.94 ± 0.10 | 1 | 0 |

Per-curve noise (sd of Σ²_i − law_i): 0.0021 (ℓ = 0.5), 0.0039 (1), 0.0047 (2), 0.0048 (4), 0.0042 (6), 0.0054 (8),
0.0053 (10). Degree 1 (rounds 134–138): 0.003 at ℓ = 8 with ≈ 880 zeros. CUE with 840 eigenvalues: 0.0030, 0.0093,
0.0138, 0.0202, 0.0269, 0.0313, 0.0357 — the zeros are 2–7 times less noisy than CUE, because most of their
long-range count fluctuation is the few small-prime orbits, whose variance converges fast and is predicted per curve.

## 8. Doubts

1. The law is the ratios recipe (conjecture), with error O(t^(−1/2+ε)) claimed; at t = 20–400 the size of the
   finite-height terms is known only empirically (development level within 0.0015).
2. The rank-r central zero and the low-lying repulsion: the band starts at t = 20 (≈ 30–40 mean spacings up); no
   model of the effect; the development curves of rank 0, 1, 2 show none.
3. Estimator bias: CUE with 840 eigenvalues gives +0.010 ± 0.005 at ℓ = 8 (40 samples, inconclusive); the zeros show
   none at the 0.001 level (development). Not modelled.
4. Semistable curves only (N squarefree, 300–3000); additive primes and CM curves are not tested. For CM curves every
   inert prime has (W(1), W(2)) = (0, 4) — a strong separate test, not registered here (needs PARI's Sym² local
   factors at additive primes).
5. The off-diagonal Sym² factor is weakly tested (γ_offzeta se ≈ 0.1 on 20 curves, ≈ 0.04 expected on 138).
6. The candidates come from equations with |a4|, |a6| ≤ 30, not a complete catalogue; a bias in the a_p distribution
   (e.g. towards torsion) is possible and is carried by the per-curve predictions.
7. Development curves were chosen by seed only; the fresh set is drawn by the same seed in the same run, before any
   fresh zero existed.
8. Smith's thesis and BKS state the formula for the Selberg class (with GRH and Selberg orthogonality); the elliptic-curve
   local factors at bad primes are specialised here.

## 9. What is whose

- **Known**: the formula (Bui–Keating–Smith 2016, Conjecture RCd; Smith 2016, Heuristic 1), from the ratios conjecture
  (Conrey–Farmer–Zirnbauer) by the Conrey–Snaith method; Bogomolny–Keating for ζ and Dirichlet L-functions.
- **Derived here** (elementary consequences): the weights W_p(k) = (α^k + β^k)² written out per prime and traversal;
  the closed form of O_p and the factorisation through L(1 ± iε, Sym²E); the bad-prime factors of E; the per-prime,
  per-traversal number-variance contrasts; "less rigid than GUE".
- **Not found in print** (search below): any test of this formula, or of any arithmetic term of the pair correlation,
  on zeros of an elliptic-curve L-function.

## 10. Fresh set, frozen predictions, power (computed)

- `fresh_set.json`: 138 semistable classes, N = 307–2929, 7 per (|a_2|, |a_3|) cell (cell (2, 3): 5), seed 20261125,
  disjoint from the development set; band [20, 400]; 114,723 zeros expected (831 per curve); mean L 13.89.
- `fresh_predictions.json`, sha256 9cadbca81e70bb95a6a3462893ccbadf10fe3cc8d08b4f5998af6deb139a4f5e (computed in
  4 chunks, 16–30 s per curve; `registry_out.txt` prints design, levels and regression-matched contrasts).
- 15 fresh curves are predicted above GUE by more than 0.02 at ℓ = 6 (largest excess +0.144 at ℓ = 5–6).
- Power (`power.py`, 5 synthetic data sets per hypothesis from the development noise model; `power/power_out.txt`):
  expected errors γ_deg1 0.002, γ_A 0.002, γ_k1 0.011, γ_k1p2 0.012, γ_k1p3 0.028, γ_w1p2 0.003, γ_w1p3 0.003,
  γ_offzeta 0.036; level 0.0002–0.0005 per ℓ. Each alternative is rejected by its own reading at 27–500 errors; under
  the law every primary decision passed in 5 of 5 draws; the level cannot see the repetitions (≤ 0.0012).
- The indicator χ² against the law is conservative (35.9 expected on 60 cells under the law): the residual of an
  indicator regression contains the curve-to-curve variation from the other primes, which the law predicts. The per-curve
  template slopes (γ) are the sharp readings.

## 11. Completeness rule, checked

Deleting a pair at 20, 50, 80, 95 % of the band of each development set (80 cases): all caught by the exclusion rule
(running mean of s over 21 heights > 1.0 or |s| > 2.0). Intact development sets: running mean ≤ 0.24, |s| ≤ 1.49.
`lfunzeros` with divz = 32 on two development curves (N = 563, 1441) gives the same zeros as divz = 8 (to the last
digit), so no close pair was missed there.

## 12. CM pilot (6 curves, retrospective, not registered)

`deg2cm.py` extends the law to additive primes (W = 0, O_p = 1, Sym² factor from PARI `lfuneuler`), so CM curves can
be treated; 57 CM classes (j = 0, 1728) with 300 ≤ N ≤ 5000 were found (`cm_candidates.json`). On 6 of them
(N = 576, 800, 1089, 1568, 2304, 3136; band [20, 400]): mean |measured − law| = 0.0027 over 120 cells, against 0.0358
for "all weights 1" and 0.0058 for "repetitions flat" (`cm_pilot_out.txt`). CM curves have a_p = 0 at every inert
prime (W(1) = 0, W(2) = 4): a strong separate test for a later round (51 classes unused).

## 13. Prior art (searched 25 Nov 2026; read = full text or TeX source, abstract = abstract only)

- Bui–Keating–Smith 2016 (read, TeX): the formula for the Selberg class (Conjecture RCd); numerics on the prime side
  only (variances for ζ, τ, a conductor-37 curve).
- Smith 2016 thesis (read, TeX): the same formula (Heuristic 1); numerics on the prime side only (27a, 37a, 37b; 11a).
- Rubinstein, "Computational methods and experiments in analytic number theory", arXiv:math/0412181 (read, TeX section
  on zero statistics): pair correlation and nearest-neighbour spacings of the first 100,000 zeros of the L-functions of
  the curves of conductor 11, 14, 15, 17, 19 (and 284,410 zeros of L(s, Δ)) against GUE; "the fit is comparable to
  ζ at the same height"; no arithmetic terms.
- Bogomolny–Leboeuf 1994 (abstract): Dirichlet L-functions and Epstein zeta, not cusp forms.
- Bogomolny–Keating 2013 (read in the degree-1 rounds): Dirichlet L-functions.
- Huynh–Keating–Snaith (arXiv:0811.2304), Miller, Young, Dueñez–Huynh–Keating–Miller–Snaith, Fiorilli–Parks–Södergren
  (abstracts): lower-order terms of the one-level density of low-lying zeros in families of elliptic curves — a
  different statistic (family, near s = ½).
- Chandee–Klinger-Logan–Li and others on GL(2) pair correlation (abstracts): families, restricted support, no
  numerics on zeros.
- Verdict. Formula: **KNOWN** (Bui–Keating–Smith; Smith), with the elliptic-curve local factors **KNOWN IN SUBSTANCE**
  (the Selberg-class formula contains them; the explicit weights per traversal, the L(Sym²E) factorisation and the bad
  primes are written out here). Numerical test on zeros: **KNOWN at leading order** (GUE: Rubinstein); test of the
  arithmetic terms (weights, doubled orbit, bad primes, Sym² factor) on zeros of elliptic-curve L-functions:
  **NOT FOUND**.
