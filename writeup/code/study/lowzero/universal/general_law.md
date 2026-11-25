# The one-point law Σ e^{2iθ_L(γ)} for all L-functions: statement, derivation, special cases, first tests

Folder: `scratchpad/universal/` (nothing in the repository was changed). Labels: **proved** (algebra from stated
formulas), **recipe** (ratios-conjecture recipe), **conjecture**, **heuristic**, **computed**, **classical** (in print;
opened only where said).

## 0. In short

**Law (recipe-level conjecture).** Let L(s) = Σ λ(n)n^{−s} be a primitive L-function of degree d (the L-function of a
cuspidal automorphic representation π of GL(d) over Q, unitary normalisation), conductor Q, gamma factor
∏_j Γ_R(s + μ_j) (μ_j real here), completed Λ(s) = Q^{s/2}∏Γ_R(s + μ_j)L(s) = ε_L·Λ̃(1 − s), Λ̃(s) = conj Λ(s̄). Put

    θ_L(t) = (t/2)·ln(Q/π^d) + Σ_j Im lnΓ((½ + μ_j + it)/2),   so that   L(½ + it) = ε_L·e^{−2iθ_L(t)}·conj L(½ + it).

Then, with zeros ½ + iγ of L counted with multiplicity (GRH assumed, or read e^{2iθ_L(α_ρ)} with ρ = ½ + iα_ρ),

    Σ_{T₁<γ≤T₂} e^{2iθ_L(γ)} = −ε_L·c_L·(T₂ − T₁)/2π + o(T),       c_L = lim_{x→0+} F_L(1 + x)/x,
    F_L(s) = Σ_n conj(λ(n))·μ_L(n)·n^{−s}   (μ_L: coefficients of 1/L),
    F_L(s) = B_L(s)/L(s, π × π̃),   so   c_L = B_L(1)/Res_{s=1}L(s, π×π̃) = B_L(1)/L(1, Ad π),

B_L(s) = ∏_p F_p(p^{−s})·L_p(s, π_p × π̃_p) an Euler product with B_p = 1 + O(p^{−2}) at unramified p (absolutely
convergent at s = 1 under Ramanujan at p). With r = a/b in lowest terms the twisted sums have lines

    Σ_{T₁<γ≤T₂} e^{2iθ_L(γ)}·r^{−iγ} = −ε_L·c_L·(ρ_L(r)/√(ab))·(T₂ − T₁)/2π + o(T),
    ρ_L multiplicative, ρ_L(p^v) = p^{|v|/2}·S_p(v)/F_p(1/p),  S_p(v) = Σ_{m−h=v} conj(λ(p^m))·μ_L(p^h)·p^{−(m+h)/2}.

For a non-primitive L = ∏ L_j^{m_j}, F_L has a zero of order Σ m_j² ≥ 2 at s = 1 and **c_L = 0**: the law is a
primitivity test. Degree 1 (ζ, Dirichlet L-functions) is a theorem (c = q/φ(q), Gauss sums included, §3); degree 2 for
elliptic curves is round 148 (c_E = B_E/L(1, Sym²E)). Every step below is the elliptic-curve derivation with d, Q, μ_j,
λ, λ̄ kept general; nothing in it depends on d.

**Tests (this folder; predictions for 33 L-functions frozen before any zero was computed, sha256
946e1e3a1436dd78e87fb53bccded7ec1b93d46a27364d2622ce3668af3bdf38).** Run on 23 (10 rejected by the completeness check:
PARI's zero finder fails on non-integer coefficients). **20 primitive L-functions in degrees 1–4: pooled slope
0.985 ± 0.021, every one within 3 se; lines r = 2, 1/2, 3, 1/3: 0.981 ± 0.015.** Degree 1: four complex characters (phase
of −ε within 1–2°). Degree 2: six level-1 forms of weights 12–26 (three with ε = −1) and two S₃ Artin representations.
Degree 3: Sym² of six elliptic curves. Degree 4: two genus-2 curves (conductors 249, 277). Products: X₁(13) = L(f)L(f̄)
gives 0 and its half-phase reading carries the sum of the f and f̄ laws (1.05 ± 0.23). The two ζ_K of cubic fields give
−62 and −25 instead of 0. The first fails the registered rule at 3.2 se. A finite-height term explains both (post hoc,
§6.3): a degree-1 factor's chirp meets the other factor's lines with aligned phases. Results in §6.

## 1. Conventions (proved; checked numerically)

PARI data [N, k, Vga] (arithmetic normalisation, Λ(s) = N^{s/2}∏Γ_R(s + Vga_j)L(s), s ↔ k − s) become Q = N,
μ_j = Vga_j + (k − 1)/2, s_an = s − (k − 1)/2. The phase X_L(s) := Q^{½−s}∏_jΓ_R(1 − s + μ_j)/Γ_R(s + μ_j) satisfies
L(s) = ε_L X_L(s) L̃(1 − s), and for real μ_j, X_L(½ + it) = e^{−2iθ_L(t)} with θ_L as above (Γ_R(s) = π^{−s/2}Γ(s/2)).
Since L̃(½ − it) = conj L(½ + it): L(½ + it) = ε_L e^{−2iθ_L(t)} conj L(½ + it), so Z_L(t) = ε_L^{−½}e^{iθ_L(t)}L(½ + it)
is real. Elliptic curves: Q = N, μ = {½, 3/2}, and duplication gives θ_L = θ_E = (t/2)ln N − t ln 2π + Im lnΓ(1 + it).

*Computed* (`s0_check.py`, `s0_check_out.txt`): for all 33 objects PARI's functional-equation check `lfuncheckfeq`
≤ −47 and the identity L(½+it) = ε e^{−2iθ}conj L(½+it) holds to ≤ 6·10^{−14} at t = 7.3, 23.1, 51.7, except the order-6
nebentypus form at t = 51.7, where PARI's uninitialised `lfun` returns a wrong value (6·10^{11}); at t = 0.5, 2, 7.3 the ratio
equals ε to 10^{−15} (`dbg2.py`). Root numbers: PARI `lfunrootres`; level-1 weights 18, 22, 26 have ε = −1 (= i^k).

Stirling (proved from lnΓ(z) = (z − ½)ln z − z + ½ln 2π + O(1/z)):

    2θ_L(t) = t·ln(Q(t/2π)^d) − d·t + φ₀ + O(1/t),   φ₀ = (π/2)Σ_jμ_j − dπ/4,
    2θ_L′(t) = ln(Q(t/2π)^d) + O(t^{−2}),   2θ_L″(t) = d/t + O(t^{−3}).

(E: Σμ = 2, d = 2, φ₀ = π/2, as in round 148.) The zero density is θ_L′/π; the number of zeros in (0, T] is
M(T) = θ_L(T)/π + O(1) + S(T).

## 2. Derivation (recipe; each step as in round 148 §2.1, with the general data)

**Step 1 (contour; proved).** For a window W on (T₁, T₂], the argument principle with X_L(s)^{−1} (holomorphic and
zero-free in −½ < σ < 3/2 when all μ_j ≥ 0; in general in a strip around the critical line) gives

    Σ_γ e^{2iθ(γ)}W(γ) = (1/π)∫e^{2iθ}θ′W dt + (1/2π)∫e^{2iθ}W·(L′/L)(½ + δ + it)dt + (1/2π)∫e^{2iθ}W·(L̃′/L̃)(½ + δ − it)dt.

The first term is a boundary term O(1) (e^{2iθ}θ′ = (1/2i)(e^{2iθ})′). The third term has phase 2θ(t) + t ln n, with no
stationary point; it averages out.

**Step 2 (recipe for L′/L).** (L′/L)(s + α) = ∂_α[L(s + α)/L(s + γ)] at γ = α. The approximate functional equation
L(s + α) ≈ Σ_m λ(m)m^{−s−α} + ε_L X_L(s + α) Σ_m λ̃(m) m^{−(1−s−α)}, with λ̃(m) = conj λ(m), and 1/L(s + γ) = Σ_n μ_L(n)n^{−s−γ}.

**Step 3 (the phase cancels; proved for the swap factor).** On s = ½ + it, e^{2iθ(t)}·X_L(½ + α + it) =
X_L(½+it)^{−1}X_L(½+α+it) = exp(−α·2θ′(t) + O(α²/t)) = (Q(t/2π)^d)^{−α}(1 + O(|α|/t)). The unswapped term keeps
e^{2iθ(t)} and averages out; the swapped term loses its phase.

**Step 4 (diagonal).** The t-average of Σ_m Σ_n λ̃(m)μ_L(n) m^{−½+α+it}n^{−½−γ−it} keeps m = n:
Σ_n conj(λ(n))μ_L(n)n^{−1+α−γ} = F_L(1 + γ − α). **Which coefficient is conjugated:** the swapped L(1 − s − α) is the
dual L̃, so conj λ sits on the numerator side and μ_L of L itself on the denominator side; ε_L is the root number of L
(PARI's convention).

**Step 5 (the constant).** The integrand of Step 1 becomes (1/2π)∂_α[ε_L (Q(t/2π)^d)^{−α}F_L(1 + γ − α)] at γ = α =
(1/2π)ε_L[−ln(Q(t/2π)^d)F_L(1) − F_L′(1)]. Since F_L(1) = 0 (Step 6), this is −ε_L F_L′(1)/2π, **independent of t and of
the gamma factor** (exactly, not only to leading order: Step 3's factor equals 1 at α = 0). Integrate over (T₁, T₂]:
−ε_L c_L (T₂ − T₁)/2π with c_L = F_L′(1) = lim F_L(1 + x)/x.

**Step 6 (F_L = B_L/L(s, π×π̃); proved locally).** At an unramified p with Satake parameters α_1..α_d (|α_i| = 1 under
Ramanujan): λ(p^k) = h_k(α), μ_L(p^k) = (−1)^k e_k(α), so F_p(X) = Σ_k (−1)^k h_k(ᾱ)e_k(α)X^k = 1 − |λ(p)|²X + O(X²),
while L_p(X, π×π̃)^{−1} = ∏_{i,j}(1 − α_iᾱ_jX) = 1 − |λ(p)|²X + O(X²). Hence B_p := F_p·L_p(π×π̃) = 1 + O(X²)
(*computed* for random unitary Satake sets, d = 1..5: |B_p − 1|/X² ≤ 96, `check_Freal.py`). Since L(s, π×π̃) has a
simple pole at s = 1 for cuspidal π (Jacquet–Shalika; Rankin–Selberg theory) and L(s, π×π̃) = ζ(s)L(s, Ad π):

    F_L(1) = 0,   c_L = B_L(1)/Res_{s=1}L(s, π×π̃) = B_L(1)/L(1, Ad π).

Bad primes: B_p = F_p(1/p)·L_p(1, π_p×π̃_p) with the true local Rankin–Selberg factor; any L-function A with
A(s) = L(s, π×π̃)/ζ(s) at all p may stand for L(s, Ad π) (code: `c = ∏_{p≤P₀}F_p(1/p)A_p(1/p)/(1 − 1/p) / A(1)`).
Equivalent conditionally convergent form: c_L = ∏_p F_p(1/p)/(1 − 1/p) (limit of partial products; it converges by the
prime number theorem for Ad π).

**Reality.** For d ≤ 2 (any central character) and for every self-dual L, F_p is real and c_L real (d = 2 with
nebentypus χ: F_p = (1 − X)(1 − (|λ_p|² − 1)X), using λ̄_p = χ̄(p)λ_p). For non-self-dual L of degree ≥ 3, F_p is
complex in general (*computed*: |Im F_p(0.37)| up to 0.44 for random unitary Satake sets, `check_Freal_out.txt`); then
c_L is complex and c_L̃ = conj c_L. **The law for L̃ is −ε̄_L·conj(c_L)·(T₂ − T₁)/2π**, the complex conjugate of L's
(its zeros at positive height are the conjugates of L's zeros at negative height).

**Step 7 (lines).** With the factor r^{−it}, r = a/b, the t-average keeps m = ak, n = bk, and S_r(α, γ) =
Σ_k conj(λ(ak))μ_L(bk)(ak)^{−½+α}(bk)^{−½−γ} still contains 1/L(1 + γ − α, π×π̃) up to finitely many Euler factors;
∂_α at γ = α gives the factor ρ_L(r)/√(ab) above (as round 148, λ → conj λ in the numerator slot). Leading values:
ρ_L(p) ≈ conj λ(p), ρ_L(1/p) ≈ −λ(p).

**Status.** Steps 1, 3, 6 are proved; Steps 2, 4, 5, 7 are the recipe (the swap/diagonal heuristic of the ratios
conjecture), and the o(T) error is the recipe's T^{½+ε}. For d = 1 the result is a theorem (§3).

### 2.1 The prime side: a prime number theorem for the degree-d resonance (proved for the smoothed identity in d = 2; same proof in degree d; constant heuristic)

Explicit formula with h(t) = e^{2iθ(t)}W(t): Σ_γ h(γ) = O(1) − (1/2π)Σ_n Λ_L(n)n^{−½}∫W e^{i(2θ(t) − t ln n)}dt + (non-
stationary), −L′/L = ΣΛ_L(n)n^{−s}. Stationary point 2θ′(t_n) = ln n: n = Q(t_n/2π)^d, t_n = 2π(n/Q)^{1/d}; phase value
2θ(t_n) − t_n ln n = −d·t_n + φ₀; 2θ″ = d/t_n. Hence

    Σ_γ e^{2iθ(γ)}W(γ) ≈ −(e^{i(φ₀+π/4)}/√d)·Σ_n Λ_L(n)·n^{−½}(n/Q)^{1/(2d)}·e(−d(n/Q)^{1/d})·W(t_n).

(d = 2, Q = N, φ₀ = π/2: round 148's −(e^{3πi/4}/√2)N^{−1/4}ΣΛ_E(n)n^{−1/4}e(−2√(n/N)).) The kernel e(−d(n/Q)^{1/d}) is
the GL(d) Voronoi kernel at dual index 1. The law is equivalent (for smoothed sums) to

    Σ_{n≤x} Λ_L(n)·e(−d(n/Q)^{1/d}) ~ ε_L c_L e^{−i(φ₀+π/4)}·(2√d/(d + 1))·(x^{d+1}/Q)^{1/(2d)}       (conjecture, d ≥ 2)

(differentiate: dT/2π = (1/d)(n/Q)^{1/d}dn/n). d = 2 gives (2√2/3)e^{−3πi/4}w c_E(x³/N)^{1/4} (round 148). Over all
integers the same sum (with λ in place of Λ_L) resonates with "c = 1" by GL(d) Voronoi summation; Λ_L = μ_L ∗ (λ ln)
turns that into F_L′(1) = c_L (type I heuristic).

## 3. Degree 1: the law is a theorem, with c = q/φ(q) (proved here from PNT in progressions; checks the conventions)

χ primitive mod q, a = (1 − χ(−1))/2, μ = a, Q = q, ε_χ = τ(χ)/(i^a√q). d = 1: φ₀ = πa/2 − π/4, so e^{i(φ₀+π/4)} = i^a,
t_n = 2πn/q, kernel e(−n/q), amplitude −i^a q^{−½}. By the prime number theorem in progressions
Σ_{n≤x}Λ(n)χ(n)e(−n/q) = Σ_{b mod q}χ(b)e(−b/q)ψ(x; q, b) = χ(−1)τ(χ)x/φ(q) + o(x). With x = qT/2π:

    Σ_{0<γ≤T} e^{2iθ_χ(γ)} = −i^a q^{−½}·χ(−1)τ(χ)·(qT/2π)/φ(q) + o(T) = −ε_χ·(q/φ(q))·T/2π + o(T),

since τ(χ)/√q = i^a ε_χ and χ(−1) = (−1)^a = i^{2a}, so i^a χ(−1)τ(χ)/√q = i^{4a}ε_χ = ε_χ. **The general formula agrees:** λ = χ, μ_L = μχ, F(s) =
Σ_{(n,q)=1}μ(n)n^{−s} = 1/(ζ(s)∏_{p|q}(1 − p^{−s})), c = ∏_{p|q}(1 − 1/p)^{−1} = q/φ(q); Ad π is the zero representation
(L(1, Ad) = 1) and B_p = (1 − 1/p)^{−1} at p | q (F_p = 1 there, L_p(χ×χ̄) = ζ_p). ζ: q = 1, c = 1 (Arias de Reyna's
−ψ(T/2π)); real χ_D: −ε(k/φ(k))T/2π, measured 1.008 ± 0.007 in round 143 (STUDY.md).

## 4. Special cases with closed forms

**4.1 Level-1 holomorphic forms, weight k.** Q = 1, μ = {(k−1)/2, (k+1)/2}, ε = i^k. Ad f = Sym²f (conductor 1).
B_p = (1 − e₁X)/(1 − e₁X + e₁X² − X³) at X = 1/p, e₁ = λ_p² − 1. Petersson: **L(1, Sym²f) = (π/2)(4π)^k⟨f, f⟩/Γ(k)**
(⟨,⟩ with dμ = dx dy/y² over SL₂(Z)\H, PARI's `mfpetersson`; *computed*: the ratio equals π/2 to 10 digits for k = 12, 16,
18, 20, 22, 26, `s1b_out.txt`). So

    c_f = 2Γ(k)·B_f/(π(4π)^k⟨f, f⟩).   Δ: ⟨Δ,Δ⟩ = 1.035362·10^{−6}, L(1, Sym²Δ) = 0.6317929457, B = 1.4915511, c_Δ = 2.360823.

c_f < 0 is possible: F_2(½) = ½(1 − (λ_2² − 1)/2) < 0 iff λ_2² > 3 (|a(2)| > √3·2^{(k−1)/2}); this happens for one
weight-24 embedding (λ_2 = 1.7596, c = −0.0858): the law then predicts a **positive** sum for ε = +1. (For elliptic
curves λ_2² ≤ 2, so c_E > 0 always.)

**4.2 CM forms.** f = θ_ψ for a Hecke character ψ of an imaginary quadratic K: Ad f = χ_K ⊕ Ind(ψ/ψ̄), so
L(1, Ad f) = L(1, χ_K)·L(1, ψ/ψ̄) = (2πh_K/(w_K√|D_K|))·L(1, ψ^{1−c}); at inert p, λ_p = 0 and F_p(X) = 1 − X² (B_p =
1/((1 + X)(1 − X)·…) from the Ad factor). Not tested here.

**4.3 Nebentypus.** f with character χ (here level 13, order 6): F_p = (1 − X)(1 − (|λ_p|² − 1)X), real; L(s, Ad f) with
local factors 1/((1 − X)(1 − (|λ_p|² − 2)X + X²)) and 1/(1 − X) at the ramified prime (conductor 169; constructed and
checked by `lfuncheckfeq` −52). c real; **the law's phase is that of −ε_f** (ε_f = 0.9642 − 0.2650i), the dual f̄ has −ε̄_f.

**4.4 Sym²E (degree 3).** Sym²⊗Sym² = Sym⁴ ⊕ Sym² ⊕ 1 (SU(2)), so Ad(Sym²E) = Sym²E ⊕ Sym⁴E and
**c = B/(L(1, Sym²E)·L(1, Sym⁴E))** (both from PARI `lfunsympow`; at a multiplicative p the PARI factors are
ζ_p(s+1)ζ_p(s+2), the Steinberg St₃×St₃ factor without ζ_p(s)). ε = +1, Q = N², μ = {1, 1, 2}.

**4.5 Rankin–Selberg E₁×E₂ (degree 4).** Ad(π₁⊠π₂) = Ad π₁ ⊕ Ad π₂ ⊕ (Ad π₁ ⊠ Ad π₂): L(1, Ad) =
L(1, Sym²E₁)L(1, Sym²E₂)L(1, Sym²E₁×Sym²E₂). Not tested (PARI has `lfunmul`/`lfunrankin`-type tools; time).

**4.6 S₃ Artin representations (weight 1).** K a non-normal cubic field of discriminant d < 0, ρ the 2-dimensional
representation, ζ_K = ζ·L(ρ). ρ⊗ρ = 1 ⊕ sgn ⊕ ρ, so Ad ρ = χ_d ⊕ ρ and

    L(1, Ad ρ) = L(1, χ_d)·L(1, ρ) = (πh(d)/√|d|)·(2πh_K R_K/√|d|),   c_ρ = B_ρ·|d|/(2π²h(d)h_K R_K).

**4.7 Genus 2 (GSp(4), degree 4).** spin⊗spin = 1 ⊕ std₅ ⊕ ad₁₀, so L(1, Ad) = L(1, std₅)·L(1, ad₁₀); PARI has neither,
so c is the direct product (exact local factors p ≤ 3,000, first order to 10⁵, log-Cesàro mean, bias ν/(√P ln P), ν = 2).
Primitivity check: the mean of |λ_p|² over p ≤ 10⁵ is 0.990 (conductor 249) and 0.982 (277), against 1.997 for X₁(13).

**4.8 Products (not primitive): c = 0, and the cross-structure (recipe; proved from it).** For L = L₁L₂ (distinct
primitive factors), e^{2iθ_L} = e^{2iθ₁}e^{2iθ₂} and the zeros are the union. In Step 2 the ratio L(s+α)/L(s+γ) has four
swap terms; only the double swap carries X₁X₂ = e^{−2iθ_L}, and its diagonal is F_{L}(1 + γ − α) with F_L ~ (s − 1)²·const
(L×L̃ has a double pole). ∂_α at γ = α of A(α)F_L(1 + γ − α) is A′F_L(1) − A F_L′(1) = 0. Hence:

- Σ_{all zeros of L} e^{2iθ_L(γ)} = o(T) (no main term), and the same for the zeros of each factor with the product phase
  (each such sum carries the factor e^{2iθ_other}, a chirp, against the factor's own laws);
- each factor keeps its own law with its own phase: Σ_{zeros of L_j} e^{2iθ_j(γ)} = −ε_jc_j(T₂ − T₁)/2π, and on the union
  the zeros of the other factor contribute o(T) to it, so Σ_{all zeros} e^{2iθ_j(γ)} = −ε_jc_j(T₂ − T₁)/2π;
- if the two factors have the same gamma factor and conductor (X₁(13): L(f)·L(f̄)), θ_L = 2θ_f and the union with the
  half phase gives the sum of the two laws: Σ_{all} e^{2iθ_f(γ)} = −(ε_f + ε̄_f)c_f(T₂ − T₁)/2π = −2 Re(ε_f)c_f(T₂ − T₁)/2π.

**Primitivity test.** In general F_{∏L_j^{m_j}} has a zero of order Σm_j² (Selberg's orthonormality), so the law
detects primitivity: c ≠ 0 exactly for primitive L (when B(1) ≠ 0). The direct partial product shows it: for ζ_K
∏_{p≤P}F_p(1/p)/(1 − 1/p) = 1.048, 0.841 at P = 10⁴, 10⁵ (d = −23), falling like 1/ln P.

## 5. Degree, size of the fluctuations, where the law could fail

- **No step depends on d.** The main term is −ε c T/2π in every degree; only the resonant length grows
  (x = Q(T/2π)^d) and the Voronoi kernel changes (e(−d(n/Q)^{1/d})). The density of zeros grows like (d/2π)ln t, the
  main term stays ∝ T: the bias per zero is −ε c/ln(Q(t/2π)^d), smaller in higher degree but the same total.
- **Fluctuations (heuristic).** The diagonal of the prime sum gives Var ≈ (1/d)·Σ|Λ_L(n)|²n^{−1+1/d}Q^{−1/d}|W|², and
  with Σ_p|λ_p|²ln²p·g(p) ≈ ∫g(x)ln x dx (Rankin–Selberg) and n^{−1+1/d}dn = (d/2π)Q^{1/d}dt this is
  (1/2π)∫ln(Q(t/2π)^d)|W|²dt = the number of zeros M (random-phase value, any d). Data for d = 2 (round 148) show
  Var ≈ 0.23 M at these heights (the zeros are rigid at the Heisenberg frequency); the frozen noise model uses
  sd = √(0.25 M) and the decisions use the measured block errors.
- **Signal to noise** ≈ |c|·T/(2π·√(0.25 M)) ∝ |c|√T/√(d ln T): the law is equally testable in all degrees at fixed |c|,
  but c = B(1)/L(1, Ad) spreads more as d grows (L(1, Ad) has ~d² Euler-factor terms), so individual c can be small or
  negative.
- **Where it could change form:** (i) non-primitive L (c = 0, §4.8); (ii) if Ramanujan fails at some p, B_p is still a
  finite local factor (fine) but B(1) converges absolutely only with bounds θ < 1/4 toward Ramanujan (known for GL(2);
  for GL(d), d ≥ 3, beyond the general bounds); (iii) Siegel-type zeros of L(s, Ad π) near s = 1 would make c large;
  (iv) the sharp-window error of the explicit-formula method is O(T^{d/2}) (Arias de Reyna's κ = d case), which exceeds
  the main term for d ≥ 2: the sharp form is not provable this way, the smoothed form is (degree 2 proof note).

## 6. Tests (computed; order: constants → predictions frozen → zeros → sums)

**Order.** Constants from Euler products and L-values only (`s1_constants.py`, `s1b_sym2_level1.py`); predictions written by
`s2_freeze.py` to `predictions.json` at 11:58:39 UTC, **sha256 946e1e3a1436dd78e87fb53bccded7ec1b93d46a27364d2622ce3668af3bdf38**
(`freeze_sha256.txt`); the zeros folder was empty then and every zero file is younger than `predictions.json`
(`validate_out.txt` [1]). Sums by `s4_sums.py <sha256>` (refuses a changed file). Not registered and labelled post hoc:
the control-frequency noise, the product mechanism (§6.3), the noise factor κ (`s6_noise.py`).

**Completeness** (`s3_zeros.py`): `lfunzeros` at the default subdivision, then the 30 largest unfolded gaps re-searched
with divz 64/256 (three rounds), duplicates removed; drift of n − θ(γ_n)/π between the first and last 15 %. All 24 sets used
have |drift| ≤ 0.054 and n − θ(T₂)/π equal to the expected offset (½ for ε = +1, 0 for ε = −1, about 1.5–2 for ζ_K whose
pole adds 1) within 0.6; the gap search added 2 zeros (χ mod 5) and 16 (X₁(13), close f/f̄ pairs, minimum gap 0.0003).
The cubic-field zeros split exactly into the 138 zeros of ζ below 300 and the zeros of L(ρ) computed separately (424/424,
438/438 matched to 10⁻⁶).

**Not run (zero lists rejected by the completeness check).** The six weight-2 forms with real-quadratic coefficients
(levels 23, 29, 31) and both weight-24 embeddings: PARI `lfunzeros` returns 79 zeros to 401 (W23a; about 601 exist)
and 60 to 601 (F24a; about 685), also from an explicit real coefficient vector and at 200 digits (29 to height 100
where 107 exist), although `lfuncheckfeq` passes. The order-6 nebentypus pair (N13, N13d): 247 zeros to height 100 where
about 98 exist, and `lfun` returns 10¹⁹ at t ≈ 60. Their constants are frozen in `predictions.json` (c = 0.71, 1.30, 0.040,
2.32, 1.74, 0.245; −0.0858 and 0.527; 1.0216) but not tested; the nebentypus law enters only through X₁(13) (§6.2).

### 6.1 Primitive L-functions (`results_table.md`, `results.json`)

Slope = Re(G/pred), G = Σ_{0<γ≤T₂}e^{2iθ_L(γ)}; se from 6 equal t-blocks (in brackets the frozen heuristic √(0.25M)/|pred|).

| L | d | Q | ε | c_L | T₂ | zeros | predicted | measured | slope ± se |
|---|---|---|---|---|---|---|---|---|---|
| Δ (k=12) | 2 | 1 | +1 | 2.3608 | 600 | 683 | −225.4 | −216.7 − 7.7i | 0.961 ± 0.073 (0.058) |
| k=16 | 2 | 1 | +1 | 1.1238 | 600 | 683 | −107.3 | −107.4 − 7.5i | 1.001 ± 0.146 (0.122) |
| k=18 | 2 | 1 | −1 | 0.6750 | 600 | 683 | +64.5 | +62.7 + 8.8i | 0.973 ± 0.138 (0.203) |
| k=20 | 2 | 1 | +1 | 1.1238 | 600 | 684 | −107.3 | −92.9 − 1.4i | 0.866 ± 0.126 (0.122) |
| k=22 | 2 | 1 | −1 | 1.3910 | 600 | 685 | +132.8 | +136.9 + 7.2i | 1.031 ± 0.061 (0.099) |
| k=26 | 2 | 1 | −1 | 1.8524 | 600 | 685 | +176.9 | +167.6 + 8.4i | 0.947 ± 0.024 (0.074) |
| ρ, d = −23 | 2 | 23 | +1 | 1.8534 | 500 | 787 | −147.5 | −150.4 − 10.1i | 1.020 ± 0.079 (0.095) |
| ρ, d = −31 | 2 | 31 | +1 | 1.8792 | 500 | 811 | −149.5 | −156.9 + 12.3i | 1.049 ± 0.084 (0.095) |
| Sym² 11a | 3 | 121 | +1 | 1.9447 | 200 | 388 | −61.9 | −52.7 + 17.9i | 0.852 ± 0.216 (0.159) |
| Sym² 14a | 3 | 196 | +1 | 3.3114 | 200 | 403 | −105.4 | −97.5 + 14.5i | 0.925 ± 0.139 (0.095) |
| Sym² 15a | 3 | 225 | +1 | 3.2768 | 200 | 408 | −104.3 | −97.8 − 7.7i | 0.938 ± 0.060 (0.097) |
| Sym² 17a | 3 | 289 | +1 | 1.8559 | 200 | 416 | −59.1 | −57.9 + 7.5i | 0.979 ± 0.247 (0.173) |
| Sym² 19a | 3 | 361 | +1 | 1.1688 | 200 | 423 | −37.2 | −42.6 − 15.7i | 1.144 ± 0.217 (0.276) |
| Sym² 37a | 3 | 1369 | +1 | 0.1541 | 200 | 465 | −4.9 | −2.1 − 9.3i | 0.42 ± 2.40 (2.20) |
| genus 2, N = 249 | 4 | 249 | +1 | 0.7131 | 200 | 489 | −22.7 | −33.0 + 17.9i | 1.45 ± 0.58 (0.49) |
| genus 2, N = 277 | 4 | 277 | +1 | 0.8166 | 200 | 493 | −26.0 | −31.0 − 11.6i | 1.19 ± 0.53 (0.43) |
| χ mod 5, order 4 | 1 | 5 | 0.851+0.526i | 1.25 | 1000 | 904 | −169.2 − 104.6i | −163.5 − 105.0i | 0.976 ± 0.038 (0.076) |
| its conjugate | 1 | 5 | 0.851−0.526i | 1.25 | 1000 | 904 | −169.2 + 104.6i | −167.9 + 101.5i | 0.986 ± 0.056 (0.076) |
| χ mod 7, order 3 | 1 | 7 | 0.896−0.444i | 7/6 | 1000 | 957 | −166.4 + 82.5i | −168.8 + 88.0i | 1.025 ± 0.093 (0.083) |
| its conjugate | 1 | 7 | 0.896+0.444i | 7/6 | 1000 | 957 | −166.4 − 82.5i | −170.2 − 76.1i | 1.003 ± 0.073 (0.083) |

Pooled (Σ Re(G·conj pred)/Σ|pred|²): **all 20 primitive 0.985 ± 0.021** (heuristic 0.024); GL(1) complex 0.997 ± 0.033;
GL(2) level 1 0.963 ± 0.036; S₃ Artin 1.035 ± 0.057; Sym²E (GL(3)) 0.936 ± 0.068; genus 2 (GL(4)) 1.31 ± 0.39. Every L
passes the registered rule |slope − 1| ≤ 3 se. H0 (no main term) is rejected by 46 se. Root numbers: the three ε = −1
forms give positive sums as predicted (H-noroot would give −64.5, −132.8, −176.9); the complex characters' sums lie at
152.5°, −155.9°, −147.3°, 148.8° against the law's 153.6°, −153.6°, −148.3°, 148.3° (heuristic error 4–5°), while −c
(no root number) is 180° and −ε̄c the mirror angle (`s5_extra_out.txt`). Hann-tapered slopes agree (column in
`results_table.md`; control-frequency errors). Im(G)/|pred| is noise (largest 0.79 for genus 2, 249; 0.29 for Sym² 11a).
**Lines** r = 2, 1/2, 3, 1/3 (pooled over the 20 L, heuristic se): 0.978 ± 0.031, 0.990 ± 0.024, 0.940 ± 0.040,
0.993 ± 0.035; all **0.981 ± 0.015**. The predicted line amplitudes are large (e.g. Δ at r = 2: +115.6 against the
r = 1 sum −225.4), so the line factors ρ_L(p^{±1}) = p^{|v|/2}S_p(v)/F_p(1/p) with conj λ in the numerator slot are confirmed.
**Noise** (`s6_noise_out.txt`): measured Var/M = 0.077 (along the main term) and 0.21 (across) for the primitive L
(frozen model 0.25; χ² of the slopes against the heuristic errors 6.1 on 20); degree-1 complex characters 0.015/0.021.

### 6.2 Products

| product | reading | frozen | measured |
|---|---|---|---|
| ζ_K, d = −23 (T₂ = 300) | all zeros, product phase θ_K | 0 | −62.1 − 14.2i (ĉ = 1.30 ± 0.41: **fails** the registered rule at 3.2 block se; Hann/control noise 2.1 se) |
| | ζ zeros, θ_ζ (theorem) | −47.7 | −47.3 − 0.5i |
| | ρ zeros, θ_ρ (degree-2 law) | −88.5 | −98.2 + 11.6i |
| | ζ zeros, θ_K | 0 | −4.5 − 8.5i |
| | ρ zeros, θ_K | 0 | −57.6 − 5.7i |
| ζ_K, d = −31 | all zeros, θ_K | 0 | −25.4 − 14.2i (ĉ = 0.53 ± 0.43, passes) |
| | ζ zeros θ_ζ / ρ zeros θ_ρ | −47.7 / −89.7 | −47.3 / −84.8 + 1.3i |
| | ζ zeros θ_K / ρ zeros θ_K | 0 / 0 | −3.4 + 1.2i / −22.0 − 15.4i |
| X₁(13) = L(f)L(f̄) (T₂ = 200) | all zeros, θ_J = 2θ_f | 0 | +17.2 − 2.4i (ĉ = −0.54 ± 0.27, passes) |
| | all zeros, half phase θ_f | −2Re(ε_f)c_f T₂/2π = −62.7 | −66.0 − 4.9i (slope 1.05 ± 0.23, passes) |

The cross-structure holds: each factor keeps its own law with its own phase on the same zero list (ζ: −47.3 against
−47.7 in both fields; ρ: 1.11 and 0.95 of the law), and the X₁(13) union with the half phase carries −(ε_f + ε̄_f)c_f, the
sum of the laws of the non-self-dual f and its dual (this is the only test of the nebentypus constant c_f = 1.0216: 1.05 ± 0.23).
With the product phase the leading order (c = 0) holds for X₁(13) and d = −31 but not cleanly for d = −23.

### 6.3 Post hoc: the finite-height term of a product with a degree-1 factor (derived after seeing Z23; `s5_extra.py`)

On the zeros of the factor L(ρ), the extra factor e^{2iθ_ζ(t)} is a chirp of instantaneous frequency ln(t/2π). The
one-point spectrum of the ρ-zeros has lines at r = 1/k (density −c_ρρ(1/k)/√k per dt/2π times k^{it}); the chirp meets
the line r = 1/k at t = 2πk, where the stationary phase is 2θ_ζ(2πk) − 2πk ln k = −2πk − π/4, cancelled by the Fresnel
e^{iπ/4}: **all contributions arrive with phase 1**, each of size −c_ρρ(1/k) (2θ_ζ″ = 1/t gives √(2π·2πk) = 2π√k).
Hence Σ_{ρ zeros ≤ T} e^{2iθ_K} ≈ −c_ρ Σ_{k≤T/2π} ρ_ρ(1/k), and since Σ_k ρ_ρ(1/k)k^{−s} ≈ 1/L(s, ρ) (ρ(1/p) ≈ −λ_p,
ρ(1/p²) ≈ χ_d(p)), this is a Mertens-type sum: o(T), O(T^{½+ε}) under GRH, so the recipe's leading order (0) is not
contradicted, but at T = 300 it is large when small primes do not split (x³ − x − 1: no split prime below 59, all
μ_ρ(p) ≥ 0 up to 47). The same mechanism on the Landau lines of the ρ-zero density gives Σ_{ρ zeros} e^{2iθ_ζ(γ)} ≈
−ψ_ρ(T/2π) (Arias de Reyna's −ψ(T/2π) with Λ_ρ). A degree-2 chirp (X₁(13); ζ-zeros under e^{2iθ_ρ}) meets the lines with
phases e(−2√(k/Q)), not aligned: no such term.

| | mechanism | measured |
|---|---|---|
| d = −23: ρ zeros, θ_ζ | −ψ_ρ(47.7) = +12.0 | +11.5 − 1.4i |
| d = −23: ρ zeros, θ_K | −c_ρΣρ(1/k) = −33.1 | −57.6 − 5.7i |
| d = −31: ρ zeros, θ_ζ | −0.1 | +6.7 + 0.2i |
| d = −31: ρ zeros, θ_K | −15.9 | −22.0 − 15.4i |

With a random-phase noise √M ≈ 20 for the ρ-zeros at the product frequency (κ ≈ 1: the product phase runs at about 1.4
times the ρ-zeros' own Heisenberg frequency, where they are not rigid), the residuals are 1.2 and 0.3 sd. The mechanism
is registered in `round_draft.txt` for a blind test on 30 cubic fields.

### 6.4 Validation (`validate.py`, `validate_out.txt`, `validate_out_part4.txt`; separate code, no shared library)

Hash and file order; all 23 sums recomputed with mpmath θ (|Δ| ≤ 1.8·10⁻¹¹); counts against θ(T₂)/π; constants by other
routes: the six level-1 c via ⟨f,f⟩ (`mfpetersson`) and an independent B loop (agree to 10⁻⁶); R23, R31 via the
class-number closed form (h(−23) = h(−31) = 3, h_K = 1, R_K = 0.281200, 0.382245: agree to 10⁻⁶); Sym² 11a by an
independent direct product over `ellap` to 3·10⁵ (1.9434 against 1.9447; the validation script first used a wrong
factor at the multiplicative prime and gave 1.9288, corrected); q/φ(q). Constants: exact and direct routes agree to
≤ 0.06 % on every family (`s1_out_*.txt`); the PARI-analytic Sym² of weights ≥ 18 failed its FE check (−3 … +1) and was
replaced by the arithmetic build (FE −109 … −34; `s1b_out.txt`).

## 7. Prior art (searches with mathematical terms only)

- **Degree 1 is known.** Durkan–Hughes–Pearce-Crump, arXiv:2601.18025 (2026; **opened**, html): Σ_{0<γ≤T}χ(ρ) = −T/2π +
  O(T^{½}log²T) under RH (unconditional error Te^{−a√log T}), and for integer X = o(T) the sum Σχ(ρ)X^ρ has main term −T/2π:
  this is the ζ case of the law and its lines at integer r (for ζ, ρ(k) = 1 for all k, so the line is −(T/2π)/√k). Their
  Theorem 5 is the explicit prime-sum form (−X Σ Λ(n)e^{2πiXn} over T/2πX < n ≤ T/πX), i.e. the degree-1 resonance. Built on
  Gonek's 1984 lemma (∫χ(σ+it)r^{it}dt = 2πr^{1−σ}e^{2πir}, the d = 1 stationary phase). Arias de Reyna 2014 (−ψ(T/2π), read
  in round 148). Dirichlet characters with Gauss sums: not found written in this form (it follows from the same proof).
- **Gram's law / Gram points for degree 2:** Weishäupl, "The weak Gram law for Hecke L-functions" (arXiv:2201.08898, opened
  in round 148; search result again): sums of Z_f at Gram points, ±T/π main terms from the n = 1 term; values at Gram
  points, not zeros, and no adjoint constant.
- **Landau–Gonek formulas in the Selberg class** (Σx^ρ, i.e. our lines without the phase e^{2iθ}): Murty–Perelli,
  Murty–Zaharescu, Ford–Soundararajan–Zaharescu (fractional parts {αγ}: the arithmetic term of the density of αγ mod 1 is
  Landau's −(T/2π)Λ(x)/√x; FSZ treat general L-functions), Mazhouda (arXiv:1505.01544, abstract opened), a GL(2) Landau–Gonek
  lemma in arXiv:2407.14656 (search result). These are fixed-frequency (α) statements; **the one-point law is not a case of
  them**: its weight e^{2iθ_L(γ)} is a chirp at the Heisenberg frequency, whose prime side is the GL(d) Voronoi resonance,
  not a Landau line. They do, however, give the lines' degree-1 analogue, and the mechanism of §6.3 is "Landau lines swept
  by a chirp".
- **GL(2) resonance over primes:** Zhao 2006 / Iwaniec–Luo–Sarnak (read in round 148): bound and heuristic asymptotic with an
  unnamed constant.
- Not found: any sum Σ_ρ X_L(ρ)^{−1} (= Σe^{2iθ_L(γ)}) for an L-function of degree ≥ 2; any appearance of 1/L(1, Ad π)
  (or 1/L(1, Sym²f)) as a density for a statistic of the zeros of one L-function; any c = 0 primitivity statement of this kind.
  (arXiv:2610.21952 and arXiv:2610.26575, 2026, abstracts opened: degree 1 Gram-point samples; pair/triple correlation
  surfaces: not this.) The web search was limited (one search session, abstracts mostly); MathSciNet not available.

**Verdict.** Degree 1 (ζ and, in substance, Dirichlet L-functions): **KNOWN** (theorem). The general law for degree ≥ 2 with
c_L = B_L(1)/L(1, Ad π), its lines, c = 0 for products and the product secondary term: **NOT FOUND**. The prime-side
resonance for GL(2): **KNOWN IN SUBSTANCE** as a heuristic without constant (ILS/Zhao).

## 8. Doubts

1. For d ≥ 2 the law is recipe-level; the data (20 L, 4 degrees) support it at 2 % (pooled) but degree 3 is measured to 7 %
   and degree 4 to 40 % only. A proof needs a prime number theorem for the GL(d) resonance (type II sums), beyond reach.
2. Slopes in degree 2 and 3 run slightly low (0.963 ± 0.036, 0.936 ± 0.068; round 148: 0.989 ± 0.011): within 1–2 se; the
   same sign as round 148. Lower-order terms of the recipe at T ≤ 600 are not computed.
3. The non-self-dual degree-2 test could not be run (PARI failure); the phase −ε is confirmed only in degree 1 (theorem) and,
   for its real part, through X₁(13). The claim that c_L is complex for non-self-dual L of degree ≥ 3 (F_p complex) is
   computed for random Satake sets only; no such L was tested.
4. Z23 fails the registered "0" at 3.2 block se. The explanation (§6.3) is post hoc and checked on two fields; the noise
   model for products (no rigidity) is also post hoc. A blind test is drafted.
5. The completeness rule cannot see a missed close pair in the middle of the range if the drift stays small; X₁(13) needed 16
   additions. No set used had any sign of a missing zero after the re-search (counts within 0.6 of θ/π + offset).
6. The κ = 0.25 noise model was conservative (measured 0.08 along the main term); the registered decisions use block
   errors, which are noisy with 6 blocks (e.g. weight 26: se 0.024 against the heuristic 0.074).
7. Genus-2 constants come from the direct product (two routes agree to 0.1 % elsewhere); their primitivity (End = Z) is
   inferred from mean |λ_p|² = 0.990, 0.982 (a product would give 2), not proved.
8. c_L for Sym² 37a is small (0.154) and its slope is uninformative (0.42 ± 2.4).
9. Hypotheses used: GRH (zeros on the line; PARI's lfunzeros finds only those), Ramanujan at p for the absolute convergence of
   B (true for all classes tested: Deligne, Artin, Weil), automorphy of Ad (for the direct-product route of genus 2).

## 9. Files and command lines (folder `scratchpad/universal/`)

Every job: `(ulimit -v 1500000; timeout 600 nice .venv/bin/python -I <script> <args>)`.

| step | command | output |
|---|---|---|
| objects, FE, phase identity | `s0_check.py` | `s0_check_out.txt`, `s0_check.json` |
| constants | `s1_constants.py [keys]`; `s1b_sym2_level1.py [keys]` | `constants_*.json`, `s1_out_a/b/c.txt`, `s1b_out.txt`, `s1b_out_24.txt` |
| reality of F_p, B_p − 1 | `check_Freal.py` | `check_Freal_out.txt` |
| freeze | `s2_freeze.py` | `predictions.json`, `s2_freeze_out.txt`, `freeze_sha256.txt` |
| zeros | `s3_zeros.py <keys>` | `zeros/zeros_<key>.npz`, `s3_out_1…7.txt`; rejected: `zeros_F24a_incomplete_defaultprec.npz` |
| sums | `s4_sums.py 946e1e3a…df38` | `results.json`, `results_table.md`, `s4_out.txt` |
| extra readings | `s5_extra.py 946e1e3a…df38`; `chirp_lines.py` | `s5_extra_out.txt`, `s5_extra.json`, `chirp_lines_out.txt` |
| noise | `s6_noise.py` | `s6_noise_out.txt` |
| validation | `validate.py 946e1e3a…df38 [skip2]` | `validate_out.txt`, `validate_out_part4.txt` |
| PARI diagnostics | `probe*.py`, `dbg*.py` | (console) |
| registration draft | | `round_draft.txt` |
