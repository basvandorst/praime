# Two spectra on the same orbits: the cross-correlation of the zeros of two Dirichlet L-functions

Folder: this scratch folder (nothing in the repo was changed; the fresh zeros r138_Fn were not read). Scripts and outputs:

| script | output | what it does |
|---|---|---|
| `crosslib.py` | | Kronecker tables, per-prime coefficients A_p, B_p, θ_χ, window counts δn, pair predictions |
| `theory.py` | `theory_out.txt` | (T1) per-prime numbers, (T2) correlation coefficients and spread over pairs, (T3) closed form with PARI L-values against prime sums, (T4) the recipe's two-point function integrated, (T5) the swapped term |
| `t4fix.py` | `t4fix_out.txt` | (T4) with a converged quadrature at w = 4.5 |
| `leak.py` | `leak_out.txt`, `leak_300.json` | theory only: the finite-range form of the law (orbit leakage on the registered grid) and its effect on every reading |
| `noise.py` | `noise_out.txt` (its cache of used-pair covariances was deleted after use; rerunning regenerates it) | NOISE ONLY on the 546 used real sets: completeness, residual dispersion, jackknife se, off-line projections |
| `swaplib.py`, `swap_run.py` | `swap_*_P10000_*.npz` (crude, used sets), `swap_N*_P50000_*.npz` | the swapped (Gram-lattice) term per near pair |
| `noise2.py` | `noise2_out.txt`, `noise2_se.json` | NOISE ONLY: jackknife se of every registered reading (readings.py), incl. the swap reading |
| `predict.py` | `predict_out.txt`, `predictions_cross.json`, `predictions_cross_pairs.npz` | frozen predictions (discriminants only) and the expected readings under each hypothesis |
| `assemble_swap.py` | `assemble_swap_out.txt`, `swap_pred_NN.npz`, `swap_pred_NO.npz` | frozen swapped-term predictions |
| `swapcheck.py` | `swapcheck_out.txt` | scatter of the swap predictions from the Euler-product cutoff |
| `check.py` | `check_out.txt` | independent recomputation of the frozen pair predictions (sympy, mpmath, plain loops, closed geometric sums) |
| `readings.py`, `rcross.py`, `synth.py` | `dryrun_out.txt`, `rcross_NN_dryrun_synthetic.json` | the registered readings, the analysis script, a dry run on synthetic Poisson zeros |
| `round_draft.txt` | | the pre-registration draft |

Labels: **proved** (algebra from stated formulas), **conjecture** (the ratios recipe), **heuristic** (semiclassical or
approximate steps, or constants taken over by analogy), **computed**, **classical** (in print), **noise** (computed on
used zero sets without printing any coefficient).

## Summary

- **The law (conjecture; the ratios recipe and the semiclassical diagonal sum give the same term).** For distinct
  primitive χ₁ mod k₁, χ₂ mod k₂ put ψ = χ₁χ̄₂ (ψ(n) = 0 when (n, k₁k₂) > 1). With ε = t₂ − t₁,

      R¹²(t₁, t₂) = d̄₁(t₁) d̄₂(t₂) + (1/2π²) Re D_ψ(−iε) + R¹²_swap,
      D_ψ(z) = Σ_n Λ(n)² ψ(n) n^(−1−z) = (L′/L)′(1+z, ψ) − Σ_p (ln p)² ψ(p)² p^(−2−2z)/(1 − ψ(p)p^(−1−z))².

  L(s, ψ) is entire (ψ is never principal), so the term is finite at ε = 0 and needs no Heisenberg cutoff. The recipe
  adds nothing else to the non-swapped term. The swapped term exists for every pair; its phase e^(2i(θ₂(t₂) − θ₁(t₁)))
  drifts like (k₂/k₁)^(it), so the recipe keeps it only for k₁ = k₂ (§3).
- **The covariance of window counts (proved from the law).** With δn_j(t) = #{γ_j ∈ (t, t+w]} − (θ_j(t+w) − θ_j(t))/π,

      Cov(δn₁, δn₂)(w) = (2/π²) Σ_p Σ_m Re[ψ(p)^m] sin²(m w ln p/2)/(m² p^m)
                       = (1/π²) ln( |L(1,ψ)| / √(|L(1+iw,ψ)| |L(1−iw,ψ)|) ) − (1/π²) Σ_p Σ_(m≥2) (1/m − 1/m²) Re ψ(p)^m (1 − cos(m w ln p))/p^m,

  in the window width w (height units) only: independent of the height, and of the conductors except through ψ.
- **Per prime (proved).** For real characters p adds A_p(w)·χ₁(p)χ₂(p) + B_p(w)·1[p ∤ k₁k₂] with A_p, B_p ≥ 0: counts
  of two L-functions whose characters agree at p rise and fall together, counts of characters that disagree move
  oppositely. Peaks of A_p: 0.104 (p = 2, w ≈ 4.5 ≈ 7.4 mean spacings), 0.068 (p = 3, w ≈ 2.9), 0.040 (p = 5, w ≈ 2.0),
  0.029 (p = 7, w ≈ 1.6). Against a number variance of 0.41-0.43 this is a correlation coefficient of up to 0.25, 0.16,
  0.10, 0.07 per prime. Over the 148,785 used pairs the whole law has sd 0.078 over pairs at w = 3.9.
- **The swapped (Gram-lattice) term.** For real pairs with |D₁| = |D₂| it is exactly odd and adds nothing to equal
  windows; for nearly equal conductors it survives at finite height (up to 0.046 per pair in the fresh design; rms
  0.002-0.005 over its 3,080 near pairs). It is the cross version of the Bogomolny–Keating off-diagonal term and of the
  one-point Gram term (Σ_γ e^(2iθ(γ)) = −T/2π + …, in print for ζ): two spectra locked to nearly the same Gram lattice.
- **What the data can see (noise only, used sets).** After the regression on the predictors the pair covariances scatter
  by 0.003 (0.0039 at w = 0.3), against a predicted spread of 0.008-0.079: the per-pair noise is small because the count
  fluctuations are made of the deterministic orbit terms. Jackknife se on a 186 × 186 design: γ₂ 0.0007-0.0008, γ₃
  0.0011, γ₅ 0.0018-0.0020, γ₇ 0.0020-0.0021, γ_full 0.0006-0.0007, γ_sw ≈ 0.025 (scaled). So z against "independent
  spectra" is about 1,300 (p = 2), 900, 500, 500, 1,600 (full law), 40 (Gram-lattice term). Any design can see it; the
  question is no longer detection but whether the law holds to 0.1 %.
- **Prior art.** The leading-order independence of distinct L-functions is in print (heuristically Bogomolny–Leboeuf 1994,
  who set the cross term to zero; rigorously for restricted support Rudnick–Sarnak 1996, Liu–Ye 2002; conjecturally
  Murty–Perelli 1999). A ratios-conjecture pair correlation of two different L-functions exists in an unpublished preprint
  (Conrey–Simons–Smolyarenko) that I could not open. Formula: KNOWN IN SUBSTANCE (very likely in that preprint; the
  window-covariance closed form and the per-prime sign rule not found). Numerical test on zeros: NOT FOUND.

## 1. Setting

χ₁ mod k₁ and χ₂ mod k₂ are primitive and distinct. Zeros of L(s, χ_j) are 1/2 + iγ (GRH assumed), d_j(t) = Σ δ(t − γ).
θ_j(t) = (t/2) ln(k_j/π) + Im lnΓ((1/2 + a_j + it)/2) with a_j the parity, so N_j(t) = θ_j(t)/π + const + S_j(t) and
d̄_j = θ_j′/π = (1/2π) ln(k_j t/2π) + O(t⁻²).

ψ = χ₁χ̄₂ is a character mod lcm(k₁, k₂), zero at every p | k₁k₂. It is induced by a primitive ψ* of conductor
f | lcm(k₁, k₂), and L(s, ψ) = L(s, ψ*) ∏_(p | k₁k₂, p ∤ f) (1 − ψ*(p) p^(−s)). ψ* is principal only if χ₁ = χ₂ as
primitive characters, so for distinct primitive characters L(s, ψ) is entire and L(1 + it, ψ) ≠ 0 (**classical**). For
real characters χ_D₁χ_D₂ = (D₁D₂/·) and ψ* = χ_d with d the discriminant of Q(√(D₁D₂)).

## 2. Route (b): periodic orbits (heuristic)

1. S_j(t) = (1/π) Im log L(1/2 + it, χ_j) = (1/π) Im Σ_n χ_j(n) Λ(n)(ln n)⁻¹ n^(−1/2−it) (formal). Differentiating,
   d̃_j(t) = −(1/π) Re Σ_n χ_j(n) Λ(n) n^(−1/2) e^(−it ln n): an orbit n = p^m of period ln n, amplitude Λ(n)/√n and
   phase arg χ_j(n). The orbits are the same for every L-function; only the phases differ, and an orbit p | k_j is
   missing from spectrum j.
2. With a, b the two sums, ⟨d̃₁(t) d̃₂(t + ε)⟩_t = (1/π²)⟨Re a Re b⟩ = (1/2π²) Re⟨a b̄⟩ (⟨ab⟩ oscillates like e^(−2it)).
   The diagonal n = m of ⟨a b̄⟩ is Σ_n χ₁(n)χ̄₂(n) Λ(n)²/n · n^(iε). Hence
   **R¹²_diag(ε) = (1/2π²) Re Σ_n Λ(n)² ψ(n) n^(−1+iε) = (1/2π²) Re D_ψ(−iε).**
3. Each orbit contributes (1/2π²)(Λ(n)²/n) cos(ε ln n + arg ψ(n)): the spectra interfere through each prime power with
   the phase difference m(arg χ₁(p) − arg χ₂(p)).
4. Convergence. For one L-function the diagonal sum Σ Λ²/n diverges; the Heisenberg cutoff (Berry) or the off-diagonal
   term (Bogomolny–Keating) is essential. Here ψ is non-principal: under GRH Σ_(n>x) Λ(n)²ψ(n) n^(−1−iε) =
   O(x^(−1/2) log³(x k₁k₂)), so orbits beyond the Heisenberg time hardly contribute and the diagonal approximation is
   insensitive to where it fails (checked: the prime sum to 2·10⁵ is within 1·10⁻⁴ of the exact value, (T3)).
5. Off-diagonal orbit pairs n ≠ m (Hardy–Littlewood in progressions): the average of χ₁(n)χ̄₂(n + h) has a Fourier
   transform in h that is a product of Gauss sums of the two characters; carried through, this is the swapped term of
   the recipe (§3.3), with the functional-equation phase. I did not carry the Hardy–Littlewood computation through;
   route (a) is the complete one.

## 3. Route (a): the ratios recipe (conjecture; Conrey–Farmer–Zirnbauer, as applied in `amp/derivation.md` §1.1)

### 3.1 The ratio and its swaps

    I(α, β, γ, δ) = ∫ L(s+α, χ₁) L(1−s+β, χ̄₂) / ( L(s+γ, χ₁) L(1−s+δ, χ̄₂) ) dt,  s = 1/2 + it.

With L(s, χ) = ε_χ X_χ(s) L(1 − s, χ̄), X_χ(1/2 + iτ) = e^(−2iθ_χ(τ)), ε_χ ε_χ̄ = 1, |ε_χ| = 1:

- no swap: kept;
- one swap: oscillates like (kt/2πe)^(∓it); dropped (its average over 100 ≤ t ≤ 600 is O(1/(T ln kT)) ≈ 2·10⁻⁴ relative);
- both swapped: ε₁ ε̄₂ X_χ₁(s+α) X_χ̄₂(1−s+β), and 2(θ₂ − θ₁)(t) = t ln(k₂/k₁) + (a₂ − a₁)π/2 + O(e^(−πt/2)) (the
  parity constant is exact to all orders of Stirling's series: only even-index Bernoulli polynomials enter the imaginary
  part, and B_n(1/4) = B_n(3/4) for even n). **Non-oscillating if and only if k₁ = k₂** (either parity). For k₁ ≠ k₂ the
  recipe drops it; it then drifts at the frequency ln(k₂/k₁), slowly for close conductors (§4.5).

So the answer to "no swap term when χ₁ ≠ χ₂?" is: the swap term is there for every pair and is non-oscillating exactly
when the two conductors are equal; its coefficient is the root-number ratio ε₁ε̄₂ (= 1 for real characters).

### 3.2 The non-swapped term (proved from the recipe)

The diagonal (mh = nk) of χ₁(m)·μχ₁(h)·χ̄₂(n)·μχ̄₂(k) has, with P_xy := ψ(p) p^(−1−x−y), the Euler factor

    G_p = Σ_(a+c=b+d; c,d∈{0,1}) (−1)^(c+d) ψ(p)^(a+c) x^a y^c u^b v^d = (1 − P_βγ − P_αδ + P_γδ)/(1 − P_αβ),

because χ₁(p)^(a+c) χ̄₂(p)^(b+d) = ψ(p)^(a+c) when a + c = b + d. The cases (c, d) = (0,0), (1,1), (1,0), (0,1) give 1,
yv, −yu, −xv over 1 − xu. For p | k₁k₂, G_p = 1. At γ = α, δ = β: G_p = 1 and ∂_α log G = ∂_β log G = 0, so
∂_α∂_β G = ∂_α∂_β log G; with X = P_αβ,

    ∂_α∂_β log(1 − P_βγ − P_αδ + P_γδ) = −(ln p)² X²/(1 − X)²,   ∂_α∂_β(−log(1 − X)) = (ln p)² X/(1 − X)²,
    sum = (ln p)² X/(1 − X) = Σ_m (ln p)² ψ(p)^m p^(−m(1+α+β)).

So **∫ (L′/L)(s+α, χ₁)(L′/L)(1−s+β, χ̄₂) dt ⊃ D_ψ(α+β)**, exactly the semiclassical diagonal sum (Λ(p^m)² = (ln p)²).
With L-functions: D_ψ(z) = (L′/L)′(1+z, ψ) − B_ψ(z), B_ψ(z) = Σ_p (ln p)² ψ(p)² p^(−2−2z)/(1 − ψ(p)p^(−1−z))². For
χ₁ = χ₂ it is the repo's (L′/L)′(1+z, χ₀) − B_k(z).

### 3.3 The swapped term (proved from the recipe)

After both swaps the numerator is L(1−s−α, χ̄₁) L(s−β, χ₂); the diagonal Euler factor is G′_p = (1 + YU + XV + YV)/(1 − XU)
with XU = ψ̄(p)p^(−1+α+β), YU = −χ₀⁽¹⁾(p) p^(−1+α−γ), XV = −χ₀⁽²⁾(p) p^(−1+β−δ), YV = ψ(p) p^(−1−γ−δ) (χ₀⁽ʲ⁾ principal mod
k_j). Its polar part is L(1−α−β, ψ̄) L(1+γ+δ, ψ) / ( L(1−α+γ, χ₀⁽¹⁾) L(1−β+δ, χ₀⁽²⁾) ). The two principal factors have
poles at 1, so the term vanishes at γ = α, δ = β and ∂_α∂_β survives only through 1/L(1+z, χ₀⁽ʲ⁾) = z k_j/φ(k_j) + O(z²):

    swapped part = (k₁/φ(k₁))(k₂/φ(k₂)) ε₁ε̄₂ X_χ₁(s+α) X_χ̄₂(1−s+β) L(1+α+β, ψ) L(1−α−β, ψ̄) A₁₂(α+β),
    A₁₂(z) = ∏_p (1 − ψ(p)p^(−1−z)) (1 − χ₀⁽¹⁾(p)/p − χ₀⁽²⁾(p)/p + ψ(p)p^(−1−z)) / ((1 − χ₀⁽¹⁾(p)/p)(1 − χ₀⁽²⁾(p)/p)).

For χ₁ = χ₂ this is the repo's (k/φ(k))² L(1+z, χ₀)L(1−z, χ₀) A_k(z). At z = 0 the factor 1 − 2/p + ψ(p)/p vanishes for
p = 3, ψ(3) = −1: then the swapped term is 0 at ε = 0 ((T5): H(0) = 0 for D = ±104, ±296, 1004/−1007).

### 3.4 The two-point function and its cases

By the contour argument of Conrey–Snaith §7 (as in the self case), with ε = t₂ − t₁ and α + β = −iε:

    R¹²(t₁, t₂) = d̄₁ d̄₂ + (1/2π²) Re D_ψ(−iε)
                + (1/4π²) { (k₁/φ(k₁))(k₂/φ(k₂)) ε₁ε̄₂ e^(2i(θ₂(t₂) − θ₁(t₁))) L(1−iε, ψ) L(1+iε, ψ̄) A₁₂(−iε) + c.c. }.

The constants of the last line are taken over from the self case, where they reproduce Bogomolny–Keating's
(1/4π²)e^(iεL)|ζ|²Φ_off (**heuristic in its constants**). Partial check: with H → 1 the last line is the product of two
one-point Gram terms (§3.5), whose ζ case is a theorem. The diagonal line is checked numerically against the covariance
formula in (T4).

- **k₁ ≠ k₂** (most real pairs): the recipe keeps only the first two terms; at finite height the third survives for close
  conductors (§4.5).
- **|D₁| = |D₂|, real** (D and −D with 8 | D): the phase is e^(±iπ/2) e^(2i(θ(t₂) − θ(t₁))) with 2θ(t+y₂) − 2θ(t+y₁)
  exactly antisymmetric in y₁ ↔ y₂, and H(−ε) = conj H(ε) (real ψ; (T5)), so the integrand is odd: **0** in equal
  windows at every width and height (computed: |C_swap| < 10⁻⁵ for D = ±104, ±296). It would show in lagged windows.
- **complex χ, same modulus and parity** (χ₂ = χ̄₁ included: ψ = χ₁², ε₁ε̄₂ = ε_χ₁²): non-oscillating, relative size
  2(k/φ(k))²|L(1,ψ)|²A₁₂(0)/L², phase set by the root numbers: two spectra on the same Gram lattice, offset by
  arg(ε₁ε̄₂).
- **χ₂ = χ̄₁**: the zeros of L(s, χ) at +t and at −t are correlated through the orbits with weight Re χ(p)^(2m), plus the
  root-number term ε_χ².
- **ψ imprimitive**: through L(s, ψ) above; a prime dividing either conductor is absent (an orbit missing from one
  spectrum cannot correlate).
- **Real against complex**: the formula holds for both; for complex characters the m-th repetition of orbit p enters with
  cos(m(φ₁ − φ₂)), φ_j = arg χ_j(p).

### 3.5 The one-point Gram term (derived; known for ζ)

The recipe with one swap and the weight e^(2iθ_χ(t)), which cancels the oscillation of X_χ, gives
⟨e^(2iθ_χ(t)) (L′/L)(1/2 + δ + it, χ)⟩ = −ε_χ k/φ(k) + o(1), so the zero density of one L-function carries
d̃ ⊃ −(1/π)(k/φ(k)) Re[ε_χ e^(−2iθ_χ(t))]: a modulation of relative size 2(k/φ(k))/L at period one in N̄, pushing zeros
away from the Gram points of Z_χ (Gram's law in density form). For ζ (k = 1, ε = 1) it says
Σ_(0<γ≤T) e^(2iθ(γ)) = −T/2π + o(T), which is the theorem Σ_(0<γ≤T) e^(2iθ(γ)) = −ψ(T/2π) + O(T^(1/2) log² T) of
Arias de Reyna (arXiv:1403.7126, Theorem 3.1 at κ = 1; opened). The χ-version (with ε_χ and k/φ(k)) I did not find.
The product of two one-point terms is the swapped term with H → 1.

## 4. The statistic and the covariance

### 4.1 Window counts in the common variable t

The two spectra have different densities, so the windows are defined in the height t, common to both:
δn_j(t; w) = #{γ_j ∈ (t, t + w]} − (θ_j(t + w) − θ_j(t))/π with the exact smooth count, and C₁₂(w) = ⟨δn₁ δn₂⟩_t over
t ∈ [T₁, T₂ − w]. Subtracting the exact smooth count matters: w·d̄(t) grows by w ln 6/(2π) ≈ 1.3 over 100 ≤ t ≤ 600 at
w = 4.5 for both functions, and a raw-count covariance would pick up that common trend (≈ 0.15).

### 4.2 The covariance (proved from the law)

R¹²_diag depends on t₁, t₂ only through ε, so C₁₂(w) = (1/2π²) ∫_(−w)^(w) (w − |ε|) Re D_ψ(−iε) dε. Term by term,
∫_(−w)^(w) (w − |ε|) cos(εx + φ) dε = cos φ · 4 sin²(wx/2)/x², and Λ(p^m)²/ln²(p^m) = 1/m², so

    C₁₂(w) = (2/π²) Σ_p Σ_m Re[ψ(p)^m] sin²(m w ln p/2)/(m² p^m).

With sin²(x/2) = (1 − cos x)/2 and ln L(s, ψ) = Σ ψ(p)^m p^(−ms)/m this is (1/π²){ln|L(1,ψ)| − ½ ln|L(1+iw,ψ)| −
½ ln|L(1−iw,ψ)|} minus the absolutely convergent m ≥ 2 correction. Only Re ψ enters (equal windows are symmetric in
the two functions). Large w: C₁₂ → (1/π²) ln|L(1,ψ)/L(1+iw,ψ)| + small, the log-ratio of two values on the 1-line.

**Finite range (computed, `leak.py`).** On the actual grid the orbit cosines of different prime powers are not exactly
orthogonal and the mean of a squared cosine is ½(1 + ⟨cos 2m(t + w/2) ln p⟩). The frozen predictions use the exact grid
averages for all orbits n ≤ 300 (79 prime powers) and the T → ∞ form beyond. This finite-range part has sd 0.0007-0.0016
over pairs (zero mean) and moves the pooled amplitudes by at most 6·10⁻⁴ (γ_p), 0.006 at small w, 6 % for γ_B7.

### 4.3 Per prime and sign (proved)

For real characters ψ(p)^m = (χ₁(p)χ₂(p))^m: odd m give χ₁(p)χ₂(p), even m give 1[p ∤ k₁k₂]. So p adds
A_p(w)χ₁(p)χ₂(p) + B_p(w)1[p ∤ k₁k₂], A_p = (2/π²) Σ_(m odd) sin²(m w ln p/2)/(m²p^m), B_p = the even-m part; both ≥ 0
at every w (checked for the registered grid in `check_out.txt`). The covariance is raised when the characters agree at
p and lowered when they disagree; the second traversal (m = 2) has phase χ(p)² = 1 in both spectra and always raises
it. For complex characters the p-term is (2/π²) Σ_m cos(m(φ₁ − φ₂)) sin²(m w ln p/2)/(m²p^m): characters at right
angles at p (ψ(p) = ±i) have no m = 1 correlation from p.

### 4.4 Numbers (computed, `theory_out.txt`)

Raw windows w = 0.3, 0.6, …, 4.8 (0.49-7.9 mean spacings at L = 10.33, conductor 550 at t = 350; 0.53-8.5 at
L = 11.07, the fresh design). A_p(w):

| w | 0.3 | 0.9 | 1.5 | 2.1 | 2.7 | 3.3 | 3.9 | 4.5 | 4.8 |
|---|---|---|---|---|---|---|---|---|---|
| p = 2 | 0.0014 | 0.0116 | 0.0279 | 0.0468 | 0.0662 | 0.0843 | 0.0983 | 0.1044 | 0.1033 |
| p = 3 | 0.0020 | 0.0160 | 0.0367 | 0.0565 | 0.0678 | 0.0641 | 0.0478 | 0.0267 | 0.0166 |
| p = 5 | 0.0024 | 0.0179 | 0.0354 | 0.0401 | 0.0275 | 0.0090 | 0.0000 | 0.0088 | 0.0178 |
| p = 7 | 0.0024 | 0.0171 | 0.0287 | 0.0230 | 0.0071 | 0.0001 | 0.0107 | 0.0259 | 0.0290 |

B_2 reaches 0.0127 at w = 2.4 and is 0 at w = π/ln 2 = 4.53; B_3 ≤ 0.0056, B_5 ≤ 0.0020, B_7 ≤ 0.0010. The m = 1 peaks
are (2/π²)/p at w = π/ln p. At a fixed number of spacings the numbers move with height: for conductor 550, p = 2 at
l = 8 is 0.090 at t = 100 and 0.104 at t = 600; p = 3 at l = 4 is 0.067 at t = 100, 0.064 at t = 350.

Correlation coefficient Cov/Σ²: the BK number variance at L = 10.33 (no small bad prime) is 0.275, 0.359, 0.394, 0.414,
0.423, 0.426, 0.423, 0.412 at w = 0.3, 0.9, …, 4.5. Hence A_p/Σ² peaks at 0.25 (p = 2, w = 4.5), 0.16 (p = 3, w ≈ 2.8),
0.10 (p = 5), 0.07 (p = 7); for two even conductors (Σ² lower, p = 2 absent) p = 3, 5, 7 reach 0.20, 0.11, 0.10. The
measured mean variance of δn over the 546 used sets is 0.269-0.371 (all conductors mixed). Over the 148,785 used pairs
the full law has mean +0.002 to +0.009 and sd 0.008 (w = 0.3), 0.056 (1.5), 0.073 (2.7), 0.079 (3.9), 0.077 (4.8); p ≤ 7
carries sd 0.051-0.076 for w ≥ 1.5, the primes above 7 0.015-0.029 (at w = 0.3 the larger part, 0.0072 against 0.0035).

### 4.5 The swapped term in equal windows (computed, `swaplib.py`, (T5), `assemble_swap_out.txt`)

C_swap(w) = mean_t (1/2π²) Re[ ∫∫_[0,w]² (k₁/φ)(k₂/φ) e^(2i(θ₂(t+y₂) − θ₁(t+y₁))) H(y₂ − y₁) dy ], H(u) =
L(1−iu,ψ)L(1+iu,ψ)A₁₂(−iu) (real ψ), computed exactly in (u, v) = (y₂ − y₁, y₁ + y₂) with θ_j linearised in the window
(curvature ≤ 0.03 rad dropped) and the window integral tabulated in ln(√(k₁k₂)t/2π). With H ≡ 1 it reduces to
(1/2π²)(2 sin(Lw/2)/L)² cos(phase) ≤ 0.0019. Real values: up to 0.017 for (−1112, 1117), 0.008-0.009 for (1001, 1004),
0.003-0.006 for (401, −403), 0 for equal |D| and for far conductors (401, 929: < 10⁻⁵). Pairs with averaging factor
F = |sinc(250 ln(k₂/k₁))| > 0.1: 3,080 of the 17,205 fresh pairs (rms 0.002-0.005 per width, max 0.046), 134 of the
101,556 fresh-old pairs, 4,291 of the used pairs. The Euler product for H is truncated (primes ≤ 5·10⁴ frozen): against
2·10⁵ the predictions differ by 15 % rms with slope 0.956, and 2·10⁵ is itself not converged (1·10⁴ → 2·10⁵: slope 0.92)
(`swapcheck_out.txt`). So the swap predictions carry a scale uncertainty of order 5-10 % and a per-pair scatter of
order 15 %.

## 5. Checks (computed)

- (T3) Closed form with PARI L-values against the prime sums, 5 pairs × 3 widths: primes ≤ 2·10⁵ within 9.8·10⁻⁵,
  ≤ 10⁶ within 8.3·10⁻⁵ (the conditionally convergent tail). Frozen pair predictions use primes ≤ 2·10⁵.
- (T4) The recipe's two-point function (1/2π²) Re D_ψ(−iε), with D_ψ from PARI's L, L′, L″ at 1 − iε, the removed Euler
  factors and B_ψ, integrated over the window square: equal to the closed form within 2.8·10⁻⁷ at w = 1.2 and w = 4.5
  for (401, −403) and (−995, 997) (`t4fix_out.txt`; the 20-node rule of `theory.py` was not converged at 4.5).
- (T5) H(−ε) = conj H(ε) for real ψ; H(0) = 0 when ψ(3) = −1; C_swap = 0 for equal |D|.
- `check.py`: the per-prime tables recomputed with mpmath (max difference 1.5·10⁻¹⁵); 20 random pairs of the frozen file
  recomputed with sympy Kronecker symbols, plain loops over primes ≤ 2·10⁵ and the finite-range part with closed
  geometric sums: largest difference 6.3·10⁻⁹ (the float32 storage) (`check_out.txt`); sha256 of the pairs file checked.
- `predict.py` dry run: the registered readings computed on the law itself plus Gaussian noise of sd 0.003 per pair
  return every amplitude within 0.005 of its expected value. `rcross.py NN` on synthetic Poisson zeros (independent,
  no orbits; `dryrun_out.txt`): every amplitude within 1.2 se of 0 (se 0.05-0.15 for Poisson sets, 100 times the real
  ones: real zeros are rigid and orbit-dominated); the law is rejected, independence is not.

## 6. What a cross test measures

1. **The projection picture.** Write δn_i(t) = Σ_(p small) χ_i(p) g_p(t) + η_i(t), g_p the orbit cosine. The pooled
   per-prime reading is ≈ (Σ_i χ_i(p)⟨δn_i g_p⟩)²: each L-function's own line at frequency ln p, whose amplitude χ_i(p)
   is Landau's formula Σ_(0<γ≤T) x^(iγ) ≈ −(T/2π) χ(x) Λ(x)/√x for L(s, χ) (**classical**: Landau 1911, Gonek 1993
   uniform). So the per-prime amplitudes at fixed small p are close to a theorem; a failure would point first at the data
   (missing zeros, wrong smooth count). What they add over one-function statistics is the phase: rounds 117-136 could
   only see |χ(p)|².
2. **The conjectural content** is (i) that nothing else correlates two spectra: no off-diagonal (n ≠ m) cross term for
   k₁ ≠ k₂, measured by γ_full, γ_rest and small_w (at w = 0.3, 0.6 the orbits up to the Heisenberg range carry most of
   the spread, and the finite-range form moves small_w by 0.6 %), and by the pair residual; (ii) the swapped
   Gram-lattice term for close conductors (X6), the functional-equation part, which no one-function statistic contains
   (the self case's off-diagonal term is its |·|² analogue).
3. **Why the noise is small.** The count fluctuations are dominated by deterministic orbit terms, so a pair covariance has
   little genuinely random part; after the regression the residual is 0.003 per pair at every width, against
   0.005-0.028 for a stationary Gaussian model with the same variances (estimate √(V²·(2w/3)/500), not computed). The finite-range part (0.0007-0.0016) is about a
   third of it and is in the predictions. The per-function noise (each function's realised line amplitude) is what
   limits the pooled amplitude: the delete-one-L-function jackknife carries it.
4. Off-line single-function projections (frequencies ln y, y not a prime power) were computed as a second noise
   estimate. With the flat window they are dominated by sidelobe leakage of the p = 2 line (rms 1.1-3 % of a line,
   consistent with 2/(TΔx) at Δx ≈ 0.2) and do not measure the noise at a line; they were not used.

## 7. Noise and power (noise only, used sets; `noise_out.txt`, `noise2_out.txt`)

Completeness of the 546 used sets by round 138's rule: 0 flagged. Residual sd of the pair covariances after the
regression: 0.0039 (w = 0.3), 0.0029-0.0031 (w ≥ 0.6), identical in the subsets. Jackknife (delete one L-function with
all its pairs) standard errors of the registered readings; NN and NO from random 186-function subsets of the used sets
(two subsets each; NO with 360 old partners, the real NO has 546, so slightly smaller):

| reading | OO (148,785 pairs) | NN-like (17,205) | NO-like (66,960) | z vs H_indep on NN | z vs H_flip on NN |
|---|---|---|---|---|---|
| γ₂ | 0.0004 | 0.0007-0.0008 | 0.0004 | ~1,300 | ~2,600 |
| γ₃ | 0.0006 | 0.0011 | 0.0007 | ~900 | ~1,800 |
| γ₅ | 0.0010 | 0.0018-0.0020 | 0.0011 | ~530 | ~1,050 |
| γ₇ | 0.0010 | 0.0020-0.0021 | 0.0011-0.0012 | ~490 | ~980 |
| γ_B2 | 0.0039 | 0.008-0.009 | 0.0045 | ~120 | |
| γ_B3 | 0.0064 | 0.015-0.016 | 0.008 | ~65 | |
| γ_rest (p > 7) | 0.0008 | 0.0017 | 0.0009-0.0010 | ~590 | |
| γ_full | 0.0003 | 0.0006-0.0007 | 0.0003-0.0004 | ~1,600 | ~3,200 |
| small_w (w = 0.3, 0.6) | 0.0008 | 0.0017-0.0018 | 0.0010 | ~560 | |
| γ_sw (Gram lattice) | 0.0146 (crude predictor, Σx² 1.25) | 0.036-0.056 (subsets, Σx² 0.13-0.17) | 0.02 | ~40 (real NN: Σx² 0.56, se ≈ 0.022-0.027 scaled) | |

H_abs (orbits shared without phases) gives γ_p = 0 and γ_B2 = 6.4: rejected at several hundred se. The per-width
coefficients s_p(w) have se 2-9·10⁻⁵ (NN-like), against A_p(w) up to 0.104, so the sin²(w ln p/2) profile of each prime
(its orbit period) is measured at every width (X3).

Answer to "can the data see it": yes, overwhelmingly, in every design, including the 17,205 fresh pairs alone. The
power against independence is not the issue; the readings test the law's size to 0.1 % (γ₂, γ_full), 0.2 % (γ₅, γ₇),
1 % (γ_B2) and about 2.5 % plus the prediction's own 5-10 % for the Gram-lattice term.

## 8. The test design (summary of `round_draft.txt`)

Primary NN (Fn × Fn, fresh), secondary NO (Fn × used sets), optional OO. Readings X1 per-prime amplitudes, X2 full law,
X3 shape in w for p = 2, 3, X4 even-m terms, X5 rest and small w, X6 Gram-lattice term, X7 intercepts. Frozen expected
values (finite-range law) are within 6·10⁻⁴ of 1 for γ_p, γ_full, small_w; 1.004/1.006 for γ_B2/γ_B3 in NN. Pass bands
4 se from the data's jackknife; a computed theory allowance (D4) is stated separately.

## 9. Prior art (searched with mathematical queries; what was opened is marked)

- Bogomolny–Leboeuf, "Statistical properties of the zeros of zeta functions – beyond the Riemann case", Nonlinearity 7
  (1994), chao-dyn/9409004 (**opened**, section 3): for products of L-functions they write
  R₂ = d̄² + Σ_i R₂^(i) + Σ_(i≠j)⟨d_osc^(i)(t) d_osc^(j)(t+ε)⟩ and set the cross term to 0 because "the mean value of any
  non-principal character is zero" (their eq. 3.11-3.12): the leading-order superposition. Their numerical check (the
  three non-principal characters mod 5 near t = 10⁵) compares with that leading order only. Our D_ψ term is exactly what
  this step discards (it is non-universal, of the same order as the per-prime terms of rounds 117-136).
- Rudnick–Sarnak, Duke 1996, and Liu–Ye, "Superposition of zeros of distinct L-functions", Forum Math. 14 (2002)
  (abstracts only): distinct primitive L-functions are uncorrelated for restricted test functions (rigorous, leading
  order). Murty–Perelli, IMRN 1999 (not opened; the PDF host redirects): the pair-correlation conjecture for two
  Selberg-class functions, main term only for F = G. Murty–Zaharescu, "Explicit formulas for the pair correlation of zeros
  of functions in the Selberg class", Forum Math. 14 (2002) (not opened; no accessible copy found): its explicit formula
  in the Montgomery range presumably contains Σ Λ_F(n)Λ̄_G(n)/n-type terms; whether a lower-order cross term is stated I
  could not check.
- Conrey–Simons–Smolyarenko, "Pair correlation of different L-functions" (unpublished; listed on Conrey's AIM publication
  page; the PDF returned 404). From the indexed first page: it starts from a ratios conjecture for a ratio of two
  L-functions and derives the correlation between their zeros, assuming the same level "for simplicity". This is
  presumably §3 of this note (D_ψ and the same-level swapped term). Not opened, so not confirmed.
- Kandhil–Languasco–Moree, Math. Ann. 394 (2026), arXiv:2411.19762 (**HTML opened, first 100k characters**): define
  G_(χ₁,χ₂)(x, T) = Σ x^(i(γ₁−γ₂)) W(γ₁ − γ₂) for two characters mod q and ask how far the zeros of L(s, χ₁) and L(s, χ₂)
  are correlated; in the part read no lower-order cross term and no numerics on zeros.
- de Laat–Rolen–Tripp–Wagner, arXiv:1908.04876 (**opened**): pair correlation of Dedekind zeta functions of abelian fields;
  cross terms between characters are not isolated.
- Arias de Reyna, arXiv:1403.7126 (**opened**): Σ_(0<γ≤T) e^(2iθ(γ)) = −ψ(T/2π) + O(T^(1/2) log² T): the one-point Gram
  term for ζ (§3.5).
- Landau's formula (classical) gives each L-function's own orbit line; the cross covariance per prime is its quadratic
  form.

Verdicts. **Formula** (cross two-point function D_ψ and swapped term): KNOWN IN SUBSTANCE (the ratios recipe gives it
mechanically; very likely written out in the unpublished CSS preprint; Bogomolny–Leboeuf name the term and set it to 0).
**Window covariance, its closed form (1/π²) ln|L(1,ψ)/L(1+iw,ψ)|, the per-prime sign rule, the equal-|D| oddness, the
Gram-lattice term for close conductors**: NOT FOUND. **One-point Gram term**: KNOWN for ζ, NOT FOUND for L(s, χ).
**A numerical test of the cross term on zeros**: NOT FOUND.

## 10. Doubts

1. The swapped term's constants (1/4π², + c.c.) are taken over from the self case; only its H → 1 limit is checked
   against a theorem (Arias de Reyna, ζ). Its sign convention in the drifting phase (e^(+iΘ₀) with Θ₀ = 2θ₂ − 2θ₁ paired with
   H(t₂ − t₁)) follows the identification L₁ at t₁, L₂ at t₂ of §3.4; I did not test the opposite pairing on data (it
   would be a different, also parameter-free predictor).
2. The swap predictions are limited by the Euler product of H (5-10 % scale, 15 % per-pair scatter). Before freezing,
   recompute with SWAP_PMAX = 200000 (≈ 30 min) or with PARI L-values on a coarse u-grid; the registered allowance for X6 is
   0.2.
3. Disclosure: the noise computation on the used sets already shows that the used pair covariances are linear in the
   predictors (residual variance 25 % of the predicted variance at w = 0.3, below 1.5 % at w ≥ 0.6). It does not fix any
   coefficient, but test OO is no longer blind in that sense.
4. The per-prime amplitudes at small p are close to the explicit formula (Landau): passing them shows the orbit phases and
   the data quality, not a deep conjecture; the conjecture is in X2, X5 (small w) and X6.
5. At 0.1 % precision, unknown lower-order terms of the recipe could show. Computed: finite-range ≤ 6·10⁻⁴ (in the
   expected values), prime cutoff ≤ 10⁻⁴ per pair, one-swap terms ≈ 2·10⁻⁴. Hence the D4 allowance 0.003 (about 4 se on
   NN for γ₂).
6. The NN power comes from 186-function subsets of the used sets (L 9-11); Fn has L 10.1-11.9 and denser conductors
   (more near pairs). The residual sd was flat in w and identical in all subsets.
7. The delete-one-L-function jackknife for dyadic data is standard but not proved here to be unbiased; a bootstrap over
   L-functions would duplicate functions and create i = i pairs, so it was not used.
8. Complex characters (the phase law cos(m(φ₁ − φ₂)) with angles other than 0, π, and the root-number term of equal
   modulus) have no fresh data; r121 (88 complex characters, both half-lines, 100-500) could test them as "zeros known,
   statistic new" after predictions are frozen (not done).
9. GRH is assumed throughout; the exact smooth count uses θ_χ without the O(e^(−πt)) terms.
10. Completeness: a missing zero biases every pair of its set; the round-138 rule is required. All 546 used sets pass it.

## 11. What is whose

- **Classical / in print**: the explicit formula and Landau's formula; the ratios recipe (CFZ 2008, Conrey–Snaith 2007);
  BK's one-function law; Bogomolny–Leboeuf's superposition (with the cross term set to 0); leading-order independence
  (Rudnick–Sarnak, Liu–Ye; Murty–Perelli conjecture); the ζ one-point Gram term (Arias de Reyna); very likely the ratios
  cross pair correlation (Conrey–Simons–Smolyarenko, unpublished, unopened).
- **Derived here**: the recipe for two characters written out (G_p with ψ, the swapped term with (k₁/φ(k₁))(k₂/φ(k₂))
  and A₁₂, the case analysis including the exact oddness for |D₁| = |D₂|); the window covariance, its closed form with
  L(1 ± iw, ψ), the per-prime decomposition and sign rule; the finite-range form; the Gram-lattice term in equal windows
  for close conductors; the χ-version of the one-point Gram term.
- **Computed**: every number above; the noise and power on used sets (no coefficient looked at).
