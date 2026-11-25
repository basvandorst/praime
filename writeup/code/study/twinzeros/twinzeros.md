# Which zeros do twin primes hear? A linearised explicit formula for prime patterns

> **In short.** Proposed: twin primes hear the zeros of ζ twice as loudly as the primes do, and also hear the
> zeros of Dirichlet L-functions, including the character mod 3 at about the loudness of ζ.
> Derivation, done two ways (local densities; major arcs of the circle method):
> (1) the count of a pattern is S times the product of the local densities of its members;
> (2) linearised, each member adds its prime error ψ(a x, χ) times the mean of χ̄ over the admissible classes;
> (3) for pairs (n, n + h) the two members' terms cancel for odd χ.
> Result: Σ_{n≤x} Λ(n)Λ(n+h) − S(h)x ≈ 2S(h)[(ψ(x) − x) + Σ_χ μ(q)χ̄(h)/∏_{p|q}(p − 2) · ψ(x, χ)], where the sum runs over the
> even primitive χ of odd squarefree conductor q coprime to h.
> So the factor 2 holds, but the character mod 3 (odd) and the character mod 4 are absent for twins. The character mod 5 enters at 1/3
> of the loudness of ζ, and the characters mod 7 at 1/5.
> Discovery range (10⁶ to 2·10⁹, 10 patterns): every ζ coefficient and every L-function loudness agrees within 1.7σ,
> and the mod-3 coefficient for twins is 0.010 ± 0.043 (the guess was ±2.64).
> A frozen test on (2·10⁹, 2·10¹⁰] is ready.
> Novelty: the averaged ζ part is known (Korevaar–te Riele 2010). The linear-in-zeros structure is known for
> exceptional zeros near σ = 1 (Pintz 2018). The individual coefficient and the character selection rule were not found.

Status: **open** (heuristic conjecture; discovery-range agreement; holdout not yet run).

Everything below is in this folder (called "the folder"); the count files disc.npz (85 MB) and hold.npz are scratch files, not stored, with their sha256 in the registration.

---------------------------------------------------------------------------------------------------------------------

## 1. The conjecture

**Pattern.** Linear forms L_i(n) = a_i n + b_i (i = 1..k), Hardy–Littlewood constant S = ∏_p (1 − ν_p/p)(1 − 1/p)^(−k)
(ν_p = number of n mod p with p | ∏ L_i(n)). Count R(x) = Σ_{n≤x} ∏_i Λ(L_i(n)), error E(x) = R(x) − S x.

**Statement (heuristic; linear terms in the zeros only).**

    E(x) ≈ S · Σ_i (1/a_i) · [ (ψ(a_i x) − a_i x)  +  Σ_{χ primitive, χ ≠ 1} c_i(χ) · ψ(a_i x, χ) ]  + (rest)

    c_i(χ) = mean of conj χ(a_i n + b_i) over the residues n mod q (q = conductor of χ) with gcd(∏_j L_j(n), q) = 1.

Since ψ(y, χ) = −Σ_ρ y^ρ/ρ (+ smooth terms), each zero ρ of L(s, χ) enters E(x) as −A(ρ) x^ρ/ρ with

    A(ρ) = S · Σ_i c_i(χ) · a_i^(ρ − 1)            (for ζ: c_i = 1, so A = S · Σ_i a_i^(ρ − 1)).

Relative to the primes themselves (ψ(x) − x = −Σ x^ρ/ρ: amplitude 1), the pattern hears ζ with amplitude S·Σ a_i^(ρ−1)
and L(s, χ) with amplitude S·Σ c_i(χ) a_i^(ρ−1). For translates (all a_i = 1) only the sum Σ_i c_i(χ) matters,
because ψ(x + b_i, χ) − ψ(x, χ) = O(b_i ln x).

**Pairs (n, n + h), h even (closed form, proved from the definition in §2.3).** With S(h) = 2C₂ ∏_{p|h, p>2} (p−1)/(p−2):

    E_h(x) ≈ 2S(h) · (ψ(x) − x)  +  Σ_χ β_h(χ) ψ(x, χ),
    β_h(χ) = 2 S(h) μ(q) χ̄(h) / ∏_{p|q} (p − 2)   if χ is EVEN, q odd squarefree, gcd(q, h) = 1;   β_h(χ) = 0 otherwise.

The mechanism, mod p: the admissible n are the p − 2 classes with n ≢ 0, −h. The member n contributes
Σ_{n adm} χ̄(n) = −χ̄(−h), and the member n + h contributes −χ̄(h). The total is −χ̄(h)(1 + χ(−1)), divided by p − 2.
For odd χ the two members cancel exactly. For example, a twin pair p > 3 has p ≡ 2 and p + 2 ≡ 1 mod 3. The excess of primes in
class 2 mod 3 is exactly the deficit in class 1, so the two cancel. The character mod 4 cancels in the same way: n and n + 2
always lie in opposite classes mod 4. Higher prime powers in q cancel too, because the average over a full lift class is 0.

**Predicted coefficients, conductor ≤ 30** (full table for 8 patterns: `amp_table.md`; machine table: `theory_table.txt`).
"Loudness" means |β|/(coefficient of ζ):

| pattern | ζ (coefficient on ψ(ax) − ax) | characters present (loudness) | absent (examples) |
|---|---|---|---|
| twins h = 2 | 2S = 2.6406 | 5.1 real even, +2S/3 (1/3); 7: cubic pair (1/5); 11: four quintic (1/9); 13: real + cubic + sextic (1/11); 15: quartic pair (1/3); 17 (1/15); 19 (1/17); 21: real (−2S/5) + sextic pair (1/5); 23 (1/21); 29 (1/27) | χ₋₃ (mod 3), χ₋₄, the odd characters mod 5, 7, …, all of conductor 8, 9, 12, 16, 20, 24, 25, 27, 28 |
| cousins h = 4 | 2S = 2.6406 | as h = 2 with χ̄(4) in place of χ̄(2); 5.1 has coefficient −2S/3 (**opposite sign to twins**) | same |
| h = 6, 12 | 2S(6) = 5.2813 | conductors coprime to 3 only: 5, 7, 11, 13, 17, 19, 23, 29 | in addition 15, 21 |
| h = 10, 30 | 2S(h) | conductors coprime to 5 (and 3) | 5.1 absent for h = 10 |
| Sophie Germain (n, 2n+1) | S at scale x, S/2 at scale 2x (so A(ρ) = S(1 + 2^(ρ−1))) | **χ₋₃ at full loudness on both members** (β = −S, −S/2); **χ₋₄ only through 2n+1** (β₂ = −S/2, full loudness); χ₁₂ through 2n+1 (+S/2); odd characters are present too (no cancellation between different scales); conductors 5, 7, 11, 13, 15, 20, 21, … | conductors 8, 9, 16, 24, 25, 27 |
| triple (0, 2, 6) | 3S₃ = 8.5747 (S₃ = 2.858249) | χ₋₃: β = −S₃ (1/3); 5.1: −S₃; 15 (quartic pair 0.47, real 1/3); 21; 11; 13 | 4.1, 7.x, 8.x, … |
| triple (0, 4, 6) | 3S₃ | χ₋₃: β = **+S₃** (opposite sign to (0, 2, 6)); others as (0, 2, 6) | same |

Elementary reasons: Sophie Germain pairs (p, 2p + 1) with p > 3 have both members ≡ 2 mod 3, and 2p + 1 ≡ 3 mod 4.
The triple (p, p + 2, p + 6) lies in classes (2, 1, 2) mod 3, and (p, p + 4, p + 6) in (1, 2, 1).

**Amplitude at each zero.** In the spectrum F(γ) = ∫ E(e^u) e^(−u/2) W(u) e^(−iγu) du, a zero ρ_j = 1/2 + iγ_j of L(s, χ)
gives the line −A(ρ_j)/ρ_j · Ŵ(γ − γ_j), with Ŵ the transform of the taper (a Dirichlet kernel for a rectangular window).
For twins, at the first zero of ζ (14.1347) |A/ρ| = 2S/|ρ| = 0.187. At the first zero of L(s, χ₅) (6.6485) it is (2S/3)/|ρ| = 0.132.
At 8.0397 (χ₋₃) and 6.0209 (χ₋₄) it is 0.
Each single line is weaker than the noise (SNR ≈ 0.4–0.8 per line in the discovery window). Only the whole comb is measured, in the
least-squares fit of §3.

**Against the original guess.** The guess had A(q, χ, h) = 2μ(q)χ̄(−h)χ(−1)q/φ(q)² × (correction): one cross term, doubled.
There are two cross terms, main(q) × conj(fluct_χ) and fluct_χ × conj(main(q)). They carry χ̄(−h) and χ̄(h) respectively, so
their sum is χ̄(h)(1 + χ(−1)). The factor 2 for ζ survives (χ₀ is even). The mod-3 line does not.

## 2. Derivation

### 2.1 Local densities (heuristic)
The density of primes near y in the class c mod q is (1/φ(q))·[ψ′(y) + Σ_{χ≠χ₀} χ̄(c) ψ′(y, χ)]. This is an exact identity
for the counts, ψ(y; q, c) = (1/φ(q)) Σ_χ χ̄(c) ψ(y, χ), up to the primes dividing q.

The Hardy–Littlewood conjecture restricted to n ≡ b (mod q), b admissible, gives every admissible class the same share of S x
(the local factor at p | q does not depend on b).

Assume the pattern density in class b is S times the product of the class densities of its members, and expand to first order.
The relative fluctuation is then Σ_i [(ψ′(a_i n) − 1) + Σ_χ χ̄(L_i(b)) ψ′(a_i n, χ)]. Average it over the admissible b
and integrate dn, using ∫_0^x ψ′(a t, χ) dt = ψ(a x, χ)/a. This gives the statement of §1.

### 2.2 Circle method (major arcs; heuristic)
S(α) = Σ_{n≤N} Λ(n)e(nα). For (a, q) = 1 and (n, q) = 1, e(an/q) = Σ_{χ mod q} c_χ(a) χ(n) with
c_χ(a) = φ(q)^(−1) Σ_b χ̄(b) e(ab/q) = φ(q)^(−1) χ(a) μ(q/q*) χ̄*(q/q*) τ(χ̄*), where χ is induced by χ* mod q*.

So S(a/q + β) ≈ Σ_χ c_χ(a) Ψ_χ*(β) with Ψ_χ(β) = Σ χ(n)Λ(n)e(nβ) = δ_χ T(β) − Σ_ρ V_ρ(β).
The principal character carries μ(q)/φ(q) and therefore both the main term and the ζ-zeros.

R(N; h) = ∫|S(α)|² e(−hα) dα. On the arc at a/q the terms linear in the zeros are
main × conj(fluct) + fluct × conj(main). With ∫T conj(V_ρ) e(−hβ) dβ ≈ N^ρ̄/ρ̄ and ∫V_ρ conj(T) e(−hβ) dβ ≈ N^ρ/ρ, the
sum over a of χ̄*(a)e(−ha/q) is a Gauss sum of the induced character. For q = q*m (m squarefree, coprime to q*, gcd(h, q*) = 1)
it equals χ*(−h) χ̄*(m) τ(χ̄*) c_m(h).

Collecting terms, using τ(χ)τ(χ̄) = χ(−1)q and conj τ(χ̄) = χ(−1)τ(χ), the coefficient of N^ρ/ρ is

    −(μ(q*) q*/φ(q*)²) · Σ_{m} c_m(h)/φ(m)² · χ̄*(h)(1 + χ*(−1))
    = −S(h) μ(q*) χ̄*(h)(1 + χ*(−1)) / ∏_{p|q*}(p − 2),

because Σ_{m sqfree, (m,q*)=1} c_m(h)/φ(m)² = S(h)/∏_{p|q*}(1 − 1/(p−1)²) and q*/φ(q*)² · ∏ (p−1)²/(p(p−2)) = 1/∏(p − 2).

Linear terms need the main term at the same q, so μ(q) ≠ 0. Conductors 4 and 8 never enter, and since there is no primitive
character mod 2, q* is odd. If p | gcd(q*, h), the Gauss sum vanishes.

The result is the same as §2.1. Pintz (2018, Theorem 1, the "generalized singular series" S′(χ₀, χ, m)) has the same structure,
rigorously, for exceptional zeros near σ = 1.

### 2.3 Checks of the bookkeeping
[`check_formula.py`](check_formula.py) → `check_formula.out`. It runs three independent computations:
(1) enumeration of admissible classes (theory.local_coeffs);
(2) the closed form;
(3) the circle-method sum over q ≤ 2500 with explicit Gauss sums of the induced characters.
Results:
* (1) and (2) agree for all 167 primitive characters of conductor ≤ 30 and gaps 2, 4, 6, 8, 10, 12, 30 (max difference 2·10⁻¹²).
* (3) agrees with (1) to the truncation error (≤ 0.4 % at q ≤ 900). Examples for h = 2: χ₅: −0.8796 against −0.8802;
  the cubic character 7.2: −0.2638 − 0.4569i against −0.2641 − 0.4574i; the quartic 15.2: 0.8839i against 0.8802i;
  21.1: +0.5312 against +0.5281.
* The circle method gives exactly 0 for every odd character tested (3.1, 4.1, 5.2, 7.1, 7.4, 15.1) and for conductors 8 and 12.
  So the two cross terms of the major arcs cancel for odd χ, independently of the local-density argument.

**Model where everything is a theorem.** Take the pattern {q m + b} (one form, a = q). The rule gives S = q/φ(q), and
c(χ) = χ̄(b) if cond(χ) | q, else 0. The statement becomes E(x) = ψ(qx; q, b) − qx/φ(q) ≈ (1/φ(q))[ψ(qx) − qx + Σ_{cond χ | q, χ ≠ 1} χ̄(b) ψ(qx, χ)],
which is the classical explicit formula for primes in progressions (exact, up to the primes dividing q). The same fitting
pipeline was then checked against the theorem with explicit zeros (`disc1.py` → `disc1.out`; band [2, 1000], zeros up to 1500,
window 10⁶ to 2·10⁹):

| series | regressor(s) | fitted | theorem |
|---|---|---|---|
| (θ(x) − x)/√x | comb of ζ zeros | 0.9993 ± 0.0011 | 1 |
| θ(x, χ₋₃)/√x | comb of L(s, χ₋₃) zeros | 0.9999 ± 0.0009 | 1 |
| θ(x, χ₋₄)/√x | comb of L(s, χ₋₄) zeros | 1.0002 ± 0.0008 | 1 |
| θ(x, χ₅)/√x | comb of L(s, χ₅) zeros | 1.0003 ± 0.0009 | 1 |
| (θ(x; 3, 2) − x/2)/√x | ζ comb, χ₋₃ comb | 0.4997 ± 0.0005, −0.5000 ± 0.0004 | ½, −½ |

Here comb(u) = −2 Re Σ_{0<γ≤T} e^(iγu)/ρ. The conventions (signs, 1/ρ, phases, the √x normalisation) are right.

### 2.4 What the linear formula leaves out
* **Terms quadratic in the zeros.** From |S|² on the arcs: ∫V_ρ conj(V_ρ′) e(−hβ) dβ ≈ x^(ρ+ρ̄′−1)/(ρ+ρ̄′−1). On GRH this is
  x^(i(γ−γ′)): size x⁰, frequency γ − γ′. In Y = E/√x these terms decay like e^(−u/2). They produce no fixed-frequency line of
  order x^(1/2), at γ + γ′ or anywhere else.
  The diagonal γ = γ′ gives non-oscillating ln x-type terms; with the truncation of the major arcs they belong to the
  lower-order mean (compare the −(1/2)·ln-type corrections of Montgomery–Soundararajan).
  For patterns with three or more members, products of two unconjugated factors give x^(i(γ+γ′)), again of size x⁰.
* **Minor arcs** are not modelled. Empirically they are the dominant "noise": a random walk in u whose increments have about
  0.5–0.6 of the Poisson variance S·(ln x)² per unit u. This is measured, not derived; it matches the four-point
  Hardy–Littlewood rigidity of twin counts (STUDY.md round 31).
* **Prime powers.** The statistic uses primes only (weights ln p·ln(p+h), regressors θ(x, χ)). The difference from Λ-weights is
  a smooth c·x^(1/2) term plus O(x^(1/4)) fluctuations: frequency ≈ 0, absorbed by the nuisance regressors.
* **Secondary main terms** (S(h)(x − h), x⁰ constants): smooth, absorbed the same way.
* **Characters of conductor > 30** are not in the fit. Their lines stay in the noise, which adds noise but no bias, since they are
  nearly orthogonal to the fitted combs. Heuristically, a prime conductor p carries about p/2 even characters, each with loudness 1/(p − 2) and a zero density growing
  like ln p. Their summed power behaves like Σ_{p≤Q} ln p/(2p), which grows like ½ ln Q: each line is tiny, but all of them
  together are not negligible. This was not estimated separately; it inflates the noise, not the coefficients.

### 2.5 Interpretation (honest)
The factor 2 is the linear response of a product of two densities to a fluctuation of one of them. It appears in any model in which
pairs are products of prime events.

[`cramer_ctrl.py`](cramer_ctrl.py) → `cramer_ctrl.out` shows this on a Cramér model with no zeros (odd n "prime" with probability
2/ln n, S_model = 2). The same regression gives **4.013 ± 0.052** = 2S_model (band [2, 6000]; 3.58 ± 0.25 in band [2, 300]).

So the content of the conjecture is not "twins know about zeros". It is that the pair counts respond to the real fluctuations of
the primes like a product of local densities. Those fluctuations are rigid (zero-driven, far below Poisson at long scales).

The selection rule for the L-function lines is the same local statement made class by class. Under Hardy–Littlewood for triples,
averaged over the third member, it is the covariance identity Cov(pairs, primes) = 2S·Var(primes) at every scale.
What is new here is the explicit form and the measurement, not a new mechanism.

## 3. Statistic and fit

* **Counts** ([`tz_counts.py`](tz_counts.py), segmented numpy sieve, [`sieve_core.py`](sieve_core.py)). The grid has edges
  X_k = ⌊exp(k/2048)⌋ (integers, computed with mpmath).
  Per cell it stores: Σ ln p by residue class mod 3, 4, 5, 7, 8, 9, 11, 12, 13, 15, 17, 19, 20, 21, 23, 24, 28, 29; the same on the
  doubled grid 2X_k (for the 2n + 1 member); pair weights ln p·ln(p+h) for h = 2, 4, 6, 8, 10, 12, 30; the triples (0,2,6) and
  (0,4,6); and Sophie Germain weights ln p·ln(2p+1).
  Verified against sympy on small ranges (`check_counts.py`, all statistics at 4 grid points, two runs, both grids), and a paused
  and resumed run gives the identical file hash. Discovery file `disc.npz`: π(4·10⁹) = 189,961,812 (correct), 6,388,041 twin pairs
  and 6,178,398 Sophie Germain pairs with p ≤ 2·10⁹.
* **Series** (`tzlib.py`, `models.py`): Y(u) = E(eᵘ)/√x, regressors (θ(ax) − ax)/√x and θ(ax, χ)/√x.
* **Fit** (approach A, the comb from the primes):
  1. Difference the series in u; a random walk becomes white noise.
  2. Weight by (ln x)^(−k/2), the noise level of a k-term pattern.
  3. FFT with a rectangular window and keep the bins with γ_lo ≤ γ ≤ γ_hi.
  4. Fit weighted least squares on the real and imaginary parts, with nuisance regressors (u, u², e^(−u/2), e^(−u), constant drift).
  5. Estimate the noise variance per bin from the smoothed squared residuals (41 bins, 3 passes); errors from (XᵀWX)⁻¹.

  The prime-based regressor θ(x, χ) equals the full zero comb of L(s, χ) by the explicit formula (a theorem), so this is a
  least-squares fit of the whole predicted comb, all zeros of all included L-functions, with no zero list needed.
  Model "all": ζ free, one loudness factor for all predicted characters together, and free coefficients for a null set of absent
  characters (3.1, 4.1, 5.2, 7.1, 7.4, 8.1, 8.2, 11.1, 12.1, 9.1).
  Model "cond": one loudness per conductor.
* **Approach B** (explicit zeros): the same fit on combs built from the zero lists of ζ, L(s, χ₅), L(s, χ₋₃) and L(s, χ₋₄)
  (`/home/bas/praime/writeup/data/riemann/zeros_40000_flint.txt`, `lzeros/D5.txt`, `D-3.txt`, `D-4.txt`).
* **Error checks.**
  * Shifted regressors: twins against the prime series shifted by 1000–8000 cells, 21 shifts. The coefficients scatter with
    sd 0.151, 0.135, 0.136, 0.100 against mean reported errors 0.161, 0.143, 0.139, 0.137 (z-score sd 0.94, 0.94, 0.98, 0.73),
    so the errors are right or slightly conservative.
  * Scrambled control: Poisson pseudo-twins with the same mean (5 seeds). ζ coefficients 0.034, 0.043, −0.055, −0.044, −0.040
    (± 0.065). Loudness 0.020, −0.037, −0.054, −0.065, 0.007 (± 0.024).

## 4. Discovery-range results (x from 10⁶ to 2·10⁹; nothing above 2·10⁹ was counted for any pair, triple or SG statistic)

**Twins hear ζ with coefficient 2.675 ± 0.046 (predicted 2S(2) = 2.641).** In plain words: twin primes hear the zeros of ζ at
2.03 ± 0.04 times their own density factor S, so the factor 2 is confirmed (factor 1 is excluded at 29σ).

**They hear L(s, χ₅) at 0.833 ± 0.046 (predicted 2S/3 = 0.880).** Taking all predicted characters of conductor ≤ 30 together,
the loudness is 1.005 ± 0.017 (predicted 1; absent would be 0, excluded at 59σ).

**They do not hear χ₋₃ (0.010 ± 0.043) or χ₋₄ (0.014 ± 0.042).** The 13 null coefficients give χ² = 8.6 on 13 dof.

Band [2, 6000], model "all", `disc3.out`:

| pattern | ζ coefficient (pred.) | z | all-χ loudness (pred. 1) | z | null set χ² / dof |
|---|---|---|---|---|---|
| h = 2 | 2.675 ± 0.046 (2.641) | +0.74 | 1.005 ± 0.017 | +0.31 | 8.6 / 13 |
| h = 4 | 2.689 ± 0.047 (2.641) | +1.01 | 1.008 ± 0.017 | +0.48 | 8.9 / 13 |
| h = 6 | 5.252 ± 0.064 (5.281) | −0.46 | 0.994 ± 0.017 | −0.33 | 6.6 / 13 |
| h = 8 | 2.656 ± 0.047 (2.641) | +0.33 | 1.011 ± 0.018 | +0.62 | 8.0 / 13 |
| h = 10 | 3.449 ± 0.053 (3.521) | −1.34 | 0.968 ± 0.020 | −1.60 | 17.4 / 13 |
| h = 12 | 5.337 ± 0.064 (5.281) | +0.87 | 0.995 ± 0.017 | −0.27 | 9.7 / 13 |
| h = 30 | 6.961 ± 0.071 (7.042) | −1.13 | 0.993 ± 0.016 | −0.43 | 11.4 / 13 |
| SG, scale x | 1.297 ± 0.050 (1.320) | −0.47 | 0.995 ± 0.015 | −0.32 | 1.8 / 4 |
| SG, scale 2x | 0.692 ± 0.035 (0.660) | +0.90 | | | |
| (0, 2, 6) | 8.888 ± 0.364 (8.575) | +0.86 | 1.003 ± 0.032 | +0.09 | 16.0 / 12 |
| (0, 4, 6) | 9.000 ± 0.364 (8.575) | +1.17 | 1.033 ± 0.032 | +1.03 | 13.4 / 12 |

The 21 z-values have χ² = 13.7.

**Per conductor** (model "cond", `disc3_percond.out`; loudness, predicted 1):

| pattern | conductor: loudness |
|---|---|
| twins | 5: 0.943 ± 0.047 · 7: 0.990 ± 0.055 · 11: 0.996 ± 0.067 · 13: 1.011 ± 0.073 · 15: 1.062 ± 0.031 · 17: 1.069 ± 0.082 · 19: 0.910 ± 0.087 · 21: 0.967 ± 0.042 · 23: 1.036 ± 0.095 · 29: 0.973 ± 0.106 |
| SG | 3: 0.989 ± 0.027 · 4: 1.002 ± 0.048 · 5: 0.880 ± 0.047 (−2.6σ) · 7: 0.977 ± 0.059 · 12: 1.031 ± 0.045 · 15: 1.064 ± 0.045 · 20: 1.085 ± 0.075 · 21: 0.995 ± 0.057 |
| (0, 2, 6) | 3: 1.028 ± 0.117 (raw −2.94 ± 0.33, predicted −S₃ = −2.858) · 5: 0.999 ± 0.114 · 15: 0.975 ± 0.048 · 21: 1.003 ± 0.074 · 23: 1.49 ± 0.21 |
| (0, 4, 6) | 3: 1.035 ± 0.117 (raw +2.96 ± 0.33, predicted +S₃) · 5: 1.108 ± 0.113 · 15: 1.024 ± 0.047 · 21: 0.968 ± 0.073 · 29: 1.52 ± 0.23 |

In plain words:
* Twin primes hear the zeros of ζ at 2S(2), and the zeros of L(s, χ) for the even characters mod 5, 7, 11, 13, 15, 17, 19,
  21, 23, 29 at 2S/∏(p − 2). They hear neither the mod-3 nor the mod-4 character.
* Cousin primes hear χ₅ with the opposite sign (fitted inside the joint loudness 1.008 ± 0.017).
* Sophie Germain pairs hear ζ at both scales, x (1.297) and 2x (0.692). They hear the mod-3 character as loudly as ζ
  (0.989 ± 0.027), and the mod-4 character through 2p + 1 (1.002 ± 0.048).
* The prime triples hear the mod-3 character with opposite signs, (0, 2, 6) negative and (0, 4, 6) positive, as predicted.

**Dependence on band and window** (twins; `disc2.out`). The ζ coefficient in bands [2, 50] … [2, 6000]:
−0.16 ± 1.01, 1.79 ± 0.66, 2.24 ± 0.42, 2.95 ± 0.23, 2.87 ± 0.15, 2.72 ± 0.10, 2.70 ± 0.06, 2.695 ± 0.051.
In disjoint sub-bands: [200, 500] 3.24 ± 0.27; [500, 1000] 2.81 ± 0.20; [1000, 2000] 2.60 ± 0.13; [2000, 4000] 2.68 ± 0.09;
[4000, 6000] 2.69 ± 0.08.
For χ₅ in the same sub-bands: 0.40 ± 0.23, 0.87 ± 0.17, 0.92 ± 0.11, 0.81 ± 0.08, 0.86 ± 0.08.
By window (band [2, 1000]): 10⁴–10⁶: 2.57 ± 0.12; 10⁶–10⁸: 2.83 ± 0.18; 10⁸–2·10⁹: 3.02 ± 0.27.

**Low band [2, 300] only** (long scales, where primes are most rigid; `disc4.out`). ζ coefficient / loudness:
h2 2.48 ± 0.30 / 1.13 ± 0.10; h4 2.83 ± 0.29 / 1.00 ± 0.10; h6 5.94 ± 0.40 / 1.02 ± 0.09; h8 2.55 ± 0.29 / 1.00 ± 0.10;
h10 3.41 ± 0.33 / 0.99 ± 0.11; h12 5.45 ± 0.40 / 0.98 ± 0.10; h30 7.12 ± 0.44 / 0.95 ± 0.09;
SG 1.21 ± 0.31 and 0.63 ± 0.22 / 0.94 ± 0.08; (0,2,6) 9.7 ± 2.3 / 0.83 ± 0.19; (0,4,6) 11.9 ± 2.3 / 1.30 ± 0.19.
All within 1.7σ.

**Transfer function primes → twins** (cross-spectrum in bands after removing the predicted χ part; `disc4.out`).
H = 1.86 ± 1.04, 2.91 ± 0.86, 2.65 ± 0.51, 2.98 ± 0.32, 3.05 ± 0.21, 2.53 ± 0.15, 2.75 ± 0.09, 2.63 ± 0.06 in the bands from
[30, 60) up to [3000, 6000). It is flat at 2S = 2.64; imaginary parts are ≈ 0, within errors. The squared coherence is 0.1–0.2:
the zero-driven part is 10–20 % of the twin error's variance in these bands.

**Explicit zeros (approach B).** Twins, band [2, 1000], combs of zeros up to 1500 (`disc1.out`): ζ 2.88 ± 0.15, χ₅ 0.75 ± 0.13,
χ₋₃ 0.13 ± 0.13, χ₋₄ 0.09 ± 0.13, against predictions 2.64, 0.88, 0, 0.

## 5. Frozen holdout test

Files (in the folder): `predictions_twinzeros.json` + `.sha256`, `tz_test.py` (refuses other hashes), `registration.txt`.

Commands (from the folder, foreground, memory cap):

    (ulimit -v 1200000; timeout 600 nice /home/bas/praime/.venv/bin/python -B tz_counts.py hold 2e9 2e10 4e9 4e10)   # repeat until DONE (≈ 3 slices; no -I: it imports sieve_core.py from the folder)
    (ulimit -v 1200000; timeout 600 nice /home/bas/praime/.venv/bin/python -I -B tz_test.py)

Holdout: pairs, triples and SG for p ∈ (2·10⁹, 2·10¹⁰]; prime sums for (4·10⁹, 4·10¹⁰]. Analysis window x ∈ [2·10⁹, 2·10¹⁰]
(u from 21.42 to 23.72, 4,714 grid points). The holdout fits in bounded resumable runs (estimated ≈ 21 min of sieving in
3 slices), so the full (2·10⁹, 2·10¹⁰] is used.

Decisions (fixed in advance; PRIMARY = P1–P4):
* **P1**: twins h=2, band [2,6000], model all: |beta_zeta - 2S(2)| <= 3 sigma AND |beta_zeta - S(2)| > 5 sigma
* **P2**: twins: |loudness - 1| <= 3 sigma AND loudness/sigma >= 5
* **P3**: twins: |beta(chi_-3)| <= 3 sigma AND min |beta(chi_-3) -+ 2S(2)| > 5 sigma
* **P4**: twins: chi^2 of the 13 null coefficients (with their covariance) <= 34.53 (99.9% quantile, 13 dof)
* **S1**: gaps 4, 6, 8, 10, 12, 30: zeta coefficient and loudness, all 12 |z| <= 3.5
* **S2**: Sophie Germain (n, 2n+1): zeta(a=1) = S, zeta(a=2) = S/2, loudness = 1, chi q=3 and q=4 loudness = 1: all |z| <= 3 and chi q=3 loudness/sigma >= 5
* **S3**: triples (0,2,6), (0,4,6): chi q=3 loudness |z| <= 3 and the raw chi_-3 coefficient has the predicted sign (-S3, +S3) at >= 2 sigma, both
* **S4**: twins, band [2,300]: zeta and loudness within 3 sigma (low power)
* **S5**: twins on the explicit zeros (combs of zeta, chi_5, chi_-3, chi_-4 zeros up to 2500), band [2,1900]: predictions 2S, 2S/3, 0, 0, all |z| <= 3
* **S6**: twins per conductor 5, 7, 15, 21: loudness 1 within 3 sigma
* **verdict**: PRIMARY CONFIRMED iff P1-P4 all pass; S1-S6 are reported separately

Falsifier: If on (2e9, 2e10] the twin (h=2) zeta coefficient misses 2S(2) = 2.6406 by more than 3 sigma, or the joint L-function loudness misses 1 by more than 3 sigma, or the mod-3 character coefficient of twins exceeds 3 sigma, the linearised explicit formula for prime pairs is falsified.

Power (expected holdout errors = dry-run errors on 2·10⁸–2·10⁹, a window of the same length, × the ratio of ln x):
expected sigma (holdout): zeta(h2) ~0.114, loudness ~0.042, chi_-3 ~0.105; separation from 'factor 1' ~11.6 sigma, from 'no L-function lines' ~23.8 sigma, from the chi_-3 = +-2S guess ~25.1 sigma; S3 sign tests have power about 0.85 each. The expected σ for every statistic is listed in the JSON.

Dry run on the discovery range (`tz_test.py --dry`, window 2·10⁸–2·10⁹, `tz_test_dry.out`): all ten decisions PASS. The twin ζ
coefficient is 2.690 ± 0.102, loudness 1.010 ± 0.038, χ₋₃ 0.18 ± 0.09, null χ² 20.4/13. Sophie Germain χ₋₃ loudness
0.99 ± 0.06, χ₋₄ 0.82 ± 0.10. Triples raw χ₋₃ −3.13 ± 0.79 / +2.96 ± 0.80. Explicit-zero combs 2.64 ± 0.22, 0.96 ± 0.20,
−0.01 ± 0.20, 0.05 ± 0.19.

## 6. Prior art
Mathematical queries only. Verdict per statement:

| statement | verdict | source (what was read) |
|---|---|---|
| The ζ-part on average over the gap: (1/N)Σ_{r≤N}(ψ_{2r}(x) − 2C_{2r}x) = −(4 + …)Σ_ρ x^ρ/ρ + … | **KNOWN** (conjectural, numerically supported). Our individual coefficient 2S(2r) averages to 4, since the mean of S(2r) over r is 2. | Korevaar & te Riele, "Average prime-pair counting formula", Math. Comp. 79 (2010) 1209–1229 (published version read, CWI repository). They write that they "have not been able to find a good approximation for the remainder ω_{2r}(x)" of an individual gap. No characters. |
| Linear-in-zero terms in the twin/Goldbach major-arc formula, with generalised singular series S′(χ₁, χ₂, m) built from Gauss sums | **KNOWN IN SUBSTANCE**, for exceptional zeros near σ = 1 only, major arcs only, rigorous | Pintz, "A new explicit formula in the additive theory of primes with applications I", arXiv:1804.05561, Theorem 1, (2.4)–(2.11) (arXiv version read; publication status not checked) |
| Twin density as a product of two prime densities with ζ-zero sums | appears informally, without the singular series and without characters | Carella, arXiv:1207.6104, eq. (20) (unrefereed; read) |
| Error of the prime-pair count summed over all shifts relates to ψ(N)² and implies RH | KNOWN (different statement) | Chou, Haag, Huryn, Ledoan, arXiv:2308.14888 (read in HTML). For a fixed gap "nothing that is nontrivial is currently known". |
| Goldbach averages with zero sums (Fujii) | KNOWN for averages over n (different statement) | Fujii, Proc. Japan Acad. 1991 (not opened; seen only through citations) |
| **Individual gap h: E_h(x) ≈ 2S(h)(ψ(x) − x) + …, i.e. ζ-zeros with coefficient 2S(h)** | **NOT FOUND** | none of the sources read has it |
| **Character lines: coefficient 2S(h)μ(q)χ̄(h)/∏(p − 2) for even χ of odd squarefree conductor coprime to h; odd characters and even conductors absent; SG pairs hear χ₋₃ at full loudness and χ₋₄ through 2n + 1; the triples (0,2,6)/(0,4,6) hear χ₋₃ with opposite signs** | **NOT FOUND** | — |
| A numerical spectrum of the twin-prime error at zeros | **NOT FOUND** | Wolf, CMST 17 (2011), arXiv:1107.2809 (sign changes; abstract only). Korevaar–te Riele test averages only. |
| Chebyshev-type bias of twins mod 4 | KNOWN numerically (different statement; no explicit formula) | Sahoo, arXiv:2111.09053 (abstract). Consistent with our prediction that χ₋₄ has no linear term for twins: a mod-4 bias of twins must come from smooth or non-linear terms. |

**Unread or only partly read; each NOT FOUND verdict is provisional until these are checked:**
* Montgomery–Vaughan, "Error terms in additive prime number theory", Quart. J. Math. 24 (1973) 207–216 (paywalled).
* Fujii's 1991 papers (Acta Arith. 58; Proc. Japan Acad. 67A).
* Bogomolny–Keating, "Random matrix theory and the Riemann zeros I/II" (Nonlinearity 1995/1996) (not opened).
* Keating–Smith, "Twin prime correlations from the pair correlation of Riemann zeros", J. Phys. A 52 (2019) (abstract only).
* Montgomery–Soundararajan, "Beyond pair correlation" (2002) (abstract only; the arXiv PDF is glyph-encoded).
* Languasco–Perelli–Zaccagnini on pair correlation and primes in short intervals (not opened).
* Friedlander–Goldston (1995/1996) (not opened).
* Brent, "Irregularities in the distribution of primes and twin primes", Math. Comp. 29 (1975) (AMS page returned 403).
* Suzuki, arXiv:1504.01967 (Goldbach mean value in progressions; abstract only).
* The published version of Pintz's paper.

## 7. Doubts
1. **What is tested.** The factor 2S(h) is the linear response of a product of two densities. A Cramér model with no zeros gives
   the same regression coefficient (2S_model, `cramer_ctrl.out`). The test confirms that pair counts follow the real, zero-driven
   and rigid prime fluctuations like a product of local densities. It cannot separate "zeros" from "any product law"; the explicit
   formula makes the two the same statement for the real primes. In covariance form the conjecture is a consequence of averaged
   Hardy–Littlewood triple asymptotics (sketched in §2.5, not derived here).
2. **Where the weight sits.** Most of the precision comes from high frequencies (γ up to 6000, scales down to 2πx/γ ≈ 10³ at
   x = 10⁶). The long-scale band [2, 300] alone agrees but has errors about 6× larger (twins: ζ 2.48 ± 0.30).
3. **Heuristic, not a theorem.** The quadratic terms and the minor arcs are treated as noise. Their random-walk spectrum,
   about 0.5–0.6 of the Poisson variance, is measured, not derived. Only zeros near σ = 1 are covered by a rigorous major-arc
   formula (Pintz).
4. **Not modelled.** Characters of conductor > 30 are left out (extra noise; their summed power grows like ½ ln Q
   heuristically), and so are prime powers (smooth, absorbed by the nuisance regressors).
5. **Discovery watch items** (chance level among about 150 per-conductor and null coefficients):
   * SG conductor 5 loudness 0.880 ± 0.047 (−2.6σ);
   * h = 10: χ₋₃ null +2.3σ, ζ −1.3σ and loudness −1.6σ;
   * triples: conductor 23 (+2.3σ) and 29 (+2.3σ).
6. **Approximations in the set-up.**
   * The expected holdout errors are scaled from a dry run (× 1.11 for pairs, the ln x ratio). This assumes the twin noise keeps
     the same sub-Poisson factor at 2·10⁹–2·10¹⁰.
   * The secondary decision S3 has power ≈ 0.74 for both triples together.
   * The approach-B combs use zeros only up to 2500 (L-functions) and only four L-functions. Other characters (7, 15, 21, …)
     stay in the noise of that fit.
7. **Bookkeeping.**
   * The regressors θ(ax, χ) for the Sophie Germain pattern need primes up to 2x. In discovery, prime sums (not pair counts) up to
     4·10⁹ were computed. These are regressors tied by the explicit formula to the zeros, not test data, but they lie above
     2·10⁹.
8. **Novelty.** NOT FOUND verdicts are provisional: Montgomery–Vaughan 1973, Fujii 1991, Bogomolny–Keating and the other sources
   listed in §6 were not read in full.
