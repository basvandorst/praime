# ecross: the swapped cross term between the zeros of two elliptic-curve L-functions

Folder: this scratch folder only (nothing in the repo was changed). Zeros used: the 20 development curves of round 144
(`scratchpad/deg2/dev_zeros/`, used data) for everything; the 138 fresh curves (`writeup/data/study/r144/fresh_zeros.npz`)
**only for the one-point sum Σ_γ e^{2iθ_E(γ)} at ξ = 0** (allowed by the brief). No X₁₂ and no one-point line at ξ ≠ 0 was
computed on a fresh curve; the fresh dry runs use synthetic zeros only.

Labels: **proved** (algebra from stated formulas), **recipe** (ratios-conjecture recipe), **conjecture**, **heuristic**,
**computed**, **retrospective** (development data, looked at), **classical** (in print, opened).

| file | output | what it does |
|---|---|---|
| `elib.py` | | θ_E; PARI curve data; c_E; local factors h_p (four-fold sum `h_table`, lines form `h_lines`); ρ(p^v); ln H for all pairs; J; stationary-phase prime sum |
| `s1_curves.py` | `s1_out.txt`, `curves.npz` | 158 curves: N, w, a_p (p ≤ 10⁵), L(1, Sym²E), B_E, c_E; c_E against the direct partial product; h_table = h_lines |
| `s2_onepoint.py` | `s2_out.txt`, `onepoint_all.npz` | the one-point sums on all 158 zero sets against −w c_E T/2π and the alternatives |
| `s3_lines.py` | `s3_out.txt` | development only: lines of the one-point spectrum at ξ = ln(a/b), 12 rationals |
| `s4_prime.py`, `s4b_primes.py` | `s4_out.txt`, `s4b_out.txt` | the one-point sum as a prime sum (zeros vs primes); the resonance from primes alone |
| `s5_freeze.py` | `s5_check_out.txt`, `s5_freeze_out.txt`, `pred_dev.npz`, `pred_fresh.npz`, `freeze_manifest.json` | numerics checks of J; frozen predictions |
| `ecross_test.py` | `dev_test_out.txt`, `run_dev_*`, `dryrun_fresh_*.txt`, `run_fresh_<synthetic>` | the registered test (self-contained: numpy, scipy) |
| `s6_dry_exact.py` | `dry_exact_dev.txt`, `dry_exact_fresh.txt` | exact expectation of the line-modulated synthetic against the frozen law (convention check) |
| `s7_power.py` | `s7_power_out.txt`, `noise_model.json` | noise model from the development pairs; expected values and se on the fresh design |
| `validate.py` | `validate_out.txt` | independent recomputation (Satake/circle-method local factors, direct pair sums, mpmath θ) |
| `round_draft.txt` | | registration draft |

## 0. Verdict in short

1. **The scout's H0 is wrong, by derivation and by data.** The one-point sum of a weight-2 newform is, by the explicit
   formula and stationary phase, the prime sum −(e^{3πi/4}/√2) N^{−1/4} Σ_n Λ_E(n) n^{−1/4} e(−2√(n/N)) w(t_n) — the
   GL(2) Voronoi (Bessel J₁) kernel, not an additive twist e(nα). That kernel *resonates* at dual index 1, and the
   recipe gives a main term: **Σ_{T₁<γ≤T₂} e^{2iθ_E(γ)} = −w c_E (T₂ − T₁)/2π + o(T), c_E = B_E/L(1, Sym²E)** (B_E an
   absolutely convergent Euler product). The bounds of Perelli / Fouvry–Ganguly concern e(nα) and do not apply.
   *Computed* on all 158 zero sets: slope 0.989 ± 0.011 (sharp), 1.002 ± 0.014 (Hann); the c_E-dependence 0.991 ± 0.023;
   no-root-number coefficient 0.005 ± 0.012; "zero mean" is off by 7,500 in χ². Zeros = primes per curve to |G − R| ≤ 0.16
   (three curves), and the main term is visible in the primes alone (no zero used).
2. **The recipe's swapped cross term is exactly the sum, over rational detunings r = a/b, of products of the two
   curves' one-point lines** (proved, an identity of the recipe's Euler product): E X₁₂ = (1/4π²) w₁w₂c₁c₂ J₁₂,
   J₁₂ = ∫H₁₂K, H₁₂(u) = ∏_p Σ_v ρ₁(p^v)ρ₂(p^v) p^{−|v|} p^{iuv} = |L(1 + iu, E₁ × E₂)|² A₁₂(iu). The Rankin–Selberg
   factor is the inner product of the two line spectra; there is no joint Hardy–Littlewood singular series for
   a_{E₁}(p) against a_{E₂}(p′). The one-point lines at r = 2^{±1}, 3^{±1}, 4^{±1}, 5^{±1}, 7^{±1}, (3/2)^{±1} are
   present on the development curves with the predicted amplitudes (pooled slope 1.039 ± 0.022).
3. So **H1 is decided by the mathematics plus the one-point data**, conditional on one conjecture: the fluctuations of
   two different curves' one-point spectra are uncorrelated (a "continuum" with mean zero: shifted convolution sums of
   Λ_{E₁}, Λ_{E₂} over near-collinear prime pairs √(n₁/N₁) ≈ √(n₂/N₂) without main term). Retrospective, development
   (190 pairs): slope 1.059 ± 0.050 (σ = 0.5), H − 1 part 1.065 ± 0.092, reversed pairing 0.03 ± 0.11, no root numbers
   0.006 ± 0.051.
4. Status of the pieces: the one-point main term for GL(2) is a **recipe-level conjecture** (a prime-number theorem for
   the Voronoi resonance; heuristic asymptotic of this kind in Iwaniec–Luo–Sarnak 2000 as quoted by Zhao 2006, upper bound
   N^{5/6} by Zhao; not a theorem); the identity "cross term = Σ line products" is proved from the recipe; the zero-mean
   continuum is a conjecture.

## 1. Conventions (proved; checked numerically)

L(s, E) = Σ λ(n) n^{−s}, λ(n) = a_n/√n; Λ(s) = N^{s/2} Γ_C(s + ½) L(s, E) = w Λ(1 − s), Γ_C(s) = 2(2π)^{−s}Γ(s).
Then L(s) = w X(s) L(1 − s), X(s) = N^{½−s} Γ_C(3/2 − s)/Γ_C(s + ½), and on the line

    X(½ + it) = N^{−it} Γ_C(1 − it)/Γ_C(1 + it) = e^{−2iθ_E(t)},  θ_E(t) = (t/2) ln N − t ln 2π + Im lnΓ(1 + it).

Z_E(t) = e^{iθ_E(t)} L(½ + it) satisfies Z = w·conj Z (real coefficients): real for w = +1, imaginary for w = −1.
2θ_E′(t) = ln N − 2 ln 2π + 2 Re ψ(1 + it) = ln(N t²/4π²) + O(t^{−2}) =: L(t); 2θ_E″(t) ≈ 2/t (precisely
−2 Im ψ′(1 + it)); Stirling: 2θ_E(t) = t ln(N t²/4π²) − 2t + π/2 + O(1/t). N(T) ≈ θ_E(T)/π: mean spacing 2π/L ≈ 0.45 at
L ≈ 14 (band 20–400, N = 300–3,000). Root numbers: PARI `ellrootno`, equal to the set files' w for all 158 curves.

Local data: good p: α_p, β_p = e^{±iθ_p}, λ(p^j) = U_j(cos θ_p); μ_E(p) = −λ_p, μ_E(p²) = 1; multiplicative p:
λ(p^j) = λ_p^j, λ_p = a_p/√p = ±p^{−1/2}, μ_E(p) = −λ_p, μ_E(p²) = 0. F_p(X) := Σ_j λ(p^j) μ_E(p^j) X^j =
(1 − X)(1 − λ(p²)X) (good), 1 − X/p (multiplicative).

## 2. The one-point sum

### 2.1 Recipe (conjecture-level; the same steps as round 143 §2)

The zero density is (1/π)[θ_E′(t) + Re (L′/L)(½ + δ + it, E)], δ → 0⁺, so

    Σ_γ e^{2iθ(γ)} = (1/π)∫ e^{2iθ}θ′dt + (1/2π)∫ e^{2iθ}(L′/L)(½ + δ + it)dt + (1/2π)∫ e^{2iθ}(L′/L)(½ + δ − it)dt.

The first term is a boundary term O(1). In the second, e^{2iθ} cancels the phase of the swapped part of L(s + α):
e^{2iθ(t)}·w X(½ + α + it) = w (N t²/4π²)^{−α}(1 + O(1/t)). The t-average of L(1 − s − α)/L(s + γ) keeps the diagonal
Σ_n λ(n) μ_E(n) n^{−1−γ+α} = F(1 + γ − α), F(s) = ∏_p F_p(p^{−s}). Since F_p(X) = (1 − X)(1 − λ(p²)X) and the
Sym² local factor is [(1 − α²X)(1 − X)(1 − β²X)]^{−1},

    F(s) = B(s)/(ζ(s) L(s, Sym²E)),  B_p = (1 − e₁X)/(1 − e₁X + e₁X² − X³) (good, e₁ = λ_p² − 1, X = p^{−s}),
                                      B_p = (1 − X/p)·(1 − X)^{−1}·(1 − X/p)^{−1} = (1 − X)^{−1} → (1 − 1/p)^{−1} (p | N),

B absolutely convergent at s = 1 (B_p = 1 + O(p^{−2})). So F(1 + x) = x c_E + O(x²) with

    c_E = B_E/L(1, Sym²E) = "∏_{good p}(1 − λ(p²)/p) ∏_{p|N}(1 + 1/p)"  (regularised),

and ∂_α[w (Nt²/4π²)^{−α} F(1 + γ − α)] at γ = α → 0 equals −w c_E. The third term carries e^{4iθ} after its swap and
averages out. **Hence Σ_{T₁<γ≤T₂} e^{2iθ_E(γ)} = −w c_E (T₂ − T₁)/2π + o(T)** (recipe).
For ζ this is Arias de Reyna's −ψ(T/2π) with c = 1; for L(s, χ) round 143's −ε κ.

**Lines.** With the extra factor r^{−it}, r = a/b in lowest terms, the t-average keeps n/m = r, i.e. n = ak, m = bk:
S_r = Σ_k λ(ak) μ_E(bk) (ak)^{−½+α}(bk)^{−½−γ}, which also vanishes at γ = α (it contains 1/ζ(1 + γ − α)). So

    Σ_γ e^{2iθ(γ)} r^{−iγ} = −w c_E ρ(r)/√(ab) · (T₂ − T₁)/2π + o(T),
    ρ multiplicative, ρ(p^v) = p^{|v|/2} S_p(v)/F_p(1/p),  S_p(v) = Σ_{m−h=v} λ(p^m) μ_E(p^h) p^{−(m+h)/2}.

Leading values: ρ(p) = (1 − 1/p)(λ_p − λ(p³)/p)/F_p(1/p) ≈ λ_p, ρ(1/p) = −λ_p(1 − 1/p)/F_p(1/p) ≈ −λ_p, ρ(p^{±2}) = O(1).
So G_E(ξ) := Σ_γ e^{2iθ(γ) − iξγ} has a line at every ξ = ln(a/b); at ξ = ±ln p its size is ≈ ∓w c_E a_p/p · T/2π.

### 2.2 Why there is a main term: the prime side (heuristic stationary phase; computed)

Explicit formula with h(t) = e^{2iθ(t)}w(t): Σ_γ h(γ) = smooth − (1/2π)Σ_n Λ_E(n) n^{−½}[ĥ(ln n) + ĥ(−ln n)],
Λ_E(p^k) = (α_p^k + β_p^k) ln p. ĥ(−ln n) has no stationary point. ĥ(ln n) = ∫e^{i(2θ(t) − t ln n)}w dt is stationary at
2θ′(t_n) = ln n, t_n ≈ 2π√(n/N), with phase value −2t_n + π/2 and 2θ″ ≈ 2/t, so

    G_w ≈ R_w := −(e^{3πi/4}/√2) N^{−1/4} Σ_n Λ_E(n) n^{−1/4} e(−2√(n/N)) w(t_n).

The kernel e(−2√(n/N)) is the asymptotic of the GL(2) Voronoi kernel J₁(4π√(n/N)) at dual index 1. For λ(n) over all
integers, sums Σ λ(n) e(−2√(rn/N)) g(n) resonate (main term ∝ λ(r) X^{3/4}; Iwaniec–Luo–Sarnak, Ren–Ye). Writing
Λ_E = −μ_E ∗ (λ log) (Vaughan type I), each inner sum Σ_m λ(m) ln m e(−2√(dm/N)) resonates at dual index d with
amplitude ∝ λ(d); the outer sum over d produces Σ_d μ_E(d)λ(d)/d^s near s = 1, i.e. exactly F(1 + x) ≈ x c_E: the same
constant as the recipe (heuristic: type I terms only; the type II terms are assumed to have no main term).
Equivalently, the recipe predicts **a prime-number theorem for the Voronoi resonance**:

    Σ_{n≤x} Λ_E(n) e(−2√(n/N)) = (2√2/3) e^{−3πi/4} w c_E (x³/N)^{1/4} (1 + o(1))   (conjecture),

(from G = −w c_E T/2π with dt = π dn/√(nN)). This is not a consequence of any prime number theorem known to us; Zhao
(2006, opened) proves Σ_{n≤N} λ(n)Λ(n)e(α√n) ≪ N^{5/6} log²¹N (level 1) and quotes a conditional heuristic asymptotic
Z N^{3/4} + O(N^{5/8+ε}) of Iwaniec–Luo–Sarnak (2000, Appendix B; not opened: the scanned PDF has no text layer).

**Computed (`s4_out.txt`, three development curves, Hann taper on [20, 400], primes to N·402²/4π² ≈ 1.8–2.4·10⁶,
exact θ):** zeros G vs primes R: 0.83 − 8.53i vs 0.83 − 8.54i (N = 431); −12.77 − 10.41i vs −12.68 − 10.28i (563);
54.18 − 11.44i vs 54.28 − 11.41i (590); primes alone (k = 1) within 0.5. **From primes alone (`s4b_out.txt`, no zero):**
S(x) at x = 2.4·10⁶ for N = 590 (w = −1, c = 1.88): 16,650 + 13,536i against the law 15,607 + 15,607i; N = 782: −10,062 −
13,895i against −9,892 − 9,892i; N = 563: −3,968 − 4,914i against −5,007 − 5,007i; at the non-resonant frequency
1.13 × 2 the sums do not grow (e.g. −1,920 − 713i for N = 590). Random-sign scale of the raw sums ≈ 6,000.

### 2.3 Data: the one-point law on all 158 zero sets (`s2_out.txt`; allowed by the brief)

Band (20, 400], T = 380, law −w c_E T/2π (median |law| 77), Hann: half of it. c_E (`s1_out.txt`): median 1.27, range
0.12–3.42; B_E/L(1, Sym²E) agrees with the direct partial product to 10⁵ to 0.1 % (rms log 0.0009) and to 10⁶ (validate).

| reading | 158 curves | fresh 138 | development 20 |
|---|---|---|---|
| slope of Re G (sharp) on the law | 0.9893 ± 0.0108 | 0.9875 ± 0.0116 | 1.0013 ± 0.0323 |
| slope, Hann taper | 1.0022 ± 0.0140 | 0.9982 ± 0.0148 | 1.0285 ± 0.0462 |
| [−wT/2π, −w(c_E − 1)T/2π]: does the amplitude follow c_E? | 0.988 ± 0.021, 0.991 ± 0.023 | 0.986, 0.989 | 1.001, 1.002 |
| [law, law without w] | 0.991, 0.005 ± 0.012 | 0.986, −0.004 ± 0.013 | 1.029, 0.073 ± 0.027 |
| Im G (sharp), mean | 0.26 ± 1.31 | 0.74 ± 1.36 | −3.1 ± 4.4 |
| rms(Re G − law) sharp / Hann | 13.9 / 8.6 | 14.1 / 8.5 | 12.8 / 9.7 |

"No main term" (H0 of the scout): Σ(Re G)²/resid. var = 7,512 on 158 curves. The zeros are rigid at this frequency:
the scatter 13.9 is half the random-phase value √(number of zeros) ≈ 29. Means by root number: w = +1: −67.4 (law −67.8);
w = −1: +95.6 (law +97.0).

### 2.4 The lines on the development curves (`s3_out.txt`; development only)

| r | 2 | 1/2 | 3 | 1/3 | 4 | 1/4 | 5 | 1/5 | 7 | 1/7 | 3/2 | 2/3 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| slope (sharp) | 1.05(6) | 0.98(7) | 1.06(8) | 1.03(4) | 0.99(8) | 1.13(6) | 0.90(13) | 1.00(6) | 0.91(11) | 1.11(14) | 1.05(14) | 1.07(10) |

Pooled 1.039 ± 0.022 (Hann 1.026 ± 0.033); without the root number 0.005 ± 0.024 extra coefficient; control frequencies
(ξ = ±0.35, ±0.55, ±0.9, ±1.3): mean Re G −2.9 ± 1.7, rms 18.4. Lines at ξ ≠ 0 were **not** computed on fresh curves (they
enter the X₁₂ law and are registered as blind readings).

## 3. The two-curve recipe (proved from the recipe)

For non-isogenous E₁, E₂ the ratio average ⟨L(½+α+it, E₁) L(½+β−it, E₂)/(L(½+γ+it, E₁) L(½+δ−it, E₂))⟩ has four
recipe terms: no swap; swap E₁ only (phase X₁(½ + α + it) ~ e^{−2iθ₁}); swap E₂ only (e^{+2iθ₂}); both (e^{−2iθ₁(t₁) +
2iθ₂(t₂)}·(N₁t²/4π²)^{−α}(N₂t²/4π²)^{−β}). In X₁₂ = (1/T)ΣΣ e^{2iθ₁(γ₁) − 2iθ₂(γ₂)}K(γ₂ − γ₁) the weight cancels the double
swap's phase exactly (as in degree 1: **yes, the swap phase is e^{2i(θ₂ − θ₁)}** and X₁₂ isolates it); single swaps leave
e^{∓2iθ} and average out; the conjugate pairing carries e^{4i(θ₁ − θ₂)}; the unswapped term (no phase) enters only through
orbits n ≈ e^{L} (Heisenberg range) with weights λ₁(p)λ₂(p)(ln p)², mean zero by Rankin–Selberg orthogonality, and with
e^{2i(θ₁ − θ₂)} oscillating in t unless N₁ = N₂ (estimate < 10⁻⁴ per pair); the smooth d̄₁d̄₂ term is killed by K
(e^{−σ²L²/2} < 10⁻⁴ at σ = 0.5, L ≥ 8; relevant only for N₁ = N₂).

Double swap, t-average: the diagonal of L(1−s−α, E₁) L(s−β, E₂)/(L(s+γ, E₁) L(1−s+δ, E₂)) is the Euler product of
Y_p = Σ_{m+k=n+h} λ₁(p^m) λ₂(p^n) μ₁(p^h) μ₂(p^k) p^{−m(½−α)−n(½−β)−h(½+γ)−k(½+δ)}. Each curve contributes a factor
1/ζ(1 + γ − α), 1/ζ(1 + δ − β) that vanishes on the diagonal, so ∂_α∂_β at γ = α, δ = β keeps only the product of
the two derivatives: with x = α + β (m − h = n − k ⇒ the exponent depends on x only),

    c₁c₂ H₁₂(x) = ∏_p Y_p(x)/(1 − 1/p)²,  Y_p(x) = Σ_v S₁(v) S₂(v) p^{xv}   (regrouping by v = m − h = n − k),

**so H₁₂(x) = ∏_p h_p(x), h_p = Σ_v ρ₁(p^v)ρ₂(p^v) p^{−|v|} p^{xv}: the cross term is Σ_{r=a/b} [line₁(r)]·[line₂(r)]**
(proved; the r-th Fourier coefficient of H is ρ₁(r)ρ₂(r)/(ab)). The constant: with lines a_j(r) = −w_j c_j ρ_j(r)/√(ab)
per unit T/2π, (σ/√2π)(1/2π)Σ_r a₁(r)a₂(r)e^{−σ²ln²r/2} = (1/4π²)w₁w₂c₁c₂ Σ_r h_r σ√2π e^{−σ²ln²r/2} = (1/4π²)w₁w₂c₁c₂J₁₂
(the one-point normalisation reproduces the recipe's constant exactly). Leading Euler factor:
h_p = 1 + λ₁λ₂(p^{x−1} + p^{−x−1}) + O(p^{−2}) (computed: deviation 1.0·10⁻⁵ at p = 2003), so
**H₁₂(iu) = |L(1 + iu, E₁ × E₂)|² A₁₂(iu)** with A₁₂ absolutely convergent; Y_p(0) is not 1 (unlike the self term),
there is no pole (E₁ × E₂ is entire and non-vanishing on Re s = 1 for non-isogenous curves).

    E X₁₂ = (1/4π²) w₁w₂ c₁c₂ J₁₂,  J₁₂ = ∫ H₁₂(u) e^{−u²/2σ²} du,  c_j = B_j/L(1, Sym²E_j).

- Root numbers: enter as w₁w₂ (the swap of each L contributes its w). The sign of the one-point law cancels in the
  product, so "H1 with the wrong sign" (−law) is not produced by any sign error in the one-point convention.
- Direction: X pairs line r of E₁ with line r of E₂ (same ξ). The alternative "reversed pairing" (r with 1/r) gives
  ρ₁(p)ρ₂(1/p) ≈ −λ₁λ₂: the p-lines with the opposite sign (frozen as `predrev`).
- Bad primes: S_p(v) with the multiplicative local data; p | N₁N₂ handled exactly (code: always exact).
- Equal conductors (fresh: two classes with N = 1003): θ₁ = θ₂; nothing changes in the law (leakage terms above).
- **Checks (computed):** h_table (the four-fold sum) = h_lines (the lines form) to 3.3·10⁻¹⁴ for p ≤ 997, all a_p values,
  good and bad; Σ_r over r = 2^a3^b5^c of the line products = ∫h₂h₃h₅K to 3.6·10⁻¹⁴ (σ = 0.5).

## 4. The prime side of X₁₂ and what conjecture it is

X₁₂ = (σ/√2π)(1/T)∫e^{−σ²ξ²/2} G₁(ξ) conj G₂(ξ) dξ (exact, as in review 47). By §2.2 each G_j(ξ) is the prime sum
−(e^{3πi/4}/√2)(e^ξ/N_j)^{1/4}Σ_n Λ_{E_j}(n) n^{−1/4} e(−2√(n e^ξ/N_j)) w. So X₁₂ is a bilinear form
Σ_{n₁,n₂} Λ_{E₁}(n₁) Λ_{E₂}(n₂)(n₁n₂N₁N₂)^{−1/4} e(−2(√(n₁/N₁) − √(n₂/N₂)))·Φ over near-collinear pairs
|2π(√(n₁/N₁) − √(n₂/N₂))| ≲ σ, i.e. |N₂n₁ − N₁n₂| ≲ σ t N₁N₂/2π² (a shift range growing with the height, unlike degree 1).
Its mean has two parts:

- **Lines**: products of the two one-point resonance means at the same rational detuning r. These carry the whole
  recipe law (§3 identity). In degree 1 they were products of prime-number-theorem-in-progressions means (theorem);
  here they are products of **GL(2) Voronoi resonances over primes** (conjecture, §2.2), measured directly by the
  one-point readings.
- **Continuum**: correlations of the fluctuations of Λ_{E₁}(n₁)e(…) and Λ_{E₂}(n₂)e(…) beyond the line products. The
  recipe says mean zero. This is a statement on shifted convolution sums of two different Hecke eigenvalue sequences over
  primes along N₂n₁ ≈ N₁n₂ r (no main term). For integers instead of primes, Σ λ₁(n)λ₂(n + h) has no main term for
  cusp forms (classical: Good, Jutila, Blomer; not re-opened here); over primes it is open. A joint singular series of a_{E₁}(p) against
  a_{E₂}(p′) would need a correlation of the Frobenius traces of two non-isogenous curves at related primes; at equal
  primes the Rankin–Selberg orthogonality Σ_p λ₁(p)λ₂(p) ln p/p = O(1) excludes it on average, and nothing in the
  recipe produces one: **the recipe's Rankin–Selberg factor is not a singular series but the inner product of the two
  curves' line spectra.**

**Decision.** H0 ("no main terms, hence mean of X₁₂ zero unless the functional equations couple the curves") rests on the
wrong exponential twist and is refuted by the one-point data (7,500 in χ²; the lines at r ≠ 1 on the development
curves). H1 follows from the recipe; given the one-point lines (data), its only remaining content is the zero-mean
continuum. The X₁₂ test therefore measures (i) the lines at r ≠ 1 on fresh curves (blind), (ii) the H − 1 part of the
law (blind: reading X5 removes the measured r = 1 product), (iii) the continuum mean (X6).

## 5. Numerical predictions (frozen)

σ in t units as in round 143: σ = 0.5 (primary; ≈ 1.1 mean spacings at L = 14) and 1.0. J from exact local factors
(lines form) for p ≤ 3,000 and all bad primes, first-order lines 2λ₁λ₂cos(u ln p)/p plus the non-oscillating second-order
term −λ₁²λ₂²/p² for 3,000 < p ≤ 10⁵; u ∈ [0, 6] step 0.05, H(−u) = conj H(u).

Accuracy (`s5_check_out.txt`, 190 development pairs, in units of the prediction pred = w₁w₂c₁c₂J/4π²): P = 3·10⁴ vs 10⁵:
max 6.4·10⁻⁵; P₀ = 1,000 vs 3,000: 9.3·10⁻⁷; without the second-order term: 3·10⁻⁵; u-step 0.025 to 8: 4·10⁻¹⁴ (σ = 0.5),
8·10⁻⁸ (σ = 1). Independent recomputation (`validate.py`, no shared code: local factors from the Satake parameters by the constant
term on the unit circle (512-point FFT), p ≤ 2,000 exact, first-order lines to P = 2·10⁵, u-step 0.02 to 8σ; `validate_out.txt`):
12 pairs (6 development, 6 fresh) × 2 σ: |ΔJ| ≤ 2.1·10⁻³ on |J| up to 45, |Δpred| ≤ 4.1·10⁻⁵ (the difference is the
truncation 10⁵ vs 2·10⁵). c_E: direct partial product to 10⁶ / frozen c_E = 0.9997–1.0001 (12 curves). X by the direct
double sum over zero pairs (exact erf window weight, mpmath θ) vs the test's grid: ≤ 1.2·10⁻⁴ (8 development pairs; the
grid's end-point term); one-point sums with mpmath θ: 1.2·10⁻¹¹. Expected truncation error beyond 10⁵: below 10⁻⁴ in pred (noise per pair 0.07).

Frozen (`s5_freeze_out.txt`; sha256 in `freeze_manifest.json` and in `ecross_test.py`):
- `pred_fresh.npz` ef45af2e11c905acbfb287dd8fd786322f1b29afe251d77e81f97369af765010: 138 curves, 9,453 pairs.
  σ = 0.5: pred mean 0.0120, sd 0.1434, range −0.574 … 0.802; H = 1 part (pred1 = w₁w₂c₁c₂σ√2π/4π²) sd 0.0796,
  corr(pred, pred1) 0.712; J/(σ√2π) median 0.98, range −72 … 344 (large J goes with small c₁c₂: F_p(1/p) small when
  |a₂| = 2, |a₃| = 3; pred stays bounded). σ = 1.0: pred sd 0.220, corr 0.802.
- `pred_dev.npz` 25833cd6c2e9334db0fb9cc98e2a8235464dc4f42005ee4e3bf298bc4030c3e2: 20 curves, 190 pairs (validation).
- per pair: J, Jrev, pred, pred1, prednoroot (c₁c₂J/4π²), predrev; per curve: w, c, N, L(1, Sym²), B, the line amplitudes
  of the 13 registered r (with and without w), and the p = 2, 3 coefficients for the line-modulated synthetic (J23).

Spread over pairs and sign (σ = 0.5): 49 % of the fresh pairs have pred < 0 and 37 % pred < −0.02 (w₁w₂ = −1 for 4,745 of 9,453 pairs); mean |pred| 0.100. Nine pairs have equal conductors (N = 418, 1003, 1259 (three), 1405, 2098, 2555, 2699).

## 6. Development (retrospective) and dry runs

`dev_test_out.txt` (190 pairs, 20 curves; used data):

| reading | σ = 0.5 | σ = 1.0 | law |
|---|---|---|---|
| X1 slope of Re X on the law | 1.0585 ± 0.0503 | 1.0320 ± 0.0497 | 1 |
| X2 [H = 1 part, H − 1 part] | 1.046 ± 0.109, 1.065 ± 0.092 | 1.016 ± 0.077, 1.051 ± 0.112 | 1, 1 |
| X3 reversed-pairing coefficient | 0.029 ± 0.109 | 0.065 ± 0.110 | 0 |
| X4 no-root coefficient | 0.006 ± 0.051 | 0.003 ± 0.062 | 0 |
| X5 H − 1 after removing the measured r = 1 product | 1.066 ± 0.091 | 1.044 ± 0.113 | 1 |
| X6 continuum mean | −0.0091 ± 0.0062 | −0.0108 ± 0.0080 | 0 |
| X7 Im X mean | 0.0001 ± 0.0072 | 0.0025 ± 0.0097 | 0 |

Dry runs (synthetic zeros; `run_dev_*`, `dryrun_fresh_*.txt`, `dry_exact_*.txt`):
- Poisson (no term): every reading ≈ 0 (fresh: X1 0.007 ± 0.008, X5 0.007 ± 0.011, O1 −0.003 ± 0.015, O2 0.005 ± 0.011).
  X6 is −0.0052 ± 0.0015 there: X6 assumes the law for the lines not measured, so its H0 value is −mean(pred − Σ_R law
  products), not 0 (registered with that expectation).
- Jittered Gram lattice (the r = 1 line only, amplitude −w c_E): fresh, seed 7: X2 = [1.076 ± 0.042, −0.002 ± 0.011], X5
  −0.002 ± 0.010, O1 1.038 ± 0.019, O2 0.011 ± 0.011; seed 11: X2 = [0.985 ± 0.036, 0.003 ± 0.011], X5 0.010 ± 0.010,
  O1 1.002 ± 0.016, O2 0.002 ± 0.011 (the 1.038 was a 2-se fluctuation): the test separates "r = 1 only" from the full law.
- Line-modulated Poisson (all lines r = 2^a3^b with the law's amplitudes, intensity clipped at 0 on 0.5 % of t): the exact
  expectation of X (smoothing of the intensity, no sampling) agrees with the frozen law restricted to p = 2, 3
  (w w c c J23/4π²): slope 0.9985, corr 0.99997 (development, σ = 0.5); fresh design: slope 1.0002, corr 0.99996
  (σ = 0.5), 0.9997 (σ = 1); the sampled synthetic X on its exact expectation: 0.995 (σ = 0.5), 1.004 (σ = 1), residual
  sd 0.090. This checks the conventions of J (pairing direction, normalisation, window) against the X computation end to end.
- X6 under H0 and H1-const is −mean(pred − Σ₁₃ law products) = −0.0055 (σ = 0.5), −0.0049 (σ = 1); Poisson dry run −0.0052 ±
  0.0015. The 13 registered lines carry the law with correlation 0.955 (σ = 0.5).

## 7. Noise and power (`s7_power_out.txt`; noise model from the development pairs only)

Residual of Re X after the law (σ = 0.5): sd 0.071 = per-curve additive 0.018, per-curve multiplicative 0.105 (× the
pair's mean), pair 0.070; σ = 1: sd 0.099. The model over-predicts the development jackknife (0.065 vs 0.050), so the
fresh se below are conservative by about 30 %. Fresh design (9,453 pairs, jackknife over 138 curves):

| reading (σ = 0.5) | se | H1 | H0 | −H1 | no root | H = 1 only | reversed | half |
|---|---|---|---|---|---|---|---|---|
| X1 slope | 0.023 | 1 | 0 | −1 | 0.178 | 0.395 | 0.148 | 0.5 |
| X2 H − 1 coefficient | 0.020 | 1 | 0 | −1 | 0.232 | 0 | −0.454 | 0.5 |
| X3 reversed coefficient | 0.015 | 0 | 0 | 0 | 0.109 | 0 | 1 | 0 |
| X4 no-root coefficient | 0.014 | 0 | 0 | 0 | 1 | −0.034 | −0.007 | 0 |
| X5 H − 1 after r = 1 removal | 0.043 | 1 | 0 | — | 0.248 | 0 | −0.442 | 0.5 |
| X8 mean of Re X | 0.003 | 0.0120 | 0 | −0.0120 | 0.0856 | 0.0024 | 0.0002 | 0.0060 |
| O2 lines r ≠ 1 pooled | 0.009 | 1 | 0 | — | −0.190 | 0 | — | — |

So H0 is rejected by X1 at ≈ 40 se, "H = 1 only" at ≈ 26 se, "no root numbers" by X4 at ≈ 70 se, the reversed pairing
by X3 at ≈ 65 se, "half" at ≈ 22 se; the blind H − 1 coefficient X5 is measured to ≈ 4 %. σ = 1: X1 se 0.023, X5 0.047.

## 8. Prior art (mathematical queries; opened sources marked)

- Weishäupl, "The weak Gram law for Hecke L-functions", arXiv:2201.08898 (**opened**, ar5iv): for a level-1 Hecke
  eigenform of weight ≥ 12, Σ_{T<t_ν≤2T} ω(t_ν)Z(t_ν) over even/odd Gram points = ±T/π + O(T^{3/4+ε}), i.e. the mean of
  (−1)^ν Z(t_ν) is 2a(1): the n = 1 bias of the approximate functional equation, theorem-level, for values of Z at Gram
  points. It does not treat zeros or Σ_γ e^{2iθ(γ)}.
- Zhao, "Oscillations of Hecke eigenvalues at shifted primes", Rev. Mat. Iberoam. 22 (2006), arXiv:math/0508129
  (**opened**, PDF text): Σ_{n≤N} λ(n)Λ(n)e(α√n) ≪ N^{5/6}log²¹(3N) (level 1, weight ≥ 12); quotes Iwaniec–Luo–Sarnak's
  conditional heuristic asymptotic S(N) = Z N^{3/4} + O(N^{5/8+ε}), Z ≠ 0. This is the prime-side form of the one-point
  main term (for the resonant α); the zero-side statement with the explicit constant −w c_E (c_E = B/L(1, Sym²E)) and
  root number was not found.
- Iwaniec–Luo–Sarnak, Publ. IHES 91 (2000) (arXiv:math/9901141 opened: Hypothesis 4 on Σ_{p≡a(c)} e(2√p/c), a GL(1) sum;
  the Appendix B quoted by Zhao is not in the arXiv v1 and the numdam scan has no text layer: **not read**).
- Perelli; Fouvry–Ganguly (Compositio 2014, strong orthogonality of μ, additive characters and cusp-form coefficients):
  linear twists e(nα); search snippets only; they do not bear on the Bessel-resonant twist e(2√(n/N)).
- Iwaniec–Luo–Sarnak, Ren–Ye (arXiv:1411.6022), Liu–Zhang (Proc. AMS 2017): resonance for Σ λ_f(n)e(α n^β) over integers
  at β = 1/2, α = ±2√(q/N) (search results; not opened).
- Rudnick–Sarnak 1996 (GUE n-level correlations of one primitive automorphic L-function; search results): no statement on
  two different L-functions found. Murty–Perelli (pair correlation in the Selberg class; search results): one function.
- Conrey–Simons–Smolyarenko, "Pair correlation of different L-functions" (aimath.org/~kaur/publications/u8.pdf):
  HTTP 404 again; the indexed snippet: derived from a ratios conjecture, equal levels assumed.
- Bui–Keating–Smith 2016 / Smith 2016 (read in round 144): one Selberg-class function.
- Liu–Williams–Zaharescu, arXiv:2407.14656 (abstract **opened**): pair correlation of one shifted newform L-function;
  no cross correlation.
- Milinovich–Ng, arXiv:1306.0854 (**opened**, ar5iv, partly): discrete moments of L′(ρ, f); no phase sums over zeros.
- Murmurations: not searched further (they concern averages of a_p over families against conductor, not correlations
  between two fixed curves).

**Verdicts.** One-point main term for GL(2) on the prime side (resonance over primes): **KNOWN IN SUBSTANCE** as a
conditional heuristic (ILS via Zhao; constant not seen). Zero-side one-point law with −w c_E, c_E = B/L(1, Sym²E), and its
lines at rational ξ: **NOT FOUND**. The recipe's swapped cross term between two elliptic curves (w₁w₂c₁c₂|L(1+iu, E₁×E₂)|²A₁₂)
and its identity with the product of line spectra: **NOT FOUND** (CSS unpublished, equal level, unopened). Any numerical
test of a lower-order cross term between zeros of two GL(2) L-functions: **NOT FOUND**.

## 9. Doubts

1. The one-point law for GL(2) is recipe-level; the data (158 curves, 1 %) support it, but its proof would need a prime
   number theorem for the Voronoi resonance (type II sums with e(2√(dm/N)) below N^{3/4}), beyond current methods.
2. The X₁₂ test is largely pre-decided: its H = 1 part is the product of the per-curve one-point sums, which were computed
   on all fresh curves before registration (allowed). Only the r ≠ 1 lines (X5, O2/O3) and the continuum (X6) are blind.
   A reviewer may judge the X₁₂ part redundant with a blind one-point line test; it costs seconds, so it is kept.
3. The development X readings run 3–7 % high (1.059 ± 0.050, lines 1.039 ± 0.022): within 1–2 se, but the same sign
   everywhere. If real, candidates are the finite-height terms of the recipe and the window/edge terms of the lines.
   A 2 % theory allowance is stated in the draft (not in the bands).
4. The noise model comes from 20 curves (190 pairs); the decisions use the data's jackknife.
5. ρ(r) and c_E depend on a_p, L(1, Sym²E) (PARI) and B (Euler product to 10⁵): checked against the direct partial
   product (0.1 %); the PARI Sym² normalisation at bad primes was checked in round 144.
6. Equal-conductor pairs and pairs with large |J| (small c) have no special treatment; X9 bins by |ln N₂/N₁|.
7. GRH and completeness of the zero lists are assumed (round 144's rule; all 138 complete after deduplication).
8. The prime-only check of §2.2 is per curve and noisy (2–4 random-sign units); the zeros–primes identity is the sharp
   check (|G − R| ≤ 0.16).
9. ILS's Appendix B was not read: the claim that it concerns exactly this resonance at α = 2 (level 1) rests on Zhao's
   quotation.
