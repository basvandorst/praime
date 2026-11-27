# Given the primes, positivity is sharp: characters, Ramanujan's Δ and Maass forms from Weil's form (27 November 2026)

**In short.** Weil's explicit formula turns the Riemann hypothesis of an L-function into the statement that a quadratic form is positive. Here the primes are given and the local data are unknown: one probability measure per prime, on the value χ(p) of a real character or on the Satake angle of a form on GL(2) or GL(3). The only other input is the archimedean factor (and, for odd forms, the root number). All findings are numerical (float64, primes below 40 at most); nothing is proved.

1. **Real characters.** The smallest conductor that keeps the form positive is 2.9999999996 (odd, support 14) and 4.9999999994 (even, support 24), and the optimum is the character mod 3, mod 5 at the primes of the support (the last prime excepted), with 0 at the ramified prime. At support 8 the optimum is not yet the character (2.9999977, 4.99989).
2. **Level 1, weight as a real number.** Nothing is feasible below weight 12; at 12 the optimum is τ(p)/p^(11/2) for p = 2, 3, 5, 7 to 3 decimals at support 16. With root number −1, weight 14 is excluded and the first feasible weight is 18.00, with the Hecke eigenvalues of the weight-18 form within 0.01 at support 16.
3. **Maass forms of SL(2,ℤ), eigenvalue as a real number.** Nothing is feasible below the first eigenvalue of each parity. At support 16 (primes up to 13), on a grid of step 0.01 and with the allowance t* ≥ −0.002: with parity and root number alone, six of the ten eigenvalues below 20 are isolated and the four odd ones above 15.9 lie in one feasible stretch; when the forms already found are used as constraints on the next one (their Rankin–Selberg products, rounds 171 and 172), the feasible set consists of ten pieces of length at most 0.05, one at each eigenvalue. Run as a chain that feeds in only its own output (no answer key), this returns the ten eigenvalues within 0.0052 (nine within 0.0016) and a₂, a₃ within 0.0014 for nine of them (0.022 for the tenth).
4. **SL(3,ℤ).** One of the five lowest forms, (13.5966, 4.7647), lies within 0.035 of the edge of the feasible set (t* = +0.0016 there) and its a₂ is returned within 0.011; the others do not.

5. **Weight 2 with whole-number a_p.** The smallest admissible level is 10.995, with a₂, a₃, a₅, a₇ = −2, −1, 1, −2 (the elliptic curve of conductor 11); with root number −1 it is 36.91, with −2, −3, −2, −1 (the curve 37a).

Eight rounds were registered before their runs (165 to 172). The first strong claims (islands at every eigenvalue from L(f) and Sym² alone, every low form of SL(3,ℤ) on the edge) failed as registered; rounds 170 and 171 passed, round 172 passed two of its three decisions and missed the third by 0.002; sections 5 to 8g give each verdict. For the odd character problem a computer-assisted bound with no hypothesis on zeros gives a smallest conductor of at least 3 − 3.81·10⁻²² (section 8f; referee's report in section 9). In print: the positivity method with local corrections and bounds that come close (r > 7.2, conductor > 10.323, weight near 11). Not found in print by two searches: that the bound is exact when the primes are given, the local data are the unknown and the test function is optimised, and that the optimum is the local data of the extremal form. Two reviews with their own codes reproduced the numbers and corrected the note (section 9).

## 1. What is in print

- **The method.** Explicit-formula positivity bounds for discriminants and conductors: Odlyzko, Serre, Poitou 1977, Diaz y Diaz 1980 (minimum over the splitting types at each prime, one fixed kernel), Mestre 1986 (conductor of an elliptic curve > 10.323 with one test function), Miller 2002 (no cusp form on SL(2,ℤ) with r ≤ 7.2; regions for SL(3,ℤ)), Chenevier–Lannes and Chenevier–Taïbi (classification of algebraic forms of level 1 and small weight).
- **How close it comes.** Chenevier–Lannes (arXiv 1409.7616, §IX.3): without prime terms their criterion changes sign between motivic weights 9 and 11, "sa précision est assez étonnante". Mestre: 10.323 against 11. Diaz y Diaz did not optimise because the quadratic minima were known. Odlyzko 1990 leaves the optimal kernel open.
- **The unstructured problem** (any positive measure on the prime side): N. Johns, Zenodo, October 2026, and its companion note on critical conductors: the optimum is 1 for ζ and lies strictly below the true value in higher degree (2.981 against 3, 4.985 against 5), with optimal measures that are not arithmetic. Section 5e of [../bootstrap.md](../bootstrap.md) found the same numbers.
- **Not found** (prior-art reports of 26–27 November, both "partly in print"): exactness of the structured bound; the local data of the extremal form as the optimum; feasible sets that break into pieces at the automorphic spectrum. Sources that could not be opened: Serre's "Minorations de discriminants", Miller 1996, Bian 2010, the Chenevier–Lannes book beyond §IX.

## 2. Set-up (derived; checked numerically)

Test functions f live on [−a, a] in u = log x; X = e^(2a) is the support; g(w) = ∫ f(v + w) f̄(v) dv. For an L-function with Gamma factor ∏ᵢ Γ_R(s + κᵢ) (the κᵢ complex), conductor q, no pole, and −L′/L = Σ c(n) n^(−s):

`Q(f) = (log q) g(0) + Σᵢ A_κᵢ(g) − Σ_n n^(−1/2) [c(n) g(log n) + c̄(n) g(−log n)]`

`A_κ(g) = −(γ + log π) g(0) + ∫₀^∞ [2e^(−2u) g(0) − e^(−(1/2+κ)u) g(u) − e^(−(1/2+κ̄)u) g(−u)] / (1 − e^(−2u)) du`

- Derivation: Gauss's formula for ψ at (1/2 + κ + it)/2 integrated against |f̂(t)|², as in section 2 of [../bootstrap.md](../bootstrap.md); the two prime sums come from the two sides of the functional equation, so for a form that is not self-dual the matrix is Hermitian and the test functions are complex.
- For a true L-function, Q(f) = Σ_γ |f̂(γ)|², which is ≥ 0 under its Riemann hypothesis. A negative eigenvalue of the truncated matrix (basis e^(iπku/a), |k| ≤ N) proves that the form is not positive; a positive matrix proves nothing (Rayleigh–Ritz).
- **The unknown.** For each prime p ≤ X one probability measure μ_p on the local parameter; the coefficient at p^m is the mean of the local weight under μ_p, the same measure for every m and for every L-function in the problem:
  - real character: x ∈ [−1, 1], weight x^m log p (x = 0: ramified);
  - GL(2), angle θ ∈ [0, π]: L(f): 2 cos mθ; Sym²: 1 + 2 cos 2mθ; Sym³: 2 cos mθ + 2 cos 3mθ; Sym⁴: 1 + 2 cos 2mθ + 2 cos 4mθ;
  - GL(3), conjugacy class of SU(3) with eigenvalues e^(iφⱼ): L(F): τ_m = Σⱼ e^(imφⱼ); adjoint: |τ_m|² − 1.
- **The optimum.** t* = the largest t such that every form of the problem is ≥ 0 with conductor e^(−kt) (k = 1 for the standard L-function, k for Sym^k, 2 for the adjoint). q* = e^(−t*) is the smallest admissible conductor; at level 1 the test is the sign of t*. All constraints are linear in the measures, so t* is the value of a convex problem (Kelley cutting planes with an LP, HiGHS). By minimax it equals the conductor bound "with local corrections" optimised over all positive combinations of test functions (the dual optimum has rank at least 2 in each parity block, so one test function does not attain it; second review) against the worst local data.
- **Root number.** For an odd Maass form or a holomorphic form with root number −1 the L-function vanishes at 1/2, so Q_f(φ) ≥ |∫φ|²; sections 3 to 7 do not use this, section 8 does.
- Gamma factors used. Maass form of SL(2,ℤ), parameter r, parity ε: L(f): Γ_R(s+ε±ir); Sym²: Γ_R(s)Γ_R(s±2ir); Sym³: Γ_R(s+ε±3ir)Γ_R(s+ε±ir); Sym⁴: Γ_R(s)Γ_R(s±2ir)Γ_R(s±4ir). Holomorphic weight κ, with Γ_C(s) = Γ_R(s)Γ_R(s+1): L(f): Γ_C(s+(κ−1)/2); Sym²: Γ_R(s+1)Γ_C(s+κ−1); Sym³: Γ_C(s+3(κ−1)/2)Γ_C(s+(κ−1)/2); Sym⁴: Γ_R(s)Γ_C(s+κ−1)Γ_C(s+2κ−2). SL(3,ℤ), parameters (λ₁, λ₂, λ₃ = −λ₁−λ₂): L(F): ∏ Γ_R(s+iλⱼ); adjoint: Γ_R(s)² ∏_d Γ_R(s±id), d = λ₁−λ₂, 2λ₁+λ₂, λ₁+2λ₂.
- **Checks of the set-up** (computed): the engine reproduces the eigenvalues of [../weil.py](../weil.py) for ζ and is positive for χ₋₄; every true form in the answer keys passes every form at supports 8, 16, 32 (sections 4, 5, 7); the SL(3) forms fail when their coefficients are conjugated.

## 3. Real characters: conductors 3 and 5 exactly, with the character values ([chars.py](chars.py), [chars_precise.py](chars_precise.py))

One form, Γ_R(s + κ), no pole, one measure on [−1, 1] per prime. Smallest admissible conductor ([chars_precise_out.txt](../../../../data/study/maass/chars_precise_out.txt)):

| support X | odd (κ = 1) | even (κ = 0) |
|---|---|---|
| 8 | 2.9999977254 | 4.9998922240 |
| 14 | 2.9999999996 | 4.9999998848 |
| 24 | (LP failed at this tolerance) | 4.9999999994 |

From support 14 (odd) and 24 (even) on, the optimal measures are point masses ([chars_out.txt](../../../../data/study/maass/chars_out.txt), support 40): odd: −1, 0, −1, +1, −1, +1, −1, +1, −1, −1, +1, +1 at p = 2, 3, …, 37, the character mod 3; even: −1, −1, 0, −1, +1, −1, −1, +1, −1, +1, +1 at p = 2, …, 31, the character mod 5 (the last prime of the support, 37, is at −0.94 against −1).

- This changes finding 2 of [../bootstrap.md](../bootstrap.md) ("without a pole positivity did not pin the points"): that search left the positions free; with the primes given and one measure per prime the Dirichlet case is pinned.
- At support 8 the optimum is not the character: a measure that is not a point mass does better by 6.3·10⁻⁷ (odd) and 2.0·10⁻⁵ (even) in log q, and at support 14 (even) by 2.0·10⁻⁸ (second review, with a dual bound; our 0.9185 at p = 7 shows the same). The second review's values with exact moments instead of a grid: 2.9999981, 3.0000000, 3.0000000 (odd) and 4.9998949, 4.9999999, 5.0000000 (even) at supports 8, 14, 24.
- **What removes the loss of the unstructured problem (2.981 and 4.985) is not the coupling of the powers of a prime.** With the weights at p, p², p³, … decoupled (each in [−1, 1] times log p) the second review still finds 3.000000 and 5.000000 at support 24 (2.999998 and 4.998633 at support 8). The unstructured problem differs in that the positions are free and the object carries a pole; with the primes given and bounded weights there is no gap. (Our first reading, that tying the powers to one measure removes the gap, was wrong.)

## 4. Level 1 with the weight as a real number: 12, and Ramanujan's τ ([scan_hol.py](scan_hol.py), [levels.py](levels.py), [truth_hol.py](truth_hol.py))

Registered as the control of round 165. Four forms, support 8 ([hol_X8.jsonl](../../../../data/study/maass/hol_X8.jsonl)): t* rises from −2.390 at κ = 2 through −0.0366 at κ = 11.75 and −0.00012 at κ = 12 to +0.0229 at κ = 12.25, and stays positive up to κ = 30. So the smallest admissible weight is 12.00 within the grid.

| p | 2 | 3 | 5 | 7 |
|---|---|---|---|---|
| optimum at κ = 12, support 8 | −0.5297 | +0.5975 | +0.6969 | −0.4095 |
| optimum at κ = 12, support 16 | −0.530 | +0.598 | +0.691 | −0.376 |
| τ(p)/p^(11/2) | −0.5303 | +0.5987 | +0.6912 | −0.3765 |

- Weight 14 carries no cusp form but is feasible in the convex problem at support 8 (t* = 0.165): the measures there are not point masses.
- **Levels above 1.** The same curve read as a conductor: smallest admissible level 10.906, 4.963, 2.967, 1.951 at weights 2, 4, 6, 8 at support 8 and 10.936, 4.969, 2.983, 1.958 at support 16 ([levels_out.txt](../../../../data/study/maass/levels_out.txt), every prime may also take the two Steinberg local factors of a prime that divides a squarefree level: a relaxation); the first newforms have levels 11, 5, 3, 2. Not exact at these supports: the optima are mixtures. The true newforms of levels 11, 5, 3, 2, 14, 15, 7 (coefficients from PARI) pass all four forms at their levels ([truth_hol_out.txt](../../../../data/study/maass/truth_hol_out.txt)).

## 5. Round 165: Maass forms of SL(2,ℤ) with one measure per prime ([round165.txt](../../../../data/study/preregister/round165.txt), [maass.py](maass.py), [scan.py](scan.py), [decide.py](decide.py))

Registered before the scan (commit 0a48eda8); pilot of four points declared. Grid r = 1 to 20 (step 0.05 at supports 8 and 16, 0.1 at 32), both parities, four forms; answer key: the ten forms with r < 20 from the Booker–Strömbergsson–Venkatesh files ([lmfdb_key.json](../../../../data/study/maass/lmfdb_key.json)).

- **Truth passes** ([truth_out.txt](../../../../data/study/maass/truth_out.txt)): with the true Hecke eigenvalues every one of the ten forms has smallest eigenvalue ≥ −2·10⁻¹⁴ in L(f) and Sym² and ≥ 0.015 in Sym³ and Sym⁴, at supports 8, 16, 32.
- **Below the first eigenvalue nothing is feasible** (D1, pass at every support): the largest t* below r₁ − 0.3 is −0.128, −0.172, −0.205 (odd) and −0.211, −0.216, −0.261 (even) at supports 8, 16, 32. In print: r > 7.2 (Miller 2002).
- **The first eigenvalue of each parity is isolated.** Odd: t* = −0.0022 at r = 9.55 (support 8), with −0.013 at 9.50 and −0.029 at 9.75, and nothing feasible up to 12.15. Even: −0.008 at 13.80, nothing else feasible below 17.70. The maximum is a cusp narrower than the grid, so the registered D2 (a grid point within 0.05 with t* ≥ −0.002) fails by 0.0002 to 0.023; at the exact eigenvalues the true forms are feasible (first item).
- **Hecke eigenvalues** (D5; support 16, grid point nearest to the eigenvalue):

| form | a₂ | a₃ | a₅ | a₇ |
|---|---|---|---|---|
| 9.5337 found | −1.065 | −0.464 | −0.281 | −0.738 |
| true | −1.068 | −0.456 | −0.291 | −0.745 |
| 12.1730 found | +0.277 | −1.202 | +0.038 | +0.436 |
| true | +0.289 | −1.202 | +0.040 | +0.448 |
| 13.7798 found | +1.542 | +0.257 | +0.725 | −0.227 |
| true | +1.549 | +0.247 | +0.737 | −0.261 |

  The third odd form (14.3585) lies inside a feasible stretch and its eigenvalues are not returned (differences up to 1.07): D5 fails for odd, passes for even.
- **Islands (D3): fail for odd forms.** Even, support 16: a cusp at 13.80 (t* = −0.011 on the grid) and the feasible pieces [17.70, 17.90] and [19.35, 19.40], against the true 13.780, 17.739, 19.423; 7 feasible grid points in all (29 at support 8, 3 at support 32). Odd: everything from 12.20 to 20 is feasible at supports 8, 16 and 32 (157, 157 and 79 feasible grid points). D4 (a local maximum at every eigenvalue) fails for odd as well.
- Only L(f) and Sym² are active (at all but a handful of grid points); Sym³ and Sym⁴ are slack, and so are the twists of both by the six real characters of conductor ≤ 8 ([twists.py](twists.py), [twists_out.txt](../../../../data/study/maass/twists_out.txt): t* unchanged to 5 decimals at four values of r).
- Reading: the feasible stretch above 12.2 consists of mixtures. At r = 13 the optimal measure at p = 2 has two atoms, at θ ≈ 0.9 and 1.7; at the edges 9.55 and 12.2 every measure is a point mass. Section 6 tests this.

VERDICT of round 165 as registered: D1 pass; D2 fail (grid); D3 fail (odd); D4 fail (odd), pass (even); D5 fail (odd), pass (even).

## 6. Round 166: one genuine angle per prime ([round166.txt](../../../../data/study/preregister/round166.txt), [pm.py](pm.py))

Registered after round 165 (commit 3c2bcdde): every prime carries one angle, so the unknown is an Euler product and the problem is no longer convex. Branch and bound over windows of angles: the upper bound of a box is the convex problem restricted to it, the lower bound the value at one point of it. Two forms (L(f), Sym²), odd parity, support 16, grid step 0.1. A grid point is feasible if its upper bound is ≥ −0.002.

**Incomplete.** The run was stopped after 34 of 80 grid points (6,735 s; the points near the threshold need the most boxes), in favour of round 168. Completed: r = 12.0 to 13.5 and 16.0 to 17.7 ([pm_eps1_X16_k2.jsonl](../../../../data/study/maass/pm_eps1_X16_k2.jsonl)).

- r = 12.0, …, 13.5 (16 points): all infeasible, upper bounds between −0.031 and −0.0024. In the convex problem all of 12.2 to 13.5 was feasible (t* up to +0.21). The true eigenvalue 12.1730 lies between the grid points 12.1 and 12.2, both infeasible: the island is narrower than the grid, so the registered E2 (the grid point just above each eigenvalue is feasible) fails there.
- r = 16.0, …, 16.6 (7 points): infeasible. r = 16.7, …, 17.7 (11 points): feasible. This component begins at the first grid point above the eigenvalue 16.6443, and the best point there has a₂ = −1.52, a₃ = +0.92 against the true −1.540, +0.977 (E1 and E4 hold for the one component seen).
- Not run: 13.6 to 15.9, 17.8 to 20, support 32, the even parity beyond three points, and the holomorphic control. E3 (fewer feasible points than in the convex problem) holds on the part that was run: 11 feasible grid points against 32.

VERDICT of round 166: not decided (incomplete). What it shows: mixtures were the reason for most of the feasible stretch of round 165.

## 7. Round 167: SL(3,ℤ) ([round167.txt](../../../../data/study/preregister/round167.txt), [gl3.py](gl3.py), [gl3_scan.py](gl3_scan.py), [gl3_decide.py](gl3_decide.py))

Registered before any solve (commit 5e4ba0eb). Two forms: L(F), which is not self-dual (Hermitian matrix, complex test functions), and the adjoint L-function L(F × F̃)/ζ. One measure per prime on the conjugacy classes of SU(3) (a 30 × 30 grid of angles modulo permutations). Support 8. Answer key: 1,050 forms from the data repository of Farmer, Koutsoliotas, Lemurell and Roberts ([gl3_key.json](../../../../data/study/maass/gl3_key.json), [gl3_sources.md](../../../../data/study/maass/gl3_sources.md)); the reader agent found by a functional-equation test that 242 entries of the file carry the coefficients of the dual form, and corrected them.

- **Check of the engine and of the key** ([gl3_truth.py](gl3_truth.py)): the Hermitian engine agrees with [engine.py](engine.py) on two self-dual data to 2.5·10⁻¹⁴. The eight forms of smallest Laplace eigenvalue pass both forms (L(F) between 1.4·10⁻¹⁰ and 2.7·10⁻⁵, adjoint between 0.27 and 1.9); with conjugated coefficients seven of the eight fail (−0.26 to −1.25), the eighth is self-dual.
- **G3, pass:** none of the 237 coarse grid points with λ₁² + λ₁λ₂ + λ₂² ≤ 180 is feasible. The smallest value on the feasible set is 222.25, at (13.5, 2.5); the lowest form has 239.31. In print: Miller's bound 87.
- **The form (13.5966, 4.7647) is close to the edge:** t* = +0.0016 at its exact parameters (so it lies inside the feasible set, by that much), an infeasible grid point at distance 0.035, and the optimum has a₂ = 1.040 − 0.382i against the true 1.048 − 0.375i (difference 0.011; a₃: 0.048).
- **The others are not** (G2 and G4 fail): the nearest infeasible point of the local grid lies at 0.185, 0.348, 0.176 from forms 1, 2, 4 (smallest Laplace eigenvalue first), and form 5 has none within its 9 × 9 grid; a₂ is returned within 0.081 and 0.058 for forms 1 and 4 and not at all for forms 2 and 5 (1.31, 0.79). FKL's region, which assumes |a_p| ≤ 3, leaves 0.13 to 0.52.
- G1 (truth feasible) passes: t* = +0.019, +0.167, +0.002, +0.030, +0.261 at the five forms.
- Both forms are active at 844 of the 861 coarse points. On the self-dual line λ₂ = 0 the feasible set begins at λ₁ = 16.5; the first form there is the symmetric square of the first Maass form of SL(2,ℤ), at 2 × 9.5337.

VERDICT of round 167 as registered: G1 pass, G2 fail, G3 pass, G4 fail. One of the five lowest forms lies within 0.035 of the edge of the convex problem at support 8 and its Hecke eigenvalue is returned within 0.011; the problem with genuine conjugacy classes was not run.

## 8. Round 168: the central zero ([round168.txt](../../../../data/study/preregister/round168.txt), [decide168.py](decide168.py))

The review of round 165 (section 9) found that the registered problem left out a piece of archimedean data, the root number. An odd Maass form has Λ(s) = −Λ(1 − s), so L(1/2, f) = 0 and, under its Riemann hypothesis, Q_f(φ) ≥ |∫φ|². Round 168 (registered before my runs, commit 90c7c6ea; the review's own numbers at support 8 declared) replaces Q_f by Q_f − |∫φ|² in the convex problem of round 165. Two forms, odd parity, r = 9.00 to 20.00 in steps of 0.01.

- **Truth passes** (H1; [truth_center_out.txt](../../../../data/study/maass/truth_center_out.txt)): the seven odd forms have smallest eigenvalue ≥ −2.6·10⁻¹⁴ with the central zero removed, at support 16.
- **Feasible set** ([decide168_X8_out.txt](../../../../data/study/maass/decide168_X8_out.txt), [decide168_X16_out.txt](../../../../data/study/maass/decide168_X16_out.txt)):

| support | feasible pieces (grid 0.01) | true odd eigenvalues |
|---|---|---|
| 8 | [9.53, 9.54], [12.17, 12.19], [14.28, 14.43], [15.82, 20.00] | 9.5337, 12.1730, 14.3585, |
| 16 | {9.54}, [12.17, 12.19], [14.35, 14.37], [15.93, 20.00] | 16.1381, 16.6443, 18.1809, 19.4847 |

- "Feasible" means t* ≥ −0.002 at a grid point. At the grid points themselves t* is slightly negative (−0.0013, −0.0006, −0.0002 at 9.54, 12.17, 14.36 at support 16): the true islands are much narrower than the grid (t* = −1.2·10⁻⁶ at r = 14.3585 by the second review), and a piece narrower than the step that lies between two grid points cannot be excluded by a scan.
- **The three lowest odd eigenvalues are islands in this sense, and each returns its Hecke eigenvalues** (support 16, nearest grid point): a₂, a₃ = −1.068, −0.455 (true −1.068, −0.456); +0.288, −1.201 (true +0.289, −1.202); −0.230, +0.694 (true −0.231, +0.696).
- **Above 15.9 the stretch remains** (H2, H3, H4 fail there, as the declared weak spot at support 8 and also at support 16): 384 of the 415 feasible grid points at support 16 are farther than 0.03 from every eigenvalue, and the Hecke eigenvalues of the four forms inside the stretch are not returned.

VERDICT of round 168 as registered: H1 pass; H2, H3, H4 fail (all failures above r = 15.9).

**Controls with exact answers** (registered before their runs, commit 2c3b6f97):
- *Level 1, weight as a real number, with the central zero* ([hol_X8_center.jsonl](../../../../data/study/maass/hol_X8_center.jsonl), [hol_X16_center.jsonl](../../../../data/study/maass/hol_X16_center.jsonl)). Cusp forms of level 1 with root number −1 exist at weights 18, 22, 26, … and not at 14. Found: t* = −0.470 at κ = 14, −0.026 at 17.75, +0.00001 at 18 and +0.026 at 18.25 (support 8); −0.027 at 17.75 and −0.0002 at 18 (support 16). So weight 14 is excluded and the smallest admissible weight is 18.00 within the grid: the prediction holds. The second review, with exact moments: smallest weight 11.9999983 and 17.99853 at support 8, 12.0000000 and 17.9999986 at support 16. The optimum at κ = 18: −1.4585, −0.3771, −1.1750 at p = 2, 3, 5 (support 8) and −1.4575, −0.3783, −1.1697, +0.2032 at p = 2, 3, 5, 7 (support 16), against the true −1.4584, −0.3770, −1.1745, +0.2115 (coefficients −528, −4284, −1025850, 3225992 from PARI). The agreement at support 8 is partly an accident of our angle grid: the true optimum there does better than the form by 1.5·10⁻⁴ and has a₅ = −1.182, a₇ = 0.186 (second review). At support 16 the optimum is the form's data within 0.01.
- *Weight 2 with the central zero* ([levels_center_out.txt](../../../../data/study/maass/levels_center_out.txt)): smallest admissible level 32.46 at support 8 and 34.76 at support 16; the first newform of weight 2 with root number −1 has level 37. The prediction "below 37" holds; "within 3 % at support 16" fails (6.0 %).

**Even forms on the fine grid** (not registered; [scan_eps0_X16_k2_J61.jsonl](../../../../data/study/maass/scan_eps0_X16_k2_J61.jsonl), [scan_eps0_X32_k2_J61.jsonl](../../../../data/study/maass/scan_eps0_X32_k2_J61.jsonl)). Convex problem of round 165, r = 13.00 to 20.00 in steps of 0.01, support 16: the feasible set is {13.78}, [17.69, 17.94], [19.32, 19.43], one piece for each of the three even eigenvalues 13.7798, 17.7386, 19.4235, and nothing else (round 165 covers r < 13). At support 32 the pieces are {13.78}, [17.70, 17.88], [19.38, 19.43]. Hecke eigenvalues at support 16: at 13.78: 1.549, 0.247, 0.736 (true 1.549, 0.247, 0.737); at 19.42: −0.698, 1.555, −0.039 (true −0.693, 1.562, −0.038); at 17.74, inside its wider piece, they are not returned.

**Together, at support 16** (primes up to 13), on the grids of step 0.01 (0.05 below r = 9 and, for even forms, below 13): of the ten eigenvalues below 20, six are isolated (9.53, 12.17, 13.78, 14.36, 17.74, 19.42), each in its own piece of {t* ≥ −0.002}, with no piece that holds no eigenvalue; the four odd ones above 15.9 lie in one remaining stretch.

## 8b. Round 170: weight 2 with whole-number a_p — conductors 11 and 37 ([round170.txt](../../../../data/study/preregister/round170.txt), [integral.py](integral.py))

Registered before its run (commit 5f673abf). The problem of section 4 at weight 2 (four forms, Steinberg columns), with the measure at each prime restricted to the angles for which √p · 2 cos θ is a whole number (the a_p of a rational newform; Hasse's bound). Smallest admissible level ([integral_out.txt](../../../../data/study/maass/integral_out.txt)):

| support | no condition at the centre | with the central zero |
|---|---|---|
| 8 | 10.947 | 36.260 |
| 16 | 10.992 | 36.710 |
| 32 | 10.995 | 36.906 |
| first rational newform | level 11 | level 37 (root number −1) |

- Optimum at support 32, weight of the largest atom in brackets: a₂, a₃, a₅, a₇ = −2 (1.00), −1 (0.99), +1 (0.98), −2 (0.98), the values of the elliptic curve of conductor 11. With the central zero: −2 (1.00), −3 (1.00), −2 (0.89), −1 (0.73), the values of the curve 37a.
- At the larger primes of the support the optimum is not the curve's (a₁₃ = 5 or 6 against 4; a₁₉ = −2 against 0): the measures there are still mixtures and the bound stays 0.05 % and 0.25 % below 11 and 37.
- Without the restriction the bounds were 10.936 and 34.76 (sections 4 and 8). In print: conductor > 10.323 with one test function and the same restriction at p = 2, 3 (Mestre 1986).

VERDICT of round 170 as registered: I1 pass (10.992 and 10.995 lie in [10.95, 11.0005]; the four a_p as predicted). I2 pass at support 32 (36.906 in [36.8, 37.005]; the four a_p as predicted), not at support 16 (36.710; a₇ = 0 or +1).

## 8c. Beyond the first object: prescribed values (not registered; [chars_cond.py](chars_cond.py), [chars_depth.py](chars_depth.py))

Real characters again, support 40, with the values at the first k primes prescribed and one measure at the others. The bound is compared with the smallest conductor of a real primitive character of the same parity that has the prescribed values.

- **Values at 2 and 3 prescribed** ([chars_cond_out.txt](../../../../data/study/maass/chars_cond_out.txt)): the bound is exact (ratio 1.00000, the whole character returned at 11 or 12 of the 12 primes) in 8 of the 18 cases, namely conductors 3, 4, 7, 8 (odd) and 5, 8, 12, 13 (even); 16.995 against 17 in one more; and far below in the other nine (4.98 against 19, 5.07 against 24, 7.28 against 21, 6.30 against 11, …, 68.9 against 73).
- **Depth** ([chars_depth_out.txt](../../../../data/study/maass/chars_depth_out.txt)): for each fundamental discriminant D with |D| ≤ 44, the smallest k at which the bound reaches 0.9999·|D|: k = 0 for −3 and 5; 1 for −4, −7, 8; 2 for −8, 12, 13; 3 for −11, −15, 24; 4 for −20, 21; 5 for −24; 6 for 17 and −23; not reached with k ≤ 6 for the thirteen others (for example −19: 18.91; 28: 27.91; 37: 36.99; −43: 19.7).
- Reading: exactness is a property of the first object. A later one becomes extremal once enough of its local data are prescribed, and then the rest of its data is returned; how much is needed varies from one character to the next. This is the same picture as the Maass forms above 15.9.

## 8d. Round 169: central zero and genuine angles together ([round169.txt](../../../../data/study/preregister/round169.txt))

Registered (commit a019323f). [pm.py](pm.py) with the central zero, odd parity, support 16, r = 15.9, 16.0, …; stopped by its time limit at 18.3 (25 of 42 points, 3,268 s; [pm_eps1_X16_k2_center.jsonl](../../../../data/study/maass/pm_eps1_X16_k2_center.jsonl)).

- 15.9 to 16.6 (8 points): infeasible (upper bounds −0.032 to −0.0025). The convex problem of round 168 was feasible from 15.93 on.
- 16.7 to 18.3 (17 points): feasible. The best point at 16.7 has a₂, a₃ = −1.53, +0.95 (the form at 16.6443 has −1.540, +0.977): the piece that begins at 16.64 runs on to the right.

VERDICT of round 169 as registered: K1 fails (feasible grid points up to 18.3 that are farther than 0.1 from every eigenvalue); incomplete above 18.3. Genuine angles alone do not close the stretch at support 16.

## 8e. Round 171: the sequential bootstrap — all ten eigenvalues below 20 ([round171.txt](../../../../data/study/preregister/round171.txt), [mixed_scan.py](mixed_scan.py), [decide171.py](decide171.py))

Registered (commit 6ae98af0, the pilot with four known forms declared; even parity added before its run, commit 470dcca9). The forms found so far are used as constraints on the next one. For a known form g and the unknown form f, L(f × g) has Gamma factor Γ_R(s + d ± i(r + r_g)) Γ_R(s + d ± i(r − r_g)) (d = 0 for equal parities, 1 otherwise), conductor 1, no pole, and the weight (2 cos mθ_p)(2 cos mθ_p(g)) log p at p^m, which is linear in the unknown measure of f. Problem at each r: L(f) (with the central zero for odd f), Sym², and L(f × g) for every true form g of either parity with r_g ≤ r − 0.05; the Hecke eigenvalues of g are taken from the answer key. Convex problem, support 16, step 0.01.

- **Truth passes** ([mixed_truth_out.txt](../../../../data/study/maass/mixed_truth_out.txt)): for the seven odd forms every form of the problem, including up to nine Rankin–Selberg forms, has smallest eigenvalue ≥ −2.6·10⁻¹⁴.
- **Odd, r = 15.90 to 20.00** ([decide171_odd_out.txt](../../../../data/study/maass/decide171_odd_out.txt)): the feasible set is [16.14, 16.18], [16.65, 16.69], [18.18, 18.23], [19.49, 19.53] (21 grid points; 408 in round 168), against the eigenvalues 16.1381, 16.6443, 18.1809, 19.4847.
- **Even, r = 13.00 to 20.00** ([decide171_even_out.txt](../../../../data/study/maass/decide171_even_out.txt)): the feasible set is {13.78}, {17.74}, {19.43} (3 grid points; 39 without the Rankin–Selberg forms), against 13.7798, 17.7386, 19.4235.
- **Hecke eigenvalues at the first feasible grid point of each piece:**

| eigenvalue | a₂ found | true | a₃ found | true |
|---|---|---|---|---|
| 16.1381 | +1.156 | +1.162 | −1.279 | −1.282 |
| 16.6443 | −1.539 | −1.540 | +0.973 | +0.977 |
| 18.1809 | +0.374 | +0.374 | +0.101 | +0.102 |
| 19.4847 | −1.681 | −1.700 | −0.610 | −0.615 |
| 13.7798 | +1.549 | +1.549 | +0.247 | +0.247 |
| 17.7386 | −0.767 | −0.765 | −0.975 | −0.978 |
| 19.4235 | −0.686 | −0.693 | +1.563 | +1.562 |

VERDICT of round 171 as registered: M1, M2, M3 pass for both parities; the added prediction (the piece at 17.74 shorter than 0.1) holds: it is one grid point.

**The spectrum below 20.** Rounds 165, 168 and 171 together, at support 16 (the primes 2, 3, 5, 7, 11, 13): on the grids, {t* ≥ −0.002} consists of ten pieces, of length at most 0.05, one at each of the ten eigenvalues of the modular surface below 20, and each piece returns a₂ and a₃ of its form within 0.02. Caveats: a grid of step 0.01 (0.05 below r = 9 and, for even forms, below 13) cannot exclude a feasible piece narrower than the step; in round 171 the Hecke eigenvalues of the forms below (for p ≤ 13) enter from the answer key (section 8g removes this); float64; and the hypotheses of section 10.

## 8g. Round 172: the self-contained chain — ten eigenvalues and their Hecke eigenvalues with no key ([round172.txt](../../../../data/study/preregister/round172.txt), [chain.py](chain.py), [decide172.py](decide172.py))

Registered (commit 649ab761). Round 171 without the answer key: every form used as a constraint is the bootstrap's own output, its estimate r̂ of the eigenvalue and its â_p for p ≤ 13, and errors are free to propagate. Input: the primes 2, 3, 5, 7, 11, 13, the Gamma factors, the parity and, for odd forms, the central zero. Grid points are processed in ascending r for both parities; when a feasible piece begins, a form is recorded (r̂ = the left crossing of t* = −0.0005, or the maximiser of t* if the piece is narrower than that; â_p = the means of the optimal measures at r̂). Only the grid points that rounds 165 and 168 (which use no key) left feasible are solved: 505 points, 1,700 s. The key is read at the end ([chain_X16_out.txt](../../../../data/study/maass/chain_X16_out.txt), [decide172_out.txt](../../../../data/study/maass/decide172_out.txt)).

| r̂ found | parity | true eigenvalue | â₂ | true a₂ | â₃ | true a₃ | largest error, p ≤ 7 | p = 11, 13 |
|---|---|---|---|---|---|---|---|---|
| 9.53443 | odd | 9.53370 | −1.0679 | −1.0683 | −0.4566 | −0.4562 | 0.0005 | 0.0047 |
| 12.17376 | odd | 12.17301 | +0.2895 | +0.2893 | −1.2017 | −1.2019 | 0.0005 | 0.0120 |
| 13.77916 | even | 13.77975 | +1.5494 | +1.5493 | +0.2471 | +0.2469 | 0.0007 | 0.0048 |
| 14.35814 | odd | 14.35851 | −0.2313 | −0.2309 | +0.6960 | +0.6956 | 0.0004 | 0.0049 |
| 16.13660 | odd | 16.13807 | +1.1613 | +1.1619 | −1.2816 | −1.2820 | 0.0022 | 0.0050 |
| 16.64339 | odd | 16.64426 | −1.5397 | −1.5402 | +0.9769 | +0.9775 | 0.0011 | 0.0114 |
| 17.73718 | even | 17.73856 | −0.7657 | −0.7655 | −0.9771 | −0.9778 | 0.0056 | 0.0138 |
| 18.17936 | odd | 18.18092 | +0.3745 | +0.3741 | +0.1012 | +0.1020 | 0.0008 | 0.0056 |
| 19.42290 | even | 19.42348 | −0.6939 | −0.6928 | +1.5610 | +1.5624 | 0.0021 | 0.0055 |
| 19.48992 | odd | 19.48471 | −1.6780 | −1.7002 | −0.6087 | −0.6146 | 0.0331 | 0.0618 |

- **N1, pass:** the chain records exactly ten forms below 20, one at each true eigenvalue, of the right parity; the largest error of r̂ is 0.0052, and 0.0016 for the first nine.
- **N2, fail by 0.002:** a₂ and a₃ are within 0.0014 for the first nine forms; for the tenth â₂ is off by 0.022 (allowed: 0.02). That piece, [19.48, 19.53], was located by its maximum, 0.005 to the right of the eigenvalue, because t* at its first grid point was below the threshold of the rule.
- **N3, pass:** the 28 feasible grid points lie in ten pieces, one per form, and none lies outside [r_j − 0.02, r_j + 0.15].
- The declared weak spot (a fed-back error pushes a true form out) did not occur at this support and grid.

VERDICT of round 172 as registered: N1 pass, N2 fail (one form, by 0.002), N3 pass.

## 8f. A certificate for conductor 3 ([theory/exact3.md](theory/exact3.md); computer-assisted, one code; referee's report pending when this was written)

For the odd character problem at infinite support (positivity required for every smooth test function of compact support) a proof agent turned the numerical exactness into a bound with no hypothesis on zeros.

- **Lemma (proved there).** Let Ξ_χ(t) = Λ(1/2 + it, χ₋₃) and ĝ = Ξ_χ² · B with B an even polynomial, B ≥ 0. Then (log 3) g(0) + A₁(g) = 2 Σ_n χ(n) Λ(n) n^(−1/2) g(log n) exactly, wherever the zeros of L(s, χ₋₃) lie: Λ·Λ′·B is entire and odd under s ↦ 1 − s, so its integral over a vertical line vanishes.
- **Consequence.** A family of local measures that is feasible at conductor q satisfies log q ≥ log 3 − 2·Defect(g)/g(0), where Defect(g) = Σ_p [P_p(χ(p)) − min over [−1, 1] of P_p] and P_p(x) = Σ_m x^m log p · p^(−m/2) g(m log p).
- **Certified** (ball arithmetic at 600 bits; the primes above 400 bounded by 10⁻²⁹⁷): the smallest admissible conductor is at least 2.7400 with B = 1, 3 − 3.1·10⁻⁶ at degree 8, 3 − 4.5·10⁻¹⁰ at degree 12, 3 − 9.1·10⁻¹⁷ at degree 16 and **3 − 3.81·10⁻²²** at degree 22.
- Under the Riemann hypothesis for L(s, χ₋₃) the character itself is feasible at 3, so the value lies in [3 − 3.81·10⁻²², 3]. That it equals 3 needs Defect → 0 along a sequence of polynomials; this is not proved (the defect stalls at degrees 18 and 20 and falls again at 22).
- Status: one code and one author; 156 checks of its own ([theory/validate_exact3.py](theory/validate_exact3.py)); not yet refereed; no literature search for this lemma. The same construction for ζ is section 5f of [../bootstrap.md](../bootstrap.md) (in print: Johns).

## 9. The review of round 165 ([review/review.md](review/review.md))

An adversarial agent recomputed the problem with its own code (Legendre basis, the form in x-space, an exact moment formulation without an angle grid, its own barrier solver).

- No bug found in [engine.py](engine.py) or [maass.py](maass.py). Its values of t* at ten points agree with ours within 7·10⁻⁴; at the first eigenvalue it finds t* = 1.7·10⁻⁹, so our −0.0003 there came from the angle grid.
- **Width of the first island at support 8** (its number, not recomputed by us): t* ≥ 0 exactly for r in [9.533695256, 9.533696035], width 7.8·10⁻⁷; the true value is 9.533695261. Four primes and the Gamma factor give the first eigenvalue to seven digits.
- **The root number** (section 8).
- **A lemma** (theirs; the one-line check is ours): Q_odd − Q_even = ½ ∫ |φ̂(t)|² Σ_ρ sech(π(t − ρ)) dt ≥ 0, with ρ = ±r for L(f). It follows from the reflection formula, Re[ψ(3/4 + ix) − ψ(1/4 + ix)] = π sech(2πx). So without the central zero the odd problem admits every even eigenvalue, which is the bump of the odd curve at 13.8.
- A larger basis (N = 24, 48, 96) does not change the sign of t* at r = 13, 15, 18 (0.1521, 0.1520, 0.1520, …): the feasible stretch of round 165 is not an artifact of a small basis.
- With the wrong parity the first odd form fails (−0.94); the first even form with odd Gamma factors does not (+1.9·10⁻⁴), as the lemma says it must not.
- Not checked by the review: the holomorphic control, supports above 8, anything beyond float64.

**Second review** ([review2/review2.md](review2/review2.md)), of this note, again with its own code (Legendre basis, exact moments per prime, a barrier solver with a dual bound for every value):
- reproduced: the characters (table of section 3), the holomorphic curves to 10⁻⁴, t* at nine values of r for the odd problem with the central zero at support 16 within 1.3·10⁻⁴, the seven odd forms passing with the central zero, and t* = +0.0017, a₂ = 1.039 − 0.383i for the SL(3) form (13.5966, 4.7647); the Gamma factors agree with PARI's lfunsympow at weight 2;
- corrected, and taken over above: exactness for the characters only from support 14 and 24 on; our explanation of the loss of the unstructured problem (wrong); "feasible pieces" are pieces of {t* ≥ −0.002} on a grid; the agreement at weight 18 and support 8 was partly an accident of the angle grid; the hypotheses in section 10; the minimax statement; "on the edge" for the SL(3) form;
- not checked by it: the fine scans as a whole, support 32, G2 and G3, rounds 166, 169, 170, 171, interval arithmetic.

## 10. Size of the result

- **Known before:** the method, local corrections with a fixed kernel, and the remark that the bounds come close (section 1).
- **New here, as far as two searches found** (numerical; nothing proved):
  1. with the primes given and one measure per prime, the bound is exact to 9 digits for the two smallest real characters (from support 14 and 24 on) and to the grid for the first cusp form of level 1 of each sign (weights 12 and 18) and the first Maass form of each parity;
  2. the optimum is the local data of that form: the values of the character, τ(p), Hecke eigenvalues to 2–4 digits;
  3. for Maass forms of SL(2,ℤ) the feasible set at support 16 separates six of the ten eigenvalues below 20 into single pieces once the parity and the root number are put in, and all ten once the forms below are used as constraints (round 171), also when only the bootstrap's own output is fed back (round 172);
  4. on SL(3,ℤ) one of the five lowest forms lies within 0.035 of the edge and its a₂ is returned within 0.011.
- **What it is not.** Not a proof of anything about L-functions: feasibility of the true objects assumes the Riemann hypothesis of every L-function used (L(χ), L(f), Sym², and Sym³, Sym⁴, the adjoint, L(f × g) where they enter) and real Satake angles (a theorem of Deligne for holomorphic forms; for the Maass forms of the key it was checked on the data, |a_p| ≤ 1.74). Not exact above level 1 at the supports reached (10.94 against 11, 34.8 against 37). Not a set of islands for the odd forms above 15.9 at support 16 without the sequential constraints. The registered claims in their strong forms failed: D2, D3, D4 (odd), D5 (odd), G2, G4, H2, H3, H4, K1; passed: D1, G1, G3, H1, I1, I2 (support 32), M1, M2, M3, N1, N3 and the two weight controls; N2 missed by 0.002.
- **Why it may matter.** A bound that is exact, with the extremal object as its optimum, is the situation of the Cohn–Elkies bound in dimensions 8 and 24 and of Johns' statement for ζ. If it holds in general for "the first form of each archimedean type", it is a converse theorem in which positivity replaces the functional equations of twists; and the method is a way to compute first eigenvalues and Hecke eigenvalues from the primes and the Gamma factor alone.

## 11. Open points

1. Exactness for conductor 3: section 8f gives 3 − 3.81·10⁻²² (computer-assisted, not refereed); open: that the defect tends to 0, a referee, a second code, and the same certificates for conductor 5, weight 12 and the first Maass form.
2. Whether the levels 11 and 37 are reached exactly as the support grows (10.995 and 36.91 at support 32 with whole-number a_p), and the levels 5, 3, 2 at weights 4, 6, 8.
3. Rounds 166 and 169 are incomplete. The chain of round 172 solved only the grid points that rounds 165 and 168 left feasible, at one support and one grid; a second code, a finer grid, eigenvalues above 20, levels above 1 and the same chain on SL(3,ℤ) were not run.
4. Interval arithmetic for the infeasibility statements; a second code for sections 3, 4, 7 and 8 (the review covered round 165 at support 8 only).
5. SL(3,ℤ) with genuine conjugacy classes, larger support, and the forms of the key beyond the five lowest.
6. Sources not opened by the prior-art searches (section 1).

## 12. Validation

[validate.py](validate.py) reads the data files and checks every number quoted in sections 3 to 8g (literature values and the reviews' own numbers excepted); for section 8f it checks that the quoted numbers are those of [theory/exact3.md](theory/exact3.md), whose own validation is [theory/validate_exact3.py](theory/validate_exact3.py).
