# Given the primes, positivity is sharp: characters, Ramanujan's Δ and Maass forms from Weil's form (27 November 2026)

**In short.** Weil's explicit formula turns the Riemann hypothesis of an L-function into the statement that a quadratic form is positive. Here the primes are given and the local data are unknown: one probability measure per prime, on the value χ(p) of a real character or on the Satake angle of a form on GL(2) or GL(3). The only other input is the archimedean factor (and, for odd forms, the root number). All findings are numerical (float64, primes below 40 at most); nothing is proved.

1. **Real characters.** The smallest conductor that keeps the form positive is 2.9999999996 (odd) and 4.9999999994 (even), and the optimum is the character mod 3, mod 5 at every prime of the support, with 0 at the ramified prime.
2. **Level 1, weight as a real number.** Nothing is feasible below weight 12; at 12 the optimum is τ(p)/p^(11/2) for p = 2, 3, 5, 7 to 3 decimals. With root number −1, weight 14 is excluded and the first feasible weight is 18.00, again with its Hecke eigenvalues.
3. **Maass forms of SL(2,ℤ), eigenvalue as a real number.** Nothing is feasible below the first eigenvalue of each parity. At support 16 six of the ten eigenvalues below 20 (9.53, 12.17, 13.78, 14.36, 17.74, 19.42) each lie in their own feasible piece, with no piece that holds no eigenvalue; five of them return a₂ and a₃ to 0.01 or better. The four odd eigenvalues above 15.9 lie in one stretch that stays feasible.
4. **SL(3,ℤ).** One of the five lowest forms, (13.5966, 4.7647), lies on the edge of the feasible set and returns its a₂ within 0.011; the others do not.

5. **Weight 2 with whole-number a_p.** The smallest admissible level is 10.995, with a₂, a₃, a₅, a₇ = −2, −1, 1, −2 (the elliptic curve of conductor 11); with root number −1 it is 36.91, with −2, −3, −2, −1 (the curve 37a).

Five rounds were registered before their runs (165 to 168 and 170). Their strong claims (islands at every eigenvalue, every low form on the edge) failed as registered; sections 5 to 8 give each verdict. In print: the positivity method with local corrections and bounds that come close (r > 7.2, conductor > 10.323, weight near 11). Not found in print by two searches: that the bound is exact when the local data are tied to one measure per prime and the test function is optimised, and that the optimum is the local data of the extremal form.

## 1. What is in print

- **The method.** Explicit-formula positivity bounds for discriminants and conductors: Odlyzko, Serre, Poitou 1977, Diaz y Diaz 1980 (minimum over the splitting types at each prime, one fixed kernel), Mestre 1986 (conductor of an elliptic curve > 10.323 with one test function), Miller 2002 (no cusp form on SL(2,ℤ) with r ≤ 7.2; regions for SL(3,ℤ)), Chenevier–Lannes and Chenevier–Taïbi (classification of algebraic forms of level 1 and small weight).
- **How close it comes.** Chenevier–Lannes (arXiv 1409.7616, §IX.3): without prime terms their criterion changes sign between motivic weights 9 and 11, "sa précision est assez étonnante". Mestre: 10.323 against 11. Diaz y Diaz did not optimise because the quadratic minima were known. Odlyzko 1990 leaves the optimal kernel open.
- **The unstructured problem** (any positive measure on the prime side): N. Johns, Zenodo, October 2026, and its companion note on critical conductors: the optimum is 1 for ζ and lies strictly below the true value in higher degree (2.981 against 3, 4.985 against 5), with optimal measures that are not arithmetic. Section 5e of [../bootstrap.md](../bootstrap.md) found the same numbers.
- **Not found** (prior-art reports of 26–27 November, both "partly in print"): exactness of the structured bound; the local data of the extremal form as the optimum; feasible sets that break into pieces at the automorphic spectrum. Sources that could not be opened: Serre's "Minorations de discriminants", Miller 1996, Bian 2010, the Chenevier–Lannes book beyond §IX.

## 2. Set-up (derived; checked numerically)

Test functions f live on [−a, a] in u = log x; X = e^(2a) is the support; g(w) = ∫ f(v + w) f̄(v) dv. For an L-function with Gamma factor ∏ᵢ Γ_R(s + κᵢ) (the κᵢ complex, closed under conjugation), conductor q, no pole, and −L′/L = Σ c(n) n^(−s):

`Q(f) = (log q) g(0) + Σᵢ A_κᵢ(g) − Σ_n n^(−1/2) [c(n) g(log n) + c̄(n) g(−log n)]`

`A_κ(g) = −(γ + log π) g(0) + ∫₀^∞ [2e^(−2u) g(0) − e^(−(1/2+κ)u) g(u) − e^(−(1/2+κ̄)u) g(−u)] / (1 − e^(−2u)) du`

- Derivation: Gauss's formula for ψ at (1/2 + κ + it)/2 integrated against |f̂(t)|², as in section 2 of [../bootstrap.md](../bootstrap.md); the two prime sums come from the two sides of the functional equation, so for a form that is not self-dual the matrix is Hermitian and the test functions are complex.
- For a true L-function, Q(f) = Σ_γ |f̂(γ)|², which is ≥ 0 under its Riemann hypothesis. A negative eigenvalue of the truncated matrix (basis e^(iπku/a), |k| ≤ N) proves that the form is not positive; a positive matrix proves nothing (Rayleigh–Ritz).
- **The unknown.** For each prime p ≤ X one probability measure μ_p on the local parameter; the coefficient at p^m is the mean of the local weight under μ_p, the same measure for every m and for every L-function in the problem:
  - real character: x ∈ [−1, 1], weight x^m log p (x = 0: ramified);
  - GL(2), angle θ ∈ [0, π]: L(f): 2 cos mθ; Sym²: 1 + 2 cos 2mθ; Sym³: 2 cos mθ + 2 cos 3mθ; Sym⁴: 1 + 2 cos 2mθ + 2 cos 4mθ;
  - GL(3), conjugacy class of SU(3) with eigenvalues e^(iφⱼ): L(F): τ_m = Σⱼ e^(imφⱼ); adjoint: |τ_m|² − 1.
- **The optimum.** t* = the largest t such that every form of the problem is ≥ 0 with conductor e^(−kt) (k = 1 for the standard L-function, k for Sym^k, 2 for the adjoint). q* = e^(−t*) is the smallest admissible conductor; at level 1 the test is the sign of t*. All constraints are linear in the measures, so t* is the value of a convex problem (Kelley cutting planes with an LP, HiGHS). By minimax it equals the conductor bound "with local corrections" optimised over all test functions against the worst local data that are consistent across the powers of a prime.
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

The optimal measures are point masses ([chars_out.txt](../../../../data/study/maass/chars_out.txt), support 40): odd: −1, 0, −1, +1, −1, +1, −1, +1, −1, −1, +1, +1 at p = 2, 3, …, 37, the character mod 3; even: −1, −1, 0, −1, +1, −1, −1, +1, −1, +1, +1 at p = 2, …, 31, the character mod 5 (the last prime of the support, 37, is at −0.94 against −1).

- This changes finding 2 of [../bootstrap.md](../bootstrap.md) ("without a pole positivity did not pin the points"): that search left the positions free; with the primes given and one measure per prime the Dirichlet case is pinned.
- It also locates the loss in the unstructured problem (2.981 and 4.985): there the weights at p and p² are independent. Tying them to one measure removes the gap.

## 4. Level 1 with the weight as a real number: 12, and Ramanujan's τ ([scan_hol.py](scan_hol.py), [levels.py](levels.py), [truth_hol.py](truth_hol.py))

Registered as the control of round 165. Four forms, support 8 ([hol_X8.jsonl](../../../../data/study/maass/hol_X8.jsonl)): t* rises from −2.390 at κ = 2 through −0.0366 at κ = 11.75 and −0.00012 at κ = 12 to +0.0229 at κ = 12.25, and stays positive up to κ = 30. So the smallest admissible weight is 12.00 within the grid.

| p | 2 | 3 | 5 | 7 |
|---|---|---|---|---|
| optimum at κ = 12, support 8 | −0.5297 | +0.5975 | +0.6969 | −0.4095 |
| optimum at κ = 12, support 16 | −0.530 | +0.598 | +0.691 | −0.376 |
| τ(p)/p^(11/2) | −0.5303 | +0.5987 | +0.6912 | −0.3765 |

- Weight 14 carries no cusp form but is feasible in the convex problem at support 8 (t* = 0.165): the measures there are not point masses.
- **Levels above 1.** The same curve read as a conductor: smallest admissible level 10.906, 4.963, 2.967, 1.951 at weights 2, 4, 6, 8 at support 8 and 10.936, 4.969, 2.983, 1.958 at support 16 ([levels_out.txt](../../../../data/study/maass/levels_out.txt), with Steinberg columns for primes dividing the level); the first newforms have levels 11, 5, 3, 2. Not exact at these supports: the optima are mixtures. The true newforms of levels 11, 5, 3, 2, 14, 15, 7 (coefficients from PARI) pass all four forms at their levels ([truth_hol_out.txt](../../../../data/study/maass/truth_hol_out.txt)).

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
- **The form (13.5966, 4.7647) is on the edge:** t* = +0.0016 at its exact parameters, an infeasible grid point at distance 0.035, and the optimum has a₂ = 1.040 − 0.382i against the true 1.048 − 0.375i (difference 0.011; a₃: 0.048).
- **The others are not** (G2 and G4 fail): the nearest infeasible point of the local grid lies at 0.185, 0.348, 0.176 from forms 1, 2, 4 (smallest Laplace eigenvalue first), and form 5 has none within its 9 × 9 grid; a₂ is returned within 0.081 and 0.058 for forms 1 and 4 and not at all for forms 2 and 5 (1.31, 0.79). FKL's region, which assumes |a_p| ≤ 3, leaves 0.13 to 0.52.
- G1 (truth feasible) passes: t* = +0.019, +0.167, +0.002, +0.030, +0.261 at the five forms.
- Both forms are active at 844 of the 861 coarse points. On the self-dual line λ₂ = 0 the feasible set begins at λ₁ = 16.5; the first form there is the symmetric square of the first Maass form of SL(2,ℤ), at 2 × 9.5337.

VERDICT of round 167 as registered: G1 pass, G2 fail, G3 pass, G4 fail. One of the five lowest forms is extremal for the convex problem at support 8 and returns its Hecke eigenvalue; the problem with genuine conjugacy classes was not run.

## 8. Round 168: the central zero ([round168.txt](../../../../data/study/preregister/round168.txt), [decide168.py](decide168.py))

The review of round 165 (section 9) found that the registered problem left out a piece of archimedean data, the root number. An odd Maass form has Λ(s) = −Λ(1 − s), so L(1/2, f) = 0 and, under its Riemann hypothesis, Q_f(φ) ≥ |∫φ|². Round 168 (registered before my runs, commit 90c7c6ea; the review's own numbers at support 8 declared) replaces Q_f by Q_f − |∫φ|² in the convex problem of round 165. Two forms, odd parity, r = 9.00 to 20.00 in steps of 0.01.

- **Truth passes** (H1; [truth_center_out.txt](../../../../data/study/maass/truth_center_out.txt)): the seven odd forms have smallest eigenvalue ≥ −2.6·10⁻¹⁴ with the central zero removed, at support 16.
- **Feasible set** ([decide168_X8_out.txt](../../../../data/study/maass/decide168_X8_out.txt), [decide168_X16_out.txt](../../../../data/study/maass/decide168_X16_out.txt)):

| support | feasible pieces (grid 0.01) | true odd eigenvalues |
|---|---|---|
| 8 | [9.53, 9.54], [12.17, 12.19], [14.28, 14.43], [15.82, 20.00] | 9.5337, 12.1730, 14.3585, |
| 16 | {9.54}, [12.17, 12.19], [14.35, 14.37], [15.93, 20.00] | 16.1381, 16.6443, 18.1809, 19.4847 |

- **The three lowest odd eigenvalues are islands, and each returns its Hecke eigenvalues** (support 16, nearest grid point): a₂, a₃ = −1.068, −0.455 (true −1.068, −0.456); +0.288, −1.201 (true +0.289, −1.202); −0.230, +0.694 (true −0.231, +0.696).
- **Above 15.9 the stretch remains** (H2, H3, H4 fail there, as the declared weak spot at support 8 and also at support 16): 384 of the 415 feasible grid points at support 16 are farther than 0.03 from every eigenvalue, and the Hecke eigenvalues of the four forms inside the stretch are not returned.

VERDICT of round 168 as registered: H1 pass; H2, H3, H4 fail (all failures above r = 15.9).

**Controls with exact answers** (registered before their runs, commit 2c3b6f97):
- *Level 1, weight as a real number, with the central zero* ([hol_X8_center.jsonl](../../../../data/study/maass/hol_X8_center.jsonl), [hol_X16_center.jsonl](../../../../data/study/maass/hol_X16_center.jsonl)). Cusp forms of level 1 with root number −1 exist at weights 18, 22, 26, … and not at 14. Found: t* = −0.470 at κ = 14, −0.026 at 17.75, +0.00001 at 18 and +0.026 at 18.25 (support 8); −0.027 at 17.75 and −0.0002 at 18 (support 16). So weight 14 is excluded and the smallest admissible weight is 18.00 within the grid: the prediction holds. The optimum at κ = 18: −1.4585, −0.3771, −1.1750 at p = 2, 3, 5 (support 8) and −1.4575, −0.3783, −1.1697, +0.2032 at p = 2, 3, 5, 7 (support 16), against the true −1.4584, −0.3770, −1.1745, +0.2115 (coefficients −528, −4284, −1025850, 3225992 from PARI).
- *Weight 2 with the central zero* ([levels_center_out.txt](../../../../data/study/maass/levels_center_out.txt)): smallest admissible level 32.46 at support 8 and 34.76 at support 16; the first newform of weight 2 with root number −1 has level 37. The prediction "below 37" holds; "within 3 % at support 16" fails (6.0 %).

**Even forms on the fine grid** (not registered; [scan_eps0_X16_k2_J61.jsonl](../../../../data/study/maass/scan_eps0_X16_k2_J61.jsonl), [scan_eps0_X32_k2_J61.jsonl](../../../../data/study/maass/scan_eps0_X32_k2_J61.jsonl)). Convex problem of round 165, r = 13.00 to 20.00 in steps of 0.01, support 16: the feasible set is {13.78}, [17.69, 17.94], [19.32, 19.43], one piece for each of the three even eigenvalues 13.7798, 17.7386, 19.4235, and nothing else (round 165 covers r < 13). At support 32 the pieces are {13.78}, [17.70, 17.88], [19.38, 19.43]. Hecke eigenvalues at support 16: at 13.78: 1.549, 0.247, 0.736 (true 1.549, 0.247, 0.737); at 19.42: −0.698, 1.555, −0.039 (true −0.693, 1.562, −0.038); at 17.74, inside its wider piece, they are not returned.

**Together, at support 16** (primes up to 13): of the ten eigenvalues below 20, six are isolated (9.53, 12.17, 13.78, 14.36, 17.74, 19.42), each in its own piece, with no piece that holds no eigenvalue; the four odd ones above 15.9 lie in one remaining stretch.

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

## 9. The review of round 165 ([review/review.md](review/review.md))

An adversarial agent recomputed the problem with its own code (Legendre basis, the form in x-space, an exact moment formulation without an angle grid, its own barrier solver).

- No bug found in [engine.py](engine.py) or [maass.py](maass.py). Its values of t* at ten points agree with ours within 7·10⁻⁴; at the first eigenvalue it finds t* = 1.7·10⁻⁹, so our −0.0003 there came from the angle grid.
- **Width of the first island at support 8** (its number, not recomputed by us): t* ≥ 0 exactly for r in [9.533695256, 9.533696035], width 7.8·10⁻⁷; the true value is 9.533695261. Four primes and the Gamma factor give the first eigenvalue to seven digits.
- **The root number** (section 8).
- **A lemma** (theirs; the one-line check is ours): Q_odd − Q_even = ½ ∫ |φ̂(t)|² Σ_ρ sech(π(t − ρ)) dt ≥ 0, with ρ = ±r for L(f). It follows from the reflection formula, Re[ψ(3/4 + ix) − ψ(1/4 + ix)] = π sech(2πx). So without the central zero the odd problem admits every even eigenvalue, which is the bump of the odd curve at 13.8.
- A larger basis (N = 24, 48, 96) does not change the sign of t* at r = 13, 15, 18 (0.1521, 0.1520, 0.1520, …): the feasible stretch of round 165 is not an artifact of a small basis.
- With the wrong parity the first odd form fails (−0.94); the first even form with odd Gamma factors does not (+1.9·10⁻⁴), as the lemma says it must not.
- Not checked by the review: the holomorphic control, supports above 8, anything beyond float64.

## 10. Size of the result

- **Known before:** the method, local corrections with a fixed kernel, and the remark that the bounds come close (section 1).
- **New here, as far as two searches found** (numerical; nothing proved):
  1. with the primes given and one measure per prime, the bound is exact to 9 digits for the two smallest real characters and to the grid for the first cusp form of level 1 of each sign (weights 12 and 18) and the first Maass form of each parity;
  2. the optimum is the local data of that form: the values of the character, τ(p), Hecke eigenvalues to 2–4 digits;
  3. for Maass forms of SL(2,ℤ) the feasible set at support 16 separates six of the ten eigenvalues below 20 into single pieces, once the parity and the root number are put in;
  4. on SL(3,ℤ) one of the five lowest forms is extremal and returns a₂.
- **What it is not.** Not a proof of anything about L-functions: feasibility of the true forms assumes their Riemann hypotheses (for L(f), Sym², the adjoint) and real Satake angles. Not exact above level 1 at the supports reached (10.94 against 11, 34.8 against 37). Not a set of islands for the odd forms above 15.9 at support 16. The registered claims in their strong forms failed: D2, D3, D4 (odd), D5 (odd), G2, G4, H2, H3, H4.
- **Why it may matter.** A bound that is exact, with the extremal object as its optimum, is the situation of the Cohn–Elkies bound in dimensions 8 and 24 and of Johns' statement for ζ. If it holds in general for "the first form of each archimedean type", it is a converse theorem in which positivity replaces the functional equations of twists; and the method is a way to compute first eigenvalues and Hecke eigenvalues from the primes and the Gamma factor alone.

## 11. Open points

1. A proof of exactness in the simplest case (odd real characters, conductor 3): a test function with non-negative transform that vanishes at the zeros of L(s, χ₋₃), whose local polynomial Σ_m g(m log p) p^(−m/2) x^m is minimal on [−1, 1] at x = χ₋₃(p) for every p.
2. Whether the levels 11 and 37 are reached exactly as the support grows (10.995 and 36.91 at support 32 with whole-number a_p), and the levels 5, 3, 2 at weights 4, 6, 8.
3. Round 166 is incomplete; round 169 (central zero and genuine angles together on the stretch above 15.9) was running when this was written.
4. Interval arithmetic for the infeasibility statements; a second code for sections 3, 4, 7 and 8 (the review covered round 165 at support 8 only).
5. SL(3,ℤ) with genuine conjugacy classes, larger support, and the forms of the key beyond the five lowest.
6. Sources not opened by the prior-art searches (section 1).

## 12. Validation

[validate.py](validate.py) reads the data files and checks every number quoted in sections 3 to 8c (literature values and the review's own numbers of section 9 excepted).
