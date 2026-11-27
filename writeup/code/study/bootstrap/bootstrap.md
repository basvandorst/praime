# The primes as the answer to a positivity problem (26 November 2026)

**In short.** Weil's explicit formula turns the Riemann hypothesis into the statement that a quadratic form is positive. We keep the Gamma factor and the pole, and make the prime side unknown: any set of real numbers x > 1 with the weights of an Euler product. Three findings, all numerical.

1. Started with no primes, the search returns 2, 3, 5, 7 to 15 digits or more, and 11 to 7 digits. Each prime, moved alone, is pinned more tightly as the support grows (2 + 10⁻³⁸ fails at support log 11).
2. The mechanism is the pole. Without a pole (three Dirichlet data tried), or when Euler factors with a minus sign are allowed, positivity did not pin the points.
3. With a pole there is a variational form: every prime system has a smallest conductor that keeps the form positive. The smallest value found over all positive systems is 1 for Γ_R(s); for Γ_R(s)Γ_R(s+1) and Γ_R(s)² it rises with the support and stays below 3 and 5, the smallest discriminants, and from support log 10 on the best systems found sit near the prime ideal norms of Q(√−3) and Q(√5). That these are the limits is a conjecture.

**Follow-up (27 November 2026): [maass/maass.md](maass/maass.md).** With the primes given and one measure per prime as the unknown, the bound is exact for the two smallest real characters (conductors 3 and 5, with the character values), for the first cusp forms of level 1 (weights 12 and 18, with their Hecke eigenvalues) and for the first Maass forms of SL(2,ℤ). This qualifies finding 2: the Dirichlet case is pinned once the primes are given.

Status of each statement: computed (best found by a search), not proved. The first step (one prime, one support) is in print: Connes and Consani 2021.

**Correction (26 November 2026, late): prior art.** Two preprints of October 2026 (N. Johns, Zenodo; produced by AI agents, not refereed) treat the same problem in its convex form: all pairs of positive measures that solve the explicit formula with the Gamma factor and pole of ζ. They contain the smallest-conductor reading, the statement that the optimum is exactly conductor 1, certified values for Q(√5) and Q(√−3) that agree with ours, and a claimed proof with an explicit extremal function. Our two prior-art searches missed them. Sections 5d and 5e are therefore not new; section 5f says what is theirs, what we checked ourselves, and what we did not find there.

## 1. What is in print

- Connes and Consani, *Spectral triples and ζ-cycles*, arXiv 2106.01715, section 2.3: at support log 3, with the one prime 2 replaced by a real p, positivity fails at p = 1.9999 and at p = 2.0005. Section 2.2: the archimedean part alone loses positivity just past support log 2. They move one prime at one support. Our matrix reproduces their numbers ([cc_check_out.txt](cc_check_out.txt)): smallest eigenvalue at support log 3: 5.68·10⁻⁸ (theirs: below 6·10⁻⁸); at support log 11: 2.93·10⁻⁴⁸ with 97 basis functions and still falling (theirs: 2.389·10⁻⁴⁸); negative at 1.9999 and at 2.0005.
- Discriminant bounds from the explicit formula (Odlyzko; Serre; Poitou, *Sur les petits discriminants*, 1977; Diaz y Diaz, *Tables minorant la racine n-ième du discriminant*, 1980). They fix one kernel F ≥ 0 and, as a "local correction", take for each small prime the splitting type with the smallest contribution. Printed bounds for |d|: imaginary quadratic 2.9657 without and 2.9921 with local corrections; real quadratic 4.9436 and 4.9640 (Diaz y Diaz, Tables 4 and 3; Odlyzko's GRH table: 2.962 and 4.951). So the minimum over splitting types, prime by prime at integer positions, is in print. Odlyzko 1990 lists the optimal kernel (Problem 2.2) and the share of the prime ideals (Problem 6.1) as open.
- N. Johns, *Positive solutions of the explicit formula for ζ(s): near-criticality and uniqueness* (Part I, Zenodo, doi 10.5281/zenodo.23197915, version 1.4, October 2026) and *Sharpness and uniqueness for positive solutions of the explicit formula for ζ(s)* (Part II, doi 10.5281/zenodo.23238100, version 1.3, October 2026). Found by the third search (26 November 2026, late). They pose the linear program over all pairs of positive measures (zero side, prime side) with ζ's archimedean data, prove the duality with the smallest conductor, and state the sharpness (optimum = conductor 1) and the uniqueness of ζ's pair, first as conjectures (Part I), then as theorems (Part II). Details and what we read: section 5f.
- Not found in print by the first two searches (about 60 pages opened; the sources in section 7 were not read), and not found in the two Johns preprints by a keyword search (they were not read in full): primes found from no input by a search on the prime side, the pinning as a function of the support, the role of the pole, the price of a finite edit, the Euler-product class with whole-number multiplicities and its link to the prime ideal norms of quadratic fields, cubic Gamma factors. In print there, and so not new here: the smallest conductor over all positive measures, its value 1 for ζ, and the fact that it is below the discriminant for Q(√5) and Q(√−3).

## 2. Set-up (derived; checked numerically)

Test functions f live on [−a, a] in the variable u = log x; L = 2a is the support and X = e^L. Put g(w) = ∫ f(v + w) f(v) dv. For a Gamma factor ∏ᵢ Γ_R(s + κᵢ) (κᵢ ∈ {0, 1}), conductor q, and a simple pole at s = 1 or none:

`Q(f) = [pole] 2∫ g(w) cosh(w/2) dw + (log q) g(0) + Σᵢ A_κᵢ(g) − Σ_x Σ_(k ≥ 1) ε_x^k · log x · x^(−k/2) · 2 g(k log x)`

`A_κ(g) = −(γ + log π) g(0) + ∫₀^∞ [2e^(−2u) g(0) − 2e^(−(1/2+κ)u) g(u)] / (1 − e^(−2u)) du`

The last sum runs over the points x > 1 of the unknown system, each with a sign ε_x = ±1 (ε = +1: Euler factor (1 − x^(−s))^(−1), a positive Beurling system).

- Derivation of A_κ: Gauss's formula ψ(z) = −γ + ∫₀^∞ (e^(−t) − e^(−zt))/(1 − e^(−t)) dt at z = (1/2 + κ)/2 + it/2, integrated against the Fourier transform of g, then t = 2u.
- For the true primes and the Riemann datum, Q(f) = Σ_γ |f̂(γ)|² (Guinand–Weil), which is ≥ 0 under RH. Check ([ef_check_out.txt](ef_check_out.txt)): with a Gaussian g the two sides agree to 27 digits.
- Basis: e^(iπku/a), |k| ≤ N. The matrix is real symmetric, with entries in closed form in terms of 4(N + 1) one-dimensional integrals ([weil.py](weil.py)). A negative eigenvalue of this finite matrix proves that the form is not positive (Rayleigh–Ritz); a positive matrix proves nothing. So "pinned within" below is an outer statement that holds up to rounding and quadrature error (no interval arithmetic), and "feasible" means "no negative eigenvalue in this basis". A larger basis can only lower every eigenvalue: for a fixed system the true smallest eigenvalue is below ours, and the true smallest conductor q*(P) is above ours.
- **Conductor shift (exact).** The conductor enters as (log q)·g(0) = (log q)·‖f‖². So λ_min(P; q) = λ_min(P; 1) + log q. A system P is feasible for exactly the conductors q ≥ q*(P, L) = exp(−λ_min(P; 1)).
- **A minus-sign point near 1 acts as a conductor (in the limit x → 1).** Exactly, (1 + x^(−s))^(−m) = 2^(−m)·x^(ms/2)·cosh(s·log x/2)^(−m); as x → 1 the last factor tends to 1, so m such points act as a conductor x^m.

## 3. Results for the Riemann datum (Γ_R(s), pole, q = 1, positive weights)

**3.1 Scan** ([scan.py](scan.py), [scan_out.txt](scan_out.txt); float64, N = 12). All multisets of one or two points on a grid of step 0.01 in (1.04, X], and of three points on a grid of step 0.03:

| support X | no points | feasible grid configurations |
|---|---|---|
| 3 | −0.416 | one point: only x = 2.00; two points: (2.00, y) with y ≥ 2.96 (weight zero at the edge) |
| 4 | −0.704 | only (2.00, 3.00); best single point −0.24; best three points −0.007 |

These are values at N = 12; in the reviewer's converged basis the empty system at X = 3 has −0.430. The signs do not change.

**3.2 Pinning** ([pins.py](pins.py), [pins_out.txt](pins_out.txt); high precision). One prime moved by ±10^(−k), the others at their true values: the largest shift exponent k at which the matrix still has a negative eigenvalue.

| support X | N | smallest eigenvalue | 2 | 3 | 5 | 7 |
|---|---|---|---|---|---|---|
| 3 | 12 | 8.5·10⁻⁸ | +4, −6 | | | |
| 4 | 14 | 1.1·10⁻¹² | +8, −11 | +7, −8 | | |
| 6 | 16 | 4.2·10⁻²² | +18, −20 | +16, −18 | +9, −10 | |
| 8 | 18 | 3.3·10⁻²⁹ | +25, −27 | +23, −25 | +18, −20 | +10, −11 |
| 11 | 30 | 2.2·10⁻⁴³ | +38, −42 | +37, −40 | +33, −34 | +27, −28 |

Reading: at support log 11, 2 + 10⁻³⁸ and 2 − 10⁻⁴² are already excluded. The eigenvalues still fall with N (section 1), so the true pinning is tighter.

**3.3 From no input** ([sequential.py](sequential.py), [sequential_out.txt](sequential_out.txt)). Start with the empty system. Enlarge the support. Add a point only when the form is negative without it (grid scan for the best position), then move all points to maximise the smallest eigenvalue (Nelder–Mead in float64, then a Gauss–Newton step in 90 digits).

| support X | points returned (distance to the nearest integer) |
|---|---|
| 2.5 | 2.0027 |
| 3.5 | 2 (1.2·10⁻⁶), 3 (2.7·10⁻⁴) |
| 4.5 | 2 (1.2·10⁻¹⁰), 3 (6.6·10⁻⁹) |
| 6 | 2 (3.8·10⁻¹⁷), 3 (3.7·10⁻¹⁵), 5 (3.2·10⁻⁹) |
| 9 | 2 (1.0·10⁻¹⁸), 3 (2.6·10⁻¹⁸), 5 (1.6·10⁻¹⁶), 7 (3.4·10⁻¹²) |
| 10.5 | 2 (3.0·10⁻¹⁷), 3 (8.8·10⁻¹⁷), 5 (1.4·10⁻¹⁶), 7 (4.0·10⁻¹⁵) |
| 13.5 | 2, 3, 5, 7 (below 5·10⁻⁹), 11 (9.7·10⁻⁸), 13.002, and four stray points at the edge |

The points are the best maximiser found, not certified feasible points: the smallest eigenvalue left at the returned points is +1.9·10⁻¹¹ at X = 4.5 and between −5·10⁻¹⁷ and +10⁻⁹⁰ from X = 6 to 10.5.

The prime powers 4, 8, 9 are not points: they come with 2 and 3. Limits: beyond support 12 the float64 stage is too coarse (the signal of a new prime near the edge is below 10⁻¹⁶), the run at X = 12 failed and X = 13.5 recovered only partly. A version with the scan in high precision is not written.

## 4. The pole is the mechanism (negative results)

- **Dirichlet data, no pole** ([datum.py](datum.py), [datum_out.txt](datum_out.txt), [characters_check_out.txt](characters_check_out.txt)). For (q, κ) = (4, 1), (3, 1), (5, 0) the same search, with signs free, returns systems such as one point 1.66 with sign − (support 9, smallest eigenvalue +0.05): nothing like the character. The true characters are feasible but sit on the edge (χ₋₄ at support 9: 7.3·10⁻⁹). A reason in the same direction: with no pole to cancel, minus-sign points can supply extra conductor (section 2; the points returned lie between 1.4 and 2.9, not near 1). And for a large conductor nothing is needed at all: for q above 9.305 (parity κ = 1) or 215.33 (κ = 0) the empty system is positive for every support, since log(q/π) + Re ψ((1/2 + κ)/2 + it/2) ≥ 0 for all t (the minimum is at t = 0). The three data tried lie below these thresholds, and there the empty system does fail.
- **Riemann datum with signs free** (same output file): at support 3 the pair (1.24 with sign −, 1.98 with sign +) is feasible with margin 0.18 (N = 12; 0.167 in the reviewer's converged basis). Positive weights are needed for section 3.

## 5. The smallest conductor (positive weights, pole)

By the conductor shift, q*(P, L) = exp(−λ_min(P; 1)) is the smallest conductor that admits P at support L, and q_min(L) = min over P is the smallest conductor that admits any positive system. q_min(L) can only grow with L. A number field with ζ_K of this Gamma factor is feasible at its discriminant under GRH, so q_min(L) ≤ the smallest discriminant.

Search ([conductor.py](conductor.py), [conductor_out.txt](conductor_out.txt)): greedy addition of points and joint maximisation, continued from one support to the next. "Found" is the value of one system in our truncated basis: a better optimiser can find a lower value, and a larger basis raises the value of a given system (the X = 6 optimum is 2.995103 in the reviewer's basis of 56 Legendre polynomials, against 2.994999 here). So the entries are estimates of q_min(L), not bounds.

| X | Γ_R(s): q found, points | Γ_R(s)Γ_R(s+1): q found, points | Γ_R(s)²: q found, points |
|---|---|---|---|
| 4 | 1.000000; 2, 3 | 2.976788; 3.335, 3.338, 3.918 | 4.965711; none |
| 6 | 1.000000; 2, 3, 5 | 2.994999; 3.024, 3.935 | 4.979974; 4.335, 4.486, 5.253 |
| 8 | 1.000000; 2, 3, 5, 7 | 2.999540; 3.002, 3.996, 6.831 | 4.991092; 4.078, 4.835 |
| 10 | | (search stuck: 3.0016) | 4.997454; 4.018, 4.963 |

True systems, and a search started 1 % away from them ([conductor_true.py](conductor_true.py), [conductor_true_out.txt](conductor_true_out.txt)):

| field | discriminant | q*(true norms) at X = 4, 6, 8, 12 | search from the perturbed true norms, X = 12 |
|---|---|---|---|
| Q(√−3) | 3 | 2.998013, 2.999988, 3.000000, 3.000000 | 3.000014 at 2.99994, 4.00013, 6.962, 7.038 |
| Q(i) | 4 | 3.960526, 3.999689, 3.999990, 4.000000 | 3.000071 at 2.99971, 4.00069, 6.901, 7.106 |
| Q(√5) | 5 | 4.965711, 4.999271, 4.999977, 5.000000 | 4.997769 at 4.016, 4.968, … |
| Q(√2) | 8 | 7.244072, 7.940275, 7.996105, 7.999969 | 5.000372 at 3.997, 5.006, 8.914, … |

Reading.
- Each true system needs a conductor that rises to its discriminant.
- The free search went to the smallest discriminant of its Gamma factor in the greedy runs and from the four perturbed starts of the table, also from the prime ideals of Q(i) or Q(√2) (two of the X = 12 runs ended slightly above the discriminant, i.e. worse than the known true system). Its points approach the prime ideal norms of that field: 3 (ramified), 4 (2 inert), 7 and 7 (split) for Q(√−3); 4 (2 inert) and 5 (ramified) for Q(√5).
- At a finite support the optimum is not arithmetic (2.9995 at X = 8 with 3.002, 3.996, 6.83) and approaches the field as the support grows. The values 3.000014 and 5.000372 lie above the discriminant: there the search is worse than the known true system.
- The rational primes are the case of degree 1: in the searches (supports up to log 8) the primes are the system of smallest conductor for Γ_R(s) with a pole, and that conductor is 1.

**Independent review with a stronger optimiser** ([review/adversary_report.md](review/adversary_report.md); own code, Legendre basis, 217 restarts of differential evolution plus Nelder–Mead, up to 7 free points):

| X | best q, Γ_R(s)Γ_R(s+1) | points | ours |
|---|---|---|---|
| 4 | 2.975898 | a triple point near 3.40–3.44 | 2.976788 |
| 6 | 2.995103 | 3.02325, 3.93619 | 2.994999 |
| 8 | 2.999174 | 3.00341, 3.99195, 7.256, 7.326, 7.348, 7.384 | 2.999540 |
| 10 | 2.999784 | 3.00086, 3.99811, 7.0141, 7.0150, 9.1138 | stuck at 3.0016 |

So our values at X = 4 and 8 were not the minima. The better values still rise with the support and stay below 3. The best systems at small supports contain repeated or clustered points; the approach to 3, 4, 7, 7 is visible only from X = 10 on. For the zeta datum the reviewer found no system with smallest eigenvalue above 10⁻⁶ at supports log 4 and log 6 (up to 8 points, 115 restarts).

**Conjecture (ours, numerical evidence only).** For Γ_R(s) with a simple pole, q_min(L) → 1 and every sequence of minimisers tends to the primes; for the two quadratic Gamma factors q_min(L) → 3 and 5 with minimisers tending to the prime ideal norms of Q(√−3) and Q(√5).

Why a minimiser should be arithmetic (heuristic). At L = ∞ a system is feasible at conductor q when its "zero measure" plus the constant density (log q)/2π is a positive measure. Lowering q until this fails leaves a measure with no slack between its atoms: a pure point measure, which is what a zeta function with all zeros on the line has. This is the shape of an extremal solution in a linear program, but the set of prime systems is not convex, so it is not a proof.

**Degree 3, one run** ([conductor_cubic_out.txt](conductor_cubic_out.txt)). Gamma factor Γ_R(s)²Γ_R(s+1) (one real, one complex place; the smallest discriminant is 23, the field of x³ − x − 1, whose first prime ideal norms are 5, 7, 8).

| X | 4 | 6 | 8 | 10 | 12 | 14 |
|---|---|---|---|---|---|---|
| q found | 21.0020 | 22.5285 | 22.5525 | 22.5935 | 22.5983 | 22.6300 |

At X = 4 and 6 the best system is empty; at X = 14 the points are 5.940, 5.940, 7.334. The value rises, stays below 23, and from X = 6 on lies above the printed bound without local corrections (22.27). The points are not yet near 5, 7, 8. So degree 3 is consistent with the limit 23 and does not show it.

## 5b. The price of an edit (proved under RH, for L = ∞; refereed once, see [theory/referee_report.md](theory/referee_report.md))

Let P come from the primes by removing a finite set R of primes and adding a finite multiset S of real points > 1. Put

`h_P(t) = Σ_(p ∈ R) c(p, t) − Σ_(x ∈ S) c(x, t)`,  `c(x, t) = Re[ log x / (x^(1/2 + it) − 1) ]`.

**Statement.** Under RH, the infimum of Q_P(f)/‖f‖² over all compactly supported f equals 2·inf_t h_P(t). So the smallest conductor that admits P is `q*(P) = exp(−2·inf_t h_P(t))`.

- Removing one prime p: `q* = p^(2/(√p + 1))` (1.776 for p = 2, 2.235 for p = 3, 2.908 for p = 7).
- Adding one point x: `q* = x^(2/(√x − 1))` (11.85 for x = 6).
- If the logarithms of the points involved are linearly independent over Q, the prices add: `log q* = Σ_R 2 log p/(√p + 1) + Σ_S 2 log x/(√x − 1)`. With dependent logarithms the price can be lower than the sum.
- Consequences.
  - Every finite edit needs a conductor above 1, so among finite edits the primes are the strict minimum. Reason: h_P is an almost periodic function with mean 0; if inf h_P were 0, then h_P ≡ 0, the two prime sides would agree, and the edit would be empty.
  - Moving a prime p to x costs the price of removing p plus the price of adding x whenever log x/log p is irrational (48.3 for 2 → 2.1), however small the move. For a rational ratio m/n with m odd the cost is a little lower (2 → 2^(11/10): 2 inf h = −3.8577 against the sum −3.8601; referee's computation). Either way, at L = ∞ the cost does not tend to 1 as the move shrinks. At a finite support it does, which is the pinning of section 3.2.
  - Edits far out are cheap: both prices tend to 1 as p or x grows.
- RH enters in steps (2) and (5) of the proof (the zero sum runs over real ordinates).

Proof. (1) Q_P(f) = Q_ζ(f) + (C_ζ − C_P)(g), where C is the prime-side sum of section 2. (2) Under RH, Q_ζ(f) = Σ_γ |f̂(γ)|² ≥ 0. (3) For a point x, Σ_k log x·x^(−k/2)·2g(k log x) = (1/π)∫ |f̂(t)|² c(x, t) dt, because 2g(u) = (1/π)∫ |f̂(t)|² cos(tu) dt and Σ_k x^(−k(1/2 + it)) = 1/(x^(1/2 + it) − 1). So (C_ζ − C_P)(g) = (1/π)∫ |f̂|² h_P dt. (4) ‖f‖² = (1/2π)∫ |f̂|² dt. Hence Q_P(f)/‖f‖² ≥ 2 inf h_P. (5) Equality in the limit: h_P is continuous; choose t₀ that is not the ordinate of a zero with h_P(t₀) < inf h_P + ε, and f(u) = φ(u/T) cos(t₀u) with φ smooth of compact support; as T → ∞, |f̂|² concentrates at ±t₀ and Σ_γ |f̂(γ)|²/‖f‖² → 0. (6) For one point, Re[1/(√x·e^(iθ) − 1)] has minimum −1/(√x + 1) at θ = π and maximum 1/(√x − 1) at θ = 0; with independent logarithms the phases are jointly dense (Kronecker), so the extremes add.

Check at finite supports ([edit_cost.py](edit_cost.py), [edit_cost_out.txt](edit_cost_out.txt)). The smallest eigenvalue must fall towards 2 inf h as the support grows:

| edit | λ_min at X = 8, 16, 32, 64 | limit 2 inf h | q* |
|---|---|---|---|
| remove 2 | −0.5409, −0.5563, −0.5629, −0.5664 | −0.5742 | 1.776 |
| remove 3 | −0.6343, −0.7310, −0.7452, −0.7643 | −0.8042 | 2.235 |
| remove 7 | −0.4030, −0.7354, −0.7355, −0.8811 | −1.0675 | 2.908 |
| add 6 | −0.6177, −0.7315, −0.7315, −1.1923 | −2.4723 | 11.849 |
| move 2 to 2.1 | −0.5879, −1.0625, −1.3723, −1.8119 | −3.8780 | 48.330 |

All values lie above their limit and fall with the support. Only "remove 2" is close to its limit at X = 64; the others need the resolution of a larger support.

What this does not settle: systems that differ from the primes in infinitely many points (the global conjecture).

## 5c. What can be proved (proofs in [theory/proofs.txt](theory/proofs.txt); refereed once: no failing step, corrections listed at the end of that file)

Call a system P *feasible* if its form (Riemann datum, q = 1) is positive for every support. By Bochner–Schwartz the form is then the Fourier transform of a positive measure μ_P, the spectrum of P (for the primes under RH: the zeta zeros).

- **Theorem A.** If P and Q are feasible and μ_P ≥ μ_Q, then P = Q. So a feasible system is determined by its spectrum, and no feasible system other than the primes has all zeta zeros in its spectrum. (Idea: the two prime sides agree near 0, so the positive measure μ_P − μ_Q has a Fourier transform that vanishes near 0, hence is 0.)
- **Theorem C (RH not assumed).** If P differs from the primes by a thin change and P is feasible, then RH holds and P is the set of primes. Thin: the change ν of the prime side is Φ′ + ρ with ∫ over [−T, T] of |Φ| + |ρ| equal to o(T). Examples: a changed set S with Σ_(x ∈ S, x ≤ X) log x/√x = o(log X), in particular any finite edit; moving every prime p by o(√p) (all moved points staying above 1).
- **Proposition.** A feasible system satisfies ψ_P(x) ~ x and its zeta function (s − 1)ζ_P(s) has no zero in Re s > 1/2. The prime ideal norms of a number field other than Q are not feasible for the Riemann datum.
- **Open.** Whether the primes are the only feasible system. Another one would have to drop at least one zeta zero from its spectrum and differ from the primes in about √X/log X points up to X. Uniqueness does not follow from RH and does not imply it; it is a converse theorem with positivity in place of the functional equation (nearest known: Hamburger 1921 for ordinary Dirichlet series).

Two finite-level facts from the same work (their code: Legendre basis, [theory/wform.py](theory/wform.py)), checked with our Fourier matrix ([empty_check.py](empty_check.py), [empty_check_out.txt](empty_check_out.txt)):

- The empty system is positive up to support X = 2.1027 (theirs). Ours: 2.1151, 2.1092, 2.1071 with N = 20, 40, 60, falling towards it. So every feasible system has a point below 2.11: positivity forces a first prime.
- The primes up to p_n without p_(n+1) fail just past p_(n+1): at 1.05·p_(n+1) the smallest eigenvalue is −0.012, −0.050, −0.086, −0.122, −0.109 (±0.001) for p_(n+1) = 3, 5, 7, 11, 13. So positivity forces the next prime, each time.

## 5d. The global statement: where the proof stands (26 November 2026, evening)

Statement wanted: the primes are the only positive Beurling system whose form (Riemann datum, q = 1) is positive for every support. Not proved by us. A proof is claimed in Johns, Part II (Theorem 2, for all positive measures; not refereed; see section 5f). Item 4 below (integer weights) is Theorem 3.11 of Johns, Part I, by the same route. What we have:

- **Convex version, computed** ([relaxed.py](relaxed.py), [relaxed_out.txt](relaxed_out.txt)). If the multiplicities may be any real numbers ≥ 0 (weights on a grid of step 0.01), the form is linear in the weights and the largest smallest eigenvalue is a concave problem; cutting planes with a linear program give the global optimum on the grid: 7.0·10⁻⁴, 3.2·10⁻⁵, 1.8·10⁻⁵, 1.6·10⁻⁶ at X = 2.5, 3, 3.5, 4 (the primes alone: 5.7·10⁻⁸ at X = 3). The weights sit at 2.00–2.01 (total 1.003 at X = 3, 1.0001 at X = 4) and 3.00–3.04 (total 1.197 at X = 3.5, 1.005 at X = 4). So at a finite support the optimum is slightly off the primes and its conductor is slightly below 1 (0.9999984 at X = 4); the values fall with the support. Uniqueness, if true, is a statement about the limit.
- **Proved so far** (full proofs in [theory/global/global.md](theory/global/global.md); refereed once inside the repo, [report](theory/global/referee_report.md): no failing step, corrections listed in section 9 of the proof file; none of these uses that the multiplicities are integers).
  1. *Positivity as a half-plane condition.* P is feasible ⟺ Ξ_P′/Ξ_P continues to Re s > 1/2 with non-negative real part there, Ξ_P = ½·s(s − 1)·π^(−s/2)·Γ(s/2)·ζ_P(s) (for ζ: Lagarias's criterion). Then Ξ_P has no zero in Re s > 1/2 and |Ξ_P(σ + it)| grows with σ.
  2. *Lemma G (one-dimensional sphere-packing argument).* Let τ ≥ 0 be an even measure with no mass in (−1, 1). If its Fourier transform equals α·δ₀ + constant on (−1, 1), then τ = α·(Dirac comb − δ₀). The test function is −sin(2π|x|)/(2π) on [−1, 1], whose transform sin²(πξ)/(π²(ξ² − 1)) is ≥ 0 for |ξ| ≥ 1 and vanishes there exactly on the integers. Consequence: a Beurling zeta function with Riemann's functional equation (either sign) is ζ. This may be Ryavec's theorem (Rocky Mountain J. Math. 5 (1975) 623–627, not read); Beurling 1951 (Ark. Mat. 1, not read) characterises ζ by growth without a functional equation.
  3. *A converse theorem without positivity (Theorem E).* If ζ_P converges in some half-plane, (s − 1)ζ_P(s) is entire of finite order and |Λ_P(1 − s)| ≤ C(1 + |s|)^N·|Λ_P(s)| on a half-plane Re s ≥ c′ with 1 < c′ < 3 (Λ_P = π^(−s/2)Γ(s/2)ζ_P), then ζ_P = ζ. No functional equation is assumed: a one-sided growth comparison is enough.
  4. *Corollary F.* If P is feasible, Ξ_P is entire of finite order and only finitely many of its zeros lie off the critical line, then P is the set of primes. In particular: feasible, entire, all zeros on the line ⟹ primes.
  5. For a feasible P: P = primes ⟺ ζ_P has Riemann's functional equation ⟺ the spectrum μ_P is a sum of point masses with integer weights.
- **The gap, exactly.** A second feasible system would have to be of one of three kinds: (G1) its zeta function has no continuation of finite order across the line (this includes a continuous spectrum and point masses of non-integer weight); (G2) Ξ_P is entire of finite order and Λ_P(1 − s)/Λ_P(s) is not polynomially bounded on any right half-plane (heuristically: infinitely many zeros to the left of the line, with positive density at bounded distance); (G3) Ξ_P is meromorphic of finite order with poles (trivial zeros missing, or poles of ζ_P left of the line, shape ζ(s)/(1 − 2^(−s))). Equivalent question on the additive side: does positivity force the Fourier transform of the measure on the Beurling integers to be α·δ₀ + constant on (−1, 1)?
- **Second attack on the gap** ([theory/global/gap.md](theory/global/gap.md); not refereed). Nothing closed. What it shows:
  - *The additive tools cannot close G2.* The function ζ(s)(1 + λL^(−s)), L > 1, 0 < λ < √L, has everything Lemma G and the converse theorem use (positive coefficients, coefficient 1 at 1, no mass in (−1, 1), entire completed function of order 1) together with an infinite row of zeros off the line. Only two things exclude it: Weil positivity fails for it, and it has no positive Euler product. So a proof must use positivity and the Euler product together. Lemma G is sharp: on any shorter interval (−a, a), a < 1, the comb is not the only solution.
  - *What positivity means for an entire Ξ_P:* the argument of Λ_P(1 − s)/Λ_P(s) is non-increasing along the critical line, and the continuous part of the spectrum has mean density lim (N − N_on)/T (all zeros minus zeros on the line). In every known example of this kind the argument is bounded, the density is almost periodic with mean 0, hence 0, and the system is the primes.
  - *Two steps would close G2, both open:* (a) positivity and the Euler product forbid a positive proportion of zeros off the line; (b) a continuous spectral density that is ≥ 0 with mean 0 vanishes (proved only when it is almost periodic).
  - *G3, a partial theorem:* with finitely many missing trivial zeros and finitely many zeros off the line, sign conditions and a narrow window for one parameter remain (for one missing zero at −2: [0.8, 0.811)); a linear-programming certificate that would close the cases of one missing zero at −2 or −4 exists only up to errors of 10⁻⁴ to 10⁻³ and is not certified.
  - *G1 (no continuation): nothing moved.* Positivity gives a function of Herglotz type, and such functions can have a singular continuous spectrum.
- **Counterexample search** (the prover's linear program over absolutely continuous prime densities capped at 50, heights up to 60 on a grid): the best smallest value of Re Ξ′/Ξ on the line stayed negative (−1.1·10⁻³ to −3.7·10⁻⁴), creeping to 0 as the family grows, with the mass concentrating in peaks. No second feasible system was found.
- **Literature** (one search; main sources on positivity criteria not read): nothing found that derives a functional equation from positivity; no discrete Beurling system with a zeta function meromorphic in the plane, beyond the arithmetic ones.

## 5e. The convex version across signatures: sharp in degree 1 only (26 November 2026, night)

Let the multiplicities be any real numbers ≥ 0 (a positive prime measure). The form is linear in the measure, so the smallest conductor that admits any such measure at support X is the optimum of a concave problem; cutting planes with a linear program give the global optimum on a grid ([relaxed_fields.py](relaxed_fields.py), [relaxed_fields_out.txt](relaxed_fields_out.txt); geometric grid of ratio 1.004 plus the integers, Fourier basis, float64). Each value is a lower bound, under GRH, for the discriminant of every field of that signature, up to the grid and the rounding (not certified).

| Gamma factor (smallest discriminant) | X = 4 | 6 | 8 | 10 | 12 | 14 | where the weight sits at X = 14 |
|---|---|---|---|---|---|---|---|
| Γ_R(s) (1) | 0.9999984 | 1.0000000 | 1.0000000 | | | | at X = 8: 2.000, 3.000, 5.000, 7.000, each with weight 1.00 |
| Γ_R(s)Γ_R(s+1) (3) | 2.97588 | 2.98099 | 2.98117 | 2.98123 | 2.98124 | 2.98126 | 3.244 (1.65), 6.259 (1.92), 8.862 (0.65) |
| Γ_R(s)² (5) | 4.96578 | 4.97863 | 4.98386 | 4.98476 | 4.98480 | 4.98485 | 4.320 (1.83), 9.224 (2.30) |
| Γ_R(s)²Γ_R(s+1) (23) | | 22.5288 | 22.5523 | 22.5664 | 22.5852 | 22.5953 | 6.206 (2.82) |
| Γ_R(s)³ (49) | | 46.5468 | 48.0437 | 48.1272 | 48.1324 | 48.1486 | 8.931 (4.21) |

For Γ_R(s) the runs at X ≥ 5 stopped at the round limit with the optimum equal to 1 within 5·10⁻⁸ (X = 3: 0.9999677).

Reading.
- **Degree 1 is sharp.** Without asking for integers, the convex program returns the primes with weight 1.00 and the conductor 1. So in the convex version too, nothing beats the primes at these supports; integrality of the multiplicities comes out, it is not put in.
- **Degree 2 and 3 are not sharp.** The optimum levels off below the smallest discriminant (2.9813 against 3; 4.9848 against 5), with non-integer weights at non-integer positions. The statement of section 5 for quadratic fields therefore depends on whole-number multiplicities: with them the searches head for 3 and 5, without them a non-arithmetic measure does better.
- **Order of the bounds for imaginary quadratic fields:** 2.966 (one positive kernel, no local correction, in print) < 2.9813 (all kernels of support log 14 against all positive measures, here) < 2.9921 (local corrections over arithmetic splitting types, in print) < 3.
- **The extremal function.** At the optimum the ground state ψ of the form has ĝ = |ψ̂|² ≥ 0 and c_g(x) = Σ_k log x·x^(−k/2)·2g(k log x) ≥ 0 for every x, with equality on the support of the optimal measure (the first-order condition of the concave problem, for a simple smallest eigenvalue). For the zeta datum this is the ground state studied by Connes and Consani, up to the small shift of the points at a finite support: it is the dual solution of a convex program whose primal solution is the primes. Global uniqueness in the convex version, together with "no positive measure is admitted below conductor 1", is equivalent to: this optimum tends to 0 as the support grows. Measured: 3.2·10⁻⁵, 1.6·10⁻⁶, below 5·10⁻⁸ at X = 3, 4, 5.

## 5f. The convex optimum tends to 0: in print since October 2026, and our own check (26 November 2026, late)

Target: the optimum λ*(L) of the convex version (section 5e; conductor e^(−λ*)) tends to 0 as the support L grows. We did not prove it. It is claimed in print, and we tested the claimed extremal function with our own code.

**Duality (proved; [theory/convex/convex_limit.md](theory/convex/convex_limit.md), section 2; not refereed).** λ*(L) is the smallest value of W(g)/g(0) over functions g on [−L, L] with non-negative Fourier transform and c_g ≥ 0, where W is the pole and Gamma part of the form (section 2 with no prime terms, q = 1). λ*(L) does not increase with L, and λ*(L) ≥ 0 under RH. The limit of λ*(L) equals the largest log(1/q) over all systems that are feasible at conductor q, and equals 2π times the optimum κ*_OPS of Johns' program in his classical cone (g ≥ 0 everywhere; his main cone asks g ≥ 0 only from log 2 on). So the target says: there are such g with W(g)/g(0) as small as we like.

**The prime-power structure is not needed** (computed; [relaxed_k1.py](relaxed_k1.py), [relaxed_k1_out.txt](relaxed_k1_out.txt)). With any positive measure on the prime side (only the k = 1 term, so c_g ≥ 0 becomes g ≥ 0), the program returns the same optimum, 3.2·10⁻⁵, 1.6·10⁻⁶, 1.0·10⁻⁷ at X = 3, 4, 5, and the weights of the explicit formula: 1.000 at 2 and 3, and 0.50 at 4 (X = 5: 0.507; X = 6: 0.500). The dual solution at X = 7.5 ([relaxed_dual.py](relaxed_dual.py), [relaxed_dual_out.txt](relaxed_dual_out.txt); float64, not converged below 10⁻⁸) falls from 1 at x = 1 to 0 at x = 2, with bumps of height 1.6·10⁻³, 4.5·10⁻⁵, 2.5·10⁻⁶ between 2–3, 3–4, 4–5, and touches 0 at 2, 3, 4, 5, 7.

**What is in the two Johns preprints** (read by us: Part I, abstract, section 1, Theorem 3.11 with the idea of its proof, Propositions 5.5 and 5.11; Part II, abstract, section 1, Lemma 3.2, Proposition 3.3, Theorems 6.1 and 6.2, Corollary 8.2. Not read: the proofs in sections 4 to 7 of Part II and most proofs of Part I. Their F is our ĝ/2π, their F̂(log x/2π) is our g(log x), their κ* is λ*/2π).

| their statement | status given there | relation to this note |
|---|---|---|
| Part I, Theorem A: positive pairs exist ⟺ κ* ≥ 0; smallest conductor q_min = e^(−2πκ*) | proved | the conductor shift and the convex version of sections 5d, 5e |
| Part I, Conjecture S: κ* ≤ 0 (q_min ≥ 1) | conjecture in Part I | the target of this section |
| Part I, Theorem C: κ* ≤ 1.43·10⁻¹⁰⁶⁰, certificate Ξ(t)²·P(t²) with zeros at the first 111 prime powers | computer-assisted, unconditional | our measured 3.2·10⁻⁵ … below 5·10⁻⁸ are values at short supports |
| Part I, Proposition 5.11: Q(√5): 1.64·10⁻⁵ ≤ κ* ≤ 4.83·10⁻⁴; Q(√−3): 1.50·10⁻⁴ ≤ κ* ≤ 9.99·10⁻⁴ (critical conductors in [4.9848, 4.9995] and [2.9812, 2.9972]) | computer-assisted | section 5e: 4.98485 and 2.98126 at support log 14, rising; "degree 2 is not sharp" is theirs, certified |
| Part I, Theorem 3.11: zero measure atomic with integer weights ⟹ q ≥ 1, and at q = 1 the pair is ζ's (exponential lift, Fejér kernel, Hamburger) | proved | Lemma G, Theorem E, Corollary F of section 5d |
| Part II, Theorem 1: an explicit F = Ξ²H with F ≥ 0, F̂ ≥ 0, F̂ = 0 at log n/2π for every integer n ≥ 2, A(F) = 0; hence κ* ≤ 0 | proof given, partly computer-assisted; "not yet checked by a human expert" | proves the target (with a truncation, below) |
| Part II, Theorem 2: every positive pair is ζ's pair | follows from Theorem 1 and a support theorem of Part I | the global statement of section 5d, in the larger class |
| Part II, Corollary 3 and 8.2: RH ⟺ the optimum is exactly 0 | same | λ*(L) ≥ 0 needs RH; λ* → 0 from above is equivalent to RH plus sharpness |

The extremal function is explicit (Part II, Theorem 4). With S(m, n; c) the Kloosterman sum and a dot for the derivative with respect to the order:

`a_n = n·Σ_c (1/c)·[S(1,n;c)·J̇₂(4π√n/c) − S(−1,n;c)·İ₂(4π√n/c)] + δ_(n,1)/4`

`G(x) = C·sin²(πx)·√x·Σ_n a_n·K₀(4π√(nx)) (x > 1), G(1) = C/(32π²), C = 880.6617`

`g(log x) = Σ_m d(m)·m^(−1/2)·G(mx) (x ≥ 1; d = number of divisors)`

**Our check of this function** (numerical, float64, own code: [magic_check.py](magic_check.py), [magic_check_out.txt](magic_check_out.txt)). We computed the a_n from the formula (Kloosterman sums to c = 400 for n ≤ 5), built g, and applied our own functional W of section 2.

| test | result |
|---|---|
| a₁ … a₅ against the printed table | agree to 5·10⁻⁹ relative (4581.370532 against 4581.370553) |
| W(g)/g(0) | pole part 1.120099144, Gamma part −1.120099143, sum 1.8·10⁻⁹ |
| cosine transform of g on 0 ≤ t ≤ 60 (step 0.1), relative to its value at 0 | smallest value −9·10⁻¹² (rounding) |
| the same at the first ten zeta zeros | at most 2.4·10⁻¹⁰ in size (14.1347: 2.4·10⁻¹⁰; one unit to either side: 1.6·10⁻³ and 5.7·10⁻⁴) |
| g ≥ 0, zero at every integer ≥ 2 | by the formula, given a_n > 0 (their Theorem 6.1, proof read: the term c = 1 dominates) |
| against the dual solution of our own program (support log 7.5) | g(log x)/g(0) at x = 1.02, 1.33, 1.50, 2.50: 0.99637, 0.46193, 0.18325, 1.561·10⁻³ against 0.9964, 0.4603, 0.1821, 1.586·10⁻³ |

So the printed function behaves as claimed to the accuracy of the test (about 10⁻⁹): it is doubly positive, its transform vanishes at the zeta zeros, and its cost is 0. This is a numerical test of the final formula, not a check of the proof. One difference: their function vanishes at every integer ≥ 2, the dual solution of our finite program does not vanish at 6 (2.3·10⁻⁷ there); the extremal function at a finite support is not unique.

**From the extremal function to compact supports** (numerical; same script). Multiply g by the window w_L = (v∗v)/(v∗v)(0), v(u) = (1 − (2u/L)²)₊. The window is ≥ 0 with transform ≥ 0, so g·w_L is doubly positive with support [−L, L], and W(g·w_L)/g(0) is an upper bound for λ*(L):

| X = e^L | 3 | 4 | 5 | 8 | 20 | 100 | 10⁴ | 10⁸ |
|---|---|---|---|---|---|---|---|---|
| W(g·w_L)/g(0) | 1.9·10⁻² | 1.1·10⁻² | 7.6·10⁻³ | 4.2·10⁻³ | 1.8·10⁻³ | 6.8·10⁻⁴ | 1.5·10⁻⁴ | 3.6·10⁻⁵ |

The bound falls like 1/L² (this window is crude; the true optimum at X = 3, 4, 5 is 3.2·10⁻⁵, 1.6·10⁻⁶, 1.0·10⁻⁷).

**Truncation lemma (proved; [theory/convex/convex_limit.md](theory/convex/convex_limit.md), section 3; the computation below rechecked by hand).** Let g ≥ 0 be even and continuous with ĝ ≥ 0 and W(g) finite, and let k(x) = 1 − (3/2)x² + (3/4)|x|³ for |x| ≤ 1, (2 − |x|)³/4 for 1 ≤ |x| ≤ 2, 0 beyond (the cubic B-spline; its transform is a fourth power of sin(t/2)/(t/2), so it is ≥ 0). Then g_M(u) = g(u)·k(2u/M) is doubly positive with support [−M, M], g_M(0) = g(0), and

`W(g_M) ≤ W(g) + (12/M²)·∫₀^∞ u²·e^(−u/2)·g(u)/(1 − e^(−2u)) du ≤ W(g) + 12·c₁·g(0)/M², c₁ = 7ζ(3) + π³/4 = 16.166, 12·c₁ = 193.99`

Proof in three lines: the window leaves g(0) unchanged, so only two terms of W move; the pole term 4∫g·cosh(u/2) can only fall, because g ≥ 0 and k ≤ 1; the Gamma term rises by 2∫e^(−u/2)·g·(1 − k)/(1 − e^(−2u)), and 1 − k(x) ≤ (3/2)x², g ≤ g(0). The constant is ∫₀^∞ u²e^(−u/2)/(1 − e^(−2u)) du = Σ_j 16/(4j + 1)³ ([riemann_kernel_cost_out.txt](riemann_kernel_cost_out.txt)). Consequences:

1. If Theorem 1 of Part II holds, then λ*(M) ≤ 194/M² for every support [−M, M], so λ*(M) → 0. (λ* ≥ 0 needs RH.) This step uses only that their function is doubly positive with cost 0; it does not use their uniqueness theorem.
2. From Part I, Theorem C alone (computer-assisted there, not rechecked by us): λ*(M) ≤ 2π·1.43·10⁻¹⁰⁶⁰ + 194/M², so lim λ*(M) ≤ 8.98·10⁻¹⁰⁶⁰.
3. With nothing from the preprints but the idea (Ξ² kills the zero side): the function with transform Ξ(t)² is doubly positive, because Riemann's kernel Φ is positive, and its cost is 0.0751971 ([riemann_kernel_cost.py](riemann_kernel_cost.py); two formulas agree, numerical). So λ*(M) ≤ 0.0752 + 194/M²: no positive system is feasible below conductor e^(−0.0752) = 0.9276. This is unconditional up to a certified evaluation of one integral.

**Referee read of Part II** (by an agent, bounded; [theory/convex/convex_limit.md](theory/convex/convex_limit.md), section 4). Sections 3 and 8 (the reduction to a function of the Gamma factor alone, and the assembly of Theorem 1) were checked step by step: no failing step. Theorem 6.1 (every a_n > 0) was checked line by line. An independent numerical rebuild of F from the formula for x > 1 gave F ≥ 0 on [0, 150], zeros at the first four zeta zeros, and cost −7·10⁻¹¹ (relative). Not checked by anyone here: section 4 (the automorphic input), Theorem 5.4(b), section 7 (positivity of H at large height), and the four computer-assisted inputs.

**Our own certificates: the best polynomial at each degree** (computed by an agent, [theory/convex/certs/README.md](theory/convex/certs/README.md); every cost recomputed with separate code, [cert_check.py](cert_check.py), [cert_check_out.txt](cert_check_out.txt)). The family is ĝ(t) = Ξ(t)²·B(t), B an even polynomial with B ≥ 0. Then g = B(i·d/du)(Φ∗Φ), the zero side vanishes, and W(g) = 2·Σ Λ(n)·n^(−1/2)·g(log n). Johns, Part I uses the member with double zeros exactly at the logs of the first J prime powers (degree 4J). Here B is optimised at each degree, in two classes: A (g ≥ 0) and E (only the Euler sums Σ_k x^(−k/2)·g(k·log x) ≥ 0).

- Parity rule (proved, two lines): for large u, g(u) = B(2πi·e^u)·(Φ∗Φ)(u)·(1 + o(1)), so g ≥ 0 needs (−1)^d·b_d ≥ 0 for the top coefficient b_d of t^(2d), and B ≥ 0 needs b_d ≥ 0. For odd d this forces b_d = 0: degree 4k + 2 gives the same optimum as degree 4k. (The divisor-weighted closed form of this family and the sign (−1)^d of its leading term are Lemma 4.6 of Part I.)
- Optimum of W(g)/g(0) at each degree (numerical optimum of a convex problem: multipliers positive, stationarity residual below 10⁻⁴⁰; not interval-certified), and the certified bound (next item):

| degree of B | class A, optimum | class A, certified bound | class E, optimum |
|---|---|---|---|
| 0 | 7.5197·10⁻² | | 7.5197·10⁻² |
| 4 | 6.1940·10⁻³ | 6.19402005455·10⁻³ | 6.1424·10⁻³ |
| 8 | 3.7072·10⁻⁴ | 3.70719233613·10⁻⁴ | 7.5310·10⁻⁵ |
| 12 | 4.6462·10⁻⁵ | 4.64624013887·10⁻⁵ | 1.2559·10⁻⁷ |
| 16 | 3.5892·10⁻⁶ | 3.58920272265·10⁻⁶ | 8.9013·10⁻¹⁴ |
| 20 | 4.4383·10⁻⁸ | 4.43833859341·10⁻⁸ | 5.9355·10⁻¹⁸ |
| 24 | 1.7045·10⁻⁹ | 1.70453120117·10⁻⁹ | 4.7813·10⁻²⁵ |
| 28 | 1.3525·10⁻¹⁰ | 1.35254770865·10⁻¹⁰ | 1.1538·10⁻²⁹ |
| 32 | 1.0769·10⁻¹² | 1.07689720773·10⁻¹² | 4.4835·10⁻³⁷ |
| 36 | 1.1216·10⁻¹⁶ | 1.12161611521·10⁻¹⁶ | |
| 40 | 6.8535·10⁻¹⁹ | 6.85355154306·10⁻¹⁹ | |
| 44 | 2.3313·10⁻²¹ | 2.33130821521·10⁻²¹ | |
| 48 | 2.3758·10⁻²⁶ | 2.37581468402·10⁻²⁶ | |

- Structure of the optimum (computed, [theory/convex/certs/decay.out](theory/convex/certs/decay.out), [report.txt](theory/convex/certs/report.txt)). At degree 4K the function g (class A) or the Euler sum (class E) has exactly K double zeros. In class A they sit just above the first K prime powers and move onto them as the degree grows: e^u = 2.0722 at degree 4; 2.0057, 3.0839 at degree 8; 2.000067, 3.00175, 4.0445, 5.0923 at degree 16; at degree 48 the first nine are 2, 3, 4, 5, 7, 8, 9, 11, 13 to 10 digits, then 16.0000005, 17.000003, 19.00015. In class E they sit at the first K primes only. No zero ever sits at 6, 10 or 12. B has no real root at any degree, so the condition ĝ ≥ 0 acts only through the parity rule. The remaining cost is carried by the first prime power without a zero (n = 23 carries 99.9 % at degree 48). So each four degrees buy one more prime power, in order, and the positions come out of the optimisation, not in.
- Certificates (computer-assisted, one code, [theory/convex/certs/certify.py](theory/convex/certs/certify.py), ball arithmetic at 800 bits; not refereed). For each degree 4, 8, …, 48 a polynomial B with rational coefficients is stored in `theory/convex/certs/certificate_A_d*.json`, with: B > 0 on the real line (exact Sturm count); g > 0 on [0, U] with U = 5 (U = 6 from degree 36 on; Taylor models on up to 3524 intervals) and beyond U (Bessel asymptotics with the remainder bound DLMF 10.40.10); W(g)/g(0) inside a ball whose upper end, rounded up to 12 digits, is the bound in the table (the strings `certified_upper_bound` in those files are rounded to nearest and must not be quoted as bounds). Assumed: the closed form of Φ∗Φ as a sum of K₀ and K₁ terms (derived there; agrees with direct quadrature to 25 digits), the ball arithmetic of arb, sympy's root count. With the truncation lemma: **λ\*(M) ≤ 2.376·10⁻²⁶ + 194/M²**, so the limit of the smallest conductor is at least 1 − 2.4·10⁻²⁶, from our own computation alone. Part I's certificate is far stronger (10⁻¹⁰⁶⁰); the point of ours is that it is independent of their code.
- Recheck with separate code ([cert_check.py](cert_check.py): nothing of the agent's code, no Bessel functions; g and W from the zero side by Gauss–Legendre quadrature of Ξ(t)²·B(t), 48000 nodes, 520 bits). For B = 1: W/g(0) = 0.075197132697534236812 from the zero side and from the prime sum, the same 20 digits as the agent. For the twelve stored polynomials: the zero side and the prime sum agree to 85 digits or better, and both agree with the centre of the certified ball to all 20 digits printed. B has no root with t² ≥ 0 (own root finder). Class E optima to degree 32: the same costs to 20 digits. For degrees 8, 24 and 48, g > 0 at all 381 grid points u = k/100, k ≤ 380 (a sample, not a proof; beyond u = 3.8 the values of g are below the working precision). Not rechecked by separate code: the interval proof that g > 0 everywhere, the optimality of the optima, and the constraint of class E.
- A row of Part I's Table 1 rebuilt from scratch (same script). The member J = 10 (double zeros at the logs of 2, 3, 4, 5, 7, 8, 9, 11, 13, 16; a 20 × 20 linear system) has ∫Ξ²P = 2.788432137 (Part I: between 2.7884318 and 2.7884322), W = 4.212556303·10⁻¹⁹ (Part I, rounded up: 4.2125564·10⁻¹⁹) and W/∫ = 1.510725776·10⁻¹⁹ (Part I: κ₁₀ ≤ 1.5107258·10⁻¹⁹). In our units this is W/g(0) = 9.4922·10⁻¹⁹; the optimum at the same degree 40 is 6.8535·10⁻¹⁹, lower, as it must be. So their functional and ours agree to 8 digits at the level 10⁻¹⁹.
- Decay (computed, twelve degrees; no proof). Class A: log₁₀ of the optimum ≈ 1.74 − 1.15·n, n the first prime power without a zero (rms 0.69). No stall up to degree 48. This is evidence that the infimum over the family is 0, as Part II claims; it is not a proof.

**Status.** The target is not our result. It is claimed in Part II, not refereed. Ours here: two independent numerical tests of the extremal function (two codes, costs 1.8·10⁻⁹ and −7·10⁻¹¹), the duality and the truncation lemma for compact supports with the explicit rate 194/M², a partial referee read, the observation that the prime-power structure is not needed (also theirs: their prime measure is arbitrary), and our own certified bound 2.376·10⁻²⁶ with the best polynomial at each degree up to 48 (computer-assisted, one code; costs rechecked by a second code).

## 5g. The double zeros and the prime powers: a trap that holds at every higher degree (26 November 2026, night)

**Status.** That the double zeros of the best polynomial converge to the prime powers as the degree grows is **not proved**. It is not in print either: Johns, Part II (sections 1.5 and 9.4) lists the convergence of the fixed-node family of Part I as open, and in the neighbouring fields (sphere packing, conformal bootstrap) the same statement about roots is numerical or conjectural in every source read ([theory/zeros/gate_report.md](theory/zeros/gate_report.md)). Proved here: a lemma that confines the zeros at **all** higher degrees once one test function of one degree is known. Computed with it: intervals around the first twelve prime powers that hold for every degree ≥ 48.

Notation: V_d = functions h with transform Ξ(t)²·B_h(t), B_h an even polynomial of degree ≤ 2d in t; cone C_d = {h ∈ V_d: h ≥ 0, B_h ≥ 0}; W(h) = 2·Σ Λ(n)·n^(−1/2)·h(log n); λ_d = min of W over C_d with h(0) = 1; G_d = a minimiser. Degree 48 in section 5f is d = 24.

**The optimum is a quadrature formula for the primes** (proved; [theory/zeros/zeros_trap.md](theory/zeros/zeros_trap.md), Theorem 2.2; finite-dimensional duality with the two points at infinity, Slater point Φ∗Φ + ε·g_d). For even d there are points u_i > 0 with weights μ_i > 0, points t_i with weights ν_i > 0 and α ≥ 0 with

`W(h) = λ_d·h(0) + Σ μ_i·h(u_i) + Σ ν_i·B_h(t_i) + α·b_d(h) for every h in V_d`

and every minimiser G_d has a double zero at every u_i. So the double zeros are the nodes of a formula that reproduces the sum over the prime powers on a (d + 1)-dimensional space, with one fixed node at 0 of weight λ_d. The shape is known from sphere packing (Cohn and Gonçalves 2019, Proposition 5.1). Computed for d ≤ 24: ν = 0, α = 0, d/2 nodes, weights μ_k = 2Λ(n_k)/√n_k to 10 digits for the older nodes.

**Trap lemma** (proved; same file, Lemma 3.1). Let w ∈ V_d′ with B_w ≤ 0 on the real line and W(w) > 0. Then for every d ≥ d′ and every minimiser G_d, the set S = {u > 0: w(u) > 0} contains a double zero of G_d, and the weights of the formula satisfy μ_d(S) ≥ (W(w) − λ_d·w(0))/sup_S w > 0.

Proof. Put w into the formula of degree d. The term λ_d·w(0) is ≤ 0, because −w has a non-negative transform, so |w| ≤ −w(0). The ν-term is ≤ 0 by B_w ≤ 0, and the α-term is ≤ 0. So Σ μ_i·w(u_i) ≥ W(w) > 0, and some node lies where w > 0. ∎

The mechanism is that of Cohn and Kumar 2009 (Proposition 4.3 and Lemma 4.4, where it localises the vectors of a lattice); a statement for the optima of every higher degree was not found in the sources read, and three named books on moment problems were not read, so no claim of novelty.

**How it is used.** Take gt in C_d′ with small cost kt = W(gt), h ∈ V_d′, and w = η·h − gt. The conditions are η·B_h ≤ B_gt and η·W(h) > kt. With gt = G_d′ and h = h_k, the function of V_d′ that is 1 at the node u_k and vanishes doubly at the other nodes (top coefficient 0), W(h_k) = μ_k + λ_d′·h_k(0), the threshold is η₀ = λ_d′/W(h_k), and for η = 2η₀ the set S is one small interval around u_k of half-width about (2η/G_d′″(u_k))^(1/2). In words: a zero of every later optimum must sit where the cost-λ function G_d′ is smaller than about 2λ; near a double zero that is a window of width about √λ.

**The traps at degree 48** (three codes for the numbers: the proof agent's [theory/zeros/trap2_d24.out](theory/zeros/trap2_d24.out) and the compute agent's [theory/zeros/traps/traps.out](theory/zeros/traps/traps.out), both with the Bessel closed form, and ours from the zero side, [trap_check.py](trap_check.py), [trap_check_d24_out.txt](trap_check_d24_out.txt); certificates: [theory/zeros/traps/trapcert.py](theory/zeros/traps/trapcert.py), [d3_summary.out](theory/zeros/traps/d3_summary.out)).

| n | u_k − log n at degree 48 | radius of S around log n, η = 2η₀ (numerical, our code) | certified radius R_n | certified weight in the trap ≥ | 2Λ(n)/√n |
|---|---|---|---|---|---|
| 2 | 2.3·10⁻²⁵ | 3.03·10⁻¹³ | 3.334·10⁻¹³ | 0.490 | 0.980 |
| 3 | 4.5·10⁻²⁴ | 1.315·10⁻¹² | 1.448·10⁻¹² | 0.634 | 1.269 |
| 4 | 1.2·10⁻²² | 6.61·10⁻¹² | 7.274·10⁻¹² | 0.346 | 0.693 |
| 5 | 3.5·10⁻²² | 1.162·10⁻¹¹ | 1.278·10⁻¹¹ | 0.719 | 1.440 |
| 7 | 9.2·10⁻²⁰ | 1.99·10⁻¹⁰ | 2.189·10⁻¹⁰ | 0.735 | 1.471 |
| 8 | 9.3·10⁻¹⁸ | 2.043·10⁻⁹ | 2.248·10⁻⁹ | 0.244 | 0.490 |
| 9 | 4.5·10⁻¹⁷ | 4.591·10⁻⁹ | 5.051·10⁻⁹ | 0.366 | 0.732 |
| 11 | 1.2·10⁻¹⁵ | 2.464·10⁻⁸ | 2.713·10⁻⁸ | 0.723 | 1.446 |
| 13 | 2.3·10⁻¹³ | 3.321·10⁻⁷ | 3.653·10⁻⁷ | 0.711 | 1.423 |
| 16 | 3.0·10⁻⁸ | 1.118·10⁻⁴ | 1.230·10⁻⁴ | 0.173 | 0.347 |
| 17 | 1.8·10⁻⁷ | 2.632·10⁻⁴ | 2.896·10⁻⁴ | 0.687 | 1.374 |
| 19 | 7.9·10⁻⁶ | 1.537·10⁻³ | 1.691·10⁻³ | 0.676 | 1.351 |

**Theorem (computer-assisted, one code for the interval part; not refereed).** For every degree ≥ 48 and every best polynomial of that degree, the function g has a double zero within R_n of log n for each of the twelve prime powers n = 2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, with R_n as in the table, and the quadrature weight of the zeros in that window is at least the value in the fifth column. What the certificate proves for each n, with gt the certified function of section 5f, h a rational rounding of h_k minus Φ∗Φ/1000, and a rational η close to 2η₀: B_gt − η·B_h has no real root (exact Sturm count); η·W(h) > W(gt) (balls); gt − η·h > 0 on [0, a] and [b, 6] (Taylor models, about 3400 intervals, 800 bits) and beyond 6 (Bessel asymptotics); R_n ≥ the distance from log n to a and to b. Assumed: the lemma (proved above), the membership of gt in the cone (the certificate of 5f), the closed Bessel form, DLMF 10.40.10, the ball arithmetic of arb, sympy's root count.

**The certificates checked with separate code** ([trapcert_check.py](trapcert_check.py), [trapcert_check_out.txt](trapcert_check_out.txt); zero side, no Bessel functions; numerical, not a second proof). For each of the twelve: our prime sum W(h) agrees with the ball in the file; η·W(h) − W(gt) > 0 with the stated margin; B_gt − η·B_h has no root with t² ≥ 0; w ≤ 0 at the stated endpoints a and b and w > 0 at log n; and w ≤ 0 at every grid point u = i/200 ≤ 3.8 outside (a, b).

**Traps from lower degrees** (numerical; [theory/zeros/traps/d2.out](theory/zeros/traps/d2.out), [theory/zeros/traps/traps.out](theory/zeros/traps/traps.out)). The same construction at each degree 4K = 4, …, 48 gives one interval for each of the 78 pairs (degree, node) at η = 2η₀. For the node at 2 the radius around log 2, valid for all higher degrees, is 0.21, 0.044, 0.014, 3.8·10⁻³, 4.2·10⁻⁴, 8.1·10⁻⁵, 2.3·10⁻⁵, 2.0·10⁻⁶, 2.1·10⁻⁸, 1.6·10⁻⁹, 9.5·10⁻¹¹, 3.0·10⁻¹³ at degrees 4, 8, …, 48. The half-width divided by (2λ/(μ_k·G″(u_k)))^(1/2) tends to √2 = 1.414214, and the offset u_k − log n_k itself is far smaller, about 9.55·λ for n = 2. In 4 of the 78 cases the interval of the newest node is wide and holds a second integer (degree 12: e^u = 4.357, holds log 5; degrees 24, 28, 40 likewise), so "the k-th zero stays at the k-th prime power" is a statement about the older nodes. Our own run at degree 8 ([trap_check_d4_out.txt](trap_check_d4_out.txt)): radius 0.0437 around log 2 and 0.156 around log 3, the same as the agent's 0.0437 and 0.156.

**What convergence still needs** ([theory/zeros/zeros_trap.md](theory/zeros/zeros_trap.md), sections 4 and 5).

- Proved: if g ≥ 0 has cost W(g) with g(0) = 1, then g(log n) ≤ √n·W(g)/(2Λ(n)) at every prime power, and if g″ ≥ m > 0 on a window around log n the minimum of g in the window is within (√n·W(g)/(Λ(n)·m))^(1/2) of log n.
- Proved under three hypotheses: (A) λ_d → 0, (B) a curvature bound G_d″ ≥ m > 0 near log n, uniform in d, (C) a zero in the window (one trap gives it). Then the zero converges to log n at rate λ_d^(1/2).
- (A) is open for the polynomial family. It follows if infinitely many members of the fixed-node family of Johns, Part I lie in the cone with the growth bound of its Conjecture 6.1(b) (proved implication; certified there for five members only). Theorem 1 of Part II does not give (A): its function is not a limit of polynomials with the positivity kept, as far as is known.
- Without (A): a trap in every window around log n would give a zero in every window at all high degrees. Whether such traps exist for every window is a uniqueness statement for positive solutions of the explicit formula on this test class. It is proved only for finite perturbations of the prime measure; the open part is about measures that grow.
- Proved: the transforms of the minimisers form a normal family in the strip |Im t| < 1/2 (because G_d ≥ 0), so subsequences of G_d converge uniformly, and if λ_d → 0 every limit vanishes at every log of a prime power.
- Computed: G_d″(u_k) settles for the early nodes (n = 2: 1.3508, 1.0929, 1.05640, 1.055894 at degrees 4, 12, 20, 48) and still drifts down for n ≥ 7; the fourth moment G_d⁗(0) tends to 945.65 and the pole term to 1.1200979 ([theory/zeros/moments.out](theory/zeros/moments.out)).

**The limit is not the function of Part II** (numerical, no proof). The pole term of the best polynomials tends to 1.1200979; that of the Part II function is 1.120099144 (section 5f). The best polynomial stays positive at log 6: g(log 6) = 3.66, 3.51, 3.40, 3.33·10⁻⁷ at degrees 36 to 48, with shrinking steps, while the Part II function vanishes at every integer. So the extremal function of this problem is not unique, if the limit exists: one vanishes at the prime powers only, one at all integers. Part II says the same of its own function and the fixed-node family (section 9.4(3)).

## 6. Size of the result

- Known before: positivity pins one prime at one support (Connes–Consani); explicit-formula lower bounds for discriminants (Odlyzko, Serre, Poitou, Diaz y Diaz) with one fixed positive kernel and the worst splitting type per prime.
- Comparison: our values at support log 8 (2.9995 and 4.9911) lie above the printed bounds with local corrections (2.9921 and 4.9640). Ours are values found by a search over a larger class (real positions) with a short support, theirs are proved bounds with a long kernel; the two are not the same quantity, and ours is an upper bound for q_min(L).
- In print before us (Johns, October 2026, not refereed; found on 26 November, late): the prime side and the zero side as unknown positive measures, the smallest-conductor reading, its value 1 for ζ with an explicit extremal function, and positive slack for Q(√5) and Q(√−3).
- Not found there by a keyword search (the preprints were not read in full): primes returned from no input by a search on the prime side; the pinning table; the role of the pole; the price of an edit; whole-number multiplicities and the prime ideal norms of quadratic fields; cubic Gamma factors.
- Proved, refereed once inside the repo (no failing step; two consequences of 5b reworded, small repairs in the proof of Theorem C): the price of a finite edit (section 5b, under RH); a feasible system is determined by its spectrum; thin changes of the primes are never feasible unless trivial (section 5c, RH not assumed).
- The convex version (section 5e): the positivity bound is sharp in degree 1, where it returns the primes with weight 1, and not in degrees 2 and 3. Degrees 1 and 2 are in Johns, Part I and II, with certified bounds; our numbers agree with theirs (section 5f).
- Global uniqueness: not proved by us; claimed in Johns, Part II. Proved here (refereed once inside the repo; the integer-weight case is also Theorem 3.11 of Johns, Part I): it holds among systems whose completed zeta function is entire of finite order with finitely many zeros off the line (section 5d), and a converse theorem for Beurling zeta functions without a functional equation.
- Our own certified bound for the convex optimum in degree 1 (section 5f): W(g)/g(0) ≤ 2.376·10⁻²⁶ with an explicit polynomial of degree 48, computer-assisted with one code, costs rechecked by a second code. Weaker than the bound in print (10⁻¹⁰⁶⁰), independent of it. Not found in Part I by a keyword search (its family has the zeros placed at the prime powers by construction): the best polynomial at each degree places its double zeros near the prime powers by itself, one more per four degrees, and they converge to the prime powers.
- The double zeros (section 5g): convergence to the prime powers is not proved, here or in print. Proved: a trap lemma (the mechanism is Cohn and Kumar's) that confines the zeros at every higher degree. Computer-assisted with it: for every degree ≥ 48 the best polynomial has a double zero within 3.4·10⁻¹³ of log 2, …, within 1.7·10⁻³ of log 19 (twelve prime powers; one code for the interval part, checked numerically by a second).
- The numerical statements are search results at small supports (X ≤ 14). Degree 3 was run once and has not converged.

## 7. Open points

1. The conjecture in degree 1 (no positive system other than the primes at conductor ≤ 1) is claimed as a theorem in Johns, Part II. Open for us: an independent check of that proof (sections 4, 5 and 7 there, and its four computer-assisted inputs); a referee for [theory/convex/convex_limit.md](theory/convex/convex_limit.md). Also open: a second interval proof that g > 0 for our own certificates (section 5f; one code so far). Convergence of the double zeros to the prime powers (section 5g): open are λ_d → 0 for the polynomial family, a curvature bound uniform in the degree, and traps in every window (a uniqueness statement on the test class Ξ²·polynomials); also a referee for [theory/zeros/zeros_trap.md](theory/zeros/zeros_trap.md), a second interval code for the twelve traps, and the three unread books named in [theory/zeros/gate_report.md](theory/zeros/gate_report.md).
2. The exact rate of q_min(L) → d, and the comparison with the published GRH discriminant bounds.
3. Cubic and quartic Gamma factors: does the minimum equal the smallest discriminant (23, 49, 117, 125, …)? The one cubic run reaches 22.63 at X = 14.
4. A high-precision scan stage, to continue section 3.3 beyond 11.
5. Sources not read by the prior-art searches: Serre 1975; Poitou's Bourbaki talk; Mestre 1986; Booker 2006; Jones–Roberts (arXiv 1610.01228, smallest Artin conductors) and Dixit (arXiv 1905.04193, general Dirichlet series with positive coefficients), abstracts only; Miller 2002 and Bober et al. 2015 (highest lowest zero), abstracts only; Bochner–Chandrasekharan 1956; Bombieri 2000; Yoshida 1992; Burnol; Connes–Consani, Selecta 2021 (arXiv 2006.13771); Connes–Consani–Moscovici, arXiv 2310.18423 and 2511.22755; Ryavec 1975; Kaczorowski–Perelli on general frequencies; Kurasov–Sarnak; Olevskii–Ulanovskii; Hilberdink–Lapidus; full texts of arXiv 2611.13129, 2608.24827, 2605.20224.

## 8. Validation

[validate.py](validate.py) reads the output files and checks every number quoted in sections 1 to 5b and in the tables and results of sections 5e, 5f and 5g (literature values excepted); 141 checks.
