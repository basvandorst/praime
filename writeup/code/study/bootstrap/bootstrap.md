# The primes as the answer to a positivity problem (26 November 2026)

**In short.** Weil's explicit formula turns the Riemann hypothesis into the statement that a quadratic form is positive. We keep the Gamma factor and the pole, and make the prime side unknown: any set of real numbers x > 1 with the weights of an Euler product. Three findings, all numerical.

1. Started with no primes, the search returns 2, 3, 5, 7 to 15 digits or more, and 11 to 7 digits. Each prime, moved alone, is pinned more tightly as the support grows (2 + 10⁻³⁸ fails at support log 11).
2. The mechanism is the pole. Without a pole (three Dirichlet data tried), or when Euler factors with a minus sign are allowed, positivity did not pin the points.
3. With a pole there is a variational form: every prime system has a smallest conductor that keeps the form positive. The smallest value found over all positive systems is 1 for Γ_R(s); for Γ_R(s)Γ_R(s+1) and Γ_R(s)² it rises with the support and stays below 3 and 5, the smallest discriminants, and from support log 10 on the best systems found sit near the prime ideal norms of Q(√−3) and Q(√5). That these are the limits is a conjecture.

Status of each statement: computed (best found by a search), not proved. The first step (one prime, one support) is in print: Connes and Consani 2021.

## 1. What is in print

- Connes and Consani, *Spectral triples and ζ-cycles*, arXiv 2106.01715, section 2.3: at support log 3, with the one prime 2 replaced by a real p, positivity fails at p = 1.9999 and at p = 2.0005. Section 2.2: the archimedean part alone loses positivity just past support log 2. They move one prime at one support. Our matrix reproduces their numbers ([cc_check_out.txt](cc_check_out.txt)): smallest eigenvalue at support log 3: 5.68·10⁻⁸ (theirs: below 6·10⁻⁸); at support log 11: 2.93·10⁻⁴⁸ with 97 basis functions and still falling (theirs: 2.389·10⁻⁴⁸); negative at 1.9999 and at 2.0005.
- Discriminant bounds from the explicit formula (Odlyzko; Serre; Poitou, *Sur les petits discriminants*, 1977; Diaz y Diaz, *Tables minorant la racine n-ième du discriminant*, 1980). They fix one kernel F ≥ 0 and, as a "local correction", take for each small prime the splitting type with the smallest contribution. Printed bounds for |d|: imaginary quadratic 2.9657 without and 2.9921 with local corrections; real quadratic 4.9436 and 4.9640 (Diaz y Diaz, Tables 4 and 3; Odlyzko's GRH table: 2.962 and 4.951). So the minimum over splitting types, prime by prime at integer positions, is in print. Odlyzko 1990 lists the optimal kernel (Problem 2.2) and the share of the prime ideals (Problem 6.1) as open.
- Not found in print (two prior-art searches of 26 November 2026, about 60 pages opened; the sources in section 7 were not read): all primes free together at real positions, primes found from no input, the pinning as a function of the support, the role of the pole, the smallest conductor over all positive systems with the full quadratic form, and the questions whether it converges to the smallest discriminant and whether the minimiser tends to the prime ideal norms.

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

Statement wanted: the primes are the only positive Beurling system whose form (Riemann datum, q = 1) is positive for every support. Not proved. What we have:

- **Convex version, computed** ([relaxed.py](relaxed.py), [relaxed_out.txt](relaxed_out.txt)). If the multiplicities may be any real numbers ≥ 0 (weights on a grid of step 0.01), the form is linear in the weights and the largest smallest eigenvalue is a concave problem; cutting planes with a linear program give the global optimum on the grid: 7.0·10⁻⁴, 3.2·10⁻⁵, 1.8·10⁻⁵, 1.6·10⁻⁶ at X = 2.5, 3, 3.5, 4 (the primes alone: 5.7·10⁻⁸ at X = 3). The weights sit at 2.00–2.01 (total 1.003 at X = 3, 1.0001 at X = 4) and 3.00–3.04 (total 1.197 at X = 3.5, 1.005 at X = 4). So at a finite support the optimum is slightly off the primes and its conductor is slightly below 1 (0.9999984 at X = 4); the values fall with the support. Uniqueness, if true, is a statement about the limit.
- **Proved so far** (full proofs in [theory/global/global.md](theory/global/global.md); referee report pending; none of these uses that the multiplicities are integers).
  1. *Positivity as a half-plane condition.* P is feasible ⟺ Ξ_P′/Ξ_P continues to Re s > 1/2 with non-negative real part there, Ξ_P = ½·s(s − 1)·π^(−s/2)·Γ(s/2)·ζ_P(s) (for ζ: Lagarias's criterion). Then Ξ_P has no zero in Re s > 1/2 and |Ξ_P(σ + it)| grows with σ.
  2. *Lemma G (one-dimensional sphere-packing argument).* Let τ ≥ 0 be an even measure with no mass in (−1, 1). If its Fourier transform equals α·δ₀ + constant on (−1, 1), then τ = α·(Dirac comb − δ₀). The test function is −sin(2π|x|)/(2π) on [−1, 1], whose transform sin²(πξ)/(π²(ξ² − 1)) is ≥ 0 for |ξ| ≥ 1 and vanishes there exactly on the integers. Consequence: a Beurling zeta function with Riemann's functional equation (either sign) is ζ. This may be Ryavec's theorem (Rocky Mountain J. Math. 5 (1975) 623–627, not read); Beurling 1951 (Ark. Mat. 1, not read) characterises ζ by growth without a functional equation.
  3. *A converse theorem without positivity (Theorem E).* If (s − 1)ζ_P(s) is entire of finite order and |Λ_P(1 − s)| ≤ C(1 + |s|)^N·|Λ_P(s)| on a half-plane Re s ≥ c′ with 1 < c′ < 3 (Λ_P = π^(−s/2)Γ(s/2)ζ_P), then ζ_P = ζ. No functional equation is assumed: a one-sided growth comparison is enough.
  4. *Corollary F.* If P is feasible, Ξ_P is entire of finite order and only finitely many of its zeros lie off the critical line, then P is the set of primes. In particular: feasible, entire, all zeros on the line ⟹ primes.
  5. For a feasible P: P = primes ⟺ ζ_P has Riemann's functional equation ⟺ the spectrum μ_P is a sum of point masses with integer weights.
- **The gap, exactly.** A second feasible system would have to be of one of three kinds: (G1) its zeta function has no continuation of finite order across the line (this includes a continuous spectrum and point masses of non-integer weight); (G2) Ξ_P is entire of finite order with infinitely many zeros to the left of the line, enough of them (a positive density at bounded distance) to make up a growth excess of Λ_P(1 − s) over Λ_P(s); (G3) trivial zeros are missing. Equivalent question on the additive side: does positivity force the Fourier transform of the measure on the Beurling integers to be α·δ₀ + constant on (−1, 1)?
- **Counterexample search** (the prover's linear program over absolutely continuous prime densities, heights up to 60): the best smallest value of Re Ξ′/Ξ on the line stayed negative (−1.1·10⁻³ to −3.7·10⁻⁴), creeping to 0 as the family grows, with the mass concentrating in peaks. No second feasible system was found.
- **Literature** (one search; main sources on positivity criteria not read): nothing found that derives a functional equation from positivity; no discrete Beurling system with a zeta function meromorphic in the plane, beyond the arithmetic ones.

## 6. Size of the result

- Known before: positivity pins one prime at one support (Connes–Consani); explicit-formula lower bounds for discriminants (Odlyzko, Serre, Poitou, Diaz y Diaz) with one fixed positive kernel and the worst splitting type per prime.
- Comparison: our values at support log 8 (2.9995 and 4.9911) lie above the printed bounds with local corrections (2.9921 and 4.9640). Ours are values found by a search over a larger class (real positions) with a short support, theirs are proved bounds with a long kernel; the two are not the same quantity, and ours is an upper bound for q_min(L).
- New here, as far as the search went: the prime side as the unknown of an optimisation; primes returned from no input; the role of the pole and of positive weights; the smallest-conductor reading with quadratic fields.
- Proved, refereed once inside the repo (no failing step; two consequences of 5b reworded, small repairs in the proof of Theorem C): the price of a finite edit (section 5b, under RH); a feasible system is determined by its spectrum; thin changes of the primes are never feasible unless trivial (section 5c, RH not assumed).
- Global uniqueness: not proved. Proved (referee pending): it holds among systems whose completed zeta function is entire of finite order with finitely many zeros off the line (section 5d), and a converse theorem for Beurling zeta functions without a functional equation.
- The numerical statements are search results at small supports (X ≤ 14). Degree 3 was run once and has not converged.

## 7. Open points

1. A proof or a counterexample for the conjecture in degree 1 (is there a positive Beurling system, other than the primes, whose form is positive at conductor ≤ 1 for every support?).
2. The exact rate of q_min(L) → d, and the comparison with the published GRH discriminant bounds.
3. Cubic and quartic Gamma factors: does the minimum equal the smallest discriminant (23, 49, 117, 125, …)? The one cubic run reaches 22.63 at X = 14.
4. A high-precision scan stage, to continue section 3.3 beyond 11.
5. Sources not read by the prior-art searches: Serre 1975; Poitou's Bourbaki talk; Mestre 1986; Booker 2006; Jones–Roberts (arXiv 1610.01228, smallest Artin conductors) and Dixit (arXiv 1905.04193, general Dirichlet series with positive coefficients), abstracts only; Miller 2002 and Bober et al. 2015 (highest lowest zero), abstracts only; Bochner–Chandrasekharan 1956; Bombieri 2000; Yoshida 1992; Burnol; Connes–Consani, Selecta 2021 (arXiv 2006.13771); Connes–Consani–Moscovici, arXiv 2310.18423 and 2511.22755; Ryavec 1975; Kaczorowski–Perelli on general frequencies; Kurasov–Sarnak; Olevskii–Ulanovskii; Hilberdink–Lapidus; full texts of arXiv 2611.13129, 2608.24827, 2605.20224.

## 8. Validation

[validate.py](validate.py) reads the output files and checks every number quoted in sections 1 to 5b (literature values excepted).
