# Seven lines, one law: what Hardy–Littlewood's singular series predict

> **In short:** Seven kinds of structure in the primes were measured on four windows (10⁷–10⁸, and 5·10⁷ or 2·10⁷ numbers at 10¹⁰, 10¹², 10¹⁵): repeated numbers, Ulam's diagonals, gap statistics, residue classes, the gap sequence as a chain, its information content, and reproducibility. The textbook formulas miss most of them by 1–46 %. All of them follow from one object, the Hardy–Littlewood correlation functions S(B)/ln^k x, once three things are done right: (1) the chance that a pattern of primes has no prime inside is expanded around independent interior numbers (a Bahadur expansion, order 4), (2) counts in intervals on a fixed grid keep the residues of the grid, (3) counts along a quadratic use the Bateman–Horn constants of polynomial pairs. With no fitted parameter, 25 of 26 statistics fit within three standard errors in every window where they were tested, including the fresh window at 10¹⁵; the one miss is the far tail of the gaps at the smallest scale (−0.9 %). The same expansion, with pairs as events, predicts the gaps between consecutive twin primes and their class bias (68 of 68 tests fit).

**Code:** [writeup/code/study/](writeup/code/study/). **Data:** [writeup/data/study/](writeup/data/study/) (results, theory values, verdicts, and the pre-registrations written before each round).

**Status:** computed. The Hardy–Littlewood and Bateman–Horn conjectures are conjectures; every prediction below assumes them. The expansions are our own derivations (heuristic: the series are asymptotic and truncated, the truncation error is stated). No parameter was fitted. Every number in the tables is printed by [tables.py](writeup/code/study/tables.py) from the verdict files.

## 1. The question and the method

The protocol asked for seven independent lines (numerical, geometric, statistical, algebraic, dynamic, information-theoretic, computational) and for every result to be checked against known theory: the prime number theorem, residue classes, analytic number theory.

- **Windows.** train = [10⁷, 10⁸) (5,096,876 primes); A = [10¹⁰, 10¹⁰ + 5·10⁷) (2,171,134); B = [10¹², 10¹² + 5·10⁷) (1,808,833); C = [10¹⁵, 10¹⁵ + 2·10⁷) (578,837; opened only in round 3). Each window is split into 10 equal blocks.
- **The test.** For a statistic s with prediction T: r_b = s_b − T_b on the ten blocks, z = mean(r)/(sd(r)/√10). The theory value is recomputed per block (λ = ln x changes inside a window). The prediction fits when |z| ≤ 3. The noise comes from the primes themselves, not from a model.
- **No sieve models.** Fake primes built from a sieve explain nothing about why the primes behave as they do, so the only fakes are coin flips with chance 1/ln n (Cramér) and the window's own gaps in random order. They answer one question: is there structure at all.
- **Pre-registration.** Before each round a file in [preregister/](writeup/data/study/preregister/) fixed the formulas, the cut-offs and the pass rule. Round 1 and round 3 were run blind; round 2 was built on the train window (single-gap shares only) and then run on A and B unchanged.

## 2. Round 1: the textbook formulas ([round1.txt](writeup/data/study/preregister/round1.txt), [verdicts_r1.json](writeup/data/study/verdicts_r1.json))

The crude Gallagher–Hardy–Littlewood law says a prime near x is followed by the gap h with chance S(0, h)·e^(−h/ln x)/ln x. Montgomery and Soundararajan give Var/mean = 1 − (ln H − 1 + γ + ln 2π)/ln x for counts in [n, n + H). Lemke Oliver–Soundararajan-type class shares follow by summing the gap law over h ≡ 0 (mod q). Bateman–Horn gives the expected count on each of 195 irreducible quadratics 4n² + bn + c (b = −4 … 4 even, c odd ≤ 79), the curves Ulam's spiral draws as diagonals.

| statistic | primes (train) | round-1 formula | rel. diff | z train | z A | z B |
|---|---|---|---|---|---|---|
| two equal gaps in a row | 0.028487 | 0.019565 | +45.6 % | 56.5 | 47.6 | 48.6 |
| P(next gap larger) − P(smaller) | -4.5126e-05 | 0 |  | -0.1 | -1.8 | 1.8 |
| gap 2 | 0.074817 | 0.066831 | +11.9 % | 46.3 | 49.0 | 28.8 |
| gap 4 | 0.074876 | 0.059659 | +25.5 % | 49.3 | 62.3 | 48.7 |
| gap 6 | 0.13121 | 0.10652 | +23.2 % | 58.0 | 71.4 | 57.6 |
| gap 8 | 0.057256 | 0.047544 | +20.4 % | 44.8 | 52.9 | 32.0 |
| gap 10 | 0.073689 | 0.056592 | +30.2 % | 55.2 | 71.6 | 67.6 |
| gap 12 | 0.092776 | 0.075784 | +22.4 % | 71.6 | 52.6 | 42.0 |
| gap 18 | 0.066882 | 0.053926 | +24.0 % | 67.0 | 127.4 | 81.1 |
| gap 30 | 0.039457 | 0.036422 | +8.3 % | 21.9 | 42.8 | 24.7 |
| Var/mean, H = 50 | 0.70526 | 0.69794 | +1.0 % | 13.0 | 9.1 | 10.5 |
| Var/mean, H = 500 | 0.57012 | 0.56737 | +0.5 % | 1.8 | 0.4 | 1.5 |
| Var/mean, H = 5000 | 0.44308 | 0.43681 | +1.4 % | 1.1 | 0.4 | 0.4 |
| next prime same class mod 3 | 0.43998 | 0.44269 | -0.6 % | -12.8 | -7.0 | -6.5 |
| same class mod 4 | 0.44457 | 0.44797 | -0.8 % | -20.8 | -5.9 | -5.1 |
| same class mod 5 | 0.17296 | 0.17687 | -2.2 % | -25.7 | -8.5 | -7.0 |
| same class mod 7 | 0.094808 | 0.096624 | -1.9 % | -11.5 | -2.6 | -2.3 |
| same class mod 8 | 0.17967 | 0.1842 | -2.5 % | -24.7 | -15.2 | -4.6 |
| χ²/curve on Ulam diagonals | 0.38614 | 1 | -61.4 % | -16.3 | -3.7 | -2.5 |
| Bateman–Horn slope on Ulam diagonals | 0.99985 | 1 | -0.0 % | 0.1 | 0.2 | -0.6 |

- Every line shows structure far beyond coin flips: against the coins |z| > 8 for every statistic except the up/down symmetry, the races (too noisy, see below), the information between gap classes 5 and 10 steps apart (and 3 steps apart at B), which is zero for the primes too and, at B only, the scatter on the quadratics ([verdicts_r1.json](writeup/data/study/verdicts_r1.json)).
- What fits: the reflection symmetry of consecutive gaps (a theorem about the singular series: S is invariant under x → span − x), the Bateman–Horn means (slope 0.99985), and Montgomery–Soundararajan for H = 500 and 5000.
- What misses: every gap share (by 8–30 %), the equal gaps (46 %), the class shares (0.6–2.5 %), the variance at H = 50, and the scatter on the quadratics, which is 0.39 of Poisson instead of 1.
- The race statistics (slope of π(x;q,N) − π(x;q,R) on √x/ln x) came out 1.14, 1.32, 1.31 against the Rubinstein–Sarnak mean 1 at train; inside a window of 5·10⁷ at 10¹⁰ the race moves by about 11 primes against a noise of about 1,500, so the holdouts cannot test it. The race itself is explained by the zeros in [FINDINGS.txt](FINDINGS.txt) (line on E(x) up to 2·10⁹) and is not repeated here.

## 3. Round 2: one law for consecutive primes ([theory2.py](writeup/code/study/theory2.py), [round2.txt](writeup/data/study/preregister/round2.txt))

**Step 1 (conjecture: Hardy–Littlewood).** For a finite set B, the numbers n + b (b ∈ B) are all prime with density S(B)/λ^|B|, λ = ln n, S(B) = ∏_p (1 − ν_p(B)/p)/(1 − 1/p)^|B|, ν_p = number of residues of B mod p. Computation: for p larger than the span of B every residue is distinct, so S(B) = ∏_(p ≤ span) (…) · G_k(span) with G_k(P) = ∏_(p > P) (1 − k/p)/(1 − 1/p)^k computed once (primes to 2·10⁷ plus the tail k(k−1)/(2X ln X)). Check: S({0, 2}) = 1.3203236, twice the twin prime constant 0.6601618 (computed).

**Step 2 (identity).** B is a run of consecutive primes when its interior I = {0 < m < span, m ∉ B} has no prime. By inclusion–exclusion, P(B) = Σ_(A ⊆ I) (−1)^|A| S(B ∪ A)/λ^(|B|+|A|). Divided by S(B)/λ^|B|, this is the chance that a 0/1 process X_m on I with correlation functions E∏_(m∈A) X_m = S(B ∪ A)/(S(B) λ^|A|) avoids every point.

**Step 3 (a wrong turn, kept because it explains the size of the misses).** The textbook resummation writes the avoidance chance as exp(−Σκ₁ + Σκ₂ − …) with factorial cumulants κ. On the train window this converged, but to the wrong values: 1.6 times too high at gap 90. The exponential formula for set functions assumes the blocks of a partition can overlap; for 0/1 variables X² = X, and the diagonal is lost. At gap 6 the interior is {2, 4} and S(0, 2, 4, 6) = 0, so the exact chance is 1 − x − y (x = ρ₂/λ, y = ρ₄/λ), while the exponential gives e^(−x−y−xy).

**Step 4 (derived; heuristic truncation).** Expand around independent interior numbers instead. With p_m = ρ_m/λ, ρ_m = S(B ∪ {m})/S(B), and Z_m = (p_m − X_m)/(1 − p_m), which has mean 0, exactly

  E∏(1 − X_m) = ∏(1 − p_m) · E∏(1 + Z_m) = ∏(1 − p_m) · Σ_A E[∏_(m∈A) Z_m].

The single terms vanish. The connected parts of the pair, triple and quadruple terms are c_ab/((1 − p_a)(1 − p_b)), −κ₃/∏(1 − p) and +κ₄/∏(1 − p), where c, κ₃, κ₄ are the joint cumulants of the X's (set-partition formulas in `Pattern.log_avoid`). Exponentiating the connected parts:

  ln P(B) = ln S(B) − |B| ln λ + Σ ln(1 − p_m) + Σ c̃ − Σ κ̃₃ + Σ κ̃₄.

At gap 6 order 2 is exact. The share per prime is P(B)·λ, averaged over the window with weight 1/ln t (integrated in u = ln t).

**Step 5 (computed, train only, the development data).** The prediction divided by the actual share at train: gap 30 gives 1.165 / 1.017 / 1.005 / 1.003 at orders 1 to 4; gap 90 gives 2.06 / 1.09 / 1.01 / 0.99. Order 1 alone (independent interior numbers with their Hardy–Littlewood densities) is the form already used in [FINDINGS.txt](FINDINGS.txt) for Euler–Rodríguez chains; it is 1–16 % off at gaps 6–30. Order 4 was frozen for spans ≤ 120 and order 3 above, before the holdouts were opened.

**Everything else is a sum over patterns.** The share of gap h is P({0, h})·λ; two equal gaps sum P({0, h, 2h}) over h; the class share mod q sums the gap law over q | h and divides by the total (the law sums to 1.0010, 1.0005, 1.0003 and 1.0002 at train, A, B and C); the joint law of two consecutive gaps is P({0, h₁, h₁ + h₂}); mutual information, entropies and the Markov gain are functions of the joint law.

## 4. Counts on a grid ([theory2.py](writeup/code/study/theory2.py), `var_ratio_aligned`)

The variance of the count in [n, n + H) needs only pairs. With the exact Hardy–Littlewood pair sum, Var/mean = 1 − 1/λ + (2/(Hλ)) Σ_(h<H) (H − h)(S(0, h) − 1). This still missed H = 50 (0.7005 against 0.7053). The reason is the grid: the intervals start at multiples of H, so position a inside every interval is the same residue mod 2 and mod 5, and only the pairs (a, a + h) inside one interval count. With u(a) = 1 when n + a is prime to 10, R(h) the pair factor without the primes 2 and 5 (R(h) = R₀ ∏_(p|h, p≠2,5) (p − 1)/(p − 2), R₀ = C₂/(0.6/0.64)) and c = 1/(½·⅘)²:

  Var/mean = 1 + (2c Σ_h R(h) N_h − H²)/(Hλ),  N_h = Σ_a u(a) u(a + h).

(Derivation: E c² = Σ_a p_a + 2 Σ_(a<b) E[X_a X_b] with p_a = c^(1/2) u(a)/λ and E[X_a X_b] = c u(a) u(b) R(b − a)/λ²; Σ_a p_a = H/λ.) For H = 50 this gives 0.70556 / 0.77451 / 0.81207 at train / A / B against 0.70526 / 0.77551 / 0.81187. The H = 50 values of A and B had been seen in round 1 when this was derived; the formula has nothing to tune, and C (round 3) was not seen.

## 5. Ulam's diagonals scatter less than Poisson ([theory_geo.py](writeup/code/study/theory_geo.py))

On each quadratic f the count of prime values has mean Σ_n C(f)/ln f(n) (Bateman–Horn), and round 1 found that the counts scatter with χ²/curve = 0.386 instead of 1. The same pair argument as in section 4, with polynomials:

  Var/E = 1 − C/λ + (2/(C·N·λ)) Σ_d (N − d)(C_d − C²)  (λ-weights done exactly per n in the code),

where C_d is the Bateman–Horn constant of the pair {f(x), f(x + d)}: C_d/C² = ∏_p (1 − ω_d(p)/p)/(1 − w(p)/p)², w(p) the number of roots of f mod p and ω_d(p) = |R ∪ (R − d)|. For large p the factor is generic (ω = 2w) except at p | d (ω = w: the shifted curve has the same roots) and at p | 16d² − D, D = b² − 64c (ω = 3: one root of f(x + d) is a root of f(x)). For 4n² + 1 the neighbours repel: C₁/C² = 0.5175 (brute-force product to 3,000: 0.51751). The predicted χ²/curve is 0.406 at train, 0.694 at A and 0.831 at B, against 0.386, 0.715 and 0.721 (block z −0.0, 1.1, −1.0). This is a Montgomery–Soundararajan-type law for a curve instead of an interval: the counts are rigid because values a short step apart along the curve share prime factors in a structured way.

## 6. Rounds 2 and 3 against all four windows ([round3.txt](writeup/data/study/preregister/round3.txt), [verdicts_r3.json](writeup/data/study/verdicts_r3.json))

| statistic | primes (train) | cluster law | rel. diff train | z train | z A | z B | z C (10¹⁵) |
|---|---|---|---|---|---|---|---|
| gap 2 | 0.074817 | 0.074866 | -0.07 % | -1.0 | 1.7 | 0.2 | -0.8 |
| gap 4 | 0.074876 | 0.074866 | +0.01 % | 0.2 | -0.1 | 1.2 | -0.5 |
| gap 6 | 0.13121 | 0.13133 | -0.09 % | -1.2 | -0.2 | -0.8 | 0.2 |
| gap 8 | 0.057256 | 0.057224 | +0.06 % | 0.3 | 0.6 | 0.0 | -0.7 |
| gap 10 | 0.073689 | 0.073738 | -0.07 % | -0.7 | -0.2 | -0.1 | 1.0 |
| gap 12 | 0.092776 | 0.092726 | +0.05 % | 0.8 | 0.2 | 0.2 | -1.1 |
| gap 18 | 0.066882 | 0.067023 | -0.21 % | -1.0 | 0.4 | 1.1 | 2.5 |
| gap 30 | 0.039457 | 0.039559 | -0.26 % | -1.2 | -1.6 | -0.7 | 0.5 |
| two equal gaps in a row | 0.028487 | 0.02865 | -0.57 % | -2.5 | 1.6 | -0.4 | -0.4 |
| next prime same class mod 3 | 0.43998 | 0.44006 | -0.02 % | -0.4 | -0.8 | -1.1 | 3.0 |
| same class mod 4 | 0.44457 | 0.44424 | +0.07 % | 2.1 | 0.6 | 1.9 | -1.2 |
| same class mod 5 | 0.17296 | 0.1732 | -0.14 % | -1.5 | -0.7 | 0.4 | 0.4 |
| same class mod 7 | 0.094808 | 0.094956 | -0.16 % | -1.1 | 0.3 | -0.0 | 0.3 |
| same class mod 8 | 0.17967 | 0.17964 | +0.02 % | 0.2 | 1.3 | 1.3 | -0.2 |
| Var/mean, H = 50 | 0.70526 | 0.70556 | -0.04 % | -0.4 | 1.3 | -0.5 | -1.0 |
| Var/mean, H = 500 | 0.57012 | 0.5683 | +0.32 % | 1.2 | 0.2 | 1.3 | 0.3 |
| Var/mean, H = 5000 | 0.44308 | 0.43695 | +1.40 % | 1.1 | 0.4 | 0.3 | -0.6 |
| mutual information of consecutive gaps (bits) | 0.26881 | 0.26859 | +0.08 % | 0.8 | -0.0 | -0.8 | -1.6 |
| entropy of the gap class (bits) | 4.0345 | 4.034 | +0.01 % | 1.3 | -0.7 | -1.3 | -2.0 |
| conditional entropy H(next \| gap) (bits) | 3.7656 | 3.7654 | +0.00 % | 0.4 | -0.7 | -0.7 | -1.5 |
| four primes in a row in arithmetic progression | 0.0016987 | 0.0017314 | -1.89 % | -2.7 | 1.8 | -2.2 | -0.3 |
| χ²/curve on Ulam diagonals | 0.3861 | 0.406 | -4.9 % | -0.0 | 1.1 | -1.0 | |

- Gap classes are c = min(g/2, 20). The information quantities are estimators; their z is computed after subtracting the estimator bias measured by drawing samples of the block size from the law's own joint distribution ([bias_r3.py](writeup/code/study/bias_r3.py)). Before that correction they gave z = −3.9 (mutual information) and −6.0 (conditional entropy) at C, where a block holds only 58,000 primes for 400 cells.
- The Markov gain of the dynamic line (a cross-validated log-loss gain) is the mutual information minus a finite-sample penalty of about K(K − 1)/(2N ln 2) bits and minus the drift between the two halves of a window; it is not tested separately.

**The tail** (gaps of at least 3, 4 and 5 times ln x):

| window | ln x | gaps ≥ 3 ln x: primes / law, z | ≥ 4 ln x | ≥ 5 ln x |
|---|---|---|---|---|
| train | 17.8 | h ≥ 54: 0.03549 / 0.03583, -4.1 | h ≥ 72: 0.009629 / 0.009754, -2.3 | h ≥ 90: 0.002803 / 0.002843, -1.4 |
| A | 23.0 | h ≥ 70: 0.03806 / 0.03809, -0.3 | h ≥ 92: 0.0117 / 0.01166, 0.6 | h ≥ 116: 0.00347 / 0.003502, -0.8 |
| B | 27.6 | h ≥ 82: 0.04254 / 0.0425, 0.3 | h ≥ 110: 0.01345 / 0.01334, 1.4 | h ≥ 138: 0.004216 / 0.004241, -0.7 |
| C | 34.5 | h ≥ 104: 0.04137 / 0.04138, -0.1 | h ≥ 138: 0.01405 / 0.01427, -1.4 | h ≥ 172: 0.004492 / 0.004612, -1.3 |

The one miss is the tail at train: the law is 0.9–1.4 % too high for gaps of 54 and more. There λ ≈ 17.8 and h/λ ≥ 3, where the series is clearly asymptotic: at gap 90 the step from order 3 to order 4 is still 2.4 %. At 10¹⁰ and above, the same h/λ fits.

**The joint law cell by cell** (Pearson χ² over the 19 × 19 pairs of consecutive gaps ≤ 38, cells with expectation below 5 merged):

| window | gap pairs in the cells | χ² | cells | χ²/cell | z |
|---|---|---|---|---|---|
| train | 4203489 | 234.3 | 276 | 0.849 | -1.8 |
| A | 1480206 | 261.1 | 276 | 0.946 | -0.6 |
| B | 1040754 | 250.0 | 276 | 0.906 | -1.1 |
| C | 260307 | 286.4 | 276 | 1.038 | 0.4 |

## 7. The same law for twin primes ([theory_twins.py](writeup/code/study/theory_twins.py), [round4.txt](writeup/data/study/preregister/round4.txt), [verdicts_twins.json](writeup/data/study/verdicts_twins.json))

"The next twin pair comes h later" is the same avoidance problem with composite events: B = {0, 2, h, h + 2}, and the interior events are the pairs e_m = {m, m + 2}, 0 < m < h. Their moments are again singular series, of the union of the events, divided by λ to the number of new elements: E∏_(e∈A) X_e = S(B ∪ ∪A)/(S(B) λ^|∪A ∖ B|). The Bahadur form of section 3 applies unchanged (order 3; the steps at h = 300 and λ = 17.8 are −1.155, −1.199 and −1.2015 for orders 1 to 3). The share per twin pair is P(B)/(S(0, 2)/λ²), averaged with weight 1/ln² t. The pre-registration was written before any twin statistic was computed; the windows hold 381,332 / 124,843 / 86,479 / 22,013 twin pairs.

| statistic | twins (train) | law, order 1 / 3 | z train | z A | z B | z C |
|---|---|---|---|---|---|---|
| twin gap 6 | 0.010146 | 0.010146 / 0.010146 | -0.0 | 2.8 | -2.2 | -0.4 |
| twin gap 12 | 0.027097 | 0.027055 / 0.027055 | 0.1 | -0.3 | -0.2 | -0.9 |
| twin gap 18 | 0.019833 | 0.019746 / 0.019743 | 0.5 | 0.2 | 1.1 | -0.1 |
| twin gap 24 | 0.012265 | 0.012319 / 0.012319 | -0.4 | -0.2 | 1.3 | -0.1 |
| twin gap 30 | 0.038088 | 0.038235 / 0.038217 | -0.3 | 0.4 | 0.6 | -0.6 |
| twin gap 36 | 0.010998 | 0.010867 / 0.010865 | 1.6 | -1.6 | -0.4 | 2.1 |
| twin gap 42 | 0.034584 | 0.034844 / 0.034807 | -0.8 | 0.5 | 1.2 | 0.7 |
| twin gap 48 | 0.018695 | 0.018627 / 0.018585 | 0.6 | -1.4 | -0.5 | -0.9 |
| twin gap 60 | 0.027417 | 0.027358 / 0.027285 | 0.5 | 0.2 | 1.2 | 0.0 |
| twin gap 90 | 0.026958 | 0.027141 / 0.026927 | 0.2 | -1.0 | -1.0 | -0.7 |
| twin gap 120 | 0.020667 | 0.021143 / 0.020902 | -1.0 | 1.2 | 1.2 | 0.2 |
| twin gap 150 | 0.019914 | 0.020274 / 0.019958 | -0.3 | -1.1 | 0.4 | 1.2 |
| twin gap 180 | 0.02325 | 0.023922 / 0.023389 | -0.6 | -1.1 | 1.3 | 0.0 |
| twin gap 210 | 0.025673 | 0.026214 / 0.025402 | 1.2 | -1.1 | -0.2 | -0.4 |
| twin gap ≤ 300 | 0.72637 | 0.73687 / 0.72634 | 0.0 | 0.4 | 0.8 | -0.4 |
| same class mod 5 (gaps ≤ 300) | 0.30457 | 0.30527 / 0.30434 | 0.5 | -0.7 | 0.6 | -0.1 |
| same class mod 7 (gaps ≤ 300) | 0.16918 | 0.16929 / 0.16883 | 0.9 | -0.1 | -0.4 | 1.7 |

All 68 tests fit (|z| ≤ 3). The share of twin gaps up to 300 is 0.72637 against 0.72634 at train; order 1 (independent interior pairs) gives 0.73687. The twin pairs also avoid repeating their class: two consecutive twin pairs with the same p mod 5 occur at a share of 0.3046 among twin gaps ≤ 300 against 0.3043 predicted, a Lemke Oliver–Soundararajan bias for twin primes that follows from the same singular series.

## 8. What the seven lines turn out to be

- **Numerical, statistical, algebraic, dynamic, information-theoretic**: five views of one object, the law of runs of consecutive primes. Its only input is the singular series; its only scale is ln x. The textbook versions miss because they truncate it: the Poisson factor e^(−h/ln x) ignores that the numbers between two primes are correlated with them (round 1: up to 30 %), the exponential resummation mistreats X² = X (step 3), and the Montgomery–Soundararajan asymptotics ignore the grid (H = 50).
- **Geometric**: Ulam's diagonals are the same law for polynomials: the means are Bateman–Horn and the scatter is set by the Bateman–Horn constants of polynomial pairs.
- **Computational**: the pre-registered windows and the block test. Four windows over eight orders of magnitude, one of them opened only after the law was frozen.
- **What is not new**: the Hardy–Littlewood and Bateman–Horn conjectures, Gallagher's Poisson law, the Lemke Oliver–Soundararajan bias and the Montgomery–Soundararajan variance are literature; that a parameter-free Hardy–Littlewood count fits gap counts has been used before (Brent 1974 for small gaps; Odlyzko, Rubinstein and Wolf 1999 for jumping champions). What this study adds, as far as the files in this archive show, is the Bahadur form and its accuracy (0.1–0.3 % for gaps up to 30, at four scales), the grid term in the variance, and the variance law along Ulam's diagonals; the literature for the last two was not searched beyond this archive.

## Validation

- [engine.py](writeup/code/study/engine.py) and [stats_r3.py](writeup/code/study/stats_r3.py) compute the statistics of the primes from an exact segmented sieve of each window (checked: π(10⁶) = 78,498; the train count 5,096,876 = π(10⁸) − π(10⁷)).
- [theory_r2.py](writeup/code/study/theory_r2.py), [theory_r3.py](writeup/code/study/theory_r3.py) and [theory_geo.py](writeup/code/study/theory_geo.py) compute every prediction from the singular series alone; [report2.py](writeup/code/study/report2.py), [report3.py](writeup/code/study/report3.py) and [bias_r3.py](writeup/code/study/bias_r3.py) make the verdicts; [tables.py](writeup/code/study/tables.py) prints the tables above.
