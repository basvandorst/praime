# Seven lines, one law: what Hardy–Littlewood's singular series predict

> **In short:** Seven kinds of structure in the primes were measured on four windows (10⁷–10⁸, and 5·10⁷ or 2·10⁷ numbers at 10¹⁰, 10¹², 10¹⁵): repeated numbers, Ulam's diagonals, gap statistics, residue classes, the gap sequence as a chain, its information content, and reproducibility. The textbook formulas miss most of them by 1–46 %. All of them follow from one object, the Hardy–Littlewood correlation functions S(B)/ln^k x, once three things are done right: (1) the chance that a pattern of primes has no prime inside is expanded around independent interior numbers (a Bahadur expansion, order 4), (2) counts in intervals on a fixed grid keep the residues of the grid, (3) counts along a quadratic use the Bateman–Horn constants of polynomial pairs. With no fitted parameter the gap statistics (N(30)/N(6) at 10⁸: 0.2903 against 0.2899 counted; 30 overtakes 6 as the most common gap at 1.742·10³⁵, against 1.743·10³⁵ from the Hardy–Littlewood constants in the literature), class shares and short-interval variances fit within the block test at all four windows, including the windows at 10¹⁵ and 10¹⁸ that were opened only after the law was frozen. Three pieces were added after seeing data (the grid term, the curve law, an estimator-bias correction) and are marked as such; the curve law then passed on fresh data. What does not fit: the gap tail at ≥ 3 ln x at the smallest scale (−0.9 %, partly the truncation of the series); a small excess of the variance over the pair law for 500 ≤ H ≤ 10⁴ (open); and, with the order-4 expansion, the full distribution of the number of primes in an interval of 30 or 60 at the level of 10⁻³ (30–80 times closer than Gallagher's Poisson limit); the exact Hardy–Littlewood sums for H = 30 remove those misses and match the data to the fourth decimal. The same expansion, with pairs as events, predicts the gaps between consecutive twin primes and their class bias (68 of 68 tests fit); with Bateman–Horn's singular series it predicts the gaps between consecutive n with n² + 1 prime (3 of 45 block tests miss, 0.7 expected by chance); for the linear forms x and 2x + 1 it predicts the gaps between consecutive Sophie Germain primes (55 of 56 fit), and along n² + n + 41 the gaps between prime values (37 of 37); and in two dimensions it predicts the nearest-neighbour distances of Gaussian and of Eisenstein primes (71 of 72 tests fit in eight boxes).

**Code:** [writeup/code/study/](writeup/code/study/). **Data:** [writeup/data/study/](writeup/data/study/) (results, theory values, verdicts, and the pre-registrations written before each round).

**Status:** computed. The Hardy–Littlewood and Bateman–Horn conjectures are conjectures; every prediction below assumes them. The expansions are our own derivations (heuristic: the series are asymptotic and truncated; section 3 gives the error of order 4 against the exact inclusion–exclusion sum for gaps up to 42). No parameter was fitted. Every number in the tables is printed by [tables.py](writeup/code/study/tables.py) from the verdict files.

## 1. The question and the method

The protocol asked for seven independent lines (numerical, geometric, statistical, algebraic, dynamic, information-theoretic, computational) and for every result to be checked against known theory: the prime number theorem, residue classes, analytic number theory.

- **Windows.** train = [10⁷, 10⁸) (5,096,876 primes); A = [10¹⁰, 10¹⁰ + 5·10⁷) (2,171,134); B = [10¹², 10¹² + 5·10⁷) (1,808,833); C = [10¹⁵, 10¹⁵ + 2·10⁷) (578,837; opened only in round 3). Each window is split into 10 equal blocks.
- **The test.** For a statistic s with prediction T: r_b = s_b − T_b on the ten blocks, z = mean(r)/(sd(r)/√10). The theory value is recomputed per block (λ = ln x changes inside a window). The prediction fits when |z| ≤ 3. The noise comes from the primes themselves, not from a model.
- **No sieve models.** Fake primes built from a sieve explain nothing about why the primes behave as they do, so the only fakes are coin flips with chance 1/ln n (Cramér) and the window's own gaps in random order. They answer one question: is there structure at all.
- **Pre-registration.** Before each round a file in [preregister/](writeup/data/study/preregister/) fixed the formulas, the cut-offs and the pass rule; the commit that adds each file comes before the commit with its results (rounds 1–11: f802663, 3f6d791, ee457c0, f9d6fee, 4a15168, 6de4fb1, 204ca7e, 8cd9279, 5cbecb5, 520ab35, 92d808c). The clock times written inside the files are estimates and run 10–25 minutes later than these commits. The commit order shows the order of committing; since some result commits followed their pre-registration within a minute, the repository alone cannot prove that each plan was fixed before its computation, and the pre-registration should be read as stated, not as verified. Round 2 was built on the train window (single-gap shares only) and then run on A and B unchanged.
- **Added after seeing data** (marked where they appear): the grid term of the variance (section 4; round 2 had registered the unaligned pair sum, which missed H = 50 at train: 0.70049 against 0.70526), the variance law along quadratics (section 5; tested afterwards on fresh data in round 9), and the estimator-bias correction of the information statistics (section 6; the registered test missed at C before it).

## 2. Round 1: the textbook formulas ([round1.txt](writeup/data/study/preregister/round1.txt), [verdicts_r1.json](writeup/data/study/verdicts_r1.json))

The crude Gallagher–Hardy–Littlewood law says a prime near x is followed by the gap h with chance S(0, h)·e^(−h/ln x)/ln x. Montgomery and Soundararajan give Var/mean = 1 − (ln H − 1 + γ + ln 2π)/ln x for counts in [n, n + H). Lemke Oliver–Soundararajan-type class shares follow by summing the gap law over h ≡ 0 (mod q). Bateman–Horn gives the expected count on each of 198 irreducible quadratics 4n² + bn + c (b = −4 … 4 even, c odd ≤ 79; only 4n² ± 4n + 1 = (2n ± 1)² is left out), the curves Ulam's spiral draws as diagonals. (The pre-registration of round 1 says 195; that was a miscount in the text, the code used all 198 from the start.)

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

**Step 3 (a wrong turn, kept because it explains the size of the misses).** A first version exponentiated the Ursell functions of distinct interior points, exp(−Σκ₁ + Σκ₂ − …), with κ the connected parts of the correlation functions over distinct points. On the train window this converged, but to the wrong values: 1.6 times too high at gap 90. The exponential formula for set functions assumes the blocks of a partition can overlap; for 0/1 variables X² = X, and the diagonal terms are lost. (The true factorial-cumulant series of the number of interior primes, ln P(N = 0) = Σ(−1)^k κ_[k]/k!, is exact and contains these diagonal terms; version 1 was not that series.) At gap 6 the interior is {2, 4} and S(0, 2, 4, 6) = 0, so the exact chance is 1 − x − y (x = ρ₂/λ, y = ρ₄/λ), while the exponential gives e^(−x−y−xy).

**Step 4 (derived; heuristic truncation).** Expand around independent interior numbers instead. With p_m = ρ_m/λ, ρ_m = S(B ∪ {m})/S(B), and Z_m = (p_m − X_m)/(1 − p_m), which has mean 0, exactly

  E∏(1 − X_m) = ∏(1 − p_m) · E∏(1 + Z_m) = ∏(1 − p_m) · Σ_A E[∏_(m∈A) Z_m].

The single terms vanish. The connected parts of the pair, triple and quadruple terms are c_ab/((1 − p_a)(1 − p_b)), −κ₃/∏(1 − p) and +κ₄/∏(1 − p), where c, κ₃, κ₄ are the joint cumulants of the X's (set-partition formulas in `Pattern.log_avoid`). Exponentiating the connected parts:

  ln P(B) = ln S(B) − |B| ln λ + Σ ln(1 − p_m) + Σ c̃ − Σ κ̃₃ + Σ κ̃₄.

At gap 6 the exact value is (1 − x)(1 − y)(1 + t) with t = −xy/((1 − x)(1 − y)), while the exponentiated form gives (1 − x)(1 − y)e^t: a relative error of 9·10⁻⁶ at λ = 17.8, so the overlap problem is pushed to higher order, not removed. The share per prime is P(B)·λ, averaged over the window with weight 1/ln t (integrated in u = ln t).

**Step 5 (computed, train only, the development data).** The prediction divided by the actual share at train: gap 30 gives 1.165 / 1.017 / 1.005 / 1.003 at orders 1 to 4; gap 90 gives 2.06 / 1.09 / 1.01 / 0.99. Order 1 alone (independent interior numbers with their Hardy–Littlewood densities) is the form already used in [FINDINGS.txt](FINDINGS.txt) for Euler–Rodríguez chains; it is 1–16 % off at gaps 6–30. Order 4 was frozen for spans ≤ 120 and order 3 above, before the holdouts were opened.

**The truncation error (from the review; independent code in [review/c3_bahadur.py](writeup/code/study/review/c3_bahadur.py) → [c3_bahadur.out](writeup/data/study/review/c3_bahadur.out), [c3_bahadur_big.out](writeup/data/study/review/c3_bahadur_big.out)).** For gaps up to 42 the live interior has at most 20 points, so the exact Hardy–Littlewood inclusion–exclusion over all its subsets can be summed. Order 4 lies above the exact value by 0.017 % (gap 12), 0.062 % (18), 0.185 % (30), 0.254 % (36) and 0.322 % (42) at λ = 17.8, and by 0.012 % at gap 42 at λ = 34.5. So the agreement of 0.1–0.3 % with the primes at train is at the size of the truncation error itself, and the error grows with h/λ: it is the natural suspect for the tail miss of section 6.

**Everything else is a sum over patterns.** The share of gap h is P({0, h})·λ; two equal gaps sum P({0, h, 2h}) over h; the class share mod q sums the gap law over q | h and divides by the total (the law sums to 1.0010, 1.0005, 1.0003 and 1.0002 at train, A, B and C); the joint law of two consecutive gaps is P({0, h₁, h₁ + h₂}); mutual information, entropies and the Markov gain are functions of the joint law.

## 4. Counts on a grid ([theory2.py](writeup/code/study/theory2.py), `var_ratio_aligned`)

The variance of the count in [n, n + H) needs only pairs. With the exact Hardy–Littlewood pair sum, Var/mean = 1 − 1/λ + (2/(Hλ)) Σ_(h<H) (H − h)(S(0, h) − 1). This still missed H = 50 (0.7005 against 0.7053). The reason is the grid: the intervals start at multiples of H, so position a inside every interval is the same residue mod 2 and mod 5, and only the pairs (a, a + h) inside one interval count. With u(a) = 1 when n + a is prime to 10, R(h) the pair factor without the primes 2 and 5 (R(h) = R₀ ∏_(p|h, p≠2,5) (p − 1)/(p − 2), R₀ = C₂/(0.6/0.64)) and c = 1/(½·⅘)²:

  Var/mean = 1 + (2c Σ_h R(h) N_h − H²)/(Hλ),  N_h = Σ_a u(a) u(a + h).

(Derivation: E c² = Σ_a p_a + 2 Σ_(a<b) E[X_a X_b] with p_a = c^(1/2) u(a)/λ and E[X_a X_b] = c u(a) u(b) R(b − a)/λ²; Σ_a p_a = H/λ.) For H = 50 this gives 0.70556 / 0.77451 / 0.81207 at train / A / B against 0.70526 / 0.77551 / 0.81187. The H = 50 values of A and B had been seen in round 1 when this was derived; the formula has nothing to tune, and C (round 3) was not seen.

## 5. Ulam's diagonals scatter less than Poisson ([theory_geo.py](writeup/code/study/theory_geo.py))

On each quadratic f the count of prime values has mean Σ_n C(f)/ln f(n) (Bateman–Horn), and round 1 found that the counts scatter with χ²/curve = 0.386 instead of 1. The same pair argument as in section 4, with polynomials:

  Var/E = 1 − C/λ + (2/(C·N·λ)) Σ_d (N − d)(C_d − C²)  (λ-weights done exactly per n in the code),

where C_d is the Bateman–Horn constant of the pair {f(x), f(x + d)}: C_d/C² = ∏_p (1 − ω_d(p)/p)/(1 − w(p)/p)², w(p) the number of roots of f mod p and ω_d(p) = |R ∪ (R − d)|. For large p the factor is generic (ω = 2w) except at p | d (ω = w: the shifted curve has the same roots) and at p | 16d² − D, D = b² − 16c (ω = 3: one root of f(x + d) is a root of f(x)). For 4n² + 1 the neighbours repel: C₁/C² = 0.5175 (brute-force product to 3,000: 0.51751). The predicted χ²/curve is 0.406 at train, 0.694 at A and 0.831 at B, against 0.386, 0.715 and 0.721 (block z −0.0, 1.1, −1.0). This is a Montgomery–Soundararajan-type law for a curve instead of an interval: the counts are rigid because values a short step apart along the curve share prime factors in a structured way. The law was derived after round 1 had shown the low scatter, so the three windows above are not an independent test; round 9 below is. The collision primes p | 16d² − D are included up to 3·10⁴; including all of them moves the train value by about +0.001 (review, c5_ulam.out).

**Round 9: fresh data and a second family** ([r9.py](writeup/code/study/r9.py), [round9.txt](writeup/data/study/preregister/round9.txt), [r9.json](writeup/data/study/r9.json)). For a general f = ax² + bx + c the coincidence rule reads a²d² ≡ D (mod p). On n ∈ [10⁶, 1.05·10⁶), every value tested by gmpy2 (no sieve), the pre-registered test compares the mean of (obs − E)²/E over the curves with the mean predicted Var/E, the noise taken from the scatter over the curves:

| family | curves | prime values | χ²/curve | law | z | Bateman–Horn: mean (obs − E)/E, z |
|---|---|---|---|---|---|---|
| Ulam diagonals 4n² + bn + c (all 198) | 198 | 681,550 | 0.499 | 0.542 | -0.9 | -0.00015, -0.2 |
| Ulam diagonals, the 159 distinct curves | 159 | 547,321 | 0.531 | 0.542 | -0.2 | +0.00039, 0.3 |
| Euler's n² + n + c, c odd ≤ 199 | 100 | 361,921 | 0.494 | 0.522 | -0.4 | +0.00074, 0.5 |

Both fit. Since 4n² + 4n + c = f₋₄,c(n + 1), 39 of the 198 Ulam curves are shifted copies of others; the middle row keeps only the 159 distinct curves, so that the scatter over curves is not counted twice (review, point 10). The counts along a quadratic scatter at about half the Poisson value at this size, as the pair constants predict. The Bateman–Horn constants run from 0.62 to 6.64 (n² + n + 41) in both families; the test is on the mean over the curves, not curve by curve.

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
- The Markov gain of the dynamic line (a cross-validated log-loss gain) was compared with the mutual information in round 2 and missed in all three windows (z = −3.2, −5.4, −8.7, [verdicts_theory_r2.json](writeup/data/study/verdicts_theory_r2.json)). It is not an estimate of the mutual information: it equals the mutual information minus a finite-sample penalty of about K(K − 1)/(2N ln 2) bits (3.4·10⁻⁴ at B, where the shortfall is 3.9·10⁻⁴) and minus the drift between the two halves of a window, and it is not used further.
- The bias correction of the information statistics was added after the registered test had missed at C (t₉ = −3.9 and −6.0) and at B (conditional entropy, −3.0); it uses only the law, not the data.

**The tail** (gaps of at least 3, 4 and 5 times ln x):

| window | ln x | gaps ≥ 3 ln x: primes / law, z | ≥ 4 ln x | ≥ 5 ln x |
|---|---|---|---|---|
| train | 17.8 | h ≥ 54: 0.03549 / 0.03583, -4.1 | h ≥ 72: 0.009629 / 0.009754, -2.3 | h ≥ 90: 0.002803 / 0.002843, -1.4 |
| A | 23.0 | h ≥ 70: 0.03806 / 0.03809, -0.3 | h ≥ 92: 0.0117 / 0.01166, 0.6 | h ≥ 116: 0.00347 / 0.003502, -0.8 |
| B | 27.6 | h ≥ 82: 0.04254 / 0.0425, 0.3 | h ≥ 110: 0.01345 / 0.01334, 1.4 | h ≥ 138: 0.004216 / 0.004241, -0.7 |
| C | 34.5 | h ≥ 104: 0.04137 / 0.04138, -0.1 | h ≥ 138: 0.01405 / 0.01427, -1.4 | h ≥ 172: 0.004492 / 0.004612, -1.3 |

The one miss is the tail at train beyond 3 ln x: the law is 0.9–1.4 % too high for gaps of 54 and more (the tails beyond 4 and 5 ln x pass the test there). There λ ≈ 17.8 and h/λ ≥ 3, where the series is clearly asymptotic: at gap 90 the step from order 3 to order 4 is still 2.4 %. At 10¹⁰ and above, the same h/λ fits.

**Round 8: the next order** ([r8.py](writeup/code/study/r8.py), [round8.txt](writeup/data/study/preregister/round8.txt), [r8_tail.json](writeup/data/study/r8_tail.json), [r8_blocks.json](writeup/data/study/r8_blocks.json)). The fifth cumulant, κ₅ = Σ_π (−1)^(|π|−1) (|π|−1)! ∏_(blocks) μ(block) over the 52 set partitions of five interior points, enters with the sign (−1)⁵; every partition carries λ⁻⁵, so κ₅ = K₅/λ⁵ with K₅ built from singular series alone. Added for spans 50–120 (and order 4 in place of order 3 for spans 122–200), it is about half the order-4 term at every gap (gap 90 at λ = 17.8: −0.037 at order 4, −0.018 at order 5, in ln P), and it lowers the train tail by 0.29 % (gaps ≥ 54: 0.035827 → 0.035722; primes 0.035494). The block t₉ goes from −4.09 to −2.80 (≥ 54), −2.30 to −1.19 (≥ 72) and −1.41 to −0.37 (≥ 90). The pre-registration expected a shift of about 1 %; it is 0.29 %, enough to bring the test inside |t₉| ≤ 3 but leaving 0.64 % (2.8 binomial standard errors) to the higher orders, which keep the same sign, or to chance. With the exact-sum comparison of section 3 (order 4 is 0.32 % above the exact law already at gap 42 and λ = 17.8, growing with h), truncation is a plausible cause of most of the rest, but that is not shown: the miss is partly truncation.

**Round 17: a window at 10¹⁸** ([r17.py](writeup/code/study/r17.py) → [r17.json](writeup/data/study/r17.json), [round17.txt](writeup/data/study/preregister/round17.txt)). [10¹⁸, 10¹⁸ + 2·10⁶), 48,427 primes listed one by one with gmpy2.next_prime (no sieve; below 2⁶⁴ the test has no known exception), λ = 41.4; the frozen law with order 2 for spans above 300 (the law sums to 1.0001):

| statistic | primes near 10¹⁸ | cluster law | rel. diff | t₉ | z (binomial) |
|---|---|---|---|---|---|
| gap 2 | 0.033535 | 0.031856 | +5.27 % | 2.0 | 2.1 |
| gap 4 | 0.030768 | 0.031856 | -3.42 % | -2.0 | -1.4 |
| gap 6 | 0.061474 | 0.060384 | +1.80 % | 0.8 | 1.0 |
| gap 8 | 0.029261 | 0.028587 | +2.36 % | 1.1 | 0.9 |
| gap 10 | 0.039152 | 0.0376 | +4.13 % | 3.1 | 1.8 |
| gap 12 | 0.051521 | 0.052641 | -2.13 % | -1.6 | -1.1 |
| gap 18 | 0.047783 | 0.046447 | +2.88 % | 1.6 | 1.4 |
| gap 30 | 0.045904 | 0.045367 | +1.18 % | 0.4 | 0.6 |
| two equal gaps in a row | 0.01398 | 0.013565 | +3.06 % | 0.6 | 0.8 |
| next prime same class mod 3 | 0.46976 | 0.47069 | -0.20 % | -0.4 | -0.4 |
| same class mod 4 | 0.47269 | 0.47254 | +0.03 % | 0.1 | 0.1 |
| same class mod 5 | 0.21166 | 0.21144 | +0.10 % | 0.1 | 0.1 |
| same class mod 7 | 0.12491 | 0.1285 | -2.79 % | -2.2 | -2.4 |
| same class mod 8 | 0.2163 | 0.21408 | +1.04 % | 1.4 | 1.2 |
| entropy of the gap class (bits) | 3.471 | 3.4602 | +0.31 % | 0.9 |  |
| gaps ≥ 124 | 0.044232 | 0.044079 | +0.35 % | 0.2 | 0.2 |
| gaps ≥ 166 | 0.015673 | 0.014576 | +7.53 % | 1.6 | 2.0 |
| gaps ≥ 208 | 0.0050592 | 0.004805 | +5.29 % | 0.9 | 0.8 |

17 of 18 pass the block test (gap 10 at t₉ = 3.1 is the exception; its binomial z is 1.8); every binomial z lies within ±2.4. The law frozen on [10⁷, 10⁸) holds ten orders of magnitude higher.

**The joint law cell by cell** (Pearson χ² over the 19 × 19 pairs of consecutive gaps ≤ 38, cells with expectation below 5 merged):

| window | gap pairs in the cells | χ² | cells | χ²/cell | z |
|---|---|---|---|---|---|
| train | 4203489 | 234.3 | 276 | 0.849 | -1.8 |
| A | 1480206 | 261.1 | 276 | 0.946 | -0.6 |
| B | 1040754 | 250.0 | 276 | 0.906 | -1.1 |
| C | 260307 | 286.4 | 276 | 1.038 | 0.4 |

**Round 19: jumping champions** ([r19.py](writeup/code/study/r19.py) → [r19.json](writeup/data/study/r19.json), [round19.txt](writeup/data/study/preregister/round19.txt)). The most common gap between consecutive primes up to x is the jumping champion. N(d, x), the number of consecutive pairs with gap d up to x, is the integral of the share of gap d over the primes, so the frozen law gives it with no data. At 10⁸, MERSENNE.md section 52 counted N(30)/N(6) = 222,847/768,752 = 0.2899, where its model S(d)·e^(−d/ln x) gave 0.333 to 0.362; the cluster law gives 0.2903. Then without data: 30 overtakes 6 at x = 1.7418·10³⁵, and 210 overtakes 30 at x = 6.418·10⁴²⁵; for the gap that follows a prime near x (the local share) the crossings come earlier, at 10^34.79 and 10^425.4. Oliveira e Silva's prime-gap page ([references/web/sweet.ua.pt_tos_gaps.html-b328ab1e.txt](references/web/sweet.ua.pt_tos_gaps.html-b328ab1e.txt), with Hardy–Littlewood constants computed with Herzog) puts the two changes of jumping champion at 1.74274·10³⁵ and 6.42869·10⁴²⁵ (Odlyzko, Rubinstein and Wolf, Experimental Math. 8, 1999, proposed the primorial champions); the cluster law agrees with them to 0.05 % and 0.2 % in x. Both rest on the Hardy–Littlewood constants, so this is a consistency check of the expansion, not an independent result. (The crude model alone puts the local crossings at 10^36.2 and 10^428.8.)

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

## 8. The same law on a parabola: Landau's n² + 1 ([theory_landau.py](writeup/code/study/theory_landau.py), [round5.txt](writeup/data/study/preregister/round5.txt), [verdicts_landau.json](writeup/data/study/verdicts_landau.json))

Whether n² + 1 is prime infinitely often is Landau's fourth problem. Here the question is how the n with n² + 1 prime follow each other. The singular series of a set of shifts B is now Bateman–Horn's: S_f(B) = ∏_p (1 − ω_B(p)/p)/(1 − 1/p)^|B| with ω_B(p) = |{(r − b) mod p}|, r running over the roots of x² + 1 mod p (two roots for p ≡ 1 mod 4, none for p ≡ 3 mod 4, the root 1 for p = 2). For large p the factor is generic, (1 − |B|(1 + χ₄(p))/p)/(1 − 1/p)^|B|, except where two shifted roots coincide: r − b = r′ − b′ with r − r′ = ±2i, that is p | (b − b′)² + 4. Check: S_f({0}) = 1.37279, the Bateman–Horn constant of n² + 1 being 1.37281 (the product converges conditionally; it runs to 2·10⁷ here). The gap law of section 3 then applies with λ = ln(n² + 1): order 4 for gaps ≤ 60, order 3 above (g = 60: −2.645, −2.793, −2.802, −2.805 for orders 1 to 4 at λ = 30).

The data are the even n in F1 = [10⁶, 3·10⁶), F2 = [10⁸, 10⁸ + 10⁶) and F3 = [10¹⁰, 10¹⁰ + 2·10⁵) with n² + 1 prime, tested one by one with gmpy2 (no sieve; the values in F3 are near 10²⁰, so these are BPSW probable primes): 94,545, 37,081 and 5,974 of them.

| statistic | F1: data / law | t₉, z_bin | F2: data / law | t₉, z_bin | F3: data / law | t₉, z_bin |
|---|---|---|---|---|---|---|
| n-gap 2 | 0.04946 / 0.04916 | 0.4, 0.4 | 0.0377 / 0.03855 | -1.2, -0.9 | 0.03432 / 0.03085 | 1.5, 1.6 |
| n-gap 4 | 0.09853 / 0.09832 | 0.4, 0.2 | 0.07621 / 0.07711 | -1.0, -0.6 | 0.05809 / 0.0617 | -1.0, -1.2 |
| n-gap 6 | 0.08579 / 0.08514 | 1.0, 0.7 | 0.06823 / 0.06901 | -0.8, -0.6 | 0.05525 / 0.05652 | -0.3, -0.4 |
| n-gap 8 | 0.04585 / 0.04575 | 0.4, 0.1 | 0.03608 / 0.0371 | -1.0, -1.0 | 0.03064 / 0.0304 | 0.0, 0.1 |
| n-gap 10 | 0.1183 / 0.1203 | -2.0, -1.9 | 0.1013 / 0.1013 | -0.0, -0.0 | 0.08572 / 0.08524 | 0.2, 0.1 |
| n-gap 12 | 0.03018 / 0.0313 | -1.5, -2.0 | 0.02648 / 0.02744 | -1.2, -1.1 | 0.02093 / 0.02374 | -1.5, -1.4 |
| n-gap 20 | 0.06531 / 0.06548 | -0.2, -0.2 | 0.06559 / 0.0622 | 6.3, 2.7 | 0.05441 / 0.05689 | -0.6, -0.8 |
| n-gap 30 | 0.03823 / 0.03761 | 3.9, 1.0 | 0.04086 / 0.04078 | 0.1, 0.1 | 0.04068 / 0.04086 | -0.1, -0.1 |
| n-gap 50 | 0.01389 / 0.01352 | 1.2, 1.0 | 0.01837 / 0.01863 | -0.3, -0.4 | 0.02076 / 0.02205 | -0.9, -0.7 |
| n-gap ≥ 60 | 0.05191 / 0.05171 | 0.4, 0.3 | 0.1051 / 0.1041 | 0.6, 0.6 | 0.1671 / 0.1692 | -0.7, -0.4 |
| n-gap ≥ 100 | 0.005786 / 0.005889 | -0.4, -0.4 | 0.02033 / 0.02003 | 0.5, 0.4 | 0.0447 / 0.04675 | -0.8, -0.7 |
| same n mod 3 | 0.282 / 0.2793 | 2.5, 1.8 | 0.2893 / 0.2894 | -0.0, -0.0 | 0.2913 / 0.297 | -0.7, -1.0 |
| same n mod 5 | 0.2788 / 0.279 | -0.1, -0.2 | 0.292 / 0.289 | 1.4, 1.3 | 0.2937 / 0.2963 | -0.5, -0.5 |
| same n mod 7 | 0.09606 / 0.09597 | 0.2, 0.1 | 0.1031 / 0.1017 | 0.9, 0.9 | 0.1157 / 0.1068 | 3.1, 2.2 |
| same n mod 13 | 0.04716 / 0.04724 | -0.1, -0.1 | 0.05221 / 0.05362 | -1.1, -1.2 | 0.05676 / 0.05833 | -1.1, -0.5 |

45 tests; 3 with |t₉| > 3 (expected by chance 0.7); largest |z_bin| = 2.7.


The block test exceeds 3 in three of the 45 cases against 0.7 expected by chance (the chance of three or more is about 0.03). The class bias appears on both kinds of moduli: consecutive n repeat their class mod 3 in 28.2 % of the cases (all 3 classes are allowed, 1/3 without bias) and mod 5 in 27.9 % (the 3 allowed classes 0, 1, 4; 1/3 without bias), and the law gets both. The binomial z, added with the results and taking its scale from the counts instead of from 10 blocks (it ignores the dependence of consecutive gaps), stays within 2.7 everywhere.

**A note on the test (computed, applies to every table).** The block z divides by a standard deviation estimated from 10 blocks, so it is a Student t with 9 degrees of freedom: |t| > 3 has chance 0.015, not 0.0027. Over the 96 tests of section 6 the block test exceeds 3 four times against 1.4 expected; three of them are the information estimators before their bias is removed, the fourth is the train tail. Over the 68 twin-prime tests it exceeds 3 zero times (1.0 expected).

## 9. Second differences, and the gap sequence in bits ([r6.py](writeup/code/study/r6.py), [r7.py](writeup/code/study/r7.py), [round6.txt](writeup/data/study/preregister/round6.txt), [round7.txt](writeup/data/study/preregister/round7.txt))

**Second differences (numerical line).** The second difference of the primes, p_(n+1) − 2p_n + p_(n−1) = g_n − g_(n−1), has P(d = k) = Σ_h P(0, h, 2h + k). The 20 tests (k = −6, −2, 0, 2, 6 in four windows) all fit ([verdicts_r6.json](writeup/data/study/verdicts_r6.json); 0 block-test misses, 0.3 expected). Two gaps that differ by 2 are 2.4 times as common as two equal gaps (0.0683 and 0.0285 at train): {0, h, 2h} survives mod 3 only when 6 | h (theorem), and the law carries the rest. P(d = k) = P(d = −k) holds exactly for the law, because (0, h, 2h + k) and (0, h + k, 2h + k) are mirror images and S is invariant under reflection; the primes agree within the noise (0.068348 and 0.068196 at train).

**Bits per prime (information line).** Code the gap class c = min(g/2, 20) of every prime with each model; the cross-entropy is the code length in bits per prime (M0: P(h) ∝ e^(−h/λ) on even h; M1: × S(0, h); M2: the cluster law; M3: the cluster law as a Markov chain from the joint law; the measured conditional entropies use the Miller–Madow correction):

| model or measurement | train | A | B | C (10¹⁵) |
|---|---|---|---|---|
| M0 coins | 4.1097 | 4.0818 | 3.9672 | 3.7414 |
| M1 crude Gallagher | 4.0460 | 4.0215 | 3.9104 | 3.6891 |
| M2 cluster law, i.i.d. | 4.0345 | 4.0156 | 3.9068 | 3.6871 |
| M3 cluster law, Markov | 3.7658 | 3.7993 | 3.7268 | 3.5487 |
| measured H(c_n+1 | c_n) | 3.7657 | 3.7993 | 3.7269 | 3.5487 |
| measured H(c_n+1 | c_n, c_n-1) | 3.6336 | 3.7022 | 3.6525 | 3.4981 |

Each layer of theory removes bits; the three-point law as a Markov chain codes the primes at their measured order-1 conditional entropy (within 1.1·10⁻⁴ bits in all four windows). But knowing two previous gaps saves another 0.13 bits per prime at train (0.05 at 10¹⁵), which no three-point law can supply.

**Where those bits live: a hidden chain on the residues (dynamic line).** Two previous gaps can reveal p_n mod 3 where one cannot (a gap ≡ 0 mod 6 hides it, the gap before may show it), and likewise mod 5 and 7. So take the residue a = p_n mod Q as a hidden state: conditioning on n ≡ a (mod Q) changes only the local factors of the singular series at the primes p | Q,

  S_a(B) = S(B) · ∏_(p|Q) 1[a + b ≢ 0 (mod p) for all b ∈ B] / (1 − ν_p(B)/p)

(identity check: the average of S_a(B) over all a mod 30 is S(B) to 10⁻¹⁶), and the same Bahadur expansion gives the gap law P_a(h) of each state. The gap sequence becomes a hidden Markov chain a → a + h whose class statistics are products of 8 × 8 (Q = 30) or 48 × 48 (Q = 210) matrices T_c[a, a′] = Σ_(h∈c, a+h≡a′) P_a(h):

| window | measured H₁ / H₂ | H₁ − H₂ measured | mod 6: H₁ / H₂, H₁ − H₂ | mod 30 | mod 210 |
|---|---|---|---|---|---|
| train | 3.7657 / 3.6336 | 0.1321 | 3.7709 / 3.6668, 0.1041 | 3.7596 / 3.6327, 0.1269 | 3.7610 / 3.6290, 0.1320 |
| A | 3.7993 / 3.7022 | 0.0971 | 3.8045 / 3.7268, 0.0777 | 3.7947 / 3.7012, 0.0935 | 3.7954 / 3.6983, 0.0971 |
| B | 3.7269 / 3.6525 | 0.0744 | 3.7313 / 3.6712, 0.0601 | 3.7229 / 3.6510, 0.0719 | 3.7232 / 3.6486, 0.0746 |
| C | 3.5487 / 3.4981 | 0.0506 | 3.5534 / 3.5119, 0.0415 | 3.5467 / 3.4973, 0.0494 | 3.5467 / 3.4956, 0.0511 |

The chain reproduces the drop H₁ − H₂ to 79 % mod 6, 96 % mod 30 and to within 10⁻⁴ bits mod 210 at train (0.1320 against 0.1321; 0.0971 against 0.0971 at 10¹⁰, 0.0746 against 0.0744 at 10¹², 0.0511 against 0.0506 at 10¹⁵). The pre-registered tolerance (0.005 bits) fails for every window mod 6 (H₂ is 0.033 too high at train), holds for every H₂ mod 30 and for every value mod 210 (by 0.0047 at most); mod 30 misses H₁ at train by 0.0061. The match of the drop mod 210 is a difference of two values that are each about 0.004 too low. Q = 2310 was not tried. The chain's levels sit 0.002–0.005 bits below the measured entropies: conditional independence of consecutive gaps given the residue over-states their dependence slightly, the remainder of the three-point coupling (the long-range part through λ) that residues do not carry. In words: beyond one step back, the extra memory of the prime gaps is, to this precision, the residue of the current prime modulo 2·3·5·7.

## 10. Counts in intervals: the third cumulant, and long intervals ([r10.py](writeup/code/study/r10.py), [round10.txt](writeup/data/study/preregister/round10.txt), [r10_skew.json](writeup/data/study/r10_skew.json), [r10_var.json](writeup/data/study/r10_var.json), [r10_var_3000000000.json](writeup/data/study/r10_var_3000000000.json))

**Third cumulant.** For 0/1 variables c³ = Σ_a X_a + 6 Σ_(a<b) X_a X_b + 6 Σ_(a<b<c) X_a X_b X_c, so on the aligned grid the moments of the count are E c = Σ p_a, E c² = Σ p_a + 2 Σ μ_ab, E c³ = Σ p_a + 6 Σ μ_ab + 6 Σ μ_abc, with μ_X = S₀(X)/λ^|X| and S₀ the singular series conditioned on n ≡ 0 (mod 10) (section 9). The third cumulant κ₃ = E c³ − 3 E c E c² + 2 (E c)³ divided by the mean is 1 for Poisson counts:

| window | H = 50: measured / law, t₉ | H = 500: measured / law, t₉ |
|---|---|---|
| train | 0.3116 / 0.3119, -0.1 | 0.1279 / 0.1248, 0.2 |
| A | 0.4395 / 0.4375, 1.1 | 0.2944 / 0.2539, 1.9 |
| B | 0.5143 / 0.5154, -0.8 | 0.3340 / 0.3441, -0.5 |
| C | 0.6008 / 0.5996, 0.2 | 0.4757 / 0.4479, 1.2 |

All eight fit (at H = 500 the test is weak: window A deviates by 16 % with t₉ = 1.9). The counts are far less skewed than Poisson counts (κ₃/mean ≈ 0.31 at H = 50 near 10⁷ instead of 1), and the triple sums of the singular series give the value.

**Round 20: the whole distribution of the count** ([r20.py](writeup/code/study/r20.py) → [r20_H30.json](writeup/data/study/r20_H30.json), r20_H60_*.json, [round20.txt](writeup/data/study/preregister/round20.txt)). Gallagher proved that the Hardy–Littlewood conjecture makes the number of primes in [n, n + H) Poisson in the limit. The law gives the distribution itself: the probability generating function G(z) = E∏(1 − (1 − z)X_m) has the same Bahadur expansion with every p_m replaced by (1 − z)p_m and every k-point cumulant multiplied by (1 − z)^k; evaluated at the 64th roots of unity, an FFT returns P(exactly k primes). (Checks: the probabilities sum to 1 and the mean is H/λ to 14 digits.) For every n of the windows (sliding intervals):

| window | H | P(0): data / law / Poisson | P(2): data / law / Poisson | distance to the data: law / Poisson | block misses (cells with P ≥ 10⁻³) |
|---|---|---|---|---|---|
| train | 30 | 0.1353 / 0.1358 / 0.1858 | 0.3173 / 0.3180 / 0.2632 | 0.0024 / 0.0822 | 7 of 7 |
| A | 30 | 0.2303 / 0.2302 / 0.2718 | 0.2645 / 0.2655 / 0.2306 | 0.0012 / 0.0638 | 3 of 6 |
| B | 30 | 0.3022 / 0.3019 / 0.3377 | 0.2201 / 0.2206 / 0.1990 | 0.0006 / 0.0551 | 2 of 6 |
| C | 30 | 0.3912 / 0.3914 / 0.4195 | 0.1672 / 0.1678 / 0.1583 | 0.0008 / 0.0429 | 0 of 5 |
| train | 60 | 0.0166 / 0.0167 / 0.0345 | 0.1895 / 0.1892 / 0.1956 | 0.0016 / 0.0931 | 3 of 9 |
| A | 60 | 0.0504 / 0.0502 / 0.0739 | 0.2718 / 0.2724 / 0.2507 | 0.0013 / 0.0705 | 2 of 8 |
| B | 60 | 0.0884 / 0.0881 / 0.1140 | 0.2971 / 0.2976 / 0.2688 | 0.0006 / 0.0540 | 0 of 8 |
| C | 60 | 0.1497 / 0.1501 / 0.1760 | 0.2927 / 0.2930 / 0.2656 | 0.0010 / 0.0447 | 0 of 7 |

The primes are far from Poisson at these sizes (total-variation distance 0.04–0.09), and the law removes 97–99 % of that distance. But the pre-registered block test fails in many cells (17 of 56), mostly at small λ: the sliding counts average millions of overlapping intervals, so deviations of 10⁻³ are significant, and on the unit circle |1 − z| reaches 2, which doubles the effective expansion parameter. At train with H = 30, the step from order 3 to order 4 is of the same size as the residuals (about 10⁻³ per cell, [r20_trunc.json](writeup/data/study/r20_trunc.json)), so a higher order is the plausible missing piece; that is not shown.

**Round 21: the exact law for H = 30** ([r21.py](writeup/code/study/r21.py) → [r21.json](writeup/data/study/r21.json), [round21.txt](writeup/data/study/preregister/round21.txt)). For an unconditioned interval every set with points of both parities has S = 0, so only subsets of the 15 even or the 15 odd positions contribute, 2 × 2¹⁵ sets. With the λ-free sums T_j = Σ_{|A|=j} S(A), the exact Hardy–Littlewood distribution is P(k) = Σ_j C(j, k)(−1)^(j−k) T_j/λ^j, no expansion. T_j vanishes from j = 9 on: no admissible 9-tuple fits in 30 consecutive integers (the narrowest prime 9-tuples span 31), and the singular series knows it. Data against the exact law (P(k), data / law):

| window | k = 0 | 1 | 2 | 3 | 4 | 5 | block misses |
|---|---|---|---|---|---|---|---|
| train | 0.1353 / 0.1354 | 0.3234 / 0.3233 | 0.3173 / 0.3172 | 0.1652 / 0.1654 | 0.0494 / 0.0494 | 0.0085 / 0.0085 | 0 of 6 |
| A | 0.2303 / 0.2300 | 0.3841 / 0.3843 | 0.2645 / 0.2649 | 0.0977 / 0.0974 | 0.0207 / 0.0207 | 0.0025 / 0.0025 | 0 of 6 |
| B | 0.3022 / 0.3018 | 0.4006 / 0.4009 | 0.2201 / 0.2203 | 0.0648 / 0.0648 | 0.0110 / 0.0110 | 0.0011 / 0.0011 | 1 of 6 |
| C | 0.3912 / 0.3914 | 0.3984 / 0.3977 | 0.1672 / 0.1676 | 0.0378 / 0.0379 | 0.0050 / 0.0050 | 0.0004 / 0.0004 | 0 of 5 |

One miss in 23 cells against 0.3 expected, where order 4 had 12 in 24: the misses of round 20 at H = 30 were the truncation of the expansion. With the exact sums, the share of the 30-number intervals near 10⁷–10⁸ holding exactly k primes agrees with the Hardy–Littlewood conjecture to the fourth decimal for every k, while Poisson is off by up to 0.05.

**Long intervals.** In a window not used before, [10⁹, 2·10⁹), the exact pair sum of section 4 (no change) against the variance for H = 10⁴, 10⁵, 10⁶; after a miss at H = 10⁴ a replication in [3·10⁹, 4·10⁹) was pre-registered and run:

| H | [10⁹, 2·10⁹): measured / law, t₉ | [3·10⁹, 4·10⁹) (replication): measured / law, t₉ |
|---|---|---|
| 10,000 | 0.5007 / 0.4967, 3.41 | 0.5189 / 0.5165, 1.40 |
| 100,000 | 0.3859 / 0.3875, -0.35 | 0.4014 / 0.4116, -1.58 |
| 1,000,000 | 0.2662 / 0.2784, -1.07 | 0.2898 / 0.3068, -1.08 |

H = 10⁵ and 10⁶ fit in both windows. H = 10⁴ misses in the first window (+0.8 %, p = 0.008) and fits in the replication (+0.46 %, t₉ = 1.40), with the same sign. **Open point:** over the variance tests with 500 ≤ H ≤ 10⁴ (sections 6 and 10) the measured variance lies above the pair-sum law in 9 of 10 cases, by up to 1.4 % (about a 1 % chance for independent fair signs; but the values for H = 500 and 5000 come from the same primes, and the range was chosen after the signs had been seen, so this figure is post hoc); at H = 10⁵ and 10⁶ all four values lie below it (−0.4, −4.4, −2.5, −5.5 %), each within its noise. A relative excess of about 10⁻⁵ in the pair counts, summed over h < H, would produce it, so the variance at long H is a very sensitive probe of the pair counts; whether this is a lower-order term of the Hardy–Littlewood pair law or chance is not settled here.

**Round 12: the excess in pair counts** (diagnostic; [r12.py](writeup/code/study/r12.py) → [r12.json](writeup/data/study/r12.json), with errors per block in [r12b.py](writeup/code/study/r12b.py) → [r12b.json](writeup/data/study/r12b.json)). For intervals on the grid, Σ_j c_j² = Σ_j c_j + 2 Σ_h P_emp(h) exactly, with P_emp(h) the number of prime pairs at distance h inside one interval (counted by FFT autocorrelation), so the variance excess is, by this identity, an excess in the pair counts. The mean count lies below H/ln x (relative −2.8·10⁻⁵ in [10⁹, 2·10⁹), the familiar π(x) < li(x)); scaling the law's pairs by (1 + δ/e)² for that (a heuristic correction: pairs scale with the density squared), the remaining relative pair excess over all distances below 10⁴ is 8.0·10⁻⁶ ± 2.5·10⁻⁶ (t₉ = 3.3, ten blocks of 10⁸) in [10⁹, 2·10⁹) and 5.3·10⁻⁶ ± 3.8·10⁻⁶ (t₉ = 1.4) in [3·10⁹, 4·10⁹). How it splits by distance is not resolved: the part up to h = 3000 is −4.1·10⁻⁵ ± 6.8·10⁻⁵. Without the density correction the pairs lie below the law. The cause is not identified.

## 11. In the plane: nearest neighbours of Gaussian primes ([r11.py](writeup/code/study/r11.py), [r11_exact.py](writeup/code/study/r11_exact.py), [round11.txt](writeup/data/study/preregister/round11.txt), [verdicts_r11_dd.json](writeup/data/study/verdicts_r11_dd.json), [verdicts_r11_rep_dd.json](writeup/data/study/verdicts_r11_rep_dd.json), [r11_exact.json](writeup/data/study/r11_exact.json))

The Gaussian primes a + bi (a² + b² prime) are a point process in the plane. Hardy–Littlewood for ℤ[i] uses the singular series S_G(B) = ∏_π (1 − ν_π(B)/N(π))/(1 − 1/N(π))^|B| over the Gaussian primes π up to units: 1 + i (residue x + y mod 2), the two primes over each p ≡ 1 (mod 4) (residues x + y·i₀ and x − y·i₀ mod p, i₀² ≡ −1), and the inert p ≡ 3 (mod 4) (residues (x mod p, y mod p), norm p²). By the prime ideal theorem a lattice point of norm N is a Gaussian prime with chance (4/π)/ln N (0.05363 predicted, 0.05366 counted in the first box). The nearest neighbour of a Gaussian prime has squared distance at least r when the disc of lattice points 0 < |w|² < r around it holds no Gaussian prime: the avoidance problem of section 3 with a disc as interior, order 4 while the live disc has at most 40 points. Checks: the pair constant of {0, 1 + i} is 1.67659 (1.67665 from a brute product over Gaussian primes of norm ≤ 3000); and for the discs with norm below 8, 10 and 16 (8, 12 and 20 live points) the exact inclusion–exclusion over all 2^20 subsets agrees with order 4 to 2.4·10⁻⁴.

Four boxes, primality by gmpy2, every Gaussian prime's nearest neighbour searched to norm 82; G3 and G4 were pre-registered as replications after G1 and G2. **Correction (second review):** the first analysis counted every Gaussian prime of the boxes G1, G2 and G3 twice. These boxes are squares placed symmetrically about the diagonal a = b, and b + ai is a Gaussian prime with the same norm and the same nearest-neighbour distance whenever a + bi is; the strips of the block test then shared their data through the mirror, which made the blocks agree too well and inflated t₉ (it had reported two misses, d = 10 in G1 at t₉ = −8.1 and in G2 at 3.9). The analysis below keeps one point of every mirror pair inside a box and forms the 10 blocks as equal-count runs along b, with λ for each block from the mean ln N of its own points (107,322, 20,224, 101,580 and 201,599 independent Gaussian primes in G1–G4, of 214,644, 40,448, 203,160 and 201,599 in the boxes):

| statistic | G1: data / law, t₉ | G2: data / law, t₉ | G3: data / law, t₉ | G4: data / law, t₉ |
|---|---|---|---|---|
| d = 2 | 0.3246 / 0.3243, 0.2 | 0.2519 / 0.2505, 0.5 | 0.3089 / 0.3083, 0.3 | 0.3066 / 0.3059, 0.7 |
| d = 4 | 0.2265 / 0.2243, 1.5 | 0.1941 / 0.1911, 0.7 | 0.2172 / 0.218, -0.6 | 0.2189 / 0.217, 2.2 |
| d = 8 | 0.1514 / 0.1516, -0.2 | 0.1377 / 0.1437, -2.0 | 0.1503 / 0.151, -0.4 | 0.1502 / 0.1509, -0.8 |
| d = 10 | 0.1965 / 0.1998, -2.8 | 0.2321 / 0.2276, 1.2 | 0.2096 / 0.2073, 1.0 | 0.2073 / 0.2084, -1.5 |
| d = 16 | 0.03426 / 0.0342, 0.1 | 0.04747 / 0.04875, -0.7 | 0.03702 / 0.03738, -0.9 | 0.03716 / 0.03786, -1.4 |
| d = 18 | 0.0247 / 0.02485, -0.4 | 0.04158 / 0.04014, 0.9 | 0.02759 / 0.02795, -0.6 | 0.029 / 0.02844, 1.2 |
| d = 20 | 0.02839 / 0.02781, 1.2 | 0.0538 / 0.05469, -0.4 | 0.03229 / 0.0327, -0.6 | 0.03286 / 0.03349, -2.5 |
| d ≥ 26 | 0.01352 / 0.01321, 0.9 | 0.04134 / 0.04346, -1.5 | 0.01708 / 0.01734, -0.6 | 0.01794 / 0.01805, -0.3 |
| mean d | 6.693 / 6.698, -0.3 | 8.584 / 8.659, -1.3 | 7.031 / 7.045, -0.6 | 7.075 / 7.1, -1.8 |

All 36 tests fit. Neither the d = 10 deviations nor any other distance stands out once each Gaussian prime is counted once. The mean nearest-neighbour norm, which depends on the whole law, fits in all four boxes.

**Round 16: the hexagonal lattice** ([r16.py](writeup/code/study/r16.py), [round16.txt](writeup/data/study/preregister/round16.txt), [verdicts_r16_dd.json](writeup/data/study/verdicts_r16_dd.json), [verdicts_r16_rep_dd.json](writeup/data/study/verdicts_r16_rep_dd.json)). The Eisenstein integers a + bω (ω = e^(2πi/3), norm a² − ab + b²) are the other lattice with a ring structure. Their primes: 1 − ω (norm 3; residue a + b mod 3), the two primes over each p ≡ 1 (mod 3) (ω ↦ a root of x² + x + 1 mod p), and the inert p ≡ 2 (mod 3) (norm p²). Six units and covolume √3/2 give the density 3√3/(π ln N) per lattice point. Here no parity rule removes neighbours: a unit neighbour is allowed but suppressed (the series of {0, 1} is 0.6112, brute product 0.6112), a neighbour at norm 3 in the same class mod 1 − ω is favoured (1.2224), and a triangle 0, 1, 1 + ω is impossible. Four boxes; E3 and E4 were pre-registered as replications after E2. The reflection σ(a + bω) = a + (a − b)ω (that is −ω² times the conjugate) fixes the line a = 2b and maps E1, E2 and E3 largely onto themselves, so as for the Gaussian boxes one point of every mirror pair is kept (205,710, 40,111, 39,319 and 54,044 independent Eisenstein primes, of 274,217, 53,468, 52,466 and 54,044 in the boxes):

| statistic | E1: data / law, t₉ | E2: data / law, t₉ | E3: data / law, t₉ | E4: data / law, t₉ |
|---|---|---|---|---|
| d = 1 | 0.235 / 0.2372, -2.1 | 0.1898 / 0.1875, 1.6 | 0.183 / 0.1843, -0.5 | 0.1902 / 0.189, 0.6 |
| d = 3 | 0.3287 / 0.3268, 1.5 | 0.2817 / 0.2816, 0.1 | 0.2792 / 0.2783, 0.3 | 0.2841 / 0.2831, 0.4 |
| d = 4 | 0.1502 / 0.1509, -1.2 | 0.1472 / 0.1462, 0.4 | 0.1434 / 0.1455, -1.5 | 0.1435 / 0.1465, -2.0 |
| d = 7 | 0.1414 / 0.1417, -0.3 | 0.1488 / 0.1559, -3.8 | 0.1609 / 0.1564, 1.8 | 0.1569 / 0.1556, 0.8 |
| d = 9 | 0.06296 / 0.06252, 0.7 | 0.08222 / 0.08045, 1.3 | 0.08139 / 0.08151, -0.1 | 0.08051 / 0.07997, 0.3 |
| d = 12 | 0.0478 / 0.04692, 1.7 | 0.07115 / 0.07135, -0.3 | 0.06994 / 0.07306, -2.3 | 0.07037 / 0.07059, -0.2 |
| d = 13 | 0.01642 / 0.01595, 1.3 | 0.03052 / 0.02931, 1.3 | 0.03113 / 0.03037, 1.0 | 0.02935 / 0.02885, 0.6 |
| d ≥ 19 | 0.0112 / 0.01148, -1.4 | 0.03353 / 0.03412, -0.5 | 0.03716 / 0.03654, 0.5 | 0.03186 / 0.03309, -1.2 |
| mean d | 4.5 / 4.489, 1.4 | 5.629 / 5.634, -0.2 | 5.73 / 5.73, -0.0 | 5.575 / 5.591, -0.5 |

35 of the 36 tests fit. The miss is d = 7 in E2 (−4.5 %, t₉ = −3.8), which was the reason for the replication: E1 and E2 lie along a = 2b, the direction of 1 − ω, a mirror axis of the lattice. In E3 (same axis) and E4 (off the axis) d = 7 fits, so by the pre-registered rule it is counted as chance. Over both lattices the block test misses once in 72 tests against 1.1 expected by chance.

## 12. The residue chain, cell by cell ([r13.py](writeup/code/study/r13.py), [r13_null.py](writeup/code/study/r13_null.py), [r13b.py](writeup/code/study/r13b.py), [r13c.py](writeup/code/study/r13c.py), [round13.txt](writeup/data/study/preregister/round13.txt))

For each residue a of p and each class c of the next gap, the conditioned law of section 9 gives P_a(c); the table of counts N[a, c] is compared with N[a, ·] P_a(c) by a Pearson χ² (cells with expectation below 5 pooled per row; the pass rule |χ²/dof − 1| ≤ 3√(2/dof)):

| window | gaps | mod 30: χ²/dof (dof 79) | mod 210: χ²/dof (dof 401) | no residue: χ²/dof | mod 2310: χ²/dof (dof 3624) |
|---|---|---|---|---|---|
| train | 5,096,875 | 0.902 | 0.691 | 7117 | 0.810 |
| A | 2,171,133 | 0.702 | 0.857 | 2680 | 0.892 |
| B | 1,808,832 | 0.947 | 0.995 | 2025 | 0.944 |
| C | 578,836 | 0.833 | 0.894 | 565 | 0.925 |

- The residue carries the structure: without it χ²/dof is in the hundreds to tens of thousands, with it about 1.
- Mod 30 passes everywhere. Mod 210 passes at A, B and C and fails at train, on the low side: χ²/dof = 0.691, where 1.000 ± 0.058 is the value for counts drawn independently from the law ([r13_null.json](writeup/data/study/r13_null.json)). All the values lie at or below 1: **the counts of the primes scatter less than independent draws**, and the law is not the cause, since it is calibrated.
- Rejected explanations: drift of the law inside the window (a check run before the addenda, not pre-registered; each block drawn from its own law: 1.001 ± 0.052); and, each pre-registered before its run, stratification by p mod 11 (each class mod 2310 drawn from its own law with the actual counts, aggregated to mod 210: 1.022 ± 0.055, [r13b.json](writeup/data/study/r13b.json)); flow conservation in the chain (the mod-210 chain itself simulated: 1.014 ± 0.056, [r13c.json](writeup/data/study/r13c.json)). The data are under-dispersed against the mod-2310 law as well (last column).
- By window length inside the train window: 0.691 for W = 9·10⁷, 0.888 ± 0.022 for W = 9·10⁶ and 0.880 ± 0.023 for W = 9·10⁵. (the sub-windows are the first 9·10⁶ or 9·10⁵ numbers of each tenth of the window). The rigidity hypothesis was stated as "χ²/dof climbs to 1 as W/x falls", but that was not derived: Montgomery–Soundararajan-type rigidity scales roughly with ln W/λ, which changes only from 16.0 to 13.7 between these two sub-window lengths, so the plateau near 0.88 does not contradict it. A law-free check by the second reviewer ([review2/s13.py](writeup/code/study/review2/s13.py) → [s13.json](writeup/data/study/review2/s13.json)) splits the train window into two halves made of interleaved chunks and compares the two tables per row: chunks of 10³, 10⁴, 10⁵ and 10⁶ give 0.98, 1.10, 0.90 and 0.84 (random halves: 1.02 ± 0.07). The cell counts behave like multinomial ones at short range; the under-dispersion needs correlations over 10⁵ numbers or more, as long-range rigidity would give.
- **Open point.** The gap classes of the primes, counted per residue class, are more regular than independent, stratified or Markov draws from a calibrated law, by about 10 % at every window length and by 30 % over the whole of [10⁷, 10⁸). A likely source is the long-range rigidity of counts that the Hardy–Littlewood correlations imply (section 4: counts of primes in long intervals scatter far less than Poisson), applied to counts of patterns; deriving its size needs four-point sums and is not done here. A third mechanism was tested (addendum 3, [r13d.py](writeup/code/study/r13d.py) → [r13d.json](writeup/data/study/r13d.json)): if the gaps after the primes of one residue class had to fill an almost fixed share of the window, every row would lose one degree of freedom (expected 0.88); but the squared projection of each row's residual on the gap-length direction averages 0.66, 0.92, 0.75 and 0.74 in the four windows (about 1 under independence, 0 under such a constraint), so a full constraint of this kind is excluded (a partial one, about 25 %, is not: each value has a noise of about ±0.2 with 48 rows); removing that component leaves χ²/dof at 0.70 at train.

**Where the low χ² sits** (diagnostics, not pre-registered; [r13e.py](writeup/code/study/r13e.py) → [r13e.json](writeup/data/study/r13e.json), [r13f.py](writeup/code/study/r13f.py) → [r13f.json](writeup/data/study/r13f.json)). By gap class, the mean squared standardised residual over the residue rows is 0.29, 0.31 and 0.31 for gaps 2, 4 and 6 at train (0.29, 0.31, 0.61 at A) and rises towards 0.9 for the long gaps (gaps ≥ 40: 0.88). For gap 2 every admissible residue has the same law, so that class compares the twin counts per residue class with the prime counts per class. Both are spread far less than Poisson counts: across the 48 classes mod 210 the prime counts have Var/mean 0.08 / 0.39 / 0.36 / 0.45 (train / A / B / C), the twin counts across their 15 classes 0.40 / 0.18 / 1.00 / 0.53, the cousin counts (gap 4) 0.25 / 0.49 / 0.50 / 0.73. For the primes this is classical: by the explicit formula (under GRH) the counts in progressions deviate by an amount of order √x/(√φ(q) ln x), times a factor that grows slowly with q, against √(x/(φ(q) ln x)) for independent counts: a variance ratio of order 1/ln x (about 0.05 near 10⁸). The same evenness of the twin and cousin counts across residue classes is what the multinomial null does not have; the Hardy–Littlewood law gives the means, not this rigidity of their fluctuations. A derivation would need the explicit-formula side for prime pairs, which does not exist in proved form.

**Round 15: three consecutive gaps** ([r15.py](writeup/code/study/r15.py) → [r15.json](writeup/data/study/r15.json), [round15.txt](writeup/data/study/preregister/round15.txt)). Four consecutive primes, P(0, h₁, h₁ + h₂, h₁ + h₂ + h₃) for all 729 triples of gaps ≤ 18 (order 4), a direct test of the four-point correlations, cell by cell:

| window | triples with all gaps ≤ 18: data / law | χ²/dof (dof) | z | control (multinomial from the law) |
|---|---|---|---|---|
| train | 0.2773 / 0.2775 | 0.992 (404) | -0.1 | 1.031 ± 0.044 |
| A | 0.1586 / 0.1582 | 1.038 (404) | 0.5 | 0.995 ± 0.082 |
| B | 0.1044 / 0.1041 | 1.089 (404) | 1.3 | 0.998 ± 0.072 |
| C | 0.0606 / 0.0606 | 0.985 (404) | -0.2 | 0.999 ± 0.054 |

All four fit, and the three-gap table shows no under-dispersion, not even at train. (The two-gap table of section 6 has χ²/cell 0.849 at train and lies below 1 in three of the four windows, so a weaker effect there is not excluded.)

## 13. Sophie Germain primes: two forms with different slopes ([r14.py](writeup/code/study/r14.py), [round14.txt](writeup/data/study/preregister/round14.txt), [verdicts_r14.json](writeup/data/study/verdicts_r14.json))

A Sophie Germain prime p has 2p + 1 prime too: two linear forms, x and 2x + 1. Write the form x + b as the token 2b and the form 2x + 2m + 1 as the token 2m + 1; for an odd prime q both have the root −t/2 (mod q), so the number of roots mod q is the number of distinct tokens mod q, and only q = 2 is special (only the forms x + b have roots there). Checks: the series of {x, 2x + 1} is 1.3203236 = 2C₂, as it must be; p and p + 6 both Sophie Germain have 5.2713 (5.2725 from a brute product to 3000); p and p + 3 or p + 15 get 0. Consecutive Sophie Germain primes at distance h are the pattern of four forms at 0 and h with interior events {x + m, 2x + 2m + 1}, and the Bahadur form of the twin primes applies with two scales, ln x and ln 2x. The windows hold 367,108 / 120,953 / 84,510 / 21,704 Sophie Germain primes (train: exactly π_SG(10⁸) − π_SG(10⁷) = 423,140 − 56,032).

| statistic | SG primes (train) | law, order 1 / 3 | t₉ train | t₉ A | t₉ B | t₉ C |
|---|---|---|---|---|---|---|
| gap 6 | 0.012432 | 0.012393 / 0.012393 | 0.2 | 2.4 | 0.5 | -0.2 |
| gap 12 | 0.020389 | 0.020547 / 0.020547 | -0.7 | -0.8 | -0.6 | 0.0 |
| gap 18 | 0.025802 | 0.026059 / 0.026054 | -1.2 | -0.4 | 3.0 | -0.1 |
| gap 24 | 0.013043 | 0.013025 / 0.013025 | 0.0 | 0.6 | 0.1 | 0.6 |
| gap 30 | 0.028368 | 0.028487 / 0.028464 | -0.3 | 1.6 | -1.2 | -0.4 |
| gap 36 | 0.0092153 | 0.0092781 / 0.0092726 | -0.4 | 0.3 | -0.1 | -1.1 |
| gap 42 | 0.032511 | 0.032539 / 0.032511 | -0.0 | -1.1 | 0.4 | -1.1 |
| gap 48 | 0.018629 | 0.018443 / 0.018398 | 1.4 | 0.3 | -1.7 | -0.7 |
| gap 60 | 0.04123 | 0.040976 / 0.040859 | 1.0 | -1.5 | 0.1 | -1.1 |
| gap 90 | 0.021452 | 0.021876 / 0.02174 | -2.0 | 0.2 | 1.3 | 0.3 |
| gap 120 | 0.019847 | 0.019817 / 0.019667 | 0.7 | -0.2 | 0.3 | 0.4 |
| gap ≤ 300 | 0.70876 | 0.71689 / 0.70876 | -0.2 | 0.5 | 1.6 | -1.7 |
| same class mod 5 (gaps ≤ 300) | 0.30345 | 0.30359 / 0.30278 | 0.5 | -0.7 | 0.6 | -0.9 |
| same class mod 7 (gaps ≤ 300) | 0.16589 | 0.16605 / 0.16571 | 0.1 | -0.7 | 3.3 | 0.3 |

55 of the 56 tests fit (one miss, same class mod 7 at B, t₉ = 3.3; 0.8 expected by chance). Consecutive Sophie Germain primes avoid repeating their class mod 5: a Sophie Germain prime above 5 is ≡ 1, 3 or 4 (mod 5) (never 2, since then 5 | 2p + 1), so without bias the next one would repeat the class a third of the time; among the gaps ≤ 300 it does so in 30.3 % (law: 30.3 %) at train, as twin primes do. For gaps 6 to 30 the interior has at most four live events, and the exact inclusion–exclusion of the second review ([review2/s14.py](writeup/code/study/review2/s14.py) → [s14.json](writeup/data/study/review2/s14.json)) agrees with order 3 to seven digits; for long gaps (up to 49 live events at h = 300) no truncation error was estimated, so the five-digit agreement of the share of gaps ≤ 300 at train is partly luck.

**Round 18: Euler's polynomial n² + n + 41** ([r18.py](writeup/code/study/r18.py), [stats_r18.py](writeup/code/study/stats_r18.py), [round18.txt](writeup/data/study/preregister/round18.txt), [verdicts_euler41.json](writeup/data/study/verdicts_euler41.json)). The most prime-rich quadratic: (−163/p) = −1 for every prime p < 41, so no value is divisible by a prime below 41, and its Bateman–Horn constant is 6.6398 (product to 2·10⁷). The machinery of section 8 carries over with w_p = 1 + (−163/p) and shifted roots coinciding exactly when p | d² + 163; two neighbouring n both giving primes have the series 44.64 (a brute product to 2000: 44.72). In F1 = [10⁶, 2·10⁶), F2 = [10⁸, 10⁸ + 5·10⁵) and F3 = [10¹⁰, 10¹⁰ + 10⁵) the values are prime for 233,707, 89,970 and 14,375 of the n:

| statistic | F1: data / law | t₉, z_bin | F2: data / law | t₉, z_bin | F3: data / law | t₉, z_bin |
|---|---|---|---|---|---|---|
| n-gap 1 | 0.2363 / 0.2368 | -0.6, -0.5 | 0.1821 / 0.1824 | -0.4, -0.3 | 0.1412 / 0.146 | -1.3, -1.6 |
| n-gap 2 | 0.1753 / 0.1763 | -1.9, -1.2 | 0.1472 / 0.1457 | 1.4, 1.3 | 0.121 / 0.1218 | -0.3, -0.3 |
| n-gap 3 | 0.1383 / 0.1384 | -0.4, -0.2 | 0.1214 / 0.1222 | -0.6, -0.7 | 0.108 / 0.1066 | 0.4, 0.6 |
| n-gap 4 | 0.1037 / 0.1034 | 0.7, 0.4 | 0.09767 / 0.09782 | -0.2, -0.2 | 0.09267 / 0.08916 | 2.1, 1.5 |
| n-gap 5 | 0.08119 / 0.08114 | 0.1, 0.1 | 0.08169 / 0.08194 | -0.3, -0.3 | 0.07945 / 0.07787 | 0.9, 0.7 |
| n-gap 6 | 0.06148 / 0.06088 | 1.2, 1.2 | 0.06451 / 0.06582 | -1.9, -1.6 | 0.06519 / 0.06533 | -0.1, -0.1 |
| n-gap 8 | 0.03675 / 0.03588 | 2.4, 2.3 | 0.04472 / 0.04434 | 0.6, 0.6 | 0.04411 / 0.04792 | -1.9, -2.1 |
| n-gap 10 | 0.02106 / 0.02118 | -0.4, -0.4 | 0.03042 / 0.02991 | 0.8, 0.9 | 0.03618 / 0.03518 | 0.6, 0.6 |
| n-gap ≥ 15 | 0.02399 / 0.02425 | -1.2, -0.8 | 0.06241 / 0.06237 | 0.0, 0.0 | 0.1135 / 0.1137 | -0.2, -0.1 |
| same n mod 3 | 0.2499 / 0.2502 | -0.4, -0.3 | 0.2681 / 0.2703 | -1.4, -1.4 | 0.2856 / 0.2834 | 0.5, 0.6 |
| same n mod 5 | 0.1098 / 0.1101 | -0.4, -0.4 | 0.1301 / 0.1297 | 0.2, 0.3 | 0.1467 / 0.1434 | 1.3, 1.1 |
| same n mod 7 | 0.05666 / 0.05637 | 0.7, 0.6 | 0.07267 / 0.07317 | -0.6, -0.6 | 0.0871 / 0.08566 | 0.7, 0.6 |
| same n mod 41 | 8.558e-06 / 4.676e-06 | 0.7, 0.9 | 7.78e-05 / 6.079e-05 | 0.7, 0.7 | 0.000487 / 0.000287 | 1.3, 1.4 |
| same n mod 43 | 0 / 2.526e-06 | -17.5, -0.8 | 4.446e-05 / 3.863e-05 | 0.3, 0.3 | 6.957e-05 / 0.0002018 | -1.9, -1.1 |
| same n mod 47 | 0 / 8.084e-07 | -15.8, -0.4 | 3.334e-05 / 1.661e-05 | 0.7, 1.2 | 0 / 0.0001045 | (degenerate), -1.2 |

Every test with an expected count of at least 5 in its window fits (37 of them). The pre-registration included repeats of the class mod 41, 43 and 47; consecutive n in the same class mod 43 or 47 need a gap that is a multiple of 43 or 47, which the law expects about once in a whole window (the data have none), so the block test is meaningless there (all-zero blocks, marked degenerate); the binomial z stays within ±1.4.

## 14. What the seven lines turn out to be

**The scoreboard** (every pre-registered test of the law, from the verdict files; a miss is |t₉| > 3 for block tests, |z| > 3 for χ² tests; the chance of a miss per test is 0.015 and 0.0027):

| round | what | tests | misses | expected by chance |
|---|---|---|---|---|
| 2-3 | gap law, classes, variance, information, tails, CPAP-4 (4 windows; registered test) | 96 | 4 | 1.4 |
| 4 | consecutive twin primes | 68 | 0 | 1.0 |
| 5 | consecutive n with n² + 1 prime | 45 | 3 | 0.7 |
| 6 | second differences | 20 | 0 | 0.3 |
| 9 | variance and means along quadratics (fresh data) | 4 | 0 | 0.0 |
| 10 | third cumulant of interval counts | 8 | 0 | 0.1 |
| 10 | variance in long intervals (2 windows) | 6 | 1 | 0.1 |
| 11 | Gaussian primes, nearest neighbours (4 boxes) | 36 | 0 | 0.5 |
| 13 | residue chain mod 30 and 210, cell by cell (χ² too low counts as a miss) | 8 | 1 | 0.0 |
| 14 | consecutive Sophie Germain primes | 56 | 1 | 0.8 |
| 15 | joint law of three gaps, 729 cells (4 windows) | 4 | 0 | 0.0 |
| 16 | Eisenstein primes, nearest neighbours (4 boxes) | 36 | 1 | 0.5 |
| 17 | the law at 10¹⁸ | 18 | 1 | 0.3 |
| 18 | consecutive n with n² + n + 41 prime (tests with expected count ≥ 5) | 37 | 0 | 0.6 |
| 20 | count distribution in short intervals (H = 30, 60; cells with P ≥ 10⁻³) | 56 | 17 | 0.8 |
| 21 | the same for H = 30 with the exact Hardy–Littlewood sums | 23 | 1 | 0.3 |
| | **all** | **521** | **30** | **7.6** |

30 misses in 521 tests against 7.6 expected. Seventeen are cells of the count distribution with the order-4 expansion (round 20); for H = 30 the exact sums of round 21 show that these were the truncation (12 misses become 1). Four more have a known cause: three are the information estimators before their bias was measured from the law (round 3), and one is the gap tail at train, which the fifth order brings inside the test (round 8). The other 9 against 7.6 expected are not a clear excess: the law is right to about 0.1 % where the counts are in the millions and to a few per cent in the smallest samples, across all these objects and ten orders of magnitude. (Before the second review the count was 15 of 405, with 4 misses in the two-dimensional boxes; 3 of those came from counting every point twice in boxes symmetric about a mirror line, see section 11.)


- **Numerical, statistical, algebraic, dynamic, information-theoretic**: five views of one object, the law of runs of consecutive primes. Its only input is the singular series; its only scale is ln x. The textbook versions miss because they truncate it: the Poisson factor e^(−h/ln x) ignores that the numbers between two primes are correlated with them (round 1: up to 30 %), exponentiating the Ursell functions of distinct points mistreats X² = X (step 3), and the Montgomery–Soundararajan asymptotics ignore the grid (H = 50).
- **Geometric**: Ulam's diagonals are the same law for polynomials: the means are Bateman–Horn and the scatter is set by the Bateman–Horn constants of polynomial pairs. In the plane, the Gaussian primes follow the same expansion with the singular series of ℤ[i] (section 11).
- **Computational**: the pre-registered windows and the block test. Four windows over eight orders of magnitude, one of them opened only after the law was frozen.
- **What is not new**: the Hardy–Littlewood and Bateman–Horn conjectures, Gallagher's Poisson law, the Lemke Oliver–Soundararajan bias and the Montgomery–Soundararajan variance are literature; that a parameter-free Hardy–Littlewood count fits gap counts has been used before (Brent 1974 for small gaps; Odlyzko, Rubinstein and Wolf 1999 for jumping champions). Hardy–Littlewood inclusion–exclusion for differences of consecutive primes, with lower-order terms, has been studied (for example by Goldston and Ledoan, and by Funkhouser, Goldston and Ledoan for large gaps); those papers were not read for this study. What this study adds, as far as the files in this archive show, is the Bahadur form and its accuracy (0.1–0.3 % for gaps up to 30 at four scales, which is the size of its own truncation error at train), the grid term in the variance, the variance law along quadratics, and the same expansion for twin primes, for n² + 1 and for Gaussian primes; the literature for these was not searched beyond this archive.

## Validation

- [engine.py](writeup/code/study/engine.py) and [stats_r3.py](writeup/code/study/stats_r3.py) compute the statistics of the primes from an exact segmented sieve of each window (checked: π(10⁶) = 78,498; the train count 5,096,876 = π(10⁸) − π(10⁷)).
- [theory_r2.py](writeup/code/study/theory_r2.py), [theory_r3.py](writeup/code/study/theory_r3.py) and [theory_geo.py](writeup/code/study/theory_geo.py) compute every prediction from the singular series alone; [report2.py](writeup/code/study/report2.py), [report3.py](writeup/code/study/report3.py) and [bias_r3.py](writeup/code/study/bias_r3.py) make the verdicts; [tables.py](writeup/code/study/tables.py) prints the tables above.
- [validate.py](writeup/code/study/validate.py) checks that every table in this file is exactly what tables.py prints from the data files, that 23 numbers quoted in the text match the data, and three identities (S(0, 2) = 2C₂; the residue-conditioned series average to S; S_f({0}) against the Bateman–Horn constant of n² + 1).
- An independent review with its own code (nothing imported from the study scripts; [review/](writeup/code/study/review/), outputs in [data/study/review/](writeup/data/study/review/)) recomputed: the block statistics of window A and C from primes listed by gmpy2.next_prime (all 32 statistics agree exactly), the singular series, the grid variance (0.70556 / 0.77451 / 0.81207), the pair constants of the quadratics (C₁/C² = 0.51751 for 4n² + 1; χ²/curve 0.407 at train), the exact inclusion–exclusion against the Bahadur orders (section 3), the twin-gap shares, the counts of n² + 1 primes and the residue identity. Its corrections are in the text: the discriminant b² − 16c, the error of order 2 at gap 6, the description of version 1, the truncation error, the steps added after seeing data, the Markov-gain misses, the block misses of the n² + 1 round, the wording of section 9, the duplicate curves and the literature.
- A second independent review ([review2/](writeup/code/study/review2/), outputs in [data/study/review2/](writeup/data/study/review2/)) recomputed sections 10–13 and found that the symmetric two-dimensional boxes counted every point twice (corrected in section 11: the Gaussian misses disappear), that the split of the round-12 pair excess by distance is not resolved (only the total with its error is now given), a post hoc probability, overstated wording on rounds 13 and 15, and that the repository cannot by itself prove the pre-registrations preceded their computations (section 1). It confirmed the third cumulant (0.31187 at H = 50 against the study's 0.3119), the long-interval variances, the Gaussian series and exact sums, the Sophie Germain series and exact sums, and the data counts.
