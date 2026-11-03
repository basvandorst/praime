# Twin primes: where the hidden part sits, and how to look at it

Status: plan (3 November 2026). Nothing below is computed yet.

## The frame

What the gap night proved for single primes (FINDINGS, 31 Oct 2026):

- Everything a sieve can see is proved arithmetic and fixes the whole gap structure: the next prime is the J-th survivor of the sieve to y = x^0.4027, J geometric.
- The one unproved piece is a fair coin: a survivor of the sieve to y is prime with chance rho, independent of its neighbours to 1e-3 (tested 1e9 to 1e18).
- That coin is exactly Selberg's parity problem (FACTS 2017): the sieve cannot tell a prime from a product of two primes above y.

Twin primes are the same statement with two coins. Hardy-Littlewood (FACTS 4004) says the two coins for n and n + 2 are independent. The twin prime conjecture is one consequence of that independence, and no sieve can prove it (FACTS 2017, 10402: the best sieve bound is 4 times the truth, now below 3.4 by Wu 2004).

So the realistic targets are not a proof. They are:

1. Say precisely how much of twin structure is proved arithmetic and how much is one independence statement.
2. Measure that independence as sharply as the box allows (1e-3 to 1e-4, up to 1e18).
3. Show where the hidden information lives in three different bases: sieve depth y, major arcs q, zeta zeros.
4. Turn the frame into checkable predictions and new facts (FACTS.md rows with a check), plus one record in the user's style (problem 44).

Everything stays inside the rules: reasoning first, bounded runs (minutes, nice, timeout, 2 cores, 3 GB), every number from a second script.

## Facts we build on

| Fact | Use |
|---|---|
| 1109, 10102 | twins p ≡ 5 (mod 6); lesser member mod 30 in {11, 17, 29}, last digit 1, 7, 9 |
| 2002, 9505 | Mertens: survivors of the sieve to y are prime with chance e^γ ln y / ln x |
| 2017, 10402 | parity problem; sieve bound 8C₂ x/ln²x against the conjectured 2C₂ |
| 3014, 4004, 4006, 2045 | twin counts to 1e11 (ours), 2e16 (Nicely), 1e19 (Chaffin 2026); 2C₂ fits to 5e-5 |
| 3026, 10303, 10406 | maximal twin gaps, consecutive twin pairs (k = 9 at 1.7e14), twin chains |
| 10403 | twins first exceed Hardy-Littlewood at 1369391 (Wolf) |
| 2010, 2004 | pair correlation of the zeros (we measured F(α) = 1.01 for 1.2 ≤ α ≤ 2.5 from 5000 zeros); explicit formula |
| 2021, 4015, 3008 | Chen (p + 2 prime or P₂), Zhang / Maynard / Polymath (gap ≤ 246; ≤ 12 under EH, ≤ 6 under GEH) |
| 10605 | problem 44: f(x) ± 1 twin for x = 0..15 (Andersen) |
| FINDINGS 31 Oct | the coin model, its code (facts_coin405*.py, facts_maxgap_is.py, facts_pattern_model.py, facts_var_crossings.py), the Chebyshev-bias fit from 289 zeros |

## The angles

| # | Angle | Question it answers | Builds on | Compute (bounded) | Cost |
|---|---|---|---|---|---|
| 1 | The twin sieve and its coin | Is the twin coin the square of the prime coin? Where is its unbiased crossing c₂? | coin model, Buchstab | exact Φ₂(x, y) counts at 1e9, 1e11; 2-dimensional Buchstab function from its delay equation; ρ₂ = π₂ / Φ₂ | 20 min |
| 2 | Independence of the two coins | Is corr(prime(n), prime(n + 2) \| both survive) zero to 1e-3, 1e-4? | facts_coin405 | twin-survivor pairs near 1e9, 1e12 (exact), 1e18 (BPSW windows) | 30 min |
| 3 | Parity, measured | How is the sieve's 4x split among (p, p), (p, P₂), (P₂, p), (P₂, P₂)? Is Chen's share what independence predicts? | 2017, 4015 | census of survivor pair types by y, with factor counts | 15 min |
| 4 | Second-order Hardy-Littlewood over h | Do the pair counts Σ Λ(n)Λ(n + h) deviate from S(h) x by pure √x noise, or by a systematic term? | 3014, 1604 | all even h ≤ 1000 at x = 1e8, 1e9, 1e10 (Λ weights); deviation / √x, correlation across h | 30 min |
| 5 | The zeros in the twin error | Is π₂(x) − 2C₂ li₂(x) correlated with ψ(x) − x (the zeros), or binomial noise? | Chebyshev fit, 2004 | twin error on [1e6, 2e9] against −2S(2) Σ x^ρ/ρ from 289 zeros; the 1369391 crossing | 30 min |
| 6 | Twin races | Do twin classes mod 30, 210 show a Chebyshev bias? | 1602-1606 | race curves of the 3 (mod 30) and 15 (mod 210) classes to 1.4e11; compare with the prime race | 15 min |
| 7 | The twin constant inside the zeros | Can 2C₂ (the singular series) be read off the pair correlation of the zeros alone? | 2010, Bogomolny-Keating 1996, Goldston-Montgomery 1987 | F(α) and R₂ from Odlyzko's 2e6 zeros; fit the lower-order term, compare with the Bogomolny-Keating / Conrey-Snaith formula | 1-2 h incl. download |
| 8 | Major and minor arcs | How much of π₂ comes from q ≤ Q, and how big is the minor-arc remainder really? | circle method | S(α) = Σ Λ(n) e(nα) to N = 1e7-1e8 by FFT; truncated singular series S_Q(2); the minor-arc integral against N/Q and √N | 20 min |
| 9 | Long-range rigidity of twin counts | Is var/mean of twin counts in windows H sub-Poisson like primes (0.226 at H = 1e6), and from which scale? | facts_var_crossings, Montgomery-Soundararajan | windows H = 1e2..1e7 on 1e8 numbers near 1e9 and 1e11; coin model vs periodic model | 20 min |
| 10 | Forecasts: records of twin structure | Where are the next maximal twin gap, the least 10 consecutive twin pairs, the least 7-twin chain? | facts_maxgap_is, facts_pattern_model | out-of-sample quantiles; calibrate on 3026 and 10303 | 30 min |
| 11 | Twin generators (problem 44) | Can class-number reasoning beat 16 values x with f(x) ± 1 both prime? | 1230 method, 5005, Rabinowitsch | f = 6g; choose g so both discriminants of 6g ± 1 are non-residues mod 5..23; bounded search in that class only | 30 min |
| 12 | Bits of a twin | How many bits does "p + 2 prime" carry after the sieve to y, and does any feature beyond p mod q# help? | facts_gaprules, bits finding | conditional entropy by y; backtest 20 predictors of gap = 2 near 1e9, 1e12, trained before, tested after | 20 min |
| 13 | The barrier map | Exactly which ingredient do twins lack that a² + b⁴, x³ + 2y³, digit-free primes and bounded gaps had? | 5006, 3008 | literature only: arXiv abstracts, no computation; one page | 1 h reading |
| 14 | Brun's constant with an error bar | What does Chaffin's 1e19 count give for B₂ with the tail from angle 4's measured deviations? | 4005, 4022, 2045, conj-71 | tail integral plus second-order term; compare Nicely, Sebah-Demichel | 10 min |

## The angles in detail

### 1. The twin sieve and its coin

- Φ₂(x, y) = #{n ≤ x : gcd(n(n + 2), P(y)) = 1}, a sieve of dimension 2. Its product is V₂(y) = ½ ∏_{2<q≤y} (1 − 2/q) = 2C₂ e^{−2γ} / ln² y (Mertens squared times C₂).
- Drill-down chain, each link checked on exact data before the next:
  1. Φ₂(x, y) / (x V₂(y)) as a function of u = ln x / ln y is a 2-dimensional Buchstab function (sieve theory, dimension κ = 2). Count it exactly at x = 1e9, 1e11 for u = 2..8 and solve the delay-differential equation numerically. Proved layer.
  2. The twin coin ρ₂(y) = (twin density) / (twin-survivor density). At y = √x the survivors are the twins themselves, so Hardy-Littlewood says the sieve to √x leaves (2e^{−γ})² = 1.261 times too many: ρ₂ = e^{2γ}/4 = 0.7930, the square of the prime coin e^γ/2 = 0.8905. For large u, ρ₂ → (e^γ/u)² = ρ². Measure ρ₂ against ρ² at every y.
  3. Find the crossing c₂ where e^{2γ} Γ(3) ω₂(1/c₂) = 1 (the twin coin is unbiased). Either c₂ = 0.4027 again or a new number. Theory picks it, no fit.
- Output: a FACTS row for Φ₂ and ω₂, the measured ρ₂ table, c₂.

### 2. Independence of the two coins

- Along the twin survivors of the sieve to y = x^0.4027, record for each survivor whether n and n + 2 are prime. Test P(both) / (P(n prime) P(n + 2 prime)) and lags 1..10 along the survivor sequence, as facts_coin405.py did for single primes.
- Hardy-Littlewood is the statement that this ratio is 1. The only known correlation is the pair factor from primes q in (y, √x] that cannot hit both n and n + 2 at once (none, since q > 2 cannot divide both). So the prediction is exactly 1 at order 1/(y ln y).
- Scale: 1e9 and 1e12 exact (2.7e6 survivors as before), 1e18 with BPSW in four windows (as facts_coin405_1e18.py).
- Output: independence to 1e-3 (1e-4 at 1e9 if the sample allows). A deviation would be the first structural fact beyond the sieve; its absence is the sharp form of "the twin information is one fair coin".

### 3. Parity, measured

- For survivors of the sieve to y, classify (n, n + 2) by (prime, prime), (prime, P₂), (P₂, prime), (P₂, P₂), other. Chen's theorem proves the union of the first three columns is infinite; no theorem separates the first.
- Independent coins predict each cell from ρ and the P₂ share of survivors (Buchstab again). Compare at 1e9, 1e12.
- Output: the exact split of the sieve's factor 4 into cells. Shows in numbers why Chen stops where he stops.

### 4. Second-order Hardy-Littlewood over h

- D(h, x) = Σ_{n≤x} Λ(n)Λ(n + h) − S(h) x for every even h ≤ 1000, where S(h) = 2C₂ ∏_{q|h, q>2} (q − 1)/(q − 2).
- Questions: is D / √x of order 1 (random) and is corr(D(h), D(h')) across h zero? Do twins (h = 2) deviate like every other h? The consecutive-prime bias (1604) is a second-order ln ln x / ln x effect that must vanish here, since these are all pairs. Check that it does.
- Output: a table of D / √x by h and x; a FACTS row if the deviations are correlated across h (that would point to the zeros, angle 5).

### 5. The zeros in the twin error

- Two models for E₂(x) = π₂(x) − 2C₂ ∫₂ˣ dt / ln² t:
  - coin: binomial noise of size √(π₂(x)) ≈ √x / ln x, uncorrelated with ψ(x) − x;
  - zeros: E₂ inherits −2 S(2) Re Σ_ρ x^ρ/ρ (with the Λ-to-π conversion), hence corr with ψ(x) − x near 1.
- We already fit Chebyshev's bias to 2e9 with corr 0.98 from 289 zeros; the same machinery tests this. Also: does the model predict the first crossing at 1369391 and the later ones (Wolf 2011)?
- Output: which model, with the correlation and the residual size. Either answer is a fact worth a row.

### 6. Twin races

- Lesser twin p ≡ 11, 17, 29 (mod 30): three classes; mod 210: 15. Hardy-Littlewood gives equal shares. The Chebyshev bias of primes comes from prime squares inside ψ(x; q, a), and twins have no square term. So the coin model predicts no bias at all.
- Measure the three races to 1.4e11 (the chains sieve) and the leader statistics; compare with the prime race mod 30 over the same range.
- Output: either the first measured bias in twin classes or the clean statement that twins are bias-free where primes are not.

### 7. The twin constant inside the zeros

- Bogomolny-Keating (1995, 1996) derived the GUE pair correlation of the zeros, including its lower-order terms, from the Hardy-Littlewood conjecture for Σ Λ(n)Λ(n + h). Goldston-Montgomery (1987) showed the pair correlation conjecture is equivalent to the variance of primes in short intervals. Conrey-Snaith (2007) give the same lower-order terms from the ratios conjecture.
- So the singular series, and with it 2C₂, is encoded in the pair correlation of the zeros. We have F(α) = 1.01 for 1.2 ≤ α ≤ 2.5 from 5000 zeros. With Odlyzko's 2e6 zeros (height ~1e6) the lower-order term is a few percent and measurable.
- Plan: compute R₂(r) and F(α) finely; fit the Bogomolny-Keating correction; read the arithmetic factor; compare with 2C₂ and with the singular series sum Σ_h S(h) e(hr) that the formula contains. Then repeat on a high block (height 1e21, Odlyzko) where the correction is smaller and the fit sharper in relative terms.
- Output: "the twin prime constant is readable from the zeros to k digits". This is the angle closest to the user's phrase: the twin information hidden in another object.

### 8. Major and minor arcs

- R(N) = Σ_{n≤N} Λ(n)Λ(n + 2) = ∫₀¹ |S(α)|² e(2α) dα with S(α) = Σ Λ(n) e(nα). The major arcs q ≤ Q give N · S_Q(2), S_Q the singular series truncated at Q; the tail Σ_{q>Q} is of order 1/Q, so the minor arcs must carry N/Q plus the true error.
- Compute S(α) on a grid of 2^24 points by FFT (N = 1e7), the major-arc contribution for Q = 2..1000, and the minor-arc integral directly. Compare the minor-arc remainder with N/Q and with √N.
- This is the sieve decomposition (angle 1) in the Fourier basis: q ≤ Q major arcs correspond to sieving by primes ≤ Q. The remainder is the coin again, now as "square-root cancellation on the minor arcs", which is what no one can prove.
- Output: a plot-free table: Q, major-arc share, remainder. States where the information hides in this basis.

### 9. Long-range rigidity of twin counts

- For primes, var/mean of counts in windows H fell from 0.71 (H = 100) to 0.226 (H = 1e6), and the periodic-deletion model explained it. Measure the same for twin pairs on 1e8 numbers near 1e9 and 1e11, H = 1e2..1e7.
- Prediction from the two layers: Poisson at H ≪ ln² x, rigidity from periodic deletion above, with a logarithmic decline. No table of this exists in the literature we hold.
- Output: a FACTS row with the table and the model comparison.

### 10. Forecasts

- Maximal twin gaps (3026): run the importance-sampling model with the two-dimensional sieve and check the known records to 1e19 (A113274) as out-of-sample tests; then quantiles for the next record.
- Consecutive twin pairs (10303): k = 9 at 1.7e14; forecast the median for k = 10 with the pattern model (exact random sieve plus coin, calibrated on k = 3..8).
- Twin chains (10406): k = 6 at 7.4e13; forecast k = 7.
- Output: three predictions with quantiles; later searches by others will test them.

### 11. Twin generators (problem 44)

- Record: 16 consecutive x with f(x) ± 1 twin (Andersen, two quadratics). Reason first: f must be ≡ 0 (mod 6) on all x, so f = 6g. Both 6g − 1 and 6g + 1 must avoid small primes for all x: that is two discriminant conditions, (D₋/q) = (D₊/q) = −1 for q = 5, 7, 11, 13, ... Bateman-Horn (5005) gives the pair constant from the root counts N(q); maximise it over the class before searching.
- Search only inside the class with both discriminants non-residues mod 5..19; expect counts from the Bateman-Horn constant and stop when the model says 17 is out of reach.
- Output: a record or a calibrated statement why 16 is near the ceiling for degree 2; both are submittable.

### 12. Bits of a twin

- Near 1e9, the next gap carries 4.66 bits and the sieve to 1000 leaves 1.65. For "gap = 2 or not", compute the conditional entropy after p mod 30, 210, 2310, 30030 and after the sieve to y; compare with the coin's h(ρ₂) / ρ₂ form.
- Backtest predictors of gap = 2: p mod q#, the last K twin gaps, digit features, Chebyshev-type features. Train below, test above, as facts_gaprules.c did. Expectation: nothing beats the sieve to y; past twin gaps leak only p mod q#.
- Output: the number of hidden bits per twin, and the null result that no non-sieve feature helps.

### 13. The barrier map

- One page, literature only, arXiv abstracts read directly: what let Friedlander-Iwaniec (a² + b⁴), Heath-Brown (x³ + 2y³), Maynard (digits) and Zhang / Maynard / Polymath (bounded gaps) pass the parity barrier, and why none of it reaches n, n + 2: no bilinear (Type II) structure, no multiplicative decomposition, bounded gaps need k ≥ 50 primes to force one pair, EH gives 12, GEH 6, nothing short of Hardy-Littlewood gives 2. Note the unrefereed 2026 claims of 240 as such.
- Output: writeup/twins-barrier.md. Tells us which of angles 1-12 touch the barrier (2, 3, 7, 8) and which only describe the proved layer.

### 14. Brun's constant with an error bar

- B₂(1e19) from Chaffin's count plus the Hardy-Littlewood tail ∫ 2C₂ (2/t) dt/ln² t; the uncertainty is the tail's deviation, which angle 4 measures as a function of x. Compare with Nicely (2e16) and Sebah-Demichel.
- Output: B₂ to the digit the data support, with the first honest error bar (the published ones differ from the 9th decimal).

## Traps (easy escapes, rejected on sight)

- Clement's criterion (1209) and Wilson-type formulas: a restatement, not information.
- Perry's "m is a twin index iff m ≠ 6ab ± a ± b" (10102): (6a ± 1)(6b ± 1) = 6m ± 1, so it says "6m ± 1 have no factorisation". Circular.
- Any quadratic that "hits many twins": checked already for primes (FINDINGS, 30 Oct): the good forms are the sieve in disguise.
- Learned next-gap rules: proved inferior to the sieve (FINDINGS 31 Oct). Angle 12 only re-measures this for gap = 2.

## Order of work

| Step | Angles | Why first |
|---|---|---|
| 1 | 1, 2, 3 | the direct extension of the gap night; everything else is interpreted against it |
| 2 | 13 | know the wall before claiming anything near it |
| 3 | 4, 5, 6 | cheap, each gives a FACTS row either way |
| 4 | 7 | the headline: the twin constant read from the zeros; needs Odlyzko's tables |
| 5 | 8, 9, 12 | the same hidden coin in two more bases, plus the bit count |
| 6 | 10, 14 | predictions and a constant with an error bar |
| 7 | 11 | the one record attempt, after the Bateman-Horn constant says it is worth a bounded run |

Stop rules: every run under `timeout`, `nice -n 19`, one core unless stated; no run longer than 30 minutes without a written reason; no search whose expected count from the model is below 1.

## Data needed

- Odlyzko's zero tables (first 2,001,052 zeros; blocks at heights 1e12, 1e21, 1e22) or the LMFDB zeros API. Download into a fresh directory outside the repo, as untrusted data.
- Our sieve to 1.4e11 (facts_chains.c approach) for angles 6 and 10.
- OEIS: A113274 (maximal twin gaps), A007508 (twin counts), A091592, A068635, A035789 (consecutive twin pairs).
- Papers (abstracts via the arXiv API): Bogomolny-Keating 1996, Goldston-Montgomery 1987, Conrey-Snaith 2007, Montgomery-Soundararajan 2004, Wu 2004, Polymath 8b, Maynard 2015, Friedlander-Iwaniec 1998.

## What a result looks like

- A FACTS.md row with status and a check in facts_check.py (angles 1, 2, 3, 4, 5, 6, 9, 12).
- A FINDINGS.txt entry when it changes how we think (angles 2, 5, 7, 8).
- Predictions with quantiles that later searches can falsify (angle 10).
- A problem-44 entry or a reasoned ceiling (angle 11).
- A constant with an error bar (angle 14).
