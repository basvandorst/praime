# Mersenne exponents: a rule that lands only on them

> **In short:** The Mersenne exponents are the primes p with 2^p − 1 prime: 2, 3, 5, 7, 13, 17, 19, 31, … Asked: a rule that steps from one to the next and touches nothing else. 1. Trial factors 2kq + 1 remove about half of the primes, but among the rest the share of exponents grows only like the logarithm of the factor bound: 1.8% → 4.1% below 10⁴, as the Wagstaff heuristic predicts (23.7 expected, 22 found). 2. In 110 tests, no arithmetic feature of the 51 exponents below 83,880,077 beats that heuristic. 3. So each step must prove its prime: the next q whose 2^q − 1 survives trial factoring and passes Lucas–Lehmer. The page has it as transformation 142, as a list without a circle, and as a test of the known exponents alone.

**Page:** [prime.html](writeup/code/riemann/prime.html): transformation 142 "Next Mersenne exponent" (indicator "Mersenne exponents" counts the hits), and the panel "Mersenne exponents, one after another" at the bottom (Start, and "Only the known exponents").

**Code:** [code/mersenne/](writeup/code/mersenne/), rerun by [run_all.sh](writeup/code/mersenne/run_all.sh).

**Status:** the algorithm is proved (Lucas–Lehmer is a theorem; trial factors are exact). Hit rates, feature tests and timings are computed. The probability model is heuristic.

## 1. The question

132 (Golomb's formula) gives the next prime, and from 1 its chain hits the exponents 2, 3, 5, 7 and then misses 11, 23 and 29: 7 of its first 10 steps are exponents. The aim is 100%: every step a Mersenne exponent, none skipped.

## 2. A filter is not enough (computed; heuristic model)

Every factor of 2^q − 1 has the form 2kq + 1 with 2kq + 1 ≡ ±1 (mod 8) (theorem). Removing the primes with such a factor for k ≤ K ([tf.py](writeup/code/mersenne/tf.py), [rate.py](writeup/code/mersenne/rate.py), output [rate.txt](writeup/code/mersenne/rate.txt); deeper K from [touch/tftable.out](writeup/code/mersenne/touch/tftable.out)):

| K | primes < 10⁴ left | share that are exponents (22) | primes < 10⁵ left | share (28) |
|---|---|---|---|---|
| 0 | 1,229 | 1.8% | 9,592 | 0.29% |
| 100 | 805 | 2.7% | 6,781 | 0.41% |
| 10⁴ | 596 | 3.7% | 5,211 | 0.54% |
| 10⁵ | 532 | 4.1% | 4,665 | 0.60% |
| 10⁷ | 428 | 5.1% | | |

The heuristic P(2^q − 1 prime | no factor below Q) ≈ e^γ·log₂(Q)/q predicts 23.7 exponents below 10⁴ (true 22) and 30.0 below 10⁵ (true 28). The share grows like log Q: reaching 50% needs no factor below about 2^(0.28·q), which is the work of a primality proof ([touch/bound.out](writeup/code/mersenne/touch/bound.out)).

## 3. No hidden structure (computed)

All 51 exponents, against all 4,882,467 primes below 83,880,077 (GIMPS: every exponent below it tested and verified, milestone page read on 15 November 2026). The null is P = e^γ·log₂(a·p)/p, a = 2 for p ≡ 3 and a = 6 for p ≡ 1 (mod 4). Output in [structure/results.txt](writeup/code/mersenne/structure/results.txt), from [analysis.py](writeup/code/mersenne/structure/analysis.py); an independent recheck (gmpy2 primes, plain integers, its own Poisson-binomial code) in [check.py](writeup/code/mersenne/structure/check.py) → [check.txt](writeup/code/mersenne/structure/check.txt) reproduces every count, expectation and p-value.

- Calibration: 51 observed, 50.3 expected (a = 2/6), 47.6 (a = 2). By mod 4: 31 with p ≡ 1 (25.2 expected), 19 with p ≡ 3 (24.1). The data lean towards a = 6 (Bayes factor 1.9) without deciding it.
- 110 tests: p mod 3, 5, 7, 8, 12, 16, 24, 30, 60, 120; the forms 2^k ± 1 and 4^k ± 3; closeness to a power of 2; binary weight and runs; 2p + 1, p + 2, p − 2 prime; the spacing of log log p. Bonferroni threshold 4.5·10⁻⁴; the smallest p-value is 0.023, and 49% of simulated null sets have a smaller one.
- That smallest p-value belongs to a control that is a theorem: p ≡ 3 (mod 4) with 2p + 1 prime (p > 3) never gives an exponent (Euler: 2p + 1 divides 2^p − 1); 0 found, 3.74 expected. A class with no exponents at all becomes visible only once it would hold more than about 8 of them.
- The one feature found in the exponents below 10⁵ (2^k ± 1 or 4^k ± 3: 10 against 7.0) lies entirely at p ≤ 127 and has nothing to test above 10⁵.

## 4. How few tests per exponent (computed)

[touch/](writeup/code/mersenne/touch/): trial factoring to 10⁷ ([tfdeep.py](writeup/code/mersenne/touch/tfdeep.py)), P−1 run for real with gmpy2 ([pm1.py](writeup/code/mersenne/touch/pm1.py), tables [pm1table_1e4.out](writeup/code/mersenne/touch/pm1table_1e4.out), [pm1table_1e5.out](writeup/code/mersenne/touch/pm1table_1e5.out)), orderings ([order.out](writeup/code/mersenne/touch/order.out)).

- Best share of exponents among the Lucas–Lehmer tests: 22 of 350 (6.3%) below 10⁴ (trial factoring to 10⁷ plus P−1 to B1 = 10⁵), 28 of 3,931 (0.71%) below 10⁵.
- P−1 (B1 = 1000, B2 = 3000) costs more than the tests it saves below about 4·10⁴, and 0.42 and 0.30 Lucas–Lehmer tests per removal in [4·10⁴, 7·10⁴) and [7·10⁴, 10⁵).
- Ordering the candidates by their posterior probability is the ascending order (0 positions differ); ordering by the prior with the mod-4 factor finds the first 25 exponents in 895 tests instead of 946 but needs 3,557 instead of 3,392 for all 28.
- The only proofs that 2^q − 1 is composite are a factor or a full Lucas–Lehmer (or PRP) test of about q squarings. A base-2 Fermat test never rejects 2^q − 1 (theorem: 2^q ≡ 1 mod 2^q − 1 and q divides 2^q − 2). So a rule that touches only exponents would need a test nobody has; a complete search tests every prime that survives the factors.

## 5. The algorithm on the page (proved; times measured)

Next exponent after p: the next prime q with no factor 2kq + 1 for k ≤ 100,000·(q/7000)² (at least 1,000; sieved by the primes below 2^14, exact double arithmetic below 2^38), and Lucas–Lehmer s₀ = 4, s ← s² − 2 mod 2^q − 1, s_(q−2) = 0. Below q = 4,200 in BigInt; above, the squaring is an irrational-base FFT (Crandall–Fagin, as in GIMPS) in Float64Array, with the largest rounding error checked at every step (a step above 0.35 is redone at twice the length). Every zero found up to 50,000 (all of them in transformation 142 and in the list) is recomputed in exact BigInt. The readable source is [engine/best.js](writeup/code/mersenne/engine/best.js); the page holds the same functions without comments.

| q | BigInt, ns per squaring | FFT | ratio |
|---|---|---|---|
| 4,001 | 5,025 | 5,236 | 0.96 |
| 10,007 | 19,799 | 10,953 | 1.81 |
| 20,011 | 59,705 | 22,642 | 2.64 |
| 30,011 | 106,020 | 34,213 | 3.10 |

On the page (headless Chrome, one worker; a computer with more cores runs several tests at once):

- **142**: every arrow lands on an exponent; destinations up to 150,000. Only the primes inside the circle are tested; above 3,000 in Web Workers, and the page redraws when they finish. Circles of 100 numbers: from 19,900 in 4.6 s, from 44,450 in 26.7 s, from 86,200 in 94 s. The chain from 1 on a circle of 1,000: 14 of 14 steps are exponents.
- **One after another**: the same rule with no circle and no cap: 9,941 after 22.5 s, 11,213 after 32 s, 19,937 after 170 s.
- **Only the known exponents**: Lucas–Lehmer on the 38 known exponents up to 6,972,593 and on nothing else: 31 of them in 2.6 minutes. One test on one core takes 13 s at 86,243 (Node, a full test) and, from the measured time per squaring, about 1.4 minutes at 216,091, 19 at 756,839, 82 at 1,398,269, 6 hours at 2,976,221 and 2 days at 6,972,593. This touches only exponents and proves each one, but the list is its input: it cannot find a new one.

## Validation

- [ll_check.py](writeup/code/mersenne/ll_check.py) → [ll_check.txt](writeup/code/mersenne/ll_check.txt): Lucas–Lehmer (gmpy2) on all 1,229 primes below 10⁴ gives exactly the 22 exponents.
- [engine/page_check.js](writeup/code/mersenne/engine/page_check.js) and [page_check.py](writeup/code/mersenne/engine/page_check.py) → [page_check.txt](writeup/code/mersenne/engine/page_check.txt), on the functions cut out of the page itself: the chain of 142 from 1 is exactly the 23 exponents up to 11,213 and every start 2 … 11,212 gives the next one; the test on every prime below 12,000 finds exactly them, and all 759 factors it reports divide 2^q − 1; the FFT residues of 18 primes from 4,201 to 30,011 equal those of gmpy2 (largest rounding error 0.0215); the zero test of the verifier is right in 46 of 46 cases.
- [engine/test_best.js](writeup/code/mersenne/engine/test_best.js) → [test_best.txt](writeup/code/mersenne/engine/test_best.txt): the readable engine, the same checks plus 36 FFT residues against BigInt.
- [structure/check.py](writeup/code/mersenne/structure/check.py) rechecks section 3; [touch/verify.py](writeup/code/mersenne/touch/verify.py) → [verify.out](writeup/code/mersenne/touch/verify.out): Lucas–Lehmer confirms the 28 exponents below 10⁵, and below 20,000 every prime is either removed by a proved factor or tested (761 tests), leaving exactly the listed exponents. [tfcheck.py](writeup/code/mersenne/touch/tfcheck.py) checks the factor tables against a plain search.
- Above 20,000 the completeness of the list (no exponent missed between the known ones) rests on GIMPS, not on these runs; the page's list and 142 recompute it as far as they run.
