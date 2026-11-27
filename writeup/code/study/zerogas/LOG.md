# Lab log: the zeros of zeta as a gas (short loops)

Each entry: hypothesis with a predicted number written BEFORE the run, the test, the result, the lesson.
Nothing here is a claim until it is promoted (registration, independent check, literature check).

## Loop 1 (theory only, no data). Constants of the sine process (bulk of GUE), by Fredholm determinants ([sine.py](sine.py)).
- Result: entropy of one spacing 0.533031 nats = 0.769001 bits; correlation of consecutive spacings -0.30837;
  information between consecutive spacings 0.051468 nats = 0.074252 bits. Checks: spacing variance 0.179994 (known 0.180).
- For comparison, the prime gaps (STUDY.md): one gap carries 0.26965 bits about the next; here 0.0743 bits.

## Loop 2. Lag correlations of spacings: theory against Odlyzko's zeros near 10^12, 10^21, 10^22 ([loop2.py](loop2.py)).
- Prediction (before the run): the correlation of spacings k apart follows the sine process: lag 1 = -0.308; lags 2..6
  from the variances of the k-th neighbour distance (computed in the same script before the data are read); every
  data value within 3 sigma (sigma about 0.01 for 10^4 zeros; the three files are independent).
- Prediction 2: the partial sums 1 + 2 sum_{j<=k} r_j tend to 0 (the gas is rigid: spectral density 0 at frequency 0).
- Result: theory r_1..r_7 = -0.30837, -0.07720, -0.03311, -0.01827, -0.01157, -0.00798, -0.00584 (variances of the
  k-th neighbour distance 0.17999, 0.24898, 0.29017: the known GUE values). Data: 20 of 21 values within 2 sigma
  (largest |z| 2.31, lag 1 at 10^12: -0.331 against -0.308). Both predictions hold; 1 + 2 sum r_j falls to 0.075 at lag 7 (tail like -0.281/k^2).
- Lesson: at these heights 10^4 zeros cannot see anything beyond the sine process in the spacings (noise 0.01). The
  primes must be looked for where they are large: in the slow part (the counting function), not in single spacings.

## Loop 3. The hidden variables of the zero gas are the phases t ln p mod 2 pi of the small primes ([loop3.py](loop3.py)).
- Reasoning: the prime gas has hidden variables n mod q (STUDY section 15: knowing them predicts the next prime). By
  the explicit formula the zero gas should have the dual ones: the phases t ln p of small primes. S(t) = N(t) - smooth
  is, heuristically, P_y(t) = -(1/pi) sum_{p^m <= y} sin(m t ln p)/(m p^(m/2)) plus a remainder.
- Prediction (before the run), zeros 10^12+1..10^12+10^4, S taken at the midpoints between zeros (exact, from the index):
  (a) regression slope of S on P_y is 1.0 within 0.1 for y = 10, 100, 1000, 10^4;
  (b) the variance removed, Var(S) - Var(S - P_y), equals Var(P_y) = (1/(2 pi^2)) sum_{p^m<=y} 1/(m^2 p^m) within 15 %;
  (c) with y = 10^4 more than half of the variance of S is explained.
- Result (zeros near 10^12): Var S at midpoints 0.1059. Explained: y = 10: 59.0 %, 100: 83.1 %, 1000: 91.3 %, 10^4: 93.5 %,
  10^5: 93.8 %. Slopes 0.980, 0.952, 0.919, 0.896, 0.891. (a) holds to y = 1000, fails narrowly at 10^4 (0.896);
  (b) holds for y = 10, 100 (4 %, 10 %), fails for y >= 1000 (18 %, 25 %); (c) holds.
- Lesson: the four primes 2, 3, 5, 7 fix 59 % of the counting fluctuation at midpoints. The slope falls with y: large
  primes enter with a weight below 1. Var(S_mid) is far below the uniform-t variance: sampling at midpoints filters.

## Loop 4. A universal response curve: how strongly a prime moves the zeros ([loop4.py](loop4.py)).
- Reasoning: let W be the prime wave. If the zeros sit where x + W(x) = n - 1/2 (x = unfolded height), then at the
  midpoint of two zeros S_mid = (W(x_n) + W(x_{n+1}))/2: the wave averaged over the two ends of the gap, not taken at the
  midpoint. For a prime with u = ln p / ln(t/2 pi) (cycles per mean spacing) the regression coefficient of S_mid on the
  wave of p AT THE MIDPOINT is then beta(u) = E[cos(pi u s)], s a spacing of the sine process: a universal curve with no
  parameter, 1 at u = 0, zero near u = 0.5, NEGATIVE beyond (the wave is sampled out of phase).
- Prediction (before the run): the binned regression coefficients at heights 10^12, 10^21, 10^22 follow beta(u) within
  3 sigma for u up to 0.6 (primes to 3 10^6 at 10^12), and the coefficient is negative in the bin u in (0.55, 0.6) at 10^12.
- Result (10^12, primes to 3 10^6, 13 bands of u): slope 0.985, 0.954, 0.898, 0.812, 0.696, 0.573, 0.429, 0.274, 0.097,
  -0.078, -0.176, -0.314, -0.308 against beta = 0.992, 0.962, 0.907, 0.825, 0.719, 0.600, 0.470, 0.334, 0.199, 0.070,
  -0.050, -0.156, -0.204. The shape and the sign change are there (data: zero near u = 0.46; model: 0.51; negative beyond,
  as predicted). The registered "within 3 sigma" FAILS: 8 of 13 bands are off by more than 3 sigma, always below the
  model, by 0.01 at small u and by 0.15 near u = 0.5.
- Lesson: the picture (the wave is read at the two ends of the gap) is right to first order; something lowers the
  response further as u grows. Candidates: (1) the zeros follow a wave of frequency u only partly (a factor c(u) < 1
  at the zeros themselves), (2) the spacing is correlated with the phase. Measure c(u) directly: next loop.

## Loop 5. Is the response curve universal in u? Same bands at 10^21 ([loop4.py](loop4.py) on the second file).
- Prediction (before the run): for u < 0.33 the slopes at 10^21 equal those at 10^12 in the same band within 3 sigma
  (0.985, 0.954, 0.898, 0.812, 0.696, 0.573, 0.429), although the primes involved are different (u = ln p / ln(t/2 pi)).

## Loop 6. The response at the zeros themselves ([loop6.py](loop6.py)).
- Exact identity: S_mid(n) = (S(x_n) + S(x_{n+1}))/2 with S(x_n) = n - 1/2 - x_n the mean of the two one-sided values.
- Prediction (before the run): the regression coefficient c(u) of S(x_n) on the wave of the band taken at the zero is 1
  within 0.02 for u < 0.2 and falls below 1 for larger u; and b(u) (loop 4) is c(u) E[cos(pi u s)] within 3 sigma.
- Result of loop 6 (10^12): c(u) = 0.992, 0.988, 0.982, 0.969, 0.943, 0.925, 0.914, 0.844, 0.838, 0.734, 0.672, 0.718 in the
  bands of loop 4. "c = 1 within 0.02 for u < 0.2" fails narrowly (0.969 +- 0.008 at u = 0.18). The product c(u) beta(u)
  reproduces b(u) for u <= 0.4 (0.985/0.985, 0.950/0.954, 0.891/0.898, 0.799/0.812, 0.678/0.696, 0.554/0.573,
  0.429/0.429, 0.282/0.274) and fails beyond (0.167/0.097, 0.051/-0.078, -0.033/-0.176, -0.112/-0.314).
- Lesson: two effects, both of the gas itself: the zeros follow a wave of frequency u only partly (c(u) < 1), and near
  u = 0.5 the spacing is tied to the phase. Both should be properties of the sine process: the wave of a band of primes
  IS the Fourier mode of the zero density at that frequency (explicit formula). So the whole curve should be reproduced
  by random unitary matrices with their own Fourier modes in place of the primes: loop 7.
- Result of loop 5 (10^21, primes to 2.5 10^6): slopes 0.991, 0.962, 0.898, 0.805, 0.712, 0.562, 0.458 (the last band only
  to u = 0.318) against 0.985, 0.954, 0.898, 0.812, 0.696, 0.573, 0.429 at 10^12. Prediction holds: same curve in u,
  with other primes (the third band has 7 primes at 10^12 and 116 at 10^21; the sixth 156 and 45,417).

## Loop 7. Random unitary matrices in place of zeta ([loop7.py](loop7.py)).
- Prediction (before the run): for eigenvalues of Haar unitary matrices (N = 400), with the Fourier modes of their own
  density in place of the primes, the two curves equal those of the zeros at 10^12 within 0.02 in every band:
  b(u) = 0.985, 0.954, 0.898, 0.812, 0.696, 0.573, 0.429, 0.274, 0.097, -0.078, -0.176, -0.314 and
  c(u) = 0.992, 0.988, 0.982, 0.969, 0.943, 0.925, 0.914, 0.844, 0.838, 0.734, 0.672, 0.718.
- Result of loop 7 (N = 400, 300 matrices): b(u) = 0.998, 0.966, 0.904, 0.810, 0.701, 0.576, 0.419, 0.262, 0.110, -0.046, -0.204,
  -0.316 and c(u) = 1.000, 0.996, 0.985, 0.966, 0.950, 0.928, 0.891, 0.858, 0.811, 0.770, 0.722, 0.655 (errors 0.001 to
  0.014). Against the zeros at 10^12: b within 0.02 in 10 of 12 bands (0.032 at u in [0.45, 0.5), 0.028 at [0.5, 0.55),
  both within 2 sigma of the zeta noise; the first bands differ by 0.012 because their mean u differs). c within the
  zeta noise except [0.55, 0.6) (0.718 +- 0.031 against 0.655). The registered "within 0.02 in every band" fails
  narrowly for b in two bands; the curves agree within the noise.
- The matrices also give the part not reachable at 10^12: b = -0.443, -0.513, -0.317, +0.120, +0.318 and
  c = 0.570, 0.410, 0.128, -0.168, -0.305 for u in [0.6, 0.7), [0.7, 0.8), [0.8, 0.9), [0.9, 1), [1, 1.2).
- Lesson: the response of the zeros to a prime is a property of the sine process alone. Beyond u = 1 the matrices still
  respond (their high modes are functions of the low ones); for zeta that region is the primes above t/2 pi.

## Loop 8. Primes above t/2 pi: low zeros, where they are small ([loop8.py](loop8.py); Odlyzko's first 100,000 zeros).
- Prediction (before the run), zeros number 60,001 to 100,000 (t from 47,000 to 75,000), all prime powers n with
  u = ln n / ln(t/2 pi) up to 1.2: b(u) and c(u) follow the matrix curves within 0.05 in every band; in particular
  b < 0 on [0.5, 0.9), b > 0 on [1, 1.2), c < 0 on [0.9, 1.2).
- Result of loop 8 (zeros 60,001..100,000, prime powers to u = 1.2; matrices in brackets):
  c = 0.980 (0.998), 0.954 (0.975), 0.906 (0.939), 0.849 (0.874), 0.781 (0.790), 0.678 (0.688), 0.567 (0.570), 0.415 (0.410),
  0.181 (0.128), -0.123 +- 0.088 (-0.168), then -0.632 +- 0.069 and -0.224 +- 0.047 on [1, 1.1), [1.1, 1.2) (matrices -0.305 on [1, 1.2));
  b = 0.947 (0.982), 0.851 (0.857), 0.641 (0.638), 0.368 (0.341), 0.054 (0.032), -0.275 (-0.260), -0.450 (-0.443), -0.501 (-0.513),
  -0.351 (-0.317), +0.063 +- 0.159 (+0.120), then +0.692 +- 0.133, +0.080 +- 0.055 (matrices +0.318 on [1, 1.2)).
  Within 0.05: holds for b in all bands to 0.9, for c in all but [0.8, 0.9) (0.053); the bands above 0.9 are noisy.
  Signs as predicted: b < 0 on [0.5, 0.9), b > 0 on [1, 1.2), c < 0 on [0.9, 1.2).
- Lesson: primes above t/2 pi do move the zeros, with the opposite sign. The low bands differ from the band averages of
  the matrices because one or three prime powers sit at one end of a band: compare at the exact u, not per band.

## Loop 9. A closed formula for c(u), derived from the three-point function of the sine process ([loop9.py](loop9.py)).
- Derivation (by hand; steps in the note to come). With s_q the Fourier coefficient of S at frequency u, sum_n S(x_n) W(x_n)
  = int S W dN = int S W dx - (1/2) int S^2 W' dx (because S dS = d(S^2)/2 at the jumps). The first term gives 1. The
  second is the bispectrum of the density: E[rho_u conj(rho_u1) conj(rho_u2)]/N = 1 - (1-|u|)+ - (1-|u1|)+ - (1-|u2|)+
  + 2 (1 - range{0, u, u1})+ for u = u1 + u2 (three-point cluster function 2 Re K12 K23 K31, K-hat the indicator of
  [-1/2, 1/2]). It vanishes when all three frequencies are below 1 and their signs agree (Diaconis-Shahshahani) and is
  u1 - 1 or u beyond. Integrating against 1/(u1 u2):
      c(u) = 1 - (1+u) ln(1+u) - (1-u) ln(1-u)                          (0 <= u <= 1),
      c(u) = 1 + 2 u^2 ln u - u(u-1) ln(u-1) - u(u+1) ln(u+1)           (1 <= u <= 2),
  zero at u = 0.90, c(1) = 1 - 2 ln 2 = -0.386.
- Prediction (before the run): matrices (N = 600, per mode, bands of width 0.05): |c - formula| < 3 sigma in every band
  for u < 1; the zeros at 10^12 (per band, formula averaged with the weights of the primes) within 3 sigma; the low
  zeros within 3 sigma for u < 0.9. The formula for u > 1 is less certain (I expect it may fail).
- Result of loop 9 (N = 600, 120 matrices, 30 bands of width 0.05 up to u = 1.5): the formula agrees in 29 bands within
  2.1 sigma (largest: z = +3.9 in [0.35, 0.40): 0.8722 +- 0.0042 against 0.8561; the errors come from 10 groups of 12
  matrices and may be too small). Also for u > 1, where I expected a failure: -0.323 against -0.325 on [1, 1.05).
  The registered "every band within 3 sigma for u < 1" fails in that one band.
- Lesson: the formula stands for random matrices on the whole range tested. Next: the zeros themselves at the exact u.

## Loop 10. The formula against the zeros ([loop10.py](loop10.py)).
- Prediction (before the run): c per band equals the formula (averaged over the prime powers of the band with their
  weights 1/(m^2 p^m)) within 3 sigma: at 10^12 (12 bands to u = 0.6), at 10^21 (7 bands to 0.33), at 10^22 (7 bands)
  and for the zeros 80,001..100,000 (bands of 0.1 to u = 1.5). Low zeros may deviate (arithmetic corrections).
- Result of loop 10. At 10^12 (12 bands): 9 bands within 2.3 sigma; the three lowest bands lie BELOW the formula by
  0.0062, 0.0149, 0.0119 (z = -3.3, -4.7, -4.0). Low zeros (80,001..100,000): bands above u = 0.5 within 3.2 sigma
  (z = -0.8, +2.9, +1.5, +3.1, +0.4, -3.2, -1.9, -0.1, +1.4, +1.1), the bands below far off: 0.9801 against 0.9944,
  0.9547 against 0.9797, 0.9055 against 0.9410, 0.8537 against 0.8777 (errors 0.001 to 0.003). The registered prediction
  fails for small u; the deficit shrinks with height.
- Lesson (reasoning): matrix modes are independent, prime powers are not: the wave of p^2 is phase-locked to the square
  of the wave of p. The derivation of loop 9 dropped exactly such triples because they vanish for matrices.

## Loop 11. The arithmetic correction, prime power by prime power ([loop11.py](loop11.py)).
- Derivation: for the wave of n = p^m (amplitude alpha_n = 1/(pi m p^(m/2)), frequency f_n = m ln p, L = ln(t/2 pi)) the
  exact resonances p^i p^(m-i) = p^m and p^(m+j)/p^j = p^m add to the numerator
      Delta = -(alpha_n f_n/8) [ sum_{i=1..m-1} alpha_{p^i} alpha_{p^(m-i)} - 2 sum_{j>=1} alpha_{p^(m+j)} alpha_{p^j} ],
  and the wave of n^2 adds 2 pi f_n alpha_{n^2}/L to the denominator:
      c_n = [ c(u_n) + 4 pi Delta/(L alpha_n^2) ] / [ 1 + 2 pi f_n alpha_{n^2}/L ].
  Leading terms: a prime loses ln p/(2 p L); a prime square loses 2 ln p / L (n = 4 at L = 9.3: 0.15).
- Prediction (before the run): for every prime power n <= 32 the measured coefficient equals c_n within 3 sigma or
  0.005, at 10^12 and for the zeros 80,001..100,000; in particular n = 4 at low height is near 0.83, not 0.98.
- Result of loop 11: all 36 coefficients (18 prime powers n <= 32, two heights) agree with c_n within 1.1 sigma.
  Low zeros (L = 9.28): n = 2: 0.9801 +- 0.0007 (c_n 0.9807; matrices alone 0.9944); n = 3: 0.9696 (0.9701); n = 4:
  0.8730 +- 0.0033 (0.8749; matrices 0.9777); n = 8: 0.7119 +- 0.0117 (0.7148); n = 9: 0.7509 (0.7577); n = 16: 0.5025 +- 0.020
  (0.5151); n = 25: 0.5662 (0.5770); n = 27: 0.4278 (0.4406); n = 32: 0.2498 +- 0.046 (0.2832).
  At 10^12 (L = 24.47): n = 2: 0.9940 +- 0.0031 (0.9938); n = 4: 0.9586 +- 0.018 (0.9572); n = 9: 0.9171 +- 0.015 (0.9206);
  n = 27: 0.8178 +- 0.033 (0.8177).
- Lesson: the response of the zeros to a prime power is the random-matrix curve c(u) plus an exact arithmetic term
  from the resonances among the powers of one prime. Two parameter-free formulas, derived by hand, confirmed in minutes.
- Reading (not yet tested): exact resonances a b = n and a/b = n give the arithmetic term; the NEAR resonances, pairs of
  prime powers with a close to n b, must give the matrix term 1 - c(u). That is a Hardy-Littlewood sum: the place where
  the law of STUDY.md (the primes as a gas) and the gas of zeros meet. Next loops: (1) derive 1 - c(u) from the
  Hardy-Littlewood pair law for (b, n b + h); (2) the same for the midpoint curve b(u); (3) literature check before any claim.

## Loop 12. How much of the position of a zero do the primes up to y fix? ([loop12.py](loop12.py))
- Reasoning: the dual of STUDY section 15 (knowing p mod q for q <= y predicts the next prime). Here: knowing the phases
  t ln p for p <= y. If the waves are orthogonal, the variance of S at the zeros explained by the prime powers n <= y is
  E(y) = sum_{n <= y} c_n^2 alpha_n^2 / 2, with c_n from loops 9 and 11 (no parameter).
- Prediction (before the run), zeros 80,001..100,000: the measured explained variance (regression on the waves in bands of
  u, cumulated) equals E(y) within 3 % at u = 0.2, 0.4, 0.6, 0.8, 1.0; and by u = 1.5 it reaches Var(S at zeros) within
  5 % (my guess: the zeros are driven linearly by the primes; this part may fail).
- Result of loop 12 (zeros 80,001..100,000, Var S at the zeros 0.10137): explained by n = 2 alone 24.9 %, by n <= 6: 53.4 %,
  n <= 16: 69.1 %, n <= 104: 86.8 %, n <= 10,833 (u < 1): 96.9 %. Data over theory: 1.036, 1.043, 1.049, 1.048, 1.047, 1.048,
  1.051, 1.059, 1.069, 1.073. Both registered readings FAIL (3 %; and the total stops near 97 %, then falls at u < 1.3).
- Lesson: two mistakes of mine. (1) The variance of a wave over the zeros is not alpha^2/2: by Landau-Gonek it is
  alpha_n^2/2 (1 + 2 pi f_n alpha_{n^2}/L); for n = 2 that factor is 1.037, the whole first-band excess. (2) Fitting 10^4
  waves to 2 10^4 zeros overfits; beyond u = 0.6 the "explained" variance is partly noise. Do not fit: predict.

## Loop 13. Predict the zeros from the primes with the derived coefficients, nothing fitted ([loop13.py](loop13.py)).
- Prediction (before the run): with S_pred = sum_{n <= y} c_n W_n (c_n from loops 9 and 11) the variance removed from S at
  the zeros equals E(y) = sum c_n^2 alpha_n^2/2 (1 + 2 pi f_n alpha_{n^2}/L) within 2 % for u <= 1, on zeros 80,001..100,000
  and, unchanged, on zeros 40,001..60,000; the residual at u = 1 is below 10 % of Var S (rms error of a predicted zero
  below 0.1 mean spacing).
- Result of loop 13 (nothing fitted). Zeros 80,001..100,000: variance removed over predicted = 0.999 (n = 2), 1.006, 1.013,
  1.014, 1.015, 1.016, 1.019, 1.025, 1.029 (u < 0.9), 1.019 (u < 1), 1.007, 1.003. Residual: 75.1 % (n = 2 only), 46.6 %
  (n <= 6), 31.0 % (n <= 16), 13.4 % (n <= 104), 6.8 % (n <= 667), 4.8 % (n <= 4,278, u < 0.9), then 5.5 %, 6.2 %, 6.4 %.
  Rms error of a predicted zero: 0.070 mean spacings at u < 0.9. Zeros 40,001..60,000: 1.022, 1.027, 1.031, 1.021 at
  u < 0.7..1.0, residual 4.6 % at u < 0.9. "Within 2 %" holds to u = 0.7 and fails narrowly at 0.8 and 0.9 (2.5 %, 2.9 %);
  "residual below 10 % at u = 1" holds (5.5 %).
- Lesson: the two formulas predict the position of every zero from the primes below (t/2 pi)^0.9 to 0.07 of a mean
  spacing, with no fitted number. Adding primes beyond u = 0.9 makes it slightly worse: there c(u) changes sign and u was
  computed with one L for the whole window (L runs from 9.2 to 9.4). Use the local L of each zero.

## Loop 14. The same with u and the correction taken at the local height of each zero ([loop14.py](loop14.py)).
- Prediction (before the run): the residual at u < 1.2 falls below the 4.8 % of loop 13 and no longer rises after u = 0.9.
- Result of loop 14: no change (4.84 % at u < 0.9, 5.46 % at u < 1, 6.39 % at u < 1.2). The prediction FAILS: the local
  height is not the cause.
- Lesson (reasoning, to test): the window is 27,000 long in t, so waves of prime powers n and n' are orthogonal over it
  only when |ln n - ln n'| > 1/27,000, that is for n below about 4,000: exactly where the residual starts to rise
  (u = 0.9 is n = 4,278). Beyond that, single-wave coefficients cannot simply be added. The 4.8 % is an upper bound set
  by the window, not by the law. Test: a longer window should move the turning point up.

## Literature check of the response law (27 November 2026, one agent; arXiv, journals, Zenodo, GitHub, ResearchGate, viXra)
- (A) c(u) = 1 - (1+u) ln(1+u) - (1-u) ln(1-u): NOT FOUND. Closest read: the power spectrum of level displacements
  (Riser-Kanzieper, arXiv 2209.04723), a two-point quantity in the level index. Bogomolny-Keating 2013 and Conrey-Snaith
  2008 average over a window, so a sum at a fixed frequency, sum_gamma S(gamma) n^(i gamma), is not among their formulas.
- The bispectrum: the vanishing for a + b <= N is in print (Hughes-Rudnick, math/0206289, Thm 7). The value
  max(0, a + b - N) was not found printed; it follows in a few lines from hook characters (exact at finite N for a, b <= N).
- (B) the arithmetic term: KNOWN IN SUBSTANCE. It is Landau's formula sum_{gamma <= T} X^rho = -(T/2 pi) Lambda(X) + O(log T)
  applied term by term to the prime sum for S; ours is the bookkeeping and the regression framing.
- (C) the variance fractions: not found as numbers; expected in substance from Selberg's variance and the splitting of
  Gonek-Hughes-Keating.
- NOT OPENED (refused by the publisher's site): three papers of Fujii (Proc. Japan Acad. 66 (1990) "On a theorem of
  Landau II", "An additive theory of the zeros"; ASPM 21 (1992)), the most likely places for sums of S(gamma) x^(i gamma).
  Until they are read, (A) is "not found", not "new".

## Loop 15. An exact formula at finite matrix size ([loop15.py](loop15.py)).
- Derivation (by hand). For Haar unitary N x N, rho_j = sum_n exp(-2 pi i j x_n/N). Exactly, S(x_n) = const + sum_{j != 0}
  rho_j exp(2 pi i j x_n/N)/(2 pi i j) (the Fourier series of the staircase at a jump), so
      sum_n S(x_n) W_k(x_n) = N |rho_k|^2/(2 pi^2 k^2) - (1/(2 pi^2 k)) sum_{j != 0, -k} B(j, k)/j,   B(a, b) = E rho_a rho_b conj(rho_{a+b}).
  With B(a, b) = min(N, a) + min(N, b) - min(N, a + b) for a, b > 0 (hook characters; for a, b <= N it is max(0, a + b - N)),
  E|rho_k|^2 = k and sum_n W_k(x_n)^2 = 2N|s|^2 + 2 Re s^2 conj(rho_2k), for 1 <= k <= N:
      c_N(k) = [ 1 - ((N + k)(H_{N+k} - H_N) - (N - k)(H_N - H_{N-k}))/N ] / [ 1 - max(0, 2k - N)/(N k) ],
  H the harmonic numbers, c defined as the ratio of the two expectations. As N grows it tends to the curve of loop 9.
- Prediction (before the run): simulation at N = 6, 10, 16 (40,000 matrices each) agrees with c_N(k) for every k <= N
  within 3 standard errors; at these sizes c_N differs visibly from the limit curve.
- Result of loop 15. N = 6 (40,000 matrices): k = 1..6: 0.97232, 0.88767, 0.74287, 0.57213, 0.22003, -0.38788 against the exact
  35/36, 559/630, 1243/1680, 3919/6930, 229/1008, -4247/11550 (= 0.97222, 0.88730, 0.73988, 0.56551, 0.22718, -0.36771);
  largest |z| 2.4; at k = 4 the limit curve (0.51483) is 16 standard errors away, the exact formula 1.8. N = 16: largest
  |z| 1.8. N = 10: one point at z = -3.3 (k = 8), so the registered "every k within 3" FAILS in that run; repeated with
  200,000 matrices ([loop15b.py](loop15b.py)): all ten k within 2.3 standard errors (k = 8: 0.28812 +- 0.00336 against
  0.28837); the denominators agree with 1 - max(0, 2k - N)/(N k) (0.967, 0.943, 0.924, 0.908, 0.904 for k = 6..10).
- Lesson: the formula is exact at finite size; the 3.3 was noise. What remains for a theorem is to write the proof of
  B(a, b) = min(N, a) + min(N, b) - min(N, a + b) in full (hook characters) and have it refereed.

## Loop 16. Connecting a dot: the 0.94 of round 147 ([loop16.py](loop16.py)).
- The repo's round 147 (STUDY.md) modelled a prime's effect on zero spacings by the "wave closure": the zeros with p are
  the zeros without p, shifted by p's wave, the rest independent of the wave. Measured amplitude of the spacing-variance
  drop against that closure: 0.940 +- 0.011 (913 L-functions, L = 10.2 to 11.1), "not 1", unexplained; the review found
  the base spectrum is not independent of the wave (+-4 %).
- Reasoning: loops 4-15 say exactly how it is not independent: the zeros follow a wave of frequency u with a coefficient
  below 1. For spacings the relevant coefficient is c'(u): the regression of 1 - s_n = S(x_{n+1}) - S(x_n) on the wave
  difference W(x_{n+1}) - W(x_n). If S = (base) + W with Cov(base differences, wave differences) = -(1 - c') Var, then
      Var(s with p) - Var(s without p) = Var(wave difference) (2 c' - 1),   closure: c' = 1.
  So the amplitude against the closure should be the average of 2 c'(u) - 1 over the powers of p, weighted by their
  share of the wave-difference variance, with c' from random matrices (the powers of one prime are handled jointly by
  the closure's own map, so no arithmetic term here).
- Prediction (before the run): c'(u) from Haar matrices (N = 400); at L = 10.5 the predicted amplitude for p = 2, 3, 5
  lies between 0.92 and 0.97 (hand estimate 0.94 to 0.96), against the measured 0.940 +- 0.011.
- Result of loop 16: c'(u) from matrices is 1.01 to 1.03 for u < 0.4 (0.99 at 0.41, 0.94 at 0.49, 0.83 at 0.57): for
  neighbouring zeros the response is full, even slightly above 1. Predicted amplitude against the closure: 1.018 to 1.023
  at L = 10.2 to 11.1. The prediction FAILS (I had said 0.92 to 0.97); the measured 0.940 +- 0.011 is not explained.
- Lesson: my hand estimate used the coefficient for the position of a zero (c < 1) where the coefficient for the
  DIFFERENCE of neighbouring positions was needed. The deficit of c comes from far pairs (modes beyond the matrix size),
  not from neighbours. Spacings follow the wave fully; positions do not. Round 147's 0.94 stays open.
