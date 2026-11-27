# Lab log: the zeros under the heat flow (de Bruijn-Newman), short loops

Each entry: hypothesis with a predicted number written BEFORE the run, the test, the result, the lesson.
Nothing here is a claim until it is promoted (independent check, literature check, proof).

## The idea (27 November 2026, reasoning only)
- The heat flow H_t(z) = int e^{t u^2} Phi(u) cos(zu) du multiplies the n-th term of the Dirichlet series by
  b_n = exp((t/4) ln^2 n) (Polymath 15). With ln n = sum of the logs of its prime factors,
  b_n = prod_i exp((t/4) l_i^2) prod_{i<j} exp((t/2) l_i l_j): a gas of primes with a pair interaction (t/2) ln p ln q.
  At t = 0 the primes are free (Euler product). RH says the zeros are real at t = 0; they are real for t >= Lambda and
  Lambda >= 0 (Rodgers-Tao). So: zeros real <=> the interaction is not repulsive, and the free gas sits on the edge.
- The logarithm of an interacting gas is a sum over connected graphs (Mayer). So for t > 0 the "von Mangoldt function"
  of H_t lives on ALL integers: ln(n) kappa_t(n), kappa_t(pq) = e^{(t/4)(l_p + l_q)^2} - e^{(t/4)(l_p^2 + l_q^2)}.
  Prediction: the zeros of H_t show lines at composite numbers.
- The same flow on the characteristic polynomial of a unitary N x N matrix multiplies coefficient j by exp(-t j (N-j))
  (the rates of the log-gas around a lattice). Partitions play the role of factorisations. In the scaling t = tau/N^2:
      K_tau(u) = u^2 e^{-2 tau u} sum_{r>=1} (1/r!) int_{sum a_i = u} prod_i (e^{2 tau a_i^2}/a_i) Phi_r(a)^2,
      Phi_r(a) = sum over connected graphs on r points of prod_edges (e^{2 tau a_i a_j} - 1),      0 < u < 1,
  first term u exp(-2 tau u (1-u)). For zeta tau = t ln^2(T/2 pi)/4. Same series from integers and from matrices.

## Loop 1. Matrices: the exact partition formula against simulation ([loop1.py](loop1.py))
- Exact claim (k <= N): E|p_k(t)|^2 = e^{-2tNk} sum_{lambda |- k} c_lambda^2 z_lambda,
  c_lambda = k prod_s e^{t i_s^2} Conn_lambda(t) / z_lambda, Conn = connected-graph sum with edge weights e^{2 t i j} - 1.
- Prediction (before the run): N = 30, tau = 1, 3, 6, 4,000 matrices, k = 1..12: every simulated value within 3 standard
  errors of the exact formula; all flowed roots on the unit circle (to rounding); at tau = 3 and k = 12 the multi-part
  terms carry more than 10 % of the value.
- Result of loop 1 (N = 30, 4,000 matrices, k = 1..12, tau = 1, 3, 6): all 36 values within 2.7 standard errors of the exact
  formula (the three tau share the matrices); roots on the circle to 1e-9. The share prediction FAILS: multi-part terms
  carry 6 % at tau = 3, k = 12 (18 % at tau = 6); I had said above 10 %.
- Lesson: the exact formula stands. The cluster terms are small for small u and tau; they matter near u = 1/2 and beyond.

## Loop 2. The limit law against large matrices, and what happens beyond u = 1 ([loop2.py](loop2.py))
- Prediction (before the run): N = 100, 1,500 matrices, tau = 1 and 4: for u = 0.1 ... 0.9 the series with clusters of up to
  four parts equals the simulated K within 3 standard errors or 2 %; at u = 1 the simulated K is ABOVE 1 (a Bragg peak
  starts to grow there, since the one-part term is not damped at u = 1 and the clusters add to it).
- Result of loop 2. First run (roots by numpy from the coefficients, N = 100): garbage (K = 0.91 at u = 0.4): np.poly and
  np.roots lose all digits at degree 100; a bug of method, not a result. Second run without root finding
  ([loop2b.py](loop2b.py): for k <= N the traces have exactly the moments of independent Gaussians, so E|p_k(t)|^2 is a
  Gaussian expectation of a power-series coefficient; 20,000 samples, N = 100):
      tau = 1:  u = 0.5: 0.3116 +- 0.0022 (series 0.3090);  u = 0.9: 0.8695 +- 0.0061 (0.8766)
      tau = 4:  u = 0.5: 0.0848 +- 0.0006 (series 0.0839);  u = 0.8: 0.4666 +- 0.0034 (0.4749);  u = 0.9: 1.1480 +- 0.0084 (1.1541)
  All 18 values within 2.5 standard errors (and within 2 %); tau = 0 gives K = u as it must. At tau = 4, u = 0.9 the
  clusters carry 62 % (terms r = 1..4: 0.438, 0.509, 0.179, 0.028). The Bragg-peak part (u >= 1) is not measured yet.
- Lesson: check a simulation at tau = 0 against the known answer before reading it. The limit law stands for u < 1.

## Reasoning between loops: a "positivity principle", stated and refuted the same hour
- Observation (elementary): the cluster weights kappa_t(n) are >= 0 for every n when t >= 0 (positive bonds), and
  kappa_t(pq) < 0 for t < 0. The zeros are real for t >= Lambda, and Lambda >= 0 (Rodgers-Tao). So RH <=> "the zeros are
  real exactly as long as all cluster weights are non-negative": two sign changes at the same t = 0.
- Tempting principle: for deformations of zeta by a universal factor, real zeros <=> non-negative log-coefficients.
- Refuted by reasoning: the universal factor 1 + c u^2 (operator 1 - c d^2/dz^2, preserves real zeros) multiplies the
  n-th term by 1 + c (L/2 - ln n)^2/4; its log-coefficient at pq is 2 eps l_p l_q - eps^2 l_p l_q (L - l_p)(L - l_q) with
  eps up to 4/L^2: negative for small p, q and large c. Real zeros do not force non-negative weights. And Dirichlet
  characters have real zeros (GRH) with weights chi(n) kappa_t(n). So the coincidence at t = 0 is special to the
  Gaussian semigroup; it is a reading of Newman's conjecture in prime language, not a route to a proof.

## Loop 3. Beyond u = 1, with a stable method ([loop3.py](loop3.py))
- Method: coefficients by FFT, the flowed real function on a grid of 256 points per zero, zeros = sign changes.
- Prediction (before the run): N = 100, 1,500 matrices: (a) tau = 0 gives K = u for u < 1 and K = 1 beyond, spacing variance
  0.180 (validation); (b) u = 0.3, 0.6, 0.9 agree with the series within 3 standard errors at tau = 1, 4; (c) K(1) > 1 for
  tau = 1 and 4, larger at 4; (d) every matrix has exactly N real zeros at tau >= 0; (e) the spacing variance falls with
  tau (my guess: below 0.10 at tau = 4).
- Result of loop 3 (N = 100, 1,500 matrices, same matrices for every tau):
      tau    real zeros   spacing variance   smallest spacing   K(0.9) (series)     K(1)     K(1.2)   K(2)
      0      all N        0.1795             0.020              0.936 (0.900)       1.010    1.002    0.980
      1      all N        0.1140             0.405              0.908 (0.877)       1.246    1.225    0.977
      4      all N        0.0395             0.618              1.194 (1.154)       3.310    1.509    0.960
      16     all N        0.0025             0.847              0.611 (0.525)      15.30     0.066    4.568
  (a), (c), (d), (e) hold. (b) holds at tau = 1, 4 (largest 1.4 standard errors); at tau = 16, u = 0.9 the series with four
  cluster terms is 14 % low (more clusters are needed there). The peak at u = 1 is also given by the series, which is
  exact up to k = N: K_tau(1) = 1.240 (tau = 1) and 3.259 (tau = 4, six terms) against 1.246 +- 0.033 and 3.310 +- 0.081.
- Lesson: the flow turns the gas into a crystal: the spacing variance falls from 0.180 to 0.040 at tau = 4, a Bragg peak
  grows at u = 1 (and at u = 2 for large tau), and the series describes the whole range 0 < u <= 1. Beyond u = 1 there is
  no formula yet (the modes above the size are functions of the lower ones).

## Loop 4. Zeta: do the zeros of H_t show lines at composite numbers? ([loop4_lines.py](loop4_lines.py), [loop4_check.py](loop4_check.py))
- Prediction (frozen in the brief of the compute agent before any run): for the zeros gamma_k(t) of H_t (in zeta heights),
      sum_{window} e^{i gamma_k(t) ln n} = -(1/2 pi) ln(n) kappa_t(n) n^(-1/2) int (gamma/2 pi)^(-(t/4) ln n) d gamma  x  phase,
  kappa_t = coefficients of log sum_m e^{(t/4) ln^2 m} m^(-s). At t = 0 this is Landau's formula (lines at prime powers
  only). For t > 0: primes damped, and lines at every composite; pooled slope for composites within 3 errors of 1.
- Data: the first 12,000 zeros moved by the flow d gamma_k/dt = (1/2) sum 1/(gamma_k - gamma_j) (agent's integrator; snapshots
  t = 0, 0.05, ..., 0.6). My check of the snapshots, with no flow equation: H_t(2 gamma) = (1/8) E_G xi(1/2 - sqrt(t/2) G + i gamma)
  by Gauss-Hermite and mpmath; six zeros (numbers 500 to 8000) at t = 0.3 and 0.6: the moved zeros are roots of H_t within
  7.4e-6 (they moved by up to 0.36).
- Result (my analysis; zeros 2,000..10,000; the phase is the rigid shift t pi/8 of all zeros, imaginary parts vanish):
      t      n = 6: measured (predicted)   n = 15            n = 21           pooled slope, 139 composites <= 200
      0        -1.3   (0)                    1.5  (0)         -2.5  (0)        -
      0.1     -26.7  (-25.6)               -51.5 (-52.5)     -61.2 (-59.2)     0.984 +- 0.016
      0.3     -48.1  (-47.1)               -82.5 (-82.9)     -91.7 (-90.4)     0.985 +- 0.009
      0.6     -46.3  (-45.5)               -64.0 (-64.0)     -67.2 (-66.5)     0.980 +- 0.009
  primes 0.996 to 0.999 +- 0.003 to 0.007, prime powers 0.992 to 1.001 +- 0.011 to 0.019 (noise from half-integers, which
  is larger than the scatter seen at integers). Separate windows at t = 0.3: zeros 2,000..6,000: 0.982 +- 0.018;
  6,000..10,000: 0.991 +- 0.019; zeros 500..2,000: 1.127 +- 0.024 (there n <= 200 reaches u = 1, outside the formula).
  The prediction holds: the lines at composites are there, with the cluster weights, at more than 50 errors from zero.
- Lesson: under the heat flow the zeros "hear" every integer. The rigid shift is t pi/8 per unit t (from the imaginary
  part of the exponent shift). The pooled composite slope sits 1.5 % below 1 at every t (same zeros, so not independent):
  either noise (1.5 to 2 errors) or a next-order term; kappa(pq) is a difference of two nearly equal numbers, so a
  0.1 % correction of b_n would do it. Open.

## Loop 5. One family for both worlds? Spacing variance of zeta's flowed zeros against flowed matrices ([loop5.py](loop5.py))
- Conjecture (ours): in the scaling tau = t ln^2(T/2 pi)/4 the zeros of H_t near height T have the local statistics of
  Haar matrices flowed to tau = t N^2 (a one-parameter family from the sine process, tau = 0, to the lattice).
- Prediction (before the run): zeros 6,000..10,000 at t = 0.05, 0.1, 0.2, 0.3, 0.45, 0.6 (tau = t <L^2>/4, about 0.6 to 7.7):
  the variance of the unfolded spacings equals that of flowed matrices (N = 100, 400 matrices) at the same tau within
  10 % at every t. (At t = 0 zeta is known to sit a little below 0.180 at this height; I expect the same few per cent.)
- Result of loop 5: FAILS. Spacing variance, zeta / matrices: 0.1565 / 0.1796 at t = 0 (ratio 0.871), then 0.871, 0.866,
  0.848, 0.824, 0.778, 0.720 at tau = 0.64 ... 7.7. Zeta at this height (L = 7) is 13 % below the sine process before any
  flow, and the ratio drifts further.
- Lesson: at L = 7 the scaling limit is not reached (the primes are few and discrete: u = 0.10 for p = 2); this run does
  not test the conjecture. The test must be made where zeta is close to the sine process: Odlyzko's zeros near 10^12
  (L = 24.5). In unfolded units both worlds obey the same equation, dx_k/dtau = (1/(2 pi^2)) sum_{j != k} 1/(x_k - x_j).

## Loop 6. The same at height 10^12, 10^21, 10^22 ([loop6.py](loop6.py))
- Method: blocks of 3,333 consecutive zeros from each of Odlyzko's three tables (9 blocks), unfolded, moved by
  dx/dtau = (1/(2 pi^2)) [sum 1/(x_k - x_j) + ln((b - x_k)/(x_k - a))] (the zeros outside the block as a uniform
  continuum), statistics on the middle 60 % of each block. Matrices: N = 200, flowed exactly (loop 3's method).
- Prediction (before the run): spacing variance of zeta at tau = 0.5, 1, 2, 4 equals the matrix value within 3 pooled
  standard errors (errors from the scatter of the 9 blocks); and the amplitude of the line of p = 2, 3, 5 in the middle
  zeros falls by exp(-tau u (1 - u)), u = ln p / L, within 3 errors at tau = 4.
- Result of loop 6 (9 blocks of 3,333 zeros, middle 2,000 of each; matrices N = 200, 150 matrices):
      tau    spacing variance: zeta        matrices            z       line of p = 2, 3, 5 relative to tau = 0 (predicted)
      0      0.1779 +- 0.0019              0.1772 +- 0.0016    +0.3
      0.5    0.1397 +- 0.0014              0.1393 +- 0.0012    +0.2    0.9918 (0.9905)   0.9854 (0.9852)   0.9786 (0.9787)
      1      0.1130 +- 0.0011              0.1129 +- 0.0010    +0.1    0.9836 (0.9811)   0.9711 (0.9705)   0.9577 (0.9579)
      2      0.0770 +- 0.0007              0.0772 +- 0.0007    -0.2    0.9678 (0.9625)   0.9432 (0.9420)   0.9174 (0.9176)
      4      0.0388 +- 0.0004              0.0393 +- 0.0004    -1.0    0.9373 +- 0.013 (0.9266)   0.8901 (0.8877)   0.8420 (0.8426)
  Both predictions hold: at every tau the zeros of zeta and the eigenvalues of Haar matrices have the same spacing
  variance (largest difference 1.0 error), and the lines of 2, 3, 5 fall by exp(-tau u (1 - u)) (all within 0.8 errors).
- Lesson: where zeta is close to the sine process (L = 24 to 49), the flow takes zeta and matrices through the same
  one-parameter family, from the sine process to the lattice. The conjecture survives its first real test. (The block
  edges are treated as a continuum; only the middle 60 % is read.)

## Reports of the two agents (27 November 2026)
- Compute agent (own flow code with a block-multipole force, all zeros below the M-th and their mirrors, M = 12,000 and
  30,000; code in [agent_zeta/](agent_zeta/)): with a smooth taper over the window the noise is 1e-3 against amplitudes of
  10 to 680, so the test is nearly exact. At t = 0.4: primes 2..13 within 1e-5 of the prediction, composites 6, 10, 12, 14,
  15, 21 within 4e-5 (errors 2e-5 to 5e-5). Pooled composite slope - 1: (-1 +- 5)e-6 at t = 0.2 and (-5 +- 3)e-6 at t = 0.4
  (M = 30,000); at t = 0.6 it is (-11 +- 3)e-6, a term of order 1/T (it halves from M = 12,000 to 30,000 and vanishes when
  L is taken at each zero's starting height; post hoc). No lines at rationals. Rigid shift per unit t: -pi/8 + 7/(8 gamma).
  So the 1.5 % deficit of my loop 4 was my sharp window, not a next-order term.
- Prior-art agent (Rodgers-Tao, Polymath 15, Dobner, Kabluchko, Hall-Ho, Hall-Ho-Jalowy-Kabluchko, Ki-Kim-Lee,
  Csordas-Smith-Varga, Newman 1976 read in full; Zenodo, GitHub searched): (A) the limit law K_tau(u): NOT FOUND; its
  first term is Rodgers-Tao's factor K_{t,n} squared (arXiv 1801.05914, p. 17). (B) The series sum b_n n^-s is Polymath's
  and Dobner's; Polymath section 8.1 notes that b_n is not multiplicative for t != 0; the lines at composites and the
  cluster weights as a law: NOT FOUND (the parts are elementary). The Curie-Weiss reading for polynomials is Kabluchko
  (arXiv 2203.05533) and Newman 1976; for primes not found. (C) an exact finite-N second moment of the flowed power
  sums and local statistics at t ~ 1/N^2 or 1/log^2 T: NOT FOUND; Rodgers-Tao p. 27 name this time scale but compute no
  statistics; Polymath Remark 1.6 expects a "liquid" below exp(C/t). NOT OPENED: the Polymath wiki and blog threads,
  the original papers of Landau, Gonek, Diaconis-Shahshahani. So "not found", with those unopened.

## Loop 7. The other side, t < 0: how many zeros leave the line? ([loop7.py](loop7.py))
- Reasoning: backwards in time two neighbours attract; alone they obey g^2(tau) = g_0^2 + 2 tau/pi^2 (g in mean spacings),
  so a pair with g_0 < sqrt(2|tau|)/pi has collided and left the line. With the small-gap law of the sine process,
  density (pi^2/3) s^2, the share of zeros off the line is f(tau) = (2/(9 pi)) (2|tau|)^(3/2) for small |tau|.
  For zeta (if RH and the pair law hold): at t < 0 the share of zeros of H_t off the line near height T would be
  (2/(9 pi)) (|t| ln^2(T/2 pi)/2)^(3/2): a count behind Newman's conjecture (Rodgers-Tao: Lambda >= 0).
- Prediction (before the run): Haar matrices N = 100, 4,000 matrices, flowed to tau = -0.03, -0.06, -0.12: share of zeros
  off the circle = f(tau) within 3 standard errors or 10 %; at tau = -0.25, -0.5 the two-body formula is off by more
  than 10 % (my guess: the true share is higher, the neighbours press the pair together).

## Loop 8. The whole form factor of zeta's flowed zeros at height 10^12 to 10^22 ([loop8.py](loop8.py))
- Prediction (before the run): for the 9 blocks of loop 6, the form factor |sum e^{2 pi i u x}|^2 / M averaged over bins of
  u of width 0.1 equals the series K_tau(u) (five cluster terms) for u in (0.2, 1.0) at tau = 1 and 4 within 3 errors
  (errors from the 9 blocks), and the bin (0.95, 1.05) shows the Bragg peak (above 1 at tau = 1, above 2 at tau = 4).
- Result of loop 7 (N = 100, 4,000 matrices; share of zeros off the circle against (2/(9 pi)) (2|tau|)^(3/2)):
      tau      -0.03          -0.06          -0.12          -0.25          -0.5           -1       -2
      share    0.00096        0.00282        0.00778        0.02265        0.05923        0.1426   0.2958
      ratio    0.92 +- 0.07   0.96 +- 0.04   0.935 +- 0.024 0.906 +- 0.013 0.837 +- 0.007 0.713    0.523
  First part holds (the three smallest |tau| within 10 % and 3 errors). Second part: the formula is indeed off by more
  than 10 % from |tau| = 0.5 on, but my guess of the direction FAILS: the true share is lower, not higher. At tau = 0 no
  zero is lost (check). At tau = -2 already 30 % of the zeros have left the circle.
- Lesson: the count behind Newman's conjecture has a simple leading law, share = (2/(9 pi)) (2|tau|)^(3/2), and the
  neighbours slow the collisions down rather than speeding them up (the next term of the small-gap law, a factor
  1 - 0.16|tau|, explains about half of the deficit). For zeta this is a prediction, not a measurement: it needs RH and
  the small-gap law of the sine process, and H_t for t < 0 was not computed here.
- Result of loop 8 (9 blocks, middle 2,000 zeros of each; form factor averaged over bins of u; series with five terms):
      u bin        tau = 0: zeta (u)        tau = 1: zeta (series)     tau = 4: zeta (series)
      0.05-0.15    0.0936 +- 0.0025 (0.1)   0.0772 +- 0.0020 (0.0827)  0.0440 +- 0.0011 (0.0473)
      0.25-0.35    0.2969 +- 0.0070 (0.3)   0.1955 +- 0.0044 (0.1976)  0.0575 +- 0.0011 (0.0582)
      0.45-0.55    0.4979 +- 0.0123 (0.5)   0.3082 +- 0.0080 (0.3097)  0.0840 +- 0.0024 (0.0848)
      0.65-0.75    0.6960 +- 0.0097 (0.7)   0.4902 +- 0.0066 (0.4926)  0.2303 +- 0.0035 (0.2326)
      0.75-0.85    0.8108 +- 0.0120 (0.8)   0.6541 +- 0.0129 (0.6473)  0.4853 +- 0.0106 (0.4906)
      0.85-0.95    0.8766 +- 0.0106 (0.9)   0.8542 +- 0.0144 (0.8817)  1.1419 +- 0.0176 (1.2056)
      0.95-1.05    0.9739 +- 0.0217         1.2237 +- 0.0192           2.9828 +- 0.0334
      1.25-1.35    0.9799 +- 0.0182         1.1717 +- 0.0239           1.0054 +- 0.0156
  The seven bins from 0.15 to 0.85 agree within 1.2 errors at every tau. The registered reading ("within 3 errors on
  (0.2, 1.0)") FAILS in one cell: bin 0.85-0.95 at tau = 4 is 3.6 errors low; the same bin is 2.2 errors low at tau = 0,
  where the comparator is Montgomery's K = u (the known rounding of zeta's form factor near u = 1 at finite height).
  The lowest bin is 2.5 to 3 errors low at every tau, also at tau = 0 (few prime powers in that bin). Bragg peak: 1.22 at
  tau = 1 and 2.98 at tau = 4, as predicted (above 1, above 2).
- Lesson: the cluster series describes the form factor of zeta's flowed zeros wherever zeta at tau = 0 follows K = u;
  the two deviating bins deviate before the flow starts. A fair reading divides by the tau = 0 value.

## Referee of theory/proof.md (27 November 2026, own code theory/referee/ref.py, 85 checks, 0 failed)
- Cluster expansion, circle property for t >= 0, exact second moment for k <= N (k = N included; Diaconis-Evans Thm 2.1
  read in the Berkeley tech report, the AMS PDF refused), limit K_tau(u) for 0 < u <= 1 and also tau < 0: NO FAILING STEP.
- WRONG: the side remark "for k > N the full Hall norm is strictly larger" fails for t < 0 (gap 0 at N = 1, k = 2,
  t = -0.346574; N = 2, k = 3, t = -0.464522 and -0.149050). Not yet corrected in proof.md.
- Third-order term derived by the referee: -4/3 u^4 (1-u)^3 + 8/15 u^7 - 2/3 u^6 (matches the series to 1e-7).
  K_2(0.7) = 0.3723234 (r <= 7); finite N = 10, 20, 40 converge as 1/N^2.
- Summary written: [heatflow.md](heatflow.md).
