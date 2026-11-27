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
