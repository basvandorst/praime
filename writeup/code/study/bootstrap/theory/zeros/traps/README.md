# Traps for the double zeros of the optimal test functions (class A, deg 4..48)

Written by a compute agent on 26 November 2026. Reuses `core.py`, `kbessel.py` and the json files of
`../../convex/certs/`. Run every script from this folder as
`(ulimit -v 1500000; timeout 3000 nice python <script> <args>)`.

## Scripts (one command each)

| command | what it does | rigorous? | output |
|---|---|---|---|
| `python basevals.py 160` | m_j = g_j(0), cost vector c_j = sum_{n<=160} 2 Lambda(n) n^{-1/2} g_j(log n), j <= 24 (arb 900 bits, own K0/K1; last term 1e-331 relative) | numerics, ~108 digits | `basevals.json` |
| `./run_polish.sh` (= `python polish_hp.py d`, d = 2..24) | re-polishes each stored KKT optimum (Newton on the stationarity residual z from the stored 30-digit touch points; |z| <= 1e-107); builds h_k, W(h_k), h_k(0), h_k''(u_j), g''(u_k) | numerics, ~100 digits | `hp_d<d>.json`, `pol_d<d>.log` |
| `python cmp_stored.py` | re-polish vs stored optimum (independent Bessel code: arb/own K vs mpmath.besselk) | check | `cmp_stored.out` |
| `python grid.py` | rows g_j(u), j <= 24, on u = i/1000 (i <= 4000) and 4 + i/100 (i <= 400) | numerics | `grid.pkl` (11 MB, deleted at the end; 8 min) |
| `python traps.py 2 4 ... 24` | D1: eta_0, eta_1, S for eta = 2 eta_0 and 10 eta_0, endpoints, trap radius, weight bound, extra components | numerics | `traps_d<d>.json`, `traps_d<d>.out` |
| `python margins.py` | D1: interior local maxima of rho = eta h/gt on the grid off the main component | numerics | `margins.out` |
| `python d2.py` | D2: radius and g''(u_k) against the degree; g(log 6), g(log 10), g(log 12), b_0 | numerics | `d2.out` |
| `python trapcert.py k` (all: `./run_trapcert.sh`) | D3: rigorous trap for node k at degree 48: (i) Sturm, (iii), S c (a_k, b_k), (ii) | ball arithmetic + exact Sturm | `trapcert_k<k>.json`, `.log` |
| `python trapweight.py k` | D3: rigorous weight lower bound with a sharp sup of h on [a_k, b_k] (Taylor models on 16 pieces) | ball arithmetic | `trapweight_k<k>.json` |
| `python collect.py` | assembles `traps.out` (D1, all 78 (degree, k)) and `d3_summary.out` (D3 + cross-check of W(h)) | | |
| modules `fam.py`, `fcert.py`, `common.py` | fast arb evaluation of the family (numerics; K inputs replaced by ball midpoints); certify.py rewritten as a class (rigorous, same bounds) | | |

## Setting (as specified)

gt = optimum at half-degree d' (g(0) = 1, kt = W(gt) = lambda), h_k in V_{d'} with h(u_j) = h'(u_j) = 0 (j != k),
h(u_k) = 1, h'(u_k) = 0, top coefficient 0. eta_0 = lambda/W(h_k), eta_1 = min_{B_h > 0} B_gt/B_h,
S = {u > 0 : eta h_k - gt > 0}, trap radius = max(|a - log n_k|, |b - log n_k|) for the component (a, b) around u_k,
weight bound = (eta W(h) - kt (1 + (eta h(0) - 1)^+)) / sup_S w.  n_k = prime power nearest e^{u_k} in log scale.

## D1 results (`traps.out`, one line per (degree, k), 78 lines; json per degree)

* KKT identity: W(h_k) - lambda h_k(0) - mu_k <= 5e-60 at all 78 (degree, k). This is an algebraic identity once
  y = A^{-T} c (h_k'(u_j) = 0 kills the z terms), so it only tests consistency. The real tests: the re-polish
  (different Bessel code) reproduces every stored lambda (25 digits), mu (10 digits), touch point (30 digits) and b
  (all 40-50 stored digits) (`cmp_stored.out`); and the D3 ball computation of W(h) (independent code path: certify.py's
  prime-sum code at log n) equals lambda h_k(0) + mu_k - DELTA W(Phi*Phi) (column Wcheck of `d3_summary.out`).
* (i): eta_1 = +infinity at all 78 (degree, k): B_{h_k}(t) < 0 at every scanned t (t = 0 and 10^{-2}..10^{4}, plus the
  negative top coefficient b_{d-1}(h_k)); also h_k(0) < 0 always (from -25.7 to -5.3e21). (iii) holds for all eta
  (b_d(h_k) = 0). So eta is limited only by the shape of S.
* gt >= 0 at all 4400 grid points at every degree (sanity of the optimum).
* eta = 2 eta_0: S is ONE interval at all 78 (degree, k) (no grid point with w > 0 outside it; limits at the other nodes
  rho_j = eta h''(u_j)/g''(u_j) < 1; all interior local maxima of rho = eta h/gt off the main component are < 0,
  `margins.out`; on (6, 8] rho <= 1.3e-4, rho(8) <= 2.3e-6, and w -> -gt at infinity since b_d(h_k) = 0). It isolates n_k (contains no other log m) in 74 cases; the 4
  exceptions are the newest node sitting between prime powers: deg 12 (e^u = 4.357, contains log 5), deg 24 (8.358,
  log 9), deg 28 (9.337, log 10), deg 40 (16.255, log 17). Weight bound / mu_k = 0.385..0.500 (= W(h_k)/(2 mu_k) -> 1/2).
* eta = 10 eta_0: one interval in 73 cases (isolating in 65); extra components in 4 cases, always around a neighbouring
  node with small g'' (rho_j > 1): deg 16 k=3 (extra near e^u = 5.09), deg 28 k=6 (9.34..10.59), deg 28 k=7 (8.08..8.09),
  deg 44 k=10 (17.06..17.17); deg 12 k=3: main component wider than 0.5 (not bracketed). Weight bound / mu -> 0.9.

### Degree 48 (eta = 2 eta_0; R = trap radius about log n; x10 = radius at 10 eta_0); all single, all isolate n_k

| k | n | u_k - log n | mu_k (= 2 Lambda(n)/sqrt n to 10 digits) | g''(u_k) | W(h_k) | h_k(0) | eta_0 | R (x2) | lb/mu | R (x10) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 2 | 2.269e-25 | 0.9802581435 | 1.0559e+00 | 0.9802581435 | -7.634e+01 | 2.4237e-26 | 3.030e-13 | 0.500 | 6.775e-13 |
| 2 | 3 | 4.532e-24 | 1.2685682012 | 4.3334e-02 | 1.2685682012 | -2.374e+03 | 1.8728e-26 | 1.315e-12 | 0.500 | 2.940e-12 |
| 3 | 4 | 1.169e-22 | 0.6931471806 | 3.1376e-03 | 0.6931471806 | -3.518e+04 | 3.4276e-26 | 6.610e-12 | 0.500 | 1.478e-11 |
| 4 | 5 | 3.519e-22 | 1.4395250311 | 4.8926e-04 | 1.4395250311 | -1.376e+05 | 1.6504e-26 | 1.162e-11 | 0.500 | 2.597e-11 |
| 5 | 7 | 9.236e-20 | 1.4709698080 | 1.6318e-06 | 1.4709698080 | -1.793e+08 | 1.6151e-26 | 1.990e-10 | 0.500 | 4.449e-10 |
| 6 | 8 | 9.279e-18 | 0.4901290717 | 4.6436e-08 | 0.4901290717 | -4.285e+09 | 4.8473e-26 | 2.043e-09 | 0.500 | 4.569e-09 |
| 7 | 9 | 4.509e-17 | 0.7324081924 | 6.1569e-09 | 0.7324081924 | -1.286e+10 | 3.2438e-26 | 4.591e-09 | 0.500 | 1.027e-08 |
| 8 | 11 | 1.249e-15 | 1.4459852557 | 1.0825e-10 | 1.4459852557 | -1.930e+12 | 1.6430e-26 | 2.464e-08 | 0.500 | 5.510e-08 |
| 9 | 13 | 2.301e-13 | 1.4227779119 | 6.0574e-13 | 1.4227779119 | -3.827e+14 | 1.6698e-26 | 3.321e-07 | 0.500 | 7.425e-07 |
| 10 | 16 | 2.980e-08 | 0.3465743901 (vs 0.3465735903) | 2.2131e-17 | 0.3465738745 | -2.170e+19 | 6.8551e-26 | 1.118e-04 | 0.500 | 2.513e-04 |
| 11 | 17 | 1.809e-07 | 1.3743173652 (vs 1.3743103385) | 1.0068e-18 | 1.3743135522 | -1.605e+20 | 1.7287e-26 | 2.632e-04 | 0.500 | 5.907e-04 |
| 12 | 19 | 7.875e-06 | 1.3513446715 (vs 1.3510012585) | 3.1388e-20 | 1.3512176298 | -5.347e+21 | 1.7583e-26 | 1.537e-03 | 0.500 | 3.518e-03 |

lambda_48 = 2.37581227957306350493855252941e-26; full digits (100) in `hp_d24.json`, endpoints a, b (40 digits) in
`traps_d24.json`.

## D2 results (`d2.out`)

* For fixed k, g''(u_k) converges as the degree grows: k = 1: 1.35078 (deg 4), 1.09286 (12), 1.05640 (20), 1.05594 (32),
  1.055894 (48), steps now ~5e-6; k = 2: 0.0963 (8) -> 0.043334 (48); k = 3..6 still drift down slowly at deg 48
  (k = 5: 1.84e-6, 1.71e-6, 1.63e-6 at deg 40, 44, 48; k = 6: 5.85e-8, 5.08e-8, 4.64e-8).
* The half-width of S is sqrt(2 eta/g''(u_k)) to leading order: radS / sqrt(2 lambda/(mu_k g''(u_k))) -> sqrt 2 =
  1.414214 (eta = 2 eta_0, W(h_k) -> mu_k), reached to 5 digits from deg 36 on for k <= 4 and about 16-20 degrees
  after a node appears for k = 5, 6. Since g''(u_k) stabilises, the trap radius tends to 0 like sqrt(lambda_{d'}):
  k = 1: 2.1e-1, 4.4e-2, 1.4e-2, 3.8e-3, 4.2e-4, 8.1e-5, 2.3e-5, 2.0e-6, 2.1e-8, 1.6e-9, 9.5e-11, 3.0e-13 (deg 4..48).
  u_k - log n_k falls like lambda (k = 1: (u_1 - log 2)/lambda = 5.72, 7.71, ..., 9.54, 9.55, 9.55 at deg 4, 8, ..., 40,
  44, 48), much faster than the radius, so the trap is centred on log n_k up to O(lambda).
* g_{d'}(log 6) = 3.20e-8 (deg 4), peak 6.00e-6 (16), then 1.80e-6, 8.64e-7, 5.29e-7, 4.02e-7, 3.66e-7, 3.51e-7,
  3.40e-7, 3.33e-7 (deg 20..48): decreasing, decrements 1.4e-8, 1.1e-8, 0.7e-8 in the last three steps: consistent with
  a positive limit near 3.2e-7, not proved. g(log 10) = 2.3e-12, g(log 12) = 1.2e-14 at deg 48, both still decreasing.
  b_0 = 2.25330413... converges (changes 1e-7, 8e-8, 5e-8 per step at the end).

## D3 results (degree 48, rigorous) (`d3_summary.out`, `trapcert_k*.json/.log`, `trapweight_k*.json`)

Set-up: gt = the cushioned rational function of `certs/certificate_A_d24.json` (eps = 1e-30, b to 60 digits);
gt(0) = g0 = [1 +/- 1e-19] is a ball, so the lemma is applied to gt/g0 with kt = W(gt)/g0 =
[2.37581468401501e-26 +/- 3.1e-41] (recomputed here; equals the certificate's ball). h = h_k (re-polish, 100 digits)
rounded to 90 significant digits, top coefficient set to 0, minus DELTA g_0 with DELTA = 1/1000 (lowers h everywhere by
DELTA Phi*Phi > 0; raises the floor of F = gt - eta' h near the other nodes from eps Phi*Phi to (eps + eta' DELTA) Phi*Phi).
eta' = 2 eta_0 rounded to 3 digits (rational), lemma eta = eta'/g0. (a_k, b_k) = numerical S widened by 10 %, rational.

Proved for every k = 1..12 (all 12 nodes):
(i) exact Sturm: q_F(s) = q_gt(s) - eta' q_h(s) has no root in [0, oo) and q_F(0) > 0;  (iii) b_24(h) = 0;
S c (a_k, b_k): F > 0 on [0, a_k] and [b_k, 6] (3377-3501 Taylor-model intervals, order 16, 800 bits, bisection depth
<= 53 below 1/64, tightest near log 19 and log 16) and on [6, oo) (asymptotic argument, N >= 2 ratio 7e-1030);
(ii) eta W(h) - kt (1 + (eta h(0) - 1)^+) = 2.37..2.38e-26 > 0 (balls; h(0) < 0 so the (.)^+ term is 0);
weight >= (eta W(h) - kt)/(eta sup_[a_k,b_k] h), sup h <= 1 + 1e-8 by Taylor models on 16 pieces (`trapweight.py`).

| k | n | a_k - log n | b_k - log n | certified radius R_k <= | certified weight >= | mu_k |
|---|---|---|---|---|---|---|
| 1 | 2 | -3.33398e-13 | 3.33398e-13 | 3.33398e-13 | 0.49039 | 0.98026 |
| 2 | 3 | -1.44713e-12 | 1.44713e-12 | 1.44713e-12 | 0.63499 | 1.26857 |
| 3 | 4 | -7.27396e-12 | 7.27396e-12 | 7.27396e-12 | 0.34679 | 0.69315 |
| 4 | 5 | -1.27760e-11 | 1.27760e-11 | 1.27760e-11 | 0.71956 | 1.43953 |
| 5 | 7 | -2.18864e-10 | 2.18864e-10 | 2.18864e-10 | 0.73540 | 1.47097 |
| 6 | 8 | -2.24719e-9 | 2.24719e-9 | 2.24719e-9 | 0.24492 | 0.49013 |
| 7 | 9 | -5.05069e-9 | 5.05069e-9 | 5.05069e-9 | 0.36631 | 0.73241 |
| 8 | 11 | -2.71205e-8 | 2.71205e-8 | 2.71205e-8 | 0.72383 | 1.44599 |
| 9 | 13 | -3.65289e-7 | 3.65294e-7 | 3.65294e-7 | 0.71143 | 1.42278 |
| 10 | 16 | -1.21849e-4 | 1.22946e-4 | 1.22946e-4 | 0.17313 | 0.34657 |
| 11 | 17 | -2.87201e-4 | 2.89565e-4 | 2.89565e-4 | 0.68764 | 1.37432 |
| 12 | 19 | -1.60459e-3 | 1.69084e-3 | 1.69084e-3 | 0.67625 | 1.35134 |

Read with the lemma: for EVERY half-degree d >= 24 (deg B >= 48), the optimal g_d has, for each of the 12 prime powers
n = 2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, a zero within R_k of log n carrying KKT weight >= the bound (about mu_k/2).
Cross-check: ball W(h) + DELTA W(Phi*Phi) agrees with lambda h_k(0) + mu_k of the re-polish to <= 4e-30 (the 30 printed
digits) for all k.

**Rigorous** (ball arithmetic with python-flint/arb, own K0/K1, exact rationals, exact Sturm counts): (i), (iii),
F > 0 on [0, a_k] u [b_k, oo) (hence S c (a_k, b_k)), (ii), the radius and weight numbers above.
**Assumed**: the trap lemma itself (specification, proved elsewhere); gt in the cone, i.e. B_gt > 0 on R and gt > 0 on
[0, oo), from the existing certificate (`certs/cert_d24.out`, certify.py on the same b_rational; not re-run here; its
cost ball was re-run and agrees); the closed Bessel form of Phi*Phi, DLMF 10.40.10, the ratio bounds in kbessel.py,
correctness of arb and of sympy's count_roots (same assumptions as certify.py).
**Not rigorous**: everything in D1/D2 (numerics; S scanned on a grid of step 1e-3 / 1e-2 plus local limits at the
nodes; KKT optimality of the stored optima; the trends).

## Pitfalls met

* The KKT check W(h_k) = lambda h_k(0) + mu_k is automatic (y = A^{-T} c and h_k'(u_j) = 0); it is informative only
  against independently computed W (stored mu, ball W of D3).
* trapcert.py's own weight bound (h(c) + |h'(c)| r + crude sup|h''| r^2/2 with coefficientwise |P|) is useless for the
  wide outer traps (k = 9..12: 1e-11..2e-5) because h_k has coefficients up to 1e22 with heavy cancellation; trapweight.py
  replaces it (Taylor models), giving ~mu_k/2. The `weight` fields inside `trapcert_k*.json` are the crude ones.
* `max depth 106` in `trapcert_k1.log` is a display bug (depth computed from the numerator of a non-dyadic width), fixed
  in fcert.py before k = 2..12 ran; the proof is unaffected.
* traps_d2/d4/d6 (deg 4, 8, 12) were produced by the first version of traps.py (sup of w by 201 samples, 700 bits, eta_1 scan with 6001
  points at 130 digits); deg >= 16 by the faster version (concavity test, 500 bits, 2401 points at 60 digits). The two
  versions were not run on the same degree; the changed parts only affect how sup_S w and eta_1 are located.
* `grid.pkl` (needed by traps.py and margins.py) was deleted to save space; regenerate with grid.py (8 min).
