# Certificates g = B(i d/du)(Phi*Phi): cost/g(0) versus deg B

All scripts run with `(ulimit -v 1500000; timeout 900 nice python <script>)` from this folder.

## The closed form used everywhere (derived, then checked numerically)

Phi(u) = 2 e^{u/2} sum_n (2y^2 - 3y) e^{-y}, y = pi n^2 e^{2u}. Each pair of theta terms in Phi*Phi gives
int z^{a-b-1} exp(-Az - C/z) dz = 2 (C/A)^{(a-b)/2} K_{a-b}(2 sqrt(AC)) with sqrt(AC) = pi m n e^u, so only N = mn enters:

    (Phi*Phi)(u) = 4 sum_{N>=1} d(N) (pi N)^{-1/2} X^{1/2} R(X),   X = pi N e^u,
    R(X) = (4X^4 + 9X^2) K0(2X) - 12 X^3 K1(2X),        d(N) = number of divisors.

d/du = X d/dX maps X^{1/2}[P0(X)K0(2X) + P1(X)K1(2X)] to the same form (rational polynomial recursions in `core.D`),
so g_j := (-1)^j (Phi*Phi)^{(2j)} (the function with transform Xi(t)^2 t^{2j}) and every g = sum b_j g_j are
finite Bessel sums with rational polynomial coefficients. Checks (`test_core.py`): int Phi = Xi(0) = 0.49712077818831...,
the closed form equals the direct convolution to 25 digits at u = 0, log 2, 1.3, is even (u = +-0.2, +-0.05), and
g_j(0) equals (1/2pi) int Xi^2 t^{2j} dt for j = 0..3 (20 digits).

Consequence (proved by this form): g(u) ~ B(2 pi i e^u) (Phi*Phi)(u) as u -> oo. B >= 0 on R needs b_d >= 0, g >= 0
at infinity needs (-1)^d b_d >= 0. So for deg B = 2d with d odd, b_d = 0: **deg 4k+2 gives exactly the deg 4k optimum**
(classes A and E). Only deg B = 0, 4, 8, ... matter.

## Scripts (one command each)

| script | does | output |
|---|---|---|
| `test_core.py`, `prec_check.py` | checks of the closed form (above); digits lost in the basis evaluation (10 of 50) | stdout |
| `b1_cost.py` | deliverable 1: B = 1 by both formulas | `b1_cost.out` |
| `arch_refine.py` | the Xi-side quadrature: 24 vs 48 Gauss-Legendre nodes per unit panel | `arch_refine.out` |
| `precompute.py 16 70` / `precompute.py 24 100 grid24.pkl` | g_j(u), j<=16 (j<=24), on u = i/400 (i<=2400) and u = log n (n<=403), 70 (100) digits, 2 processes, 2.5 (10) min | `grid.pkl`, `grid24.pkl` |
| `ximom.py J dps T file nodes` | Xi-side moments m_j, a_j | `xi_mom.pkl` (16, 50, 200, 24 nodes), `xi_mom48.pkl` (24, 80, 250, 48), `xi_mom64.pkl` (24, 110, 260, 64) |
| `lp_opt.py A|E d...` | LP on grids (HiGHS) + cutting planes + dual certificate; finds the structure (touch points) | `results_<C>_d<d>.json`, `lp*.out` |
| `polish.py A|E d [extend|re]` | exact KKT optimum by Newton (60 digits; deg >= 36 and class E deg >= 28 with `GRID=grid24.pkl DPS=90 XIMOM=xi_mom48.pkl`, `re` = re-polish a stored optimum, `MAXIT`) | `polish_<C>_d<d>.json`, `pol*.out` |
| `decay.py` | decay table and fits | `decay.out` |
| `report.py` | B coefficients, roots, touch points, g(log n)/sqrt n, cost carriers, per degree | `report.txt` |
| `archcheck.py` | every optimum: prime-sum cost recomputed from scratch vs archimedean formula (run as `XIMOM=xi_mom64.pkl DPS=110`) | `archcheck.out` |
| `kbessel.py` | own rigorous K0, K1 in arb (series / asymptotic with DLMF 10.40.10 bound) | - |
| `kbessel_test.py`, `arbtest.py`, `arbrepro.py` | K0, K1 against mpmath; arb's own bessel_k gives a useless ball at x = 153.04 | `arbtest.out` |
| `classE_check.py` | class E near u -> 0 (I0, first S values, min of g) | `classE_check.out` |
| `certify.py d eps U p` (env `SIG`, `PREC`) | rigorous certificate (ball arithmetic, 800 bits) for the class-A optimum of degree 2d | `certificate_A_d<d>.json`, `cert_d<d>.out` |
| `certify_all.sh` (env `DLIST`) | certify.py for deg 4..48, 1-5 min each | `certify_all.out` |

## Method

* LP (`lp_opt.py`): unknowns = coefficients of B in the basis of orthonormal polynomials for Xi(t)^2 dt (Cholesky of the
  Hankel matrix m_{j+k} = g_{j+k}(0)); rows scaled to max entry 1; constraints g(u) >= 0 (class A) or
  S(u) = sum_k e^{-ku/2} g(ku) >= 0 (class E) on u = i/400 (0 < u <= 6) and u = log n, B(t) >= 0 on t = i/50 (t <= 80)
  plus 400 geometric points to t = 80*1.02^399; tail sign b_d >= 0, (-1)^d b_d >= 0. Cutting planes at the local minima.
  HiGHS stops at relative violations ~1e-9 and fails (status 4) at deg 12, so the LP serves only to find the structure.
* Structure found at every degree 4K: the constraint function has exactly K double zeros u_1..u_K, B > 0 strictly
  (no real root), b_d > 0. Then b is fixed by the square linear system m.b = 1, F(u_k) = F'(u_k) = 0, and the optimum is
  the stationary point of cost(u_1..u_K); with y = A^{-T} c, d cost/d u_k = -z_k F''(u_k) where z_k is the multiplier of
  F'(u_k) = 0. `polish.py` solves z = 0 by Newton (60 digits). At the solution c = lambda m + sum mu_k F(u_k) with
  lambda = cost; all mu_k > 0 and F >= 0 on the grid, so it is the global optimum of the (convex) problem at that degree
  (LP duality), and lambda is also the lower bound for every admissible g of that degree (up to the residual |z|, 1e-40).
  Degree 4K+4 starts from the degree-4K touch points plus a scan for one new point.

## Results

### 1. B = 1 (g = Phi*Phi), `b1_cost.py`, `arch_refine.py`

g(0) = 0.319751811961621285119883398067, cost = 0.0240444194343550508370662176054, **cost/g(0) = 0.075197132697534236812**
(prime sum and archimedean formula agree to 8.6e-24 absolute with 24 Gauss-Legendre nodes per unit panel, and to all
40 digits with 48 nodes; 2 int g cosh(u/2) = 1/2 = 2 Xi(i/2)^2 exactly). n = 2 carries 98.5 % of the cost, n = 3 1.5 %.

### 2-3. Optima and decay (`decay.out`, `report.txt`)

The full tables (cost, ratio, touch points, fits, cost carriers, checks) are in `decay.out`; per degree the B coefficients,
roots of B, touch points, dual weights, g(log n)/sqrt n (n = 2..12) and S(log p) are in `report.txt`.

Main facts (all from `polish.py` / `decay.py`; every number below is in `decay.out` or `report.txt`):

* **Structure.** At deg B = 4K the optimum (both classes) has exactly K touch points: double zeros of g (class A) or of
  S (class E). They sit just ABOVE the first K prime powers 2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19 (class A) or the
  first K primes 2, 3, 5, ..., 23 (class E), and converge to log n very fast as the degree grows
  (class A, deg 48: e^{u_k} = 2, 3, ..., 13 to 10 digits, 16.0000005, 17.000003, 19.00015). The newest touch point
  sits between prime powers (e.g. e^u = 4.357 at deg 12, 16.255 at deg 40) until the next degree pins it.
  No touch point at 6, 10, 12 or any non-prime-power: g(log 6)/sqrt 6 levels off at 1.4e-7 (deg 36-48), far above
  g(log n)/sqrt n at the touched prime powers (< 1e-40 for n <= 11 at deg 48); 6 costs nothing, so it is ignored. The dual weights converge to the cost weights: mu_k -> 2 Lambda(n_k)/sqrt(n_k) (class A),
  2 log p_k (class E).
* **B >= 0 is never active**: at every degree all roots of B are complex (`report.txt`), B > 0 on R, b_d > 0. Its only
  effect is the parity rule above (deg 4k+2 = deg 4k).
* **Class E vs A**: class E lets g go slightly negative (e.g. g = -6e-5 near log 2, `classE_check.out`) because only
  S(log p) enters the cost (cost = 2 sum_p log p S(log p)); powers p^k come for free, so class E needs touch points only
  at primes and decays much faster. Near u -> 0, S(u) ~ I0/u with I0 = int_0^oo e^{-v/2} g = 0.256 > 0.
* **Who carries the remaining cost**: the first untouched prime power(s) (class A: deg 32 -> n = 13 96 %, deg 48 -> n = 23
  99.9 %); class E: the first untouched prime (> 99.9 % from deg 12 on).
* **Decay law.** No stall. Class A: log10 cost per +4 degrees drops by 1.08 1.22 0.90 1.11 1.91 1.42 1.10 2.10 3.98 2.21
  2.47 4.99 (deg 0 -> 48): irregular, set by the gaps between prime powers, and accelerating; quadratic fit in K = deg/4:
  c = -0.13 decades/K^2 (rms 0.57, linear rms 1.71). Against the first untouched prime power n_{K+1}:
  log10 cost = 1.74 - 1.15 n_{K+1} (rms 0.69), local slope -1.03 (first half) -> -1.26 (second half). Since
  n_K ~ K log K this is super-geometric in the degree, but only through the growth of n_K (and a slowly steepening
  slope). Class E: drops 1.09 1.91 2.78 6.15 4.18 7.09 4.62 7.41 (deg 0 -> 32), log10 cost ~ 3.9 - 1.69 p_{K+1}:
  clearly super-geometric. This is evidence for the conjecture inf = 0, not a proof.
* **Precision limits.** Class E deg 36 (9.3e-49) is a feasible point whose stationarity is resolved only to
  |z| = 1.4e-50 at 90 digits: an upper bound, not an established optimum; class E deg 40 was not attempted.
  Class A up to deg 48 is fine (|z|/cost < 1e-26).

### Independent check of every optimum (`archcheck.out`)

Prime-sum cost (recomputed from scratch, 110 digits) against the archimedean formula with 64-node, 110-digit Xi
moments (`xi_mom64.pkl`): relative agreement 1e-58 (deg 4) to 1e-34 (class A deg 48); class E 1e-58 to 1e-24
(deg 32). g(0) closed form against Xi side: 1e-105. (Class E deg 36: the stored 40-digit b gives a different polynomial
at the 1e-41 level, so its check confirms the identity, not the optimum.)

### 4. Rigorous certificates (`certify.py`, `certify_all.sh`, `certificate_A_d*.json`, `cert_d*.out`)

For each class-A optimum: b_j rounded to exact rationals (25 digits for deg <= 32, 60 digits for deg >= 36), eps =
10^(floor(log10 cost) - 3) added to b_0 (g + eps Phi*Phi, B + eps). Proved in ball arithmetic (python-flint arb, 800
bits; K0, K1 from `kbessel.py`, not from arb's bessel_k):

1. B > 0 on R: q(s) = B(sqrt s) has no root in [0, oo) (exact Sturm count over Q, sympy) and q(0) > 0.
2. g > 0 on [0, U] (U = 5, resp. 6 for deg >= 36): bisection from width 1/64; on each interval a Taylor polynomial of
   order 16 at the centre with ball derivatives, exact minimum of the quadratic part, terms 3..16 bounded by
   |D^k g(c)| r^k/k!, remainder by a monotone bound of |D^17 g| on the interval. Up to ~2000 intervals, depth <= 43.
3. g > 0 on [U, oo): asymptotic expansions of K0(2X), K1(2X) with 6 terms and the DLMF 10.40.10 remainder (real
   order, x > 0: remainder below the first neglected term), giving F(X)/X^M >= Flow > 0 for X >= pi e^U; the N >= 2
   terms are at most 1e-360 of the N = 1 term there.
4. cost = 2 sum_n Lambda(n) n^{-1/2} g(log n): balls for n <= 80 (X = pi N n exactly), explicit geometric tails in N
   and n; g(0) as a ball.

**Theorem (certified).** For each degree below there is an explicit even polynomial B > 0 on R with rational
coefficients (`certificate_A_d<d>.json`, key `b_rational`) such that g = B(i d/du)(Phi*Phi) is > 0 on R and
2 sum_n Lambda(n) n^{-1/2} g(log n) <= c g(0):

| deg B | certified c | optimum (KKT) |
|---|---|---|
| 4 | 6.19402005455e-3 | 6.193997990676e-3 |
| 8 | 3.70719233613e-4 | 3.707168410241e-4 |
| 12 | 4.64624013887e-5 | 4.646216109301e-5 |
| 16 | 3.58920272265e-6 | 3.589178679378e-6 |
| 20 | 4.43833859341e-8 | 4.438314549002e-8 |
| 24 | 1.70453120117e-9 | 1.704507156751e-9 |
| 28 | 1.35254770865e-10 | 1.352523664229e-10 |
| 32 | 1.07689720773e-12 | 1.076873163309e-12 |
| 36 | 1.12161611521e-16 | 1.121592070789e-16 |
| 40 | 6.85355154306e-19 | 6.853527498631e-19 |
| 44 | 2.33130821521e-21 | 2.331284170788e-21 |
| 48 | **2.37581468402e-26** | 2.375812279573e-26 |

(c = upper end of the arb ball for cost/g(0), rounded up to 12 digits. Note: the strings `certified_upper_bound` in the json files and the `certified:` column of `certify_all.out` are the same number rounded to NEAREST, so they can lie below the ball in the last digit; the rigorous object is the ball `cost_over_g0`, and the table here is the one to quote; eps per `certificate_A_d*.json`; deg 48 was run with
eps = 1e-30, the others with the rule in `certify_all.sh`. Times 23 s (deg 4) to 249 s (deg 48), 974 to 3524 intervals.)

**What is rigorous.** Steps 1-4 above, given: (a) the closed form of Phi*Phi (an exact derivation from the integral
2 (C/A)^{nu/2} K_nu(2 sqrt(AC)) = int z^{nu-1} e^{-Az-C/z} dz; it also agrees with direct quadrature of Phi*Phi to 25
digits), (b) DLMF 10.40.10 (remainder bound of the K_nu asymptotic series for real nu, x > 0) and the elementary
ratio bounds for the K0/K1 power series written in `kbessel.py`, (c) correctness of python-flint/arb ball arithmetic
(elementary functions only) and of sympy's Sturm root count. The identity cost(g) = W(g) (explicit formula) is a
theorem for this g (even, super-exponential decay, ghat = Xi^2 B of exponential decay); it is used only to read the
bound as a bound on W(g) and is checked numerically in `archcheck.out`.

**Not rigorous.** Optimality of the KKT points (numerical: multipliers > 0, stationarity residual tiny, feasibility on
the grid u = i/400 plus exact double zeros); the lower bounds lambda (same residual caveat); everything about class E
(grid feasibility only, plus the u -> 0 asymptotics); the decay law (finitely many degrees).

**Pitfall found.** python-flint 0.9.0 `arb.bessel_k(0)` at x = 6 pi e^{67/32} = 153.04 (input a 320-bit ball) returns
[+/- 5.2e-28] for a value 3.7e-68 (valid but useless, `arbrepro.py`); our `kbessel.K01` gives full precision.
In the series branch (x < 120) K01 amplifies the input radius by ~e^{2x}, so `certify.py` runs at 800 bits.
