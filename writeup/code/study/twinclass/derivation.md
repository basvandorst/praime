# Twin counts across residue classes: the linear-response term, the remainder, and a law for the spread

Folder: `writeup/code/study/twinclass/ (jobs run in a scratch folder) ` (called "the folder").
All code is new (own sieve, own singular-series code); the repository scripts were read for definitions only and
compared afterwards. Labels: **proved** (algebra, given the stated model), **computed**, **heuristic**, **conjecture**.
"Classical" marks what is in print or in the repository before tonight; "here" marks what this note derives.

## 0. The statistic and its registered prediction (rounds 41-44, classical in the repository)

For a window [L0, L0 + W), a squarefree modulus Q (210, 2310, 30030) and an even gap h (h = 2: twins), let

    X_b = #{n in [L0, L0+W): n = b mod Q, n and n + h prime},   b in A_Q(h) = {b : gcd(b(b+h), Q) = 1},  K = |A_Q(h)|.

Under Hardy-Littlewood (HL) every admissible class has the same density (the local factor at p | Q is (1-1/p)^-2 for
every admissible b), so the spread needs no law for the mean:

    u2_X = Var_b(X_b) (ddof 1) / mean_b(X_b)            (round 42: "measured")

**What the four-point prediction of rounds 42-43 contains** (r42.py `predict`, re-implemented here as
`fourpoint.decompose`, field `vX`):
1. One class only: the points n_k = b + Qk, k = 0..N-1, N = W // Q, density r2(n) = c_Q(h)/(ln n ln(n+h)) with
   c_Q(h) = prod_{p|Q}(1-1/p)^-2 prod_{p not | Q}(1 - nu_p(h)/p)/(1-1/p)^2 (nu_p(h) = #{0, h} mod p), profile at
   m_k = L0 + Q/2 + Qk.
2. Var X_b = sum_k r2(1 - r2) + 2 sum_{i>=1} w(i) (rho(i) - 1), w(i) = sum_k r2_k r2_{k+i}, rho(i) =
   S'({0,h,Qi,Qi+h})/S'({0,h})^2 with the factors at p | Q removed (they are 1 after class conditioning).
3. The sum runs over all pairs of points of the class in the window (distances Qi up to W): no truncation in distance.
4. Approximations: (a) the product over p stops at 2e7 (repo) - here it is kept in full (exponential-integral tail and
   every special prime above 2e7; changes v by -6.6e-4 relative at train mod 210, -6e-5 mod 2310, -4e-5 mod 30030,
   the "0.1-0.6 %" of the ninth review); (b) the density is the HL leading term c/(ln n ln(n+h)), no lower-order
   terms; (c) **the cross-class covariances are dropped**: the expected sample variance is
   v - mean_{b != b'} Cov(X_b, X_b')/m; (d) **the noise** is taken as chi2_{K-1}/(K-1), i.e. classes independent
   and Gaussian; (e) it is an ensemble value: no information on the realized prime errors of the window is used.
5. What it leaves out: (c), (d), (e); terms beyond HL (no explicit formula for pairs exists); prime powers (counts of
   primes only).

Reproduction (computed, `check_fp.py`, `run_all.py`; own sieve checked against sympy cell by cell, `test_sieve.py`):
every twin count of r42/r43 is reproduced exactly (381,332 / 2,984,194 / 23,988,173 / 124,843 / 86,479 / 22,013;
prime counts equal pi differences exactly), every measured u2 to 4 decimals (e.g. [1e8, 1e9): 0.1421 mod 210,
0.4214 mod 2310, 0.5996 mod 30030), and the predictions to the truncation effect: 0.43723 vs 0.43751 (train, 210),
0.54951 vs 0.54954 (train, 2310), 0.64178 vs 0.64181 ([1e8,1e9), 30030); [1e9, 1e10) mod 210: 0.4527 (review9's
chunked value 0.4492; mine 0.4492).


## 1. The linear-response term, from first principles

**1.1 The split (heuristic, round 155 in class form).** Per point n of class b let e_n = 1_{n prime} - r1(n). To first
order in the prime fluctuations, the pair indicator responds to the two members' fluctuations with the product rule:

    1_n 1_{n+h} - r2(n) = kappa(n) [e_n + e_{n+h}] + (bilinear rest),    kappa(n) = r2(n)/r1(n) = S(h) phi(Q)/(A_Q(h) ln(n+h)).

Summed over the class: X_b - E X_b = kappa Z_b + R_b, with Z_b the prime count of class b (member n) plus the prime
count of class b + h (member n + h), each prime weighted by 1/ln p (measured as sum 1/ln p, `tcsieve.py` WP), and R_b
the rest. In Lambda weights this is exactly round 155's class form E_b = (S phi/A)[Delta(x; q, b) + Delta(x; q, b+h)]
(c_Q(h) = S(h) Q/A_Q(h): **proved**, local factors p by p). The coefficient is fixed by the law, not fitted.

**1.2 The Hardy-Littlewood (HL) moments of the split (proved as algebra inside the HL ensemble; HL itself a
conjecture).** With all singular series class-conditioned (factor 1 at p | Q) and d = Qk:

    Var X      = sum r2(1-r2) + 2 sum_{k>=1} w(k) [rho(k) - 1]                              (four-point, round 42)
    Var(kZ)    = sum [r2^2/r1 (1-r1) + r2^2/r1'(1-r1') + 2 r2^2 (r2 - r1 r1')/(r1 r1')]
                 + 2 sum_k w(k) [2(s0(k)-1) + (s+(k)-1) + (s-(k)-1)]                     (pairs of single primes)
    Cov(X, kZ) = sum r2 [kappa(1-r1) + kappa'(1-r1')] + 2 sum_k w(k) [(tau1(k)-1) + (tau2(k)-1)]   (three-point)

s0 = S'({0, d}), s+- = S'({0, d +- h}) (the cross pairs (n, n'+h), (n+h, n')), tau1 = S'({0,h,d})/S'({0,h}),
tau2 = S'({0,h,d+h})/S'({0,h}). Every local factor has mean exactly 1 over k mod p (checked symbolically for each
series: e.g. four-point, p not | h: [(1-2/p) + 2(1-3/p) + (p-3)(1-4/p)]/(p(1-2/p)^2) = 1).

**1.3 Is the linear term inside the four-point sums? Yes (proved).** The four-point sum is diag + 2 sum w (rho - 1).
Write rho - 1 = lin + G with lin = 2(s0-1) + (s+-1) + (s--1). The off-diagonal of Var(kZ) is exactly 2 sum w lin
(the identity kappa kappa' r1 r1' = r2 r2'). So the four-point sums of rounds 42-43 *contain* the linear response
in full - in distance (all pairs of the window), in modulus (the class-conditioned series), with no diagonal
approximation. Locally, at each prime, rho_p = prod over the four cross pairs of their pair factors to first order in
1/p (checked case by case: p | d: 1/(1-2/p) = 1 + 2/p + ...; cross pairs (0,d), (h,d+h) give 1 + 1/p each), so the
cross-pair (linear) terms are the first-order part of rho and G is second order per prime - but G is not small in sum.

**1.4 The regression coefficient (computed).** beta_HL = Cov(X, kZ)/Var(kZ) = 1.0000 +- 0.0005 in every window and
modulus tried (18 cases, `all_decomp.json`), i.e. inside HL the three-point sums give exactly the product-rule
coefficient. (A first version had beta_HL drifting to 1.30 at N = 8e7; the cause was the tail constant of the
three-point product, -2/p^2 per prime, entered as -1/p^2: an error of 3e-9 in a log, which the window multiplies by the
mean count per class, up to 2.8e6. Fixed; every number here is after the fix.)

**1.5 Sizes (computed, the HL part per class).** [1e8, 1e9), twins:

| Q | four-point v_X | diag | off total | off linear (= prime rigidity passed on) | off rest G | P_lin = Var(kZ)/m | P_rem |
|---|---|---|---|---|---|---|---|
| 210 | 0.4438 | 0.953 | -0.510 | -0.301 | -0.209 | 0.027 | 0.417 |
| 2310 | 0.5427 | 0.943 | -0.400 | -0.273 | -0.127 | 0.081 | 0.462 |
| 30030 | 0.6418 | 0.933 | -0.291 | -0.228 | -0.063 | 0.149 | 0.493 |

Reading: about 60-78 % of the twins' rigidity is the primes' own class rigidity passed on through both members
(off linear ~ -2 kbar (1 - v_Y), v_Y the prime class law of rounds 36/135). But the linear part's *variance* is
small: P_lin = 2 kbar v_Y + small terms, because v_Y (0.065 / 0.17 / 0.29) is small: its Poisson part 2 kbar and its
rigidity -2 kbar (1 - v_Y) nearly cancel. Closed form (round 135 law for one class, **heuristic**, checked against
the exact sums: equal to 4 decimals for A, B, C; 10-15 % high for the decades): P_lin ~ 2 kbar [1 - (ln(W/Q) + S_Q - B)/lambda].

**1.6 The remainder has a ln^2 law (derived here, heuristic; computed check).** Ramanujan expansion of the twin
density of one class: f_p(n) = 1_{p not | n(n+h)}/(1 - nu_p/p) has Fourier coefficients f^_p(a) = -(1 + e(ah/p))/(p-2)
(p not | Qh), so

    rho(k) = prod_{p not | Q} [1 + (2 c_p(k) + c_p(Qk+h) + c_p(Qk-h))/(p-2)^2],   c_p = Ramanujan sum,

(checked: generic 1 - 4/(p-2)^2, p | k: p/(p-2), p | Qk+-h: 1 + (p-4)/(p-2)^2 = (1-3/p)/(1-2/p)^2 - all equal to the direct factors).
Summing over k <= J, each squarefree q coprime to Q contributes on average -1/2 times its spectral mass
prod_{p|q} m_p with m_p = sum_a |f^_p(a)|^2 = 2/(p-2) (the -1/2 is the mean of sum_{k<=J} e(k theta), theta not in Z,
as in Montgomery-Soundararajan's sum of S(d) - 1). Since m_p ~ 2/p, the Dirichlet series behaves like zeta(1+s)^2:

    F(J) = sum_{k<=J} (rho(k) - 1) = -(H_Q/4) ln^2 J + O(ln J),   H_Q = prod_{p|Q}(1-1/p)^2 prod_{p not | Q}(1+m_p)(1-1/p)^2
                                                                     = (phi(Q)/Q)^2 / prod_{p not | Q}(1 - 1/(p-1)^2).

For single primes m_p = 1/(p-1) ~ 1/p gives zeta(1+s)^1 and the familiar -(phi/2Q) ln J per pair: the twin density
has twice the spectral mass at each prime, because either member can carry the frequency. The pure assignments
(every prime of q on the same member) are the linear part (ln J); the mixed ones (some primes on n, some on n + h)
are G and carry the ln^2. Check (computed, `partial_sums.py`, `asym_fit.py`, exact sums to J = 4e6):
fitted a = 0.0154 / 0.0103 / 0.0094 against derived H_Q/4 = 0.0135 / 0.0111 / 0.0094 (Q = 210 / 2310 / 30030); the
linear part's slope 0.218 / 0.209 / 0.198 against phi/Q = 0.229 / 0.208 / 0.192.

Since c_Q(2) = (Q/phi)^2 prod_{p not | Q}(1 - 1/(p-1)^2), **c_Q H_Q = 1** (proved), so r2 H_Q = 1/lambda^2 and the
leading term of the twin class rigidity is universal:

    v_X = 1 - r2 - (L^2 - 2L + 2)/(2 lambda^2) + 2 r2 [b_Q (L - 1) + c_Q'],    L = ln(W/Q), lambda = ln x  (window short against x)

with b_Q, c_Q' the next two constants of the HL series (-0.1725, -0.1731 | -0.1669, 0.0616 | -0.1561, 0.1367, read off
the exact partial sums with a fixed; constants of the conjecture, not fitted to any count). Check against the exact
sums: A 0.6931 vs 0.6910, B 0.7869 vs 0.7854, C 0.8769 vs 0.8763 (mod 210); for the decades the same F(J) summed by
parts against the exact profile gives 0.4419 / 0.4507 / 0.4573 vs 0.4372 / 0.4438 / 0.4492 (210), i.e. +0.005..0.008
(`closed_form.py`). Leading term alone: 1 - L^2/(2 lambda^2) is far off (0.71 vs 0.44): the ln terms matter.

## 2. Measured decomposition (retrospective: every window below was seen; computed)

`run_all.py` -> `all_decomp.json`; measured Z_b uses the law's kappa (not fitted); beta = Cov_b(X, kZ)/Var_b(kZ)
across the classes is the class-level regression coefficient (law: beta_HL = 1.000). Predicted values include the two
small corrections of section 0 (c): the mean cross-class covariance of the linear part (round 135 closed form,
-Lambda(m)/(phi(m) lambda)) and the Chebyshev bias of the prime squares (-(N(c)-1) li(sqrt x)/(2 phi) per class; it
gives the measured prime spread 0.0772 against 0.0725 = 0.0650 (HL, same class) + 0.0035 (cross-class) + 0.0040 (bias) at [1e8, 1e9) mod 210). For the twin spread both are below 0.005.

| window | Q | u2X meas / pred (ratio) | u2Z meas / pred | u2R meas / pred (ratio) | beta (nominal se) |
|---|---|---|---|---|---|
| 1e7-1e8 | 210 | 0.398 / 0.444 (0.90) | 0.032 / 0.041 | 0.363 / 0.403 (0.90) | 1.06 (0.84) |
| 1e7-1e8 | 2310 | 0.404 / 0.551 (0.73) | 0.116 / 0.105 | 0.322 / 0.446 (0.72) | 0.85 (0.18) |
| 1e7-1e8 | 30030 | 0.656 / 0.660 (0.99) | 0.191 / 0.191 | 0.458 / 0.470 (0.98) | 1.02 (0.04) |
| 1e8-1e9 | 210 | 0.142 / 0.449 (0.32) | 0.027 / 0.032 | 0.113 / 0.417 (0.27) | 1.04 (0.96) |
| 1e8-1e9 | 2310 | 0.421 / 0.544 (0.77) | 0.084 / 0.082 | 0.371 / 0.462 (0.80) | 0.80 (0.20) |
| 1e8-1e9 | 30030 | 0.600 / 0.642 (0.93) | 0.147 / 0.149 | 0.449 / 0.493 (0.91) | 1.01 (0.05) |
| 1e9-1e10 | 210 | 0.199 / 0.453 (0.44) | 0.037 / 0.026 | 0.252 / 0.427 (0.59) | -0.21 (1.09) |
| 1e9-1e10 | 2310 | 0.513 / 0.539 (0.95) | 0.064 / 0.066 | 0.453 / 0.472 (0.96) | 0.97 (0.23) |
| 1e9-1e10 | 30030 | 0.605 / 0.627 (0.96) | 0.124 / 0.120 | 0.490 / 0.507 (0.97) | 0.96 (0.05) |
| A | 210 | 0.181 / 0.694 (0.26) | 0.074 / 0.118 | 0.140 / 0.575 (0.24) | 0.77 (0.59) |
| A | 2310 | 0.593 / 0.755 (0.79) | 0.189 / 0.167 | 0.524 / 0.588 (0.89) | 0.68 (0.16) |
| A | 30030 | 0.792 / 0.817 (0.97) | 0.230 / 0.227 | 0.578 / 0.590 (0.98) | 0.97 (0.04) |
| B | 210 | 1.000 / 0.787 (1.27) | 0.129 / 0.133 | 0.999 / 0.654 (1.53) | 0.50 (0.59) |
| B | 2310 | 0.870 / 0.830 (1.05) | 0.181 / 0.173 | 0.596 / 0.657 (0.91) | 1.26 (0.17) |
| B | 30030 | 0.826 / 0.873 (0.95) | 0.203 / 0.220 | 0.618 / 0.653 (0.95) | 1.01 (0.04) |
| C | 210 | 0.527 / 0.878 (0.60) | 0.097 / 0.141 | 0.431 / 0.737 (0.59) | 0.99 (0.61) |
| C | 2310 | 1.024 / 0.903 (1.13) | 0.162 / 0.172 | 0.774 / 0.731 (1.06) | 1.27 (0.18) |
| C | 30030 | 0.945 / 0.929 (1.02) | 0.215 / 0.208 | 0.762 / 0.722 (1.06) | 0.93 (0.05) |

What this shows:
1. **The linear-response part is present at its predicted size** in every window and modulus (u2Z measured/predicted
   within noise; at mod 30030 the ratio over 27 sub-windows is 1.004, sd 0.040).
2. **The regression coefficient is the law's**: mod 30030 (1485 classes) beta = 0.986 +- 0.019 over the six windows,
   0.9974 (sd 0.048) over 27 sub-windows. "Twin class counts independent of the class prime counts" (beta = 0) is
   excluded at 50 errors; "one member only" (0.5) at 25. This is round 155's class form read in fixed windows across
   classes (a projection not used before), not new structure.
3. **The linear part does not explain the evenness**: it is 6-9 % of the spread at mod 210 in the decades (16-17 % in the short windows), 12-22 % at 2310, 19-29 % at
   30030, and it is *not* low (mod 210, [1e8, 1e9): 0.027 against 0.032). The whole shortfall sits in the remainder R:
   0.113 against 0.417 at [1e8, 1e9) mod 210. Removing the linear part sharpens the fine-level deficit there
   (X: 0.934, lower-tail p 0.034; R: 0.912, p 0.007), because the linear part fits and only adds noise.

Level split (as review9; ratios to the prediction, per window):

| window | X: level 1 (210) | level 2 (2310 in 210) | level 3 (30030 in 2310) | R: level 1 | level 2 | level 3 |
|---|---|---|---|---|---|---|
| 1e7-1e8 | 0.897 | 0.718 | 1.014 | 0.901 | 0.702 | 0.999 |
| 1e8-1e9 | 0.316 | 0.818 | 0.947 | 0.271 | 0.859 | 0.922 |
| 1e9-1e10 | 0.440 | 1.002 | 0.965 | 0.590 | 0.997 | 0.966 |
| A | 0.260 | 0.841 | 0.987 | 0.244 | 0.965 | 0.990 |
| B | 1.270 | 1.023 | 0.937 | 1.527 | 0.836 | 0.950 |
| C | 0.601 | 1.194 | 1.006 | 0.585 | 1.114 | 1.056 |
| nominal se per window | 0.378 | 0.129 | 0.038 | 0.378 | 0.129 | 0.038 |

Sub-windows and noise (`levels.py`): over the 27 sub-windows (the nine equal parts of each decade, lengths 1e6..1e9
with W/x from 0.11 to 1) the remainder fits at every level, and the nominal chi2 noise is right:

    beta (mod 30030), six whole windows, inverse-variance: 0.9857 +- 0.0186
    27 sub-windows: beta 0.9974 (sd 0.0481, nominal se 0.0414); R levels 0.966 +- 0.073, 0.990 +- 0.025, 1.002 +- 0.007; X levels 1.000, 0.973, 1.002
    sd/nominal of the 27 sub-window ratios: X 1.09 1.06 1.22; R 0.94 1.04 0.99

So the decade deficit is the **cross part**: decade spread = within (sum over the nine sub-windows, fits) + cross
(products of class deviations of different sub-windows). Twins against the other pairs (n, n + h) on the same three
decades (`cross_h.py`, z = (cross measured - cross HL)/nominal sd, sd = 2 v_within sqrt((1 - sum m_s^2/m^2)/(2(K-1))),
classes and sub-windows taken as independent):

| pair gap h | cross z of R, mod 210 (3 decades; joint) | mod 2310 | mod 30030 |
|---|---|---|---|
| 2 | -0.41 -1.47 -0.86; -1.58 | -2.36 -1.26 -0.16; -2.18 | -1.06 -2.23 -0.77; -2.34 |
| 4 | -0.76 -0.27 +0.44; -0.34 | -0.95 -0.63 +1.58; -0.00 | -1.38 +0.39 +0.17; -0.47 |
| 6 | -0.46 +0.27 -0.06; -0.14 | -0.34 +0.77 -0.90; -0.27 | +0.55 +0.67 +1.47; 1.56 |
| 8 | +0.22 +0.26 +0.53; 0.58 | -1.68 +0.22 +0.74; -0.42 | +0.09 -0.13 +0.14; 0.06 |
| 10 | -0.32 -0.50 -0.42; -0.72 | -1.40 -1.45 +0.99; -1.08 | -0.00 +0.99 +2.56; 2.05 |
| 12 | -0.96 -0.96 -1.18; -1.79 | -0.25 -1.88 -0.07; -1.27 | +0.34 -0.86 +0.45; -0.04 |

h = 4..12: 15 joint z, mean -0.15, sd 0.99
twins Q=210: cross meas -0.133 -0.316 -0.212, HL -0.065 -0.064 -0.062, sd 0.167 0.171 0.174
twins Q=2310: cross meas -0.179 -0.123 -0.059, HL -0.045 -0.048 -0.049, sd 0.057 0.059 0.060
twins Q=30030: cross meas -0.043 -0.072 -0.050, HL -0.025 -0.032 -0.036, sd 0.017 0.018 0.019


Reading: for twins the class deviations of the remainder in different sub-windows of a decade anti-correlate
about 2.6 times on average (1.2 to 5 per decade) more strongly than the four-point sums say (HL cross -0.065 / -0.047 / -0.031 of Poisson at 210 / 2310 / 30030;
measured -0.22 / -0.12 / -0.055 on average). The same statistic for h = 4, 6, 8, 10, 12 (computed here for the first
time, so not selected) scatters like N(0, 1). The twins were selected (they started rounds 41-44), the three moduli
are nested and the three decades are the only ones; joint -2.2 and -2.3 at the two fine levels is a 2-3 sigma
fluctuation after selection, or a twin-specific long-range term.

Retrospective on the fifty-seventh review's counts (`rc_retro.py`, log-weighted sums mod 420 folded to 210; weights
ln p ln(p+2), normalised by sum w^2/sum w; read-only files): (1e8, 1e9] reproduces the count-based values (X 0.29,
R 0.25 of prediction); (1e9, 2e9] 1.01 / 0.99; (2e9, 1e10] 0.36 / 0.55; (1e10, 2e10] - a window where this statistic
had not been computed - 0.40 / 0.48 (lower-tail p of R 0.056); beta at 15 classes is uninformative (se ~1).

## 3. Answer to the question (does the linear term explain the evenness?)

**No.** The linear-response term is algebraically inside the four-point sums of rounds 42-43 (section 1.3), it is
measured at its predicted size and with the predicted coefficient, and it is a small part of the spread. Conditioning
on the realized prime errors of the window changes the predicted spread by kappa^2 (s2_Z - E s2_Z)/m: at most 0.045 in
the seen windows (mod 210, A and C, both downwards: A 0.694 -> 0.650 against 0.181 measured), at most 0.02 elsewhere. Nothing of the linear part was cut off: the sums run over all pairs of the window, the series are
the class-conditioned ones, and the diagonal is exact; the only approximations of rounds 42-44 that matter at all are
the dropped cross-class covariance (< 0.005) and the product cut at 2e7 (< 0.003). **The remainder X is the
non-linear (pair-pair, "mixed-frequency") part of the class counts**: its HL rigidity is the ln^2 term of section 1.6;
in sub-windows it is measured at 0.966 +- 0.073 / 0.990 +- 0.025 / 1.002 +- 0.007 of its HL value (levels 1-3), and
in the three seen decades its class deviations between sub-windows anti-correlate more than HL says (cross part
-0.22 vs -0.065 at 210, -0.12 vs -0.047 at 2310, -0.055 vs -0.031 at 30030; z joint -1.6, -2.2, -2.3), for twins only.
Size of X: in the decade windows the remainder spread is low by 0.04-0.30 (210), 0.02-0.12 (2310), 0.01-0.04 (30030)
of Poisson; summed over the three decades it is the whole of the "extra evenness".

## 4. The law for the spread (no parameter fitted to counts)

For pairs (n, n + h) in the admissible classes mod a squarefree Q, window [x, x + W), N = W/Q points per class:

    E[u2_X] = P_lin + P_rem
    P_lin = Var(kappa Z)/m  (exact two-point sums; closed form 2 kbar [1 - (ln N + S_Q - B)/lambda], round 135's
            law for each member's class, valid for W << x, 10-15 % high for decade windows)
    P_rem = [diag_R + 2 sum_k w(k) (G(k) + three-point terms)]/m, G = rho - 1 - lin    (exact)
    closed form (W << x): E[u2_X] = 1 - r2 - (L^2 - 2L + 2)/(2 lambda^2) + 2 r2 [b_Q (L - 1) + c'_Q],  L = ln N,
            r2 = c_Q(h)/lambda^2, leading coefficient 1/2 universal (c_Q H_Q = 1), b_Q, c'_Q constants of the HL series.
    noise: chi2_{K-1}/(K-1) (verified on 27 sub-windows: sd/nominal 0.94-1.22); conditional on the window's primes,
           E[u2_X | primes] = kappa^2 s2_Z(measured)/m + P_rem with the noise of R alone.

Approximations: HL k-tuple densities with the leading 1/(ln n ln(n+h)) profile (conjecture); classes treated as
exchangeable for the cross-class term (closed form of round 135, heuristic); Chebyshev bias of prime squares added for
Z (classical); prime powers ignored (counts of primes); for long windows the closed form needs the summation by parts
over the profile (accuracy +0.005..0.008 at mod 210, `closed_form.py`).

## 5. Fresh blind test (designed, frozen, not run)

Window: the decade D = [2e10, 2e11) and its nine sub-windows S_s = [2e10 (1 + s), 2e10 (2 + s)), the same shape as the
seen decades (W/x = 9; sub-windows W/x 1 to 0.11). Twin class counts mod 210/2310/30030 were never computed above
2e10 (checked: STUDY.md, LOGIC.md, MERSENNE.md mentions and preregister/*.txt; in (2e9, 2e10] the review's counts mod
420 exist). Stretches inside D used for other statistics: [2.9e10, 3.24e10) (round 139, prime class tables),
[1e11, 1.1e11) (round 64, twin gap laws), near 1.5e11 (round 140, a response stretch): other statistics, about 7 % of D.
Predictions: `freeze.py fresh` -> `predictions_twinclass.json` (sha256 in `frozen.sha256`); exact HL sums for 2310 and
30030, closed form for 210 (reported only). Test: `test_twinclass.py` (refuses other hashes of the predictions and of
tcsieve.py, measure.py, fourpoint.py; checked by a deliberate one-byte change). Dry run on the seen decade
[1e8, 1e9) (`test_dry.out`): P1 beta 1.0235 vs 0.9991 +- 0.0162 (z 1.5), P2 0.9965 +- 0.0128, P3 0.9948 +- 0.0132,
P4 Z = -2.47 (OPEN), as known for that window.

Readings, what each rejects, and power (errors from the 27 seen sub-windows):
* P1 (law): class-level response beta, mod 30030, mean over S0..S8; se 0.016; rejects "independent" (beta = 0) at about
  60 se and "one member" (0.5) at 30 se. The round-43 comparator makes no statement on beta, so this reading tests the
  linear-response part of the law against the independence model implicit in rounds 41-44's noise model.
* P2 (law): remainder at the fine level (30030 within 2310), mean over S0..S8, 1 +- 0.0128; under independence of X
  and Z the remainder would be at about 1.8: rejected at 60 se. The comparator again has no R.
* P3 (law): size of the linear part, u2Z/pZ at 30030, 1 +- 0.0132; rejects Poisson prime class counts (ratio ~1.8) at
  60 se. Known in substance (round 135's class law, passed on by 2 kappa).
* P4 (the open point): Z = (z_2310 + z_30030)/sqrt(2) of the remainder's cross part on D. HL (= this law = round 43's
  comparator, which both contain it): Z ~ N(0, 1). The seen twin excess A_twin: E[Z] = -1.85. Reading: Z <= -3 confirms
  a long-range anti-correlation beyond the four-point sums (falsifies the HL remainder for twins in decade windows;
  probability 0.0013 under HL, about 0.13 under A_twin); Z >= -0.92: chance favoured (probability 0.82 under HL, 0.18
  under A_twin); otherwise open. Power is limited by there being one fresh decade: a clean decision between HL and
  A_twin needs about four such decades (Z separation 1.85 per decade).
* Reported only: X and R ratios at the three levels for D and each S_s (both laws predict 1 for X).
Falsifier of the law: P1, P2 or P3 outside 3 errors. A reading that the comparator also passes: P4 under HL, and every
u2_X ratio; they are reported as tests of the four-point law itself, not of the linear-response decomposition.

Commands (from the folder, foreground, one at a time; timings from a seen stretch at 1e11: 2e9 numbers in 6.5 s, 95 MB):

    for s in 0 1 2 3 4 5 6 7 8; do (ulimit -v 1000000; timeout 600 nice /home/bas/praime/.venv/bin/python -I -B tcsieve.py $((20000000000*(1+s))) 20000000000 fresh_S$s.npz 1); done   # about 65 s each, 10 min total, under 100 MB
    sha256sum fresh_S*.npz          # record
    (ulimit -v 1000000; timeout 600 nice /home/bas/praime/.venv/bin/python -I -B test_twinclass.py)   # a few seconds; writes test_result.json

## 6. Prior art (mathematical queries only)

| statement | verdict | what was opened |
|---|---|---|
| class counts of single primes: variance from the HL pairs, ln(H/q) law | KNOWN (Friedlander-Goldston 1996; Goldston-Yildirim 1998; Montgomery-Soundararajan 2004; in substance Lemke Oliver-Soundararajan 2016; round 135) | Goldston-Yildirim (arXiv math/0004149) and Montgomery-Soundararajan: search snippets only; the rest as cited in STUDY.md |
| linear response of pair counts to class prime errors, factor 2 / two members, kappa = S phi/A | KNOWN IN SUBSTANCE (Bhowmik-Halupczok-Matsumoto-Suzuki 2019 for Goldbach averaged over a class; round 155 and its review) | not reopened here |
| the four-point variance of twin counts per class contains the linear term exactly; beta_HL = 1 to 5e-4; split into prime rigidity passed on (ln N) and mixed frequencies (ln^2 N) | NOT FOUND | searches: twin primes in progressions variance singular series; Barban-Davenport-Halberstam for prime pairs; Montgomery-Hooley for k-tuples; read: abstracts of Kuperberg-Rodgers-Roditty-Gershon (arXiv 2001.09513), Kowalski (arXiv 0805.4682); snippets of Meng (arXiv 0805.0052, BDH-type bound for k-free twins), Parry (arXiv 1912.03376) |
| ln^2 law: spectral mass 2/(p-2), c_Q H_Q = 1, leading -ln^2(H/q)/(2 ln^2 x) for twin counts | NOT FOUND (it is a Selberg-Delange computation on the HL series; likely implicit in Montgomery-Soundararajan's k-point machinery, which was not read) | search "variance twin primes short intervals log squared": no hit; Goldston, Illinois J. Math. 39 (1995), "Some singular series averages and the distribution of Goldbach numbers": could not open (empty fetch) |
| excess cross-window anti-correlation of twin class counts | NOT FOUND (and not established) | - |
Not opened: Granville 2007, Languasco-Zaccagnini, Fujii 1991, Montgomery-Soundararajan "Beyond pair correlation", a
paper on "pair correlation of primes in progressions" by Leung (not located by one search).

## 7. What was not checked
* The cross-part noise is nominal (independent classes and sub-windows); it could not be calibrated on independent
  decades (there are three). If it is 1.2 times larger, the twin joint z values become -1.8 and -1.9.
* b_Q, c'_Q were not derived analytically (next terms of the Selberg-Delange expansion); they are read off the exact HL
  partial sums.
* The mod-210 predictions for the fresh decade use the closed form (+0.005..0.008 high on the seen decades); they enter
  P4 not at all and level 2 with weight 14/134.
* Other pair gaps h = 4..12 were checked for the cross part only, not for P1-P3.
