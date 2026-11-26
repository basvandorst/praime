# Convex limit of the Weil-positivity problem: duality, truncation, and a referee read of Johns Part II

Folder: /tmp/claude-1000/conv_prover/ (scripts named below; all runs small, float or mpmath).
Scope change (coordinator, mid-task): the attempt to prove T2 from scratch (old P3) was stopped once
prior art (N. Johns, Parts I and II, Zenodo, Oct 2026, unrefereed, AI-produced) was found.
Delivered: P1 (duality), P2' (truncation lemma), P3' (referee read), plus the old P2 numbers.

## 0. Summary

| item | statement | label |
|---|---|---|
| P1(a) | lambda*(L) = D(L) := min{W(g)/g(0) : g p.d., supp in [-L,L], g >= 0 on [0,L]}; the min is attained; the same for the Euler class with the constraint sum_k y^{-k/2} g(k log y) >= 0 | PROVED (uses the Boas-Kac-Krein factorisation, Prokhorov, Levy, Ky Fan's minimax theorem) |
| P1(b) | lambda*(L) is non-increasing in L, both classes | PROVED |
| P1(c) | lim lambda*(L) = sup{log(1/q) : some system feasible at q}; the sup is attained when finite; Euler class too | PROVED (Helly; for the Euler class also: a limit mass m dx/x can be dropped) |
| P1(d) | lambda_inf := lim lambda*(L) = 2 pi kappa*_OPS (Johns' classical cone C_OPS) | PROVED (uses P2') |
| P2' | Truncation lemma: for g >= 0, g^ >= 0, W(g) finite, the window g_L = g k(u/L) (k = normalised cubic B-spline) is doubly positive, has support [-2L,2L], g_L(0) = g(0), and W(g_L) <= W(g) + (3/L^2) int_0^oo u^2 e^{-u/2} g(u)/(1-e^{-2u}) du <= W(g) + 3(7 zeta(3) + pi^3/4) g(0)/L^2. No decay beyond W(g) < oo is needed | PROVED |
| C1 | IF Johns II Thm 1 holds, THEN lambda*(M) <= 194.0/M^2 for all M (and <= 0.3245/M^2 with the measured constant of his g); so lambda*(M) -> 0, i.e. T2's "<= 0" holds unconditionally; ">= 0" is RH-dependent (RH => von Mangoldt feasible at q = 1) | PROVED conditional on Thm 1; 0.3245 is NUMERICAL |
| C2 | From Johns I Thm C(a) alone: lambda*(M) <= 2 pi * 1.4291572e-1060 + 194.0/M^2, so limsup lambda*(M) <= 8.98e-1060 | PROVED conditional on Johns I Thm C(a) (computer-assisted, not rechecked); the decay hypothesis holds for his certificate |
| P3' | Johns II, Sections 3 and 8: the deduction "inputs => Thm 1 => kappa*_OPS <= 0" is sound. Theorem 6.1 (a_n > 0) checked line by line and numerically. Independent numerical rebuild of F from the x > 1 formula alone: F >= 0 on [0,150], F vanishes at the first four zeta zeros, A(F)/int F = -7e-11, H_raw(0) matches to 8 digits. Sections 4, 5(b), 7 and the certificates not refereed | REFEREE READ; numerics NUMERICAL |
| old P2 | Riemann's h (g^ = Xi^2): W/g(0) = 0.0751971327 both by the cost identity and by 2 sum Lambda(n) G(n)/G(1) | NUMERICAL |

## 1. Setting and dictionary

g even real on R, g^(t) = int g(u) e^{itu} du, g(u) = (1/2pi) int g^(t) e^{-itu} dt.
Omega(t) = Re psi(1/4 + it/2) - log pi (even, increasing in |t|, Omega(0) = psi(1/4) - log pi = -5.3721, Omega = log(|t|/2pi) + O(t^-2)).
W(g) := 2 int g(u) cosh(u/2) du + A(g), A(g) := (1/2pi) int g^ Omega dt.            (1.1)
For a positive Radon measure psi on (1,oo): <psi, g> := int x^{-1/2} g(log x) dpsi(x),
Q_psi(f) := W(f*f~) - 2 <psi, f*f~>, lambda_psi(L) := inf{Q_psi(f)/||f||^2 : 0 != f in C_c^oo, supp f in an interval of length L},
lambda*(L) := sup_psi lambda_psi(L). Feasible at q: Q_psi(f) + (log q)||f||^2 >= 0 for all f in C_c^oo.
Euler class E: psi_nu = sum_{k>=1} (x -> x^k)_*(log x dnu), nu positive Radon on (1,oo); then <psi_nu, g> = int S_g dnu,
S_g(y) := log y sum_{k>=1} y^{-k/2} g(k log y).

Dictionary with Johns (Parts I, II). His F lives on the zero side: F = g^, F^(xi) = int F(t) e^{-2 pi i xi t} dt = 2 pi g(2 pi xi),
int F = 2 pi g(0), A(F) = W(g) (his Lemma I.2.2(b) is our cosh term), kappa = lambda/(2 pi).
His cones: C = {F in T : F >= 0, F^ >= 0 on [log2/2pi, oo)} (prime gap log 2) and C_OPS = {F in T : F >= 0, F^ >= 0 on R}.
Our systems may put prime mass anywhere in (1, oo), so our dual cone is g >= 0 on all of (0, oo): it is C_OPS, not C.
Hence T2 is "kappa*_OPS <= 0" (and >= 0), which is stronger than his Conjecture S (kappa* <= 0 over C).
T (test class): F even, real on R, analytic near |Im z| <= 1/2 + delta, sup (1+|z|)^2 |F(z)| < oo.

## 2. P1: duality, monotonicity, the limit

### Lemma 2.1 (integral form of A) PROVED
Let h be even, continuous, h in L^1, h^ in L^1, int |h^(t)| log(2+|t|) dt < oo. Then
A(h) = -(gamma + log pi) h(0) + int_0^oo [2 e^{-2u} h(0) - 2 e^{-u/2} h(u)]/(1 - e^{-2u}) du, the integral converging absolutely.
Proof. Gauss: psi(z) = -gamma + int_0^oo (e^{-x} - e^{-zx})/(1 - e^{-x}) dx (Re z > 0). With z = 1/4 + it/2, x = 2u and real parts:
Omega(t) = -gamma - log pi + int_0^oo 2[e^{-2u} - e^{-u/2} cos(tu)]/(1 - e^{-2u}) du. The integrand is bounded by
C e^{-u/2} for u >= 1 and by C + C min(1/u, t^2 u) for u <= 1 (use 1 - cos(tu) <= min(2, t^2u^2/2)); its u-integral is
<= C(1 + log(2+|t|)). So (t,u) -> h^(t) * integrand is absolutely integrable; Fubini and Fourier inversion
(h(u) = (1/2pi) int h^(t) cos(tu) dt) give the formula. []
For continuous p.d. h in L^1 one has h^ >= 0, h^ in L^1, int h^ = 2 pi h(0); then A(h) is finite iff int h^ log(2+|t|) < oo,
because Omega >= Omega(0) and Omega(t) >= log(|t|/2pi) - C.

### Lemma 2.2 (from test functions to p.d. functions) PROVED (cites Boas-Kac-Krein)
P_L := {g even continuous, p.d. (g^ >= 0), supp g in [-L,L], int g^ log(2+|t|) < oo}. Let psi((1, e^L]) < oo. Then
lambda_psi(L) = inf{ (W(g) - 2<psi,g>)/g(0) : 0 != g in P_L }.
Proof. ">=": f in C_c^oo with support length L gives g = f*f~ in P_L, g(0) = ||f||^2, Q_psi(f) = W(g) - 2<psi,g> (Lemma 2.1).
"<=": fix g in P_L, s in (1/2,1), 0 < eps <= (1-s)L, g_s(u) := g(u/s) (p.d.: g_s^(t) = s g^(st)). By the Boas-Kac-Krein theorem
(Krein 1940; Boas-Kac, Duke Math. J. 12 (1945); see Ehm-Gneiting-Richards, Trans. AMS 356 (2004)) a continuous p.d. function
supported in [-sL, sL] equals f0*f0~ with f0 in L^2(-sL/2, sL/2). With kappa_eps in C_c^oo(-eps/2, eps/2) even, >= 0, int = 1,
rho_eps = kappa_eps*kappa_eps: g_{s,eps} := g_s*rho_eps = f*f~, f := f0*kappa_eps in C_c^oo with support length sL+eps <= L.
So lambda_psi(L) g_{s,eps}(0) <= W(g_{s,eps}) - 2<psi, g_{s,eps}>. Let s -> 1, eps -> 0: g_{s,eps} -> g uniformly with supports in
[-L,L], so the cosh term, <psi,.> (psi finite on (1,e^L]) and g(0) converge; A(g_{s,eps}) = (1/2pi) int g^(tau) rho_eps^(tau/s)
Omega(tau/s) dtau -> A(g) by dominated convergence (0 <= rho_eps^ <= 1, |Omega(tau/s)| <= |Omega(tau)| + C for s >= 1/2). []
If psi((1, 1+delta]) = oo, then <psi, f*f~> = +oo for every f != 0 (f*f~ > 0 near 0), so lambda_psi = -oo; such psi never matter.
The same holds in E: nu((1,e^L]) = oo forces psi_nu((1,1+delta]) = oo (count the k with y^k <= 1+delta).

### Theorem 2.3 (P1) PROVED
D(L) := inf{W(g)/g(0) : g in P_L, g >= 0 on [0,L]},  D_E(L) := inf{W(g)/g(0) : g in P_L, S_g >= 0 on (1, e^L]}.
(a) lambda*(L) = D(L) and lambda*_E(L) = D_E(L); both infima are attained. D(L) >= Omega(0) = -5.3721 (cosh term >= 0).
(b) lambda*(L), lambda*_E(L) are non-increasing in L.
(c) lim_{L->oo} lambda*(L) = sup{log(1/q) : some system is feasible at q} (sup of the empty set = -oo); if finite, the sup is
    attained. The same holds in E (the sup over E, attained in E).
(d) lambda_inf := lim lambda*(L) = 2 pi kappa*_OPS.

Proof of (a). Weak duality: for g in P_L with g >= 0 on [0,L] and any psi, Lemma 2.2 gives
lambda_psi(L) g(0) <= W(g) - 2<psi,g> <= W(g). So lambda*(L) <= D(L). (In E: <psi_nu, g> = int S_g dnu >= 0.)
Strong duality. X := {g in P_L : g(0) = 1} with the topology of weak convergence of the probability measures g^ dt/2pi;
X_M := {g in X : W(g) <= M}. Since |g| <= 1 on [-L,L], |2 int g cosh| <= c_L := 4L cosh(L/2), so on X_M
int (Omega - Omega(0)) g^ dt/2pi <= M + c_L - Omega(0): the family is tight. Prokhorov gives relative compactness; a limit
probability measure mu has characteristic function g = lim g_n uniformly on compacts (Levy), so g is continuous, p.d., g(0) = 1,
supp g in [-L,L], g in L^1, hence mu = g^ dt/2pi; W is lower semicontinuous (cosh term continuous, Omega continuous and bounded
below). So X_M is compact (and metrizable); it is convex and W is affine on it.
Y_B := {psi >= 0 on (1,e^L] : psi((1,e^L]) <= B} (in E: nu with nu((1,e^L]) <= B). f(g,psi) := W(g) - 2<psi,g> is affine in each
variable and continuous in g for fixed psi (uniform convergence on [0,L]). Ky Fan's minimax theorem (Fan 1953: X compact
Hausdorff, f(.,y) lsc and convex-like, f(x,.) concave-like) gives
   min_{X_M} sup_{Y_B} f = sup_{Y_B} min_{X_M} f.
Left side: sup_{Y_B} f(g,.) = W(g) + 2B n(g), n(g) := max(0, max_{[0,L]} (-e^{-u/2} g(u))) (in E: n_E(g) := max(0, max_{[1,e^L]} (-S_g)),
where S_g extends continuously to y = 1 with S_g(1) = int_0^L e^{-u/2} g, and |S_g - S_h| <= 2||g-h||_oo).
Right side: fix g1 in X with W(g1) = w1 and take M >= w1 + 4B. For psi in Y_B and g in X \ X_M, f(g,psi) > M - 2B >= f(g1,psi),
so min_{X_M} f(.,psi) = inf_X f(.,psi) = lambda_psi(L) (Lemma 2.2); likewise min_{X_M}(W + 2Bn) = inf_X (W + 2Bn) =: m_B.
Hence lambda*(L) >= sup_{Y_B} lambda_psi(L) = m_B for every B. Let g_B attain m_B. Then W(g_B) <= m_B <= D(L), so all g_B lie in
the fixed compact X_{D(L)}, and n(g_B) <= (D(L) - Omega(0) + c_L)/(2B) -> 0. A limit point g_oo has n(g_oo) = 0 (n is continuous
under uniform convergence on [0,L]), so it is feasible for D(L), and D(L) <= W(g_oo) <= liminf m_B <= lambda*(L) <= D(L).
Equality, and g_oo attains D(L). The Euler case is word for word the same with n_E. []

Proof of (b). A function in P_L with g >= 0 on [0,L] (or S_g >= 0 on (1,e^L]) is admissible for every L' > L (g = 0 on (L,L']). []

Proof of (c). If psi is feasible at q, then lambda_psi(L) >= log(1/q) for every L, so lambda_inf >= log(1/q). Conversely let
lambda_inf > -oo and choose psi_j with lambda_{psi_j}(j) >= lambda_inf - 1/j. Local bound: for X > 1 take f_X in C_c^oo, f_X >= 0,
with g_X = f_X*f_X~ >= c_X > 0 on [0, log X] (support length log X + 2); for j >= log X + 2,
2 c_X X^{-1/2} psi_j((1,X]) <= 2<psi_j, g_X> <= W(g_X) - (lambda_inf - 1) g_X(0). So psi_j is bounded on every [1,X], and a
subsequence converges vaguely on the locally compact space [1,oo) to psi_bar = m delta_1 + psi. Every x^{-1/2} g(log x),
g = f*f~, is continuous with compact support on [1,oo), so W(g) - 2m g(0) - 2<psi,g> >= lambda_inf g(0): psi is feasible at
q = e^{-lambda_inf - 2m}. By the first part log(1/q) = lambda_inf + 2m <= lambda_inf, so m = 0 and psi attains the sup.
Euler class: psi_nu((1,b]) >= (log b/2) nu((1,b^{1/2}]) gives local bounds for nu_j on [1,oo); the vague limit is
nu_bar = m' delta_1 + nu, and int S_g dnu_j -> m' S_g(1) + int S_g dnu, with m' S_g(1) = <m' dx/x, g>. So psi_nu + m' dx/x is
feasible at e^{-lambda_inf}. Since <dx/x, f*f~> = (1/2)(1/2pi) int |f^|^2/(1/4 + t^2) dt >= 0, dropping m' dx/x keeps
feasibility: psi_nu in E is feasible at e^{-lambda_inf}. []

Proof of (d). "<=": every F in C_OPS gives g = F^(./2pi)/2pi with g >= 0, g^ = F >= 0, int g e^{|u|/2} < oo
(|F^(xi)| <= C e^{-pi(1+2 delta)|xi|}) and int g^ log < oo; Lemma 3.1 below gives lambda*(2L) <= 2 pi A(F)/int F + 3 c1/L^2.
">=": for g in P_L doubly positive, g_eps := g*rho_eps is in C_c^oo, doubly positive, so g_eps^ is in C_OPS (entire, rapidly
decreasing in strips); W(g_eps) -> W(g) (cosh(u/2) is an eigenfunction of even convolutions: the cosh term is multiplied by
int rho_eps cosh(v/2) -> 1; the A term by dominated convergence) and g_eps(0) -> g(0). So 2 pi kappa*_OPS <= D(L) for all L. []

What is used, exactly: Bochner; Boas-Kac-Krein; Gauss's integral for psi; Prokhorov and Levy's continuity theorem; Ky Fan's
minimax theorem; Helly/Banach-Alaoglu for Radon measures on [1,X]. No RH. Under RH the von Mangoldt measure is feasible at
q = 1 (Q_Lambda(f) = sum_rho |f^(gamma)|^2), so lambda*(L) >= lambda_inf >= 0. Without RH no lower bound > -oo is known to us
(Johns I Prop 5.10 claims kappa* >= -2.71e-3 at gap log 2 by a certified pair; not checked, and that is for the other cone).

## 3. P2': the truncation lemma

k(x) := 1 - (3/2)x^2 + (3/4)|x|^3 (|x| <= 1), (2-|x|)^3/4 (1 <= |x| <= 2), 0 (|x| >= 2): the centred cubic B-spline divided by its
value 2/3 at 0. Its transform is (sin(t/2)/(t/2))^4 / (2/3) >= 0, so k is p.d.; 0 <= k <= 1; and 0 <= 1 - k(x) <= (3/2) x^2.
K(u) := 4 cosh(u/2) - 2 e^{-u/2}/(1 - e^{-2u}) (negative for u < 0.28, positive after).
c1 := int_0^oo u^2 e^{-u/2}/(1 - e^{-2u}) du = sum_{k>=0} 16/(4k+1)^3 = 7 zeta(3) + pi^3/4 = 16.16597.

### Lemma 3.1 (truncation) PROVED
Let g be even, continuous, g >= 0, g p.d. (g^ >= 0), with (H1) int g(u) e^{|u|/2} du < oo and (H2) int g^ log(2+|t|) dt < oo
(equivalently: W(g) finite). For L > 0 put g_L(u) := g(u) k(u/L). Then g_L is even, continuous, g_L >= 0, g_L^ >= 0,
supp g_L in [-2L, 2L], g_L(0) = g(0), W(g_L) is finite, and
  W(g_L) = W(g) - int_0^oo K(u) g(u) (1 - k(u/L)) du                                                  (3.1)
         = W(g) - 4 int_0^oo cosh(u/2) g (1-k(u/L)) du + 2 int_0^oo e^{-u/2} g (1-k(u/L))/(1-e^{-2u}) du
        <= W(g) + (3/L^2) I_g,   I_g := int_0^oo u^2 e^{-u/2} g(u)/(1 - e^{-2u}) du <= c1 g(0).            (3.2)
Consequently, with M := 2L (support [-M, M]):
  lambda*(M) <= W(g)/g(0) + 12 I_g/(M^2 g(0)) <= W(g)/g(0) + 12 c1/M^2 = W(g)/g(0) + (84 zeta(3) + 3 pi^3)/M^2,
  12 c1 = 193.99.
Proof. g in L^1 by (H1), g^ in L^1 (p.d. and continuous); k_L in L^1 with k_L^ >= 0 in L^1; so (g k_L)^ = (g^ * k_L^)/2pi >= 0.
int (g k_L)^ log(2+|t|) < oo because log(2+|a+b|) <= log(2+|a|) + log(2+|b|) and int k_L^ log(2+|s|) ds < oo (k_L^ = O(s^-4)).
So Lemma 2.1 applies to g and to g_L. Since g_L(0) = g(0), the terms in h(0) cancel and
A(g) - A(g_L) = -2 int_0^oo e^{-u/2} g(u)(1 - k(u/L))/(1 - e^{-2u}) du, while the cosh terms differ by
4 int_0^oo cosh(u/2) g (1-k(u/L)) du >= 0 (finite by (H1)). This is (3.1). Drop the cosh difference (it is >= 0 because g >= 0),
use 1 - k(x) <= (3/2)x^2 and g <= g(0) (p.d.): (3.2). The bound for lambda* is weak duality (Lemma 2.2 / Thm 2.3(a)):
g_L is in P_{2L} and >= 0. c1 by Tonelli: int u^2 e^{-(1/2+2k)u} du = 2/(1/2+2k)^3; the sum over k = 1 mod 4 of n^-3 is
((7/8) zeta(3) + beta(3))/2 with beta(3) = pi^3/32. []
Remarks. (i) No decay of g beyond (H1) is used: the large-u part of the window error enters only through the cosh term, which
has the good sign. The rate 1/L^2 is the best possible for windows of the form k(u/L) with k p.d.: 1 - k(x) >= c x^2 near 0, and for
Johns' g the second moment int K u^2 g du = -2.29e-3 g(0) < 0 (it must be <= 0 under RH, since W(g_L) >= 0), so
W(g_L) ~ 3.4e-3 g(0)/L^2 (NUMERICAL). The triangle (1-|u|/L)_+ also works, with rate 1/L and constant 2 c0,
c0 = int_0^oo u e^{-u/2}/(1-e^{-2u}) du = pi^2/4 + 2G = 4.2993 (G Catalan's constant).
(ii) The hypotheses hold for every F in Johns' T (Lemma I.2.2(a): |F^(xi)| <= pi ||F|| e^{-pi(1+2delta)|xi|}), in particular for
his Part I certificates Xi^2 P_J(t^2) (entire, exponentially decreasing in strips by |Xi(t)| <= C(1+|Re t|)^3 e^{-pi|Re t|/4},
which I checked) and for his Part II function F = Xi^2 H (in T_delta by his Lemma 3.2(a), checked).

### Corollary 3.2 (conditional on Johns II, Theorem 1) PROVED given Theorem 1
If Theorem 1 of Part II holds (F = Xi^2 H in C_OPS, int F > 0, A(F) = 0), then with g = F^(./2pi)/2pi, W(g) = 0 and
  0 <= lambda*(M) <= 193.99/M^2   for every M > 0,
where the lower bound needs RH. With the measured constant of his g (I_g/g(0) = 0.02704, NUMERICAL from the rebuilt F of
Section 5.3): lambda*(M) <= 0.3245/M^2; the exact window values W(g_L)/g(0) are 1.49e-1, 2.97e-2, 5.51e-3, 1.12e-3, 2.47e-4, 5.77e-5
for supports [-M,M], M = 0.5, 1, 2, 4, 8, 16 (NUMERICAL, trunc_partII.py). So T2's "<= 0" holds unconditionally, given Theorem 1;
T2's ">= 0" is RH-dependent (and by Johns II Cor. 8.2, which rests on more inputs, equivalent to RH).
The deduction uses only Theorem 1(a),(c), Lemma 3.1 and weak duality; it does not use Part I's duality theorem, Theorem I.3.6,
the BRS basis, or Theorem 2.

### Corollary 3.3 (from Johns I, Theorem C(a) alone) PROVED given Theorem C(a)
Part I certifies (computer-assisted, Thm 5.1 / Prop 5.7; not rechecked here) an F = Xi^2 P_111(t^2) in C_OPS with
A(F)/int F <= 1.4291572e-1060. Lemma 3.1 applies (Remark (ii)), so for every M > 0
  lambda*(M) <= 2 pi * 1.4291572e-1060 + 193.99/M^2,   limsup lambda*(M) = lambda_inf <= 8.98e-1060.

### Old P2 numbers (kept; NUMERICAL, p2_riemann.py)
Riemann's h(y) = 2(2 pi^2 y^4 - 3 pi y^2) e^{-pi y^2}: h~(s) = pi^{-s/2} Gamma(s/2) (s/2)(s-1), so zeta h~ = xi and g0^ = Xi^2.
W(g0)/g0(0) = 0.0751971327 from the cost identity (pole part 1/2 = 2 xi(0)^2, g0(0) = 0.3197518) and 0.0751971327 from
2 sum Lambda(n) n^{-1/2} g0(log n)/g0(0) (n = 2 gives 0.07406, n = 3 0.00113): the unconditional explicit formula, checked to 1e-16.
The window lowers it: W(g0 k(./T))/g0(0) = 0.0407, 0.0640, 0.0698 for supports [-2,2], [-4,4], [-6,6]; I_{g0}/g0(0) = 0.0554.

## 4. P3': referee read of Johns Part II (Sections 3 and 8, inputs of Remarks 8.1 and 8.3, one sampled lemma)

### 4.1 Section 3 (reduction): sound
- Lemma 3.2. (a) |Xi|^2 <= C(1+|Re z|)^6 e^{-pi|Re z|/2} on the strip (Lemma 2.1, which I checked from Stirling and convexity
  bounds) times (C1) gives (1+|z|)^-2 decay: correct. (b) Stirling exponent 3 + Re s <= 4 + delta: correct. (c) The shift to
  Im w = -eta, eta in (1/2, 1/2+delta], expansion zeta(1/2+iw)^2 = sum d(m) m^{-1/2-iw}, Fubini, and the shift back are correct;
  gamma_inf^2 has a double zero at t = -i/2 that cancels the double pole of zeta^2, and gamma_inf has no poles for Im t < 5/2.
  So F^(xi) = sum_m d(m) m^{-1/2} G_H^(xi + xi_m) for every real xi: correct.
- Prop. 3.3: correct. For C_OPS it needs exactly G_H^ >= 0 for x = e^{2 pi xi} >= 1 (all arguments m x >= 1 when xi >= 0, and
  F^ is even). The vanishing at the integers propagates because d(m) >= 0 and integers are closed under multiplication; A(F) = 0
  then follows from the unconditional explicit formula (Cor. I.4.4), since Xi^2 kills every zero, on or off the line.
  This is the key idea; it replaces |zeta|^2 (pair correlations m/n, which stalled my own theta-map route) by zeta^2.

### 4.2 Section 8 (assembly): sound, given the inputs
Proof of Theorem 1 checked step by step: C = 1/H_raw(0) > 0 (Prop. 6.5: F_raw^ >= 0 as a positive combination, F_raw^(0) =
1/(32 pi^2), Fourier inversion at 0); (C2) from Cor. 7.12; (C3) from Prop. 5.11(a); (C4) on all of [0,oo) from (5.8), Thm 6.1 and
K_0 > 0, plus the value at x = 1; Prop. 3.3 then gives F in C_OPS and A(F) = 0. The C^1 claim (|F| <= C(1+|t|)^{-3+delta}, so
t F in L^1) and F^'(xi_n) = 0 are correct. Nothing wrong found. Minor: the letter C denotes both the cone and the constant
880.66. Theorem 2, Corollaries 3 and 8.2 need further inputs (Theorem I.3.6, the BRS basis decay (P5) "not re-verified" by the
author, Part I's duality) that T2's "<= 0" does not need.

### 4.3 Inputs, and what I could not check
Theorem 1 rests on: Prop. 5.9(a),(b),(c) [entire, even, growth (1.2), Gamma-only transform = G_raw^]; Thm 5.4(c),(d) [x > 1 series,
value at 1]; Prop. 5.11(a) [integer double zeros]; Thm 6.1 [a_n > 0]; Prop. 6.5 [H_raw(0) > 0]; Cor. 7.12 [H_raw > 0 on R];
the explicit formula on T (classical, BRS (1.1)); Lemma 2.1.
Computer-assisted (author's list, not rechecked by me in ball arithmetic):
 (1) Thm 6.8, moment certificate: H > 0 for |t| <= 10.355 (sixteen moments);
 (2) Prop. 7.10, inputs of the large-|t| theorem: sign of B on [0.90, 4], tail constant rho_bar(4) < 1, three integrals at t = 9;
 (3) coefficient enclosures (Kloosterman terms c <= 6e5, Weil-bound tails): used only for the value of C;
 (4) Prop. 6.4, the value of C (cosmetic).
Not checked (analytic, long): Section 4 (Thm 4.4 Fourier expansion of the odd Niebur-Poincare series; Lemma 4.6; Prop. 4.7;
Lemma 4.8; Prop. 4.9, in particular the second-difference identity Phi(tau+1) - 2Phi(tau) + Phi(tau-1) = 2 Phi~, vanishing on
iR+, absence of constant terms, growth); Thm 5.4(b) (Green/Stokes argument; structure read and plausible, estimates not
checked); (4.10) D_grow; Lemmas 5.6, 5.7; Section 7 (Thm 7.2 kernel positivity, Thm 7.3 physical-side formula, Lemmas 7.5-7.8,
Thm 7.11, Cor. 7.12); Appendix A; the Lean check (which covers only the deduction, not the inputs).
Checked: the evenness mechanism (Lemma 5.1(iv): the ODE depends on w^2 only; Lemma 5.3: both brackets of (5.5) are
S-symmetric; hence R(-w) = R(w) for the function G defined by (5.6)); Lemma 5.5 (numerically, 18 digits).
Steps that look under-justified: none found that affect Theorem 1. The weakest links are, in my judgement, Section 4 (a new
non-holomorphic construction) and the computer-assisted positivity of H, which I could only test numerically (below).

### 4.4 Sampled lemma: Theorem 6.1 (positivity of every a_n), checked properly
Proof read line by line: the split c = 1 / c = 2 / 3 <= c <= z/2 / c > z/2, the four bounds and the integral comparison
sum_{c > z/2} f(c) <= f(z/2) + int_{z/2}^oo f = kappa3 + (z/2)(kappa1 + kappa3) are correct; Lemmas B.1-B.5 re-derived
(B.1 termwise differentiation; B.2 term bounds; B.3 Schlaefli integrals, int_0^pi theta|sin 2theta| = pi; B.4 integral
representations, sqrt(1-x) >= 1 - x/2 - x^2/2; B.5(a) is the standard formula for d I_nu/d nu at nu = 2) and spot-checked on
grids (no violation). lambda_-(4 pi) = 0.85010 reproduced from the five summands (0.000250, 0.002641, 0.003145, 0.015725,
0.005683) and lambda_0(4pi) = 0.877952; lambda_- increasing (proof correct; sampled). Independent evaluation of a_n from the
Kloosterman series (c <= 1500): a_1..a_5 agree with Table 1 to relative 1e-9 ... 4e-15 (check_an.py). Verdict: correct.

### 4.5 Independent numerical rebuild of F (NUMERICAL, not a proof; check_magic.py, check_magic2.py)
Built ONLY from Theorem 4(b) for x > 1 (a_n from the Kloosterman series for n <= 12, c = 1 term beyond; agreement 1e-8 at
n = 8..12), the positive combination (3.2) for xi >= 0, and evenness of F^:
- G_raw^(x) -> 1/(32 pi^2) = 0.00316629 as x -> 1+ (0.0031661602 at x = 1.002), as Thm 5.4(d) claims;
- H_raw(0) = F(0)/Xi(0)^2 = 0.0011355098 (Prop. 6.4: 0.0011355097904);
- F(gamma_k) = -7.8e-14, -1.6e-14, -8.0e-15, -3.6e-15 at the first four zeta zeros, against F = 1e-7 ... 1e-13 nearby: F has the
  factor Xi^2. This is not built into the rebuild (the truncated Mellin integrals do not factor through zeta^2 unless the x > 1
  formula is the restriction of a self-reciprocal function), so it tests Theorem 5.4(e) together with (c);
- F >= -7e-15 on [0,150] (noise level), so H = F/Xi^2 > 0 there; H(40)/H(30) = 3.1e4, between the claimed lower rate e^{0.6435 t}
  and the allowed e^{pi t/2};
- A(F)/int F = -7.3e-11 by the cost identity (pole part 5.6445181587e-4, archimedean part -5.6445181610e-4).
Asymptotic check at x = 1 (analytic, from Bessel expansions; z^2(-I2dot(z))K0(z) = 1 - 3/(2z) + 0/z^2 - (5/16)/z^3 + ...): the
x > 1 series has leading term 1/(32 pi^4 (x-1)^2), zero residue (needed for F^'(0) = 0, i.e. evenness), and the log(x-1) and
(x-1)log(x-1) terms cancel (2 p2'(1) - p3(1) = -5/16 + 5/16 = 0), as analyticity of G at x = 1 (Thm 5.4(d)) requires.

### 4.6 Verdict on "inputs => Theorem 1 => kappa*_OPS <= 0"
Sound. The deduction in Sections 3 and 8 is correct; the inputs from Sections 4-7 are where the risk is, and four of them are
computer-assisted. With Lemma 3.1, Theorem 1 gives lambda*(M) <= 194/M^2 and hence T2's "<= 0"; nothing in this chain needs
Theorem 2, Part I's duality or the BRS basis. Numerically the object behaves exactly as claimed on every test I ran.

## 5. What is proved, what is not

Proved here: P1 (Thm 2.3 (a)-(d)); the truncation Lemma 3.1 with the explicit universal error 12(7 zeta(3) + pi^3/4)/M^2;
Corollary 3.2 (Theorem 1 => lambda*(M) <= 194/M^2 => T2's "<= 0"); Corollary 3.3 (Part I Theorem C(a) => limsup lambda* <= 8.98e-1060).
Not proved here: Theorem 1 of Johns II itself (Sections 4, 5(b), 7 and the computer-assisted certificates are not refereed);
T2's ">= 0" (needs RH; I have no unconditional lower bound). Old P3 (a proof of T2 from scratch) was stopped by the coordinator;
the theta-map route reached only Riemann's 0.0752 and its window 0.0407.

Files: convex_limit.md (this), xi_grid.py, p2_riemann.py, common.py, check_an.py, check_expansion.py, check_magic.py,
check_magic2.py, check_l55.py, trunc_partII.py; data xi_grid.npz, omega.npy, fhat.npz. lp_cos.py (theta-family LP) was stopped
at the scope change; its partial outputs (lpcos_c_*.npy) are not used.
