# Zeros of the degree-d optima: KKT with points at infinity, a trap lemma, and what convergence needs

Written by a proof agent on 26 November 2026; not refereed by a second reader. Scripts (this folder): trap.py, trap2.py (traps, weights), moments.py (moments of the
optima); outputs *.out. Inputs: ../convex/certs/polish_A_d*.json, ../convex/certs/core.py. Every number below is NUMERICAL
(mpmath, 60 to 100 digits, the stored optima taken as exact); no number is certified here.

## 0. Summary

| # | statement | label |
|---|---|---|
| 1.1 | g_j(u) = (-4 pi^2 e^{2u})^j g_0(u) (1 + O(e^{-u})) as u -> +oo | PROVED from (CF) |
| 1.3 | the minimum is attained, 0 < lambda_d, parity rule, odd d = previous even d | PROVED |
| 2.2 | T1, d even: KKT with one multiplier for both points at infinity; same multipliers for every minimiser | PROVED |
| 2.3 | T1, d odd: the multiplier of b_d is unsigned (u = oo gives -b_d, t = oo gives +b_d) | PROVED |
| 2.4 | the minimiser is unique if the dual has d/2 atoms and the Hermite matrix is non-singular | PROVED UNDER HYPOTHESIS (hypotheses NUMERICAL, d <= 24) |
| 3.1 | T2, trap lemma, in the form "w in V_d', B_w <= 0, W(w) > 0"; brief's (iii) redundant, eta h(0) <= 1 automatic | PROVED |
| 3.2 | upper bound: int v dmu_d <= W(v) - lambda_d v(0) for v in the cone; int gt dmu_d <= kt - lambda_d | PROVED |
| 3.3 | KKT weight inside the k-th trap at every higher degree = mu_k(d') up to explicit errors | PROVED UNDER HYPOTHESIS (exact KKT data at d') |
| 4.1 | T3(a): localised trap under explicit (H1)-(H4) | PROVED given (H1)-(H4); (H1)-(H3) NUMERICAL on grids for d' <= 24, (H4) only at the centre |
| 4.2 | T3(c): g(log n) <= sqrt(n) W(g)/(2 Lambda(n)); minimiser within (sqrt(n) W/(Lambda m))^{1/2} | PROVED |
| 4.3 | T3(b): zero of G_d in J converges to log n at rate lambda_d^{1/2} under (A), (B), (C) | PROVED UNDER (A),(B),(C) |
| 4.4 | (A) follows if infinitely many Hermite members F_J lie in the cone with sup a_J < pi/2 | PROVED (implication) |
| 4.5 | gt = F_J with only Props 5.6/5.7 of Johns I: mass bounds yes, one zero per prime power no | PROVED (what it gives); localisation NOT PROVED |
| 5.1 | T4: normal family, tightness, limits; lambda_d -> 0 forces the limit to vanish at every log n | PROVED (uses the explicit formula) |
| 5.2 | T4: convergence of the zeros to the prime powers | NOT PROVED; the four missing statements (M1)-(M4) are listed |

## 1. Setting

Notation as in the brief, except: the minimiser at degree d is written G_d (the brief's g_d is also the basis
element (-d^2/du^2)^d g_0). b_j(h) is the coefficient of t^{2j} in B_h; P(h) := 2 sum_n Lambda(n) n^{-1/2} h(log n)
= W(h); m := (g_j(0))_j, c := (W(g_j))_j, so h(0) = m.b and W(h) = c.b for h = sum b_j g_j.
Assumed: (CF) the closed form of g_0 = Phi*Phi (brief; Johns I Lemma 4.6; derived in certs/README.md); Phi > 0
on R (classical). (EF) the explicit formula W(h) = 2 int h(u) cosh(u/2) du + A(h) for h in V_d (classical, the zero
side vanishes because hhat = Xi^2 B_h; Johns I Cor. 4.4): used ONLY in section 5.
Correction: Phi is not entire. Its theta series converges for |Im u| < pi/4 and diverges on Im u = +-pi/4; Phi and
g_0 are real-analytic on R, which is all that is used.

**Lemma 1.1 (asymptotics). PROVED from (CF).** For every j, g_j(u) = (-4 pi^2 e^{2u})^j g_0(u) (1 + O_j(e^{-u}))
as u -> +oo. Hence for 0 != h in V_d with top non-zero coefficient b_e = b_e(h):
h(u) = b_e (-4 pi^2 e^{2u})^e g_0(u) (1 + O(e^{-u})); in particular h has the sign of (-1)^e b_e near +oo.
Proof. With D = d/du = X d/dX, the N-th term of g_j is 4 d(N) (pi N)^{-1/2} X^{1/2} [P_{0,j}(X) K0(2X) + P_{1,j}(X) K1(2X)],
X = pi N e^u, (P_{0,0}, P_{1,0}) = (4X^4 + 9X^2, -12X^3), g_j = (-1)^j D^{2j} g_0. From K0' = -K1, K1'(z) = -K0(z) - K1(z)/z:
D[X^{k+1/2} K0(2X)] = (k+1/2) X^{k+1/2} K0 - 2 X^{k+3/2} K1, D[X^{k+1/2} K1(2X)] = (k-1/2) X^{k+1/2} K1 - 2 X^{k+3/2} K0.
So each D raises the top degree by one, multiplies the top coefficient by -2 and swaps K0 and K1: P_{0,j} has top term
(-1)^j 4^{j+1} X^{2j+4}, every other term of P_{0,j}, P_{1,j} has degree <= 2j+3. Since K1(2X)/K0(2X) = 1 + O(1/X) and
K0 > 0, the N-th term equals (-4X^2)^j * [N-th term of g_0 at leading order] * (1 + O(1/X)). Terms N >= 2 are smaller
than the N = 1 term by O(N^c e^{-2 pi (N-1) e^u}), summable. Divide by the same expansion for j = 0; X = pi e^u. []
Check (NUMERICAL, asym.out): g_j/(g_0 (-4 pi^2 e^{2u})^j) for j = 1..4 is 0.929, 0.834, 0.722, 0.602 at u = 3 and
0.990, 0.977, 0.959, 0.937 at u = 5 (error of order j e^{-u}, as claimed).

**Lemma 1.2. PROVED.** Let 0 != h in V_d. (a) h has finitely many real zeros. (b) sum_n Lambda(n) n^{-1/2} |h(log n)| < oo.
(c) h(0) = (1/2pi) int Xi^2 B_h dt; if h in C_d then h(0) > 0.
Proof. (a) Every term of (CF) is analytic in e^u on Re e^u > 0 and the series converges locally uniformly there, so h is
real-analytic; by 1.1 and evenness h != 0 for |u| >= U_h; a non-zero analytic function has finitely many zeros on
[-U_h, U_h]. (b) |h(u)| <= C e^{2du} g_0(u) and g_0(u) = O(e^{5u} e^{-2 pi e^u}). (c) Fourier inversion (hhat in L^1);
B_h >= 0 and B_h != 0 (the g_j are independent: Xi^2 t^{2j} are). []

**Lemma 1.3. PROVED.** (a) K_d := {h in C_d : h(0) = 1} is compact and non-empty, so lambda_d is attained.
(b) lambda_d > 0. (c) lambda_{d+1} <= lambda_d. (d) h in C_d implies b_d(h) >= 0 and (-1)^d b_d(h) >= 0; for odd d,
C_d = C_{d-1}, lambda_d = lambda_{d-1}, with the same minimisers.
Proof. (a) K_d is closed; if unbounded, normalised elements converge to b* in C_d with |b*| = 1 and m.b* = 0, against
1.2(c); g_0/g_0(0) is in K_d. (b) W(G) = 0 with G >= 0 forces G(log n) = 0 at every prime power: infinitely many zeros,
against 1.2(a). (d) B_h >= 0 forces b_d >= 0; h >= 0 and 1.1 force (-1)^d b_d >= 0. []

## 2. T1: KKT with points at infinity

Coordinates b in R^{d+1}. Functionals e_u(b) = h_b(u) (u >= 0), f_t(b) = B_b(t) (t >= 0), tau(b) = b_d. Weights
N(u) = g_0(u)(1 + (2pi)^{2d} e^{2du}), M(t) = 1 + t^{2d}. Gamma := {e_u/N(u)} u.{f_t/M(t)} u.{(-1)^d tau, tau}.
(Evenness of h and of B_h makes u, t >= 0 enough.)

**Lemma 2.1. PROVED.** Gamma is compact. For d even there are b* and delta > 0 with y(b*) >= delta for all y in Gamma.
Hence cone(Gamma) is closed and the dual cone C_d* = {y : y(b) >= 0 for all b in C_d} equals cone(Gamma).
Proof. By 1.1, e_u/N(u) -> (-1)^d tau as u -> oo (coordinate j -> 0 for j < d, -> (-1)^d for j = d), and f_t/M(t) -> tau;
so Gamma is a continuous image of two copies of [0, oo]. Slater point (d even): b* = g_0 + eps g_d (basis element).
B_{b*} = 1 + eps t^{2d}, so f_t(b*)/M(t) >= min(1, eps) and tau(b*) = eps. g_0 > 0 (Phi > 0), g_d/g_0 is continuous on
[0, oo) and -> +oo (1.1, d even), so g_d/g_0 >= -M_d; for eps < 1/(2M_d), h_{b*}/N is continuous, positive on [0, oo) with
limit eps > 0 at oo, hence >= delta' > 0. Then conv(Gamma) is compact and lies in {y(b*) >= delta}, so 0 is not in it and
cone(Gamma) = R_+ conv(Gamma) is closed. C_d = {b : y(b) >= 0 for y in Gamma} (the two limit functionals are limits of
the others and add no constraint); the bipolar theorem gives C_d* = closed conic hull of Gamma = cone(Gamma). []

**Theorem 2.2 (T1, d even). PROVED.** There are points 0 < u_1 < ... < u_p, t_1, ..., t_q >= 0 and weights
mu_i > 0, nu_i > 0, alpha >= 0 (at most d+1 of the atoms and alpha non-zero) such that for every h in V_d
   W(h) = lambda_d h(0) + sum_i mu_i h(u_i) + sum_i nu_i B_h(t_i) + alpha b_d(h).                           (2.1)
Every minimiser G_d satisfies G_d(u_i) = G_d'(u_i) = 0, B_{G_d}(t_i) = 0, alpha b_d(G_d) = 0. Conversely ANY
non-negative (mu, nu, alpha) with (2.1) has these properties for EVERY minimiser: the multipliers are common to all
minimisers. mu has no atom at 0.
Proof. Claim c - lambda_d m in C_d*. If not, separate the point from the closed convex cone C_d* (Lemma 2.1): there is
b in C_d** = C_d with (c - lambda_d m).b < 0. If m.b > 0 then b/m.b in K_d costs less than lambda_d; if m.b = 0 then b = 0
by 1.2(c); both impossible. By Lemma 2.1 and Caratheodory, c - lambda_d m = sum_{i <= d+1} s_i y_i with y_i in Gamma,
s_i >= 0; undo the normalisations (e_u/N(u) gives an atom s/N(u) of mu at u, f_t/M(t) an atom of nu, (-1)^d tau = tau
gives alpha). Slackness: insert h = G_d into (2.1); W(G_d) = lambda_d, G_d(0) = 1, so the three non-negative sums
vanish term by term. An atom at u = 0 would contribute mu_0 G_d(0) > 0. G_d >= 0 with an interior zero has G_d' = 0
there. The converse uses only (2.1), W(G) = lambda_d and G(0) = 1. []
Numerically (polish.py data, d <= 24): nu = 0 (B has no real root), alpha = 0 (b_d > 0), p = d/2 atoms.

**Theorem 2.3 (T1, d odd). PROVED.** C_d = C_{d-1} and the minimisers coincide (1.3(d)). For every dual solution
(mu, nu, alpha') of degree d-1 there is exactly one real beta with W(h) = lambda_d h(0) + int h dmu + int B_h dnu
+ alpha' b_{d-1}(h) + beta b_d(h) for all h in V_d (the difference of the two sides vanishes on V_{d-1}, so it is a
multiple of b_d). beta has no sign: Gamma contains both tau (t -> oo) and -tau (u -> oo), 0 is in conv(Gamma) and
Lemma 2.1 fails; beta > 0 reads as mass of nu at t = oo, beta < 0 as mass of mu at u = oo.

**Proposition 2.4 (uniqueness). PROVED UNDER HYPOTHESIS.** If at even d = 2K some dual solution has K atoms
u_1..u_K and the linear map h -> (h(0), h(u_i), h'(u_i))_{i <= K} on V_d is injective, the minimiser is unique.
Proof: by 2.2 every minimiser solves the same 2K+1 linear equations G(0) = 1, G(u_i) = G'(u_i) = 0. []
polish.py solves exactly this square system with all mu_i > 0 (residual <= 1e-40), so modulo exactness the optima
for d <= 24 are unique. Nothing below needs uniqueness.

## 3. T2: the trap lemma and its upper-bound twin

**Lemma 3.1 (trap lemma). PROVED.** Let w in V_{d'} with (i') B_w <= 0 on R and (ii') W(w) > 0. Put
S := {u > 0 : w(u) > 0}. Then:
(a) w is negative definite: w(0) < 0 and |w(u)| <= -w(0) for all u; so S is bounded away from 0 and sup_S w < oo.
(b) For every even d >= d' and every dual solution (mu_d, nu_d, alpha_d) of degree d:
    int_S w dmu_d >= int w dmu_d >= W(w) - lambda_d w(0) >= W(w) > 0.
(c) mu_d(S) >= (W(w) - lambda_d w(0)) / sup_S w > 0. So S is not empty and every minimiser G_d has a double zero in S.
(d) The same holds for odd d > d' (apply (b) at d-1 >= d'), and for odd d = d' if b_{d'}(w) = 0.
Proof. (a) what = Xi^2 B_w <= 0 and in L^1, so -w is continuous positive definite (Bochner): |w(u)| <= -w(0) =
-(1/2pi) int what; w(0) = 0 would force what = 0, w = 0, against (ii'). (b) Theorem 2.2 for w in V_d:
int w dmu_d = W(w) - lambda_d w(0) - int B_w dnu_d - alpha_d b_d(w). Here -lambda_d w(0) >= 0 (lambda_d > 0, w(0) < 0),
-int B_w dnu_d >= 0 by (i'), and b_d(w) = 0 if d > d', b_d(w) <= 0 if d = d' (a polynomial <= 0 on R has top
coefficient <= 0). On (0, oo) \ S, w <= 0, and mu_d lives on (0, oo). (c) int_S w dmu_d <= mu_d(S) sup_S w, and mu_d
is carried by the zero set of G_d (2.2). []

The brief's form is the case w = eta h - gt, gt in C_{d'}, gt(0) = 1, kt = W(gt): (i) eta B_h <= B_gt is (i'), and
(ii') is eta W(h) > kt. Conversely every w of 3.1 is of that form (gt = g_0/g_0(0), h = (w + gt)/eta). Corrections:
- (iii) is implied by (i): B_gt - eta B_h >= 0 on R has a non-negative top coefficient.
- Under (i), eta h(0) - 1 = w(0) < 0 always, so the factor (1 + max(eta h(0) - 1, 0)) is 1: (ii) is just eta W(h) > kt,
  and the mass bound improves to mu_d(S) >= (eta W(h) - kt + lambda_d (1 - eta h(0)))/sup_S w.
- lambda_d <= kt is true (d >= d') but not needed.
- Where it fails to localise. (1) S may have several components; the lemma does not say which holds the zero.
  (2) S can be unbounded: by 1.1, w has the sign of (-1)^e b_e(w) at +oo (e = top index), and (i') forces b_e(w) < 0;
  so w < 0 near +oo iff e is even. With gt of exact degree d' (b_{d'}(gt) > 0, d' even) and b_{d'}(h) = 0 (the
  Hermite choice) e = d' and S is bounded. If instead eta b_{d'}(h) = b_{d'}(gt) and b_{d'-1}(w) != 0, then e = d'-1
  is odd and S contains [U, oo): the far zeros of G_d (at large prime powers) satisfy the conclusion and nothing is
  localised. (3) S is never near 0 (by (a)). (4) Under (ii') S is never empty.

**Lemma 3.2 (upper bounds). PROVED.** For even d, d'' <= d and v in C_{d''}:
int v dmu_d + int B_v dnu_d + alpha_d b_d(v) = W(v) - lambda_d v(0), with all three terms >= 0. Hence for every Borel
E in (0, oo) with inf_E v > 0: mu_d(E) <= (W(v) - lambda_d v(0))/inf_E v.
Proof: Theorem 2.2 with h = v; b_d(v) >= 0 by 1.3(d) if d'' = d, = 0 if d'' < d. []
Cases: v = gt gives the brief's int gt dmu_d <= kt - lambda_d (correct). v = g_0 gives, uniformly in d,
mu_d(E) <= (0.0240444 - 0.3197518 lambda_d)/inf_E g_0 (W(g_0) and g_0(0) from certs/README.md).

**Corollary 3.3 (weights are stable). PROVED UNDER HYPOTHESIS (K).** (K): d' = 2K even, G = G_{d'} has the exact dual
(lambda', mu' = sum_{j<=K} mu_j delta_{u_j}, nu' = 0, alpha' = 0), and the Hermite system at u_1..u_K is non-singular,
so the Hermite cardinal function h_k in V_{d'} exists: h_k(u_j) = delta_jk, h_k'(u_j) = 0, b_{d'}(h_k) = 0.
Then W(h_k) = mu_k + lambda' h_k(0) ((2.1) at d'). Let eta > eta_0 := lambda'/W(h_k) with eta B_{h_k} <= B_G,
S := {eta h_k > G}, and eps > 0 with G + eps h_k in C_{d'}. For every even d >= d', every dual solution of degree d
and every Borel E with S c E and inf_E h_k > 0:
   (mu_k + lambda' h_k(0) - lambda'/eta) / sup_S h_k  <=  mu_d(S) <= mu_d(E)  <=  (mu_k + lambda'/eps + lambda' max(h_k(0),0)) / inf_E h_k.
Proof. Lower: 3.1(c) with w = eta h_k - G, sup_S w <= eta sup_S h_k (G >= 0), dropping lambda_d(1 - eta h_k(0)) >= 0.
Upper: 3.2 with v = G + eps h_k: W(v) - lambda_d v(0) = (lambda' - lambda_d)(1 + eps h_k(0)) + eps mu_k
<= lambda'(1 + eps max(h_k(0), 0)) + eps mu_k (0 <= lambda' - lambda_d <= lambda'), inf_E v >= eps inf_E h_k. []
So the KKT mass that every higher-degree optimum puts in the k-th trap is mu_k(d') up to O(1/fac) (eta = fac eta_0),
O(lambda'/(eps mu_k)), O(lambda' |h_k(0)|/mu_k) and the variation of h_k on the trap. Numerically h_k(0) < 0 always.

## 4. T3

**Proposition 4.1 (T3(a), localised trap). PROVED given (H1)-(H4).** Let d' be even, gt in C_{d'}, gt(0) = 1,
kt = W(gt), u* > 0, rho > 0, h in V_{d'}, eta > 0 with
(H1) eta B_h <= B_gt on R;  (H2) eta W(h) > kt;  (H3) eta h <= gt on (0, oo) \ (u* - rho, u* + rho);
(H4) on I = [u* - rho, u* + rho]: gt'' >= m and h'' <= M, with kappa := (m - eta M)/2 > 0.
Let a0 = eta h(u*) - gt(u*), a1 = eta h'(u*) - gt'(u*), r := (|a1| + (a1^2 + 4 kappa a0^+)^{1/2})/(2 kappa). If r < rho,
then S c (u* - r, u* + r), and for every even d >= d' every minimiser G_d has a double zero in (u* - r, u* + r)
carrying KKT mass >= (eta W(h) - kt)/sup_S w.
Proof. For |s| <= rho, w(u* + s) <= a0 + a1 s - kappa s^2 (Taylor, w'' <= eta M - m = -2 kappa), and the right side is
< 0 for |s| > r (the roots of kappa s^2 - a1 s - a0 have modulus <= r). With (H3), S c (u* - r, u* + r); apply 3.1. []
For the exact optimum (gt = G_{d'}, u* = u_k, h = h_k): a0 = eta, a1 = 0, r = (2 eta/(m - eta M))^{1/2}, i.e. the brief's
radius (2 eta/G''(u_k))^{1/2} to first order, with the threshold eta_0 = lambda'/(mu_k + lambda' h_k(0)), not lambda'/mu_k.

Numerical check of (H1)-(H3) for gt = G_{d'}, h = h_k (NUMERICAL; trap2.py, eta = 2 eta_0, grid step 0.005 for
d' <= 8 and 0.01 for d' >= 12 on (0, u_K + 0.6 .. 1], r from the formula above with m = G''(u_k), M = h_k''(u_k)):
(H1) B_G - eta B_{h_k} has no root t^2 >= 0 in every case. The prime sum W(h_k) agrees with mu_k + lambda' h_k(0) (the
KKT identity of degree d') to the 10 digits of the stored mu_k. h_k(0) < 0 in every case (-71 at k = 1, down to
-3e10 at d' = 16, k = 8); lambda'|h_k(0)| <= 0.094 mu_k for d' >= 8 (0.23 mu_k at d' = 6, k = 3). On the grid S is one interval around u_k (or no grid point when
r < step), and the margin min G/(eta h_k) over grid points with h_k > 0 and |u - u_k| > 1.5 r + step is > 1. Beyond
the grid, w = eta h_k - G < 0 for large u by 1.1 (top index d' even, b_{d'}(w) = -b_{d'}(G) < 0), with the threshold
not computed. A grid is a sample, not a proof of (H3):
| d' | trap half-widths r_k (k = 1, 2, ...) | log n_k inside every trap | smallest margin |
|---|---|---|---|
| 8 | 3.7e-3, 1.6e-2, 7.0e-2, 8.6e-2 | yes (n_k = 2, 3, 4, 5) | 1.44 (k = 3) |
| 12 | 8.1e-5, 3.5e-4, 1.7e-3, 2.6e-3, 1.8e-2, 8.5e-2 | yes; the last trap (e^u in [7.67, 9.10]) holds 8 and 9 | 1.46 (k = 6) |
| 16 | 2.0e-6, 8.8e-6, 4.4e-5, 7.6e-5, 1.1e-3, 9.3e-3, 1.6e-2, 3.3e-2 | yes (2, 3, 4, 5, 7, 8, 9, 11) | 2.40 (k = 8) |
| 24 | 3.0e-13, 1.3e-12, 6.6e-12, 1.2e-11, 2.0e-10, 2.0e-9, 4.6e-9, 2.5e-8, 3.3e-7, 1.1e-4, 2.6e-4, 1.5e-3 | yes (2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19); u_k - log n_k from 2.3e-25 (n = 2) to 7.9e-6 (n = 19) | 75.2 (k = 12) |

Read at d' = 24 (NUMERICAL, modulo (H3) beyond the grid and modulo certification): every class-A optimum of every
degree >= 48 has a double zero within 3.0e-13 of log 2 + 2.3e-25, within 1.3e-12 of log 3, ..., within 1.5e-3 of
log 19.0001496; its KKT weight lies between W(h_k)/2 and about mu_k (3.3 with fac = 2). At d' = 24, h_k(0) reaches -5.3e21 (k = 12), but
lambda'|h_k(0)| <= 1.3e-4.
Mass bounds of 3.3 at d' = 16 (eta = 2 eta_0, eps from a grid estimate of the largest admissible eps): for every
d >= 16, 0.490 <= mu_d(S_1), mu_d(I_1) <= 0.980258; 0.634 <= mu_d(S_2), mu_d(I_2) <= 1.26857; 0.3466 <= mu_d(S_3),
mu_d(I_3) <= 0.693147; 0.720 <= mu_d(S_4), mu_d(I_4) <= 1.43953 (I_k = (u_k - r_k, u_k + r_k); the upper ends equal
2 Lambda(n_k)/sqrt(n_k) = 0.980258, 1.268568, 0.693147, 1.439525 to the digits shown; the lower ends are W(h_k)/2 because
fac = 2, and tend to W(h_k) as fac grows, at the price of traps wider by fac^{1/2}). Caveat: S_k c I_k holds only to
first order in r_k; for the wide last traps it can fail (d' = 8, k = 3: S on the grid reaches 1.485 > u_3 + r_3 = 1.4675).
At d' = 24 the upper ends are 0.980258, 1.26857, 0.693147, 1.43953, 1.47097, 0.490129, 0.732408, 1.44599, 1.42278,
0.346592, 1.3744, 1.35364 against 2 Lambda(n)/sqrt(n) = 0.980258, 1.268568, 0.693147, 1.439525, 1.470970, 0.490129,
0.732408, 1.445985, 1.422778, 0.346574, 1.374310, 1.351001 (n = 2, ..., 19).

**Theorem 4.2 (T3(c)). PROVED.** Let g be even, continuous, g >= 0, P(g) < oo (e.g. g in C_d). For every prime power
n: g(log n) <= sqrt(n) P(g)/(2 Lambda(n)). If J is an interval containing log n and g'' >= m > 0 on J, the minimiser z
of g on J satisfies |z - log n| <= (sqrt(n) P(g)/(Lambda(n) m))^{1/2}; a zero of g in J, if any, is z, and it is unique.
Proof. All terms of P(g) are >= 0. g is strictly convex on J, so g'(z)(x - z) >= 0 for x in J, and
g(log n) >= g(z) + g'(z)(log n - z) + (m/2)(log n - z)^2 >= (m/2)(log n - z)^2. A zero of g >= 0 in J is a minimiser. []
Example (NUMERICAL): d = 24, lambda = 2.3758e-26, G''(u_1) = 1.056: |u_1 - log 2| <= 2.1e-13; observed 2.3e-25. The bound
is far from sharp because n = 2 carries almost none of the cost (n = 23 carries 99.9 %): G(log 2) is ~1e-50, not ~lambda.

**Theorem 4.3 (T3(b)). PROVED UNDER (A), (B), (C).** Fix a prime power n and J = [log n - rho, log n + rho]. Assume
(A) lambda_d -> 0 (d even); (B) there are m > 0 and d_0 with G_d'' >= m on J for every minimiser of every even d >= d_0;
(C) there is d_1 such that every minimiser of every even d >= d_1 has a zero in J. Then for d >= max(d_0, d_1) G_d has
exactly one zero z_d in J, and |z_d - log n| <= (sqrt(n) lambda_d/(Lambda(n) m))^{1/2} -> 0.
Proof: Theorem 4.2 with P(G_d) = lambda_d. []
Without (A) and (B): if (C) holds for EVERY window J around log n, then dist(log n, zeros of G_d) -> 0, with no rate.
What the hypotheses need:
- (C) is implied by ONE localised trap (4.1) whose interval lies in J, at one degree d_1 (a finite check). NUMERICAL
  for the first 12 prime powers (n <= 19) from d' = 24 (table above); not certified.
- (B) is open. Reduction (PROVED, 5.1(f)): if (A) holds and sup_d G_d''''(0) < oo, then along subsequences G_d'' -> G_inf''
  uniformly, so (B) follows from G_inf''(log n) > 0 for every subsequential limit. NUMERICAL: -G_d''(0) = 16.07, 18.13,
  18.45, 18.53, ..., 18.548224 and G_d''''(0) = 680.0, 874.4, 920.3, 940.8, ..., 945.654 for d = 2, 4, 6, 8, ..., 24
  (moments.out): bounded and converging. Curvature at the zeros, G_d''(u_k(d)) for d = 8, 12, 16, 24 (trap2 outputs):
  n = 2: 1.062, 1.056, 1.056, 1.056; n = 3: 0.0468, 0.0436, 0.0434, 0.0433; n = 4: 4.14e-3, 3.41e-3, 3.18e-3, 3.14e-3;
  n = 5: 1.31e-3, 6.9e-4, 5.2e-4, 4.9e-4; n = 7: -, 1.41e-5, 2.48e-6, 1.63e-6. Decreasing in d, apparently to positive
  limits for n <= 5; for n = 7 not yet settled. So (B) is plausible for the early prime powers, with m small.
- (A) is open; see 4.4.

**Proposition 4.4 ((A) from the Hermite family). PROVED (implication).** Suppose there are infinitely many J with
(a) M_J non-singular, (b) every coefficient p_k^(J) of P_J positive and a_J := max_k ((2k)! p_k^(J))^{1/2k} <= a* < pi/2,
(c) Fhat_J >= 0 on R. Then lambda_d -> 0; more precisely lambda_{2J} <= 2 pi kappa_J <= C(a*) sum_{n > n_J} Lambda(n) n^{-1-delta}
for any fixed delta > 0.
Proof. (b) gives H_J = P_J(t^2) >= 1 on R and |H_J(z)| <= sum p_k |z|^{2k} <= cosh(a* |z|) on C. Let g_J have transform
Xi^2 H_J. By (c) and H_J >= 0, g_J/g_J(0) is in C_{2J}, and g_J(0) = (1/2pi) int Xi^2 H_J >= g_0(0) = 0.31975. For u > 0
move the line of integration of g_J(u) = (1/2pi) int Xi(t)^2 H_J(t) e^{-itu} dt to Im t = -(1/2 + delta): with
|Xi(x+iy)|^2 <= C (1+|x|)^6 e^{-pi |x|/2} for |y| <= 1/2 + delta (Stirling and the convexity bound; Johns I Lemma 2.1)
and |H_J| <= cosh(a*(|x| + 1)), |g_J(u)| <= C(a*, delta) e^{-(1/2+delta) u}, uniformly in J. g_J vanishes at the first J
prime powers (Prop. 4.5 of Johns I), so W(g_J) = 2 sum_{n > n_J} Lambda(n) n^{-1/2} g_J(log n) <= 2C sum_{n > n_J}
Lambda(n) n^{-1-delta} -> 0. []
(a), (b) are Johns I Conjecture 6.1(a), (b) restricted to those J. (c) is the finite-J cone membership certified at
five J (Prop. 5.7); it is NOT part of Conjecture 6.1, whose (c) is about the limit on [xi_2, oo) only. Johns II
Theorem 1 does not give (A): its cofactor H is entire, not a polynomial, and Part II (sec. 1.5) claims nothing on
approximating it inside the polynomial family while keeping Fhat >= 0.

**Proposition 4.5 (gt = F_J using only Props 5.6/5.7). PROVED for (a)-(c); localisation NOT PROVED.**
For J in {10, 60, 61, 110, 111} let gt := g_J/g_J(0) (transform Xi^2 P_J(t^2)/g_J(0)); by Prop. 5.7 gt is in C_{2J},
gt > 0 off the nodes +-log n_j (j <= J), double zeros there, and kt = W(gt) <= 2 pi kappa_J. For every even d >= 2J
(degree >= 4J in t), every minimiser and every dual solution:
(a) lambda_d <= 2 pi kappa_J (8.98e-1060 at J = 111).
(b) int gt dmu_d <= 2 pi kappa_J - lambda_d (3.2); so mu_d(E) <= 2 pi kappa_J/inf_E gt, and inf_E gt > 0 for every compact
E avoiding the nodes; no value of inf_E gt is stated in Part I.
(c) With h = g_0: since P_J >= 1 (Prop. 5.6), (i) holds for eta <= 1/g_J(0), and (ii) for eta > 2 pi kappa_J/W(g_0) =
261.3 kappa_J. So G_d has a zero in {eta g_0 > gt}: a set containing a neighbourhood of every node, possibly other
points, bounded (top index 2J even, coefficient < 0), and not computable from the stated facts.
What is missing for "one zero within r_k of log n_k for every k <= J and every degree >= 4J": a test function h_k
concentrated at log n_k (the Hermite cardinal function at the exact nodes) and the numbers in (H3)-(H4), i.e.
quantitative lower bounds of gt between the nodes. Prop. 5.7(b) certifies such bounds (lower bounds of
Fhat_J(xi_n + delta)/delta^2 and on Taylor pieces) but does not state them, so the localisation needs a
recomputation. Even then the radii are fixed by J: convergence needs r_k(J) -> 0, i.e. kappa_J -> 0 together with a
curvature bound at log n_k, which are (A) and (B) again.

## 5. T4: the unconditional statement

**Proposition 5.1. PROVED (uses (EF)).** Let G_d (d even) be minimisers, F_d := Ghat_d = Xi^2 B_d, lambda_inf := lim lambda_d.
(a) For |Im z| <= 1/2: |F_d(z)| <= F_d(i |Im z|) <= F_d(i/2) = B_d(i/2)/4 = pole_d/2 <= (lambda_d - Omega(0))/2 <= 2.7237,
where pole_d = 2 int G_d cosh(u/2) du, Omega(t) = Re psi(1/4 + it/2) - log pi >= Omega(0) = -5.37218.
(b) (1/2pi) int F_d (Omega - Omega(0)) dt <= lambda_d - Omega(0) (tightness).
(c) Every subsequence has a further subsequence with F_d -> F, B_d -> B_inf locally uniformly on |Im z| < 1/2, F = Xi^2 B_inf,
F_d -> F in L^1(R), int F = 2 pi, and G_d -> G_inf uniformly on R.
(d) G_inf >= 0, Ghat_inf = F >= 0, G_inf(0) = 1, F has a double zero at every non-trivial zero of zeta (on or off the
line), P(G_inf) <= lambda_inf, and pole + A at G_inf is <= lambda_inf.
(e) Every limit point of zeros of G_d (along the subsequence) is a zero of G_inf. If lambda_inf = 0, G_inf(log n) = 0
for every prime power n.
(f) If sup_d G_d''''(0) < oo, then G_d'' -> G_inf'' uniformly on R along the subsequence.
Proof. (a) G_d >= 0 gives |F_d(x+iy)| <= int G_d e^{-yu} du = F_d(iy) = int G_d cosh(yu) du, increasing in |y|;
Xi(i/2) = xi(0) = 1/2. By (EF), pole_d = lambda_d - A(G_d) and A(G_d) = (1/2pi) int F_d Omega >= Omega(0) G_d(0)
(Omega is increasing in |t|, F_d >= 0). lambda_d <= lambda_0 = 0.07520. (b) the same with A(G_d) = lambda_d - pole_d
<= lambda_d. (c) Montel on the open strip. Zeros of Xi lie in |Im t| < 1/2 (zeta has no zeros on Re s = 0, 1) and
Xi(+-i/2) = 1/2; on small circles around them |B_d| <= max|F_d|/min|Xi|^2, so B_d is locally bounded on the open strip
(maximum principle) and converges with F_d. (b) and local uniform convergence give int F = 2 pi; Scheffe gives L^1
convergence; |G_d - G_inf| <= ||F_d - F||_1/(2 pi). (d) limits of non-negative functions; Fatou for P and for the pole
term; A is lower semicontinuous (Omega continuous, bounded below). (e) uniform convergence; P(G_inf) <= 0 and every
term is >= 0. (f) int t^4 F_d = 2 pi G_d''''(0) bounded makes t^2 F_d uniformly integrable; Vitali. []
NUMERICAL: pole_d = 1.2050, 1.1301, 1.1217, 1.1204, 1.1201, ..., 1.1200979 for d = 2, 4, 6, 8, 10, ..., 24
(moments.out); the bound 2.72 in (a) is crude. (For comparison, the pole part of Johns' Part II function measured in
bootstrap.md 5f is 1.120099144; the two limits need not agree, his function vanishes at log 6 and ours does not.)

**5.2. What the unconditional statement needs (NOT PROVED).** "The double zeros of G_d converge to the prime powers"
splits into (I) for every prime power n and every window J around log n, every G_d of large degree has a zero in J,
and (II) every limit point of zeros of G_d is some log n. The missing statements:
(M1) lambda_d -> 0. Open; equivalent to "inf of W(h)/h(0) over the union of the C_d is 0"; implied by 4.4.
(M2) Non-degeneracy at log n: liminf_d min_J G_d'' > 0 (only for the rate in 4.3). Implied by (M1), bounded G_d''''(0)
and G_inf''(log n) > 0 for every limit (5.1(f)). Numerically plausible (the moments converge; curvature table in 4.3).
(M3) A trap for every window: for every n and J there are d' and w in V_{d'} with B_w <= 0, P(w) > 0, {w > 0} c J.
By 3.1, (M3) gives (I) at once, WITHOUT (M1): every G_d with d >= d' has a zero in J. (PROVED implication.)
(M4) Every subsequential limit G_inf (5.1(c), which needs no (M1)) is > 0 off {+-log n}. By uniform convergence,
(M4) gives (II). (PROVED implication.) (M4) is not implied by (M1): Johns' Part II function has cost 0 and vanishes
at log 6. NUMERICAL: G_d(log 6)/sqrt 6 levels off at 1.4e-7 (certs/README.md), so this family's limit seems not to.
Rate: (M1) + (M2) + (M3) give |z_d - log n| = O(lambda_d^{1/2}) (Theorem 4.3).
The determinacy statement behind (M1) and (M3):
(Det) If lambda >= 0, mu >= 0 on (0, oo), nu >= 0 on [0, oo) and P(h) = lambda h(0) + int h dmu + int B_h dnu for every
h in V = union of the V_d, then lambda = 0, nu = 0, mu = P.
(Det) => (M1) and (Det) => (M3): SKETCH, NOT PROVED. For (M1): a vague limit of the KKT data (lambda_d, mu_d, nu_d)
represents P on V with lambda = lim lambda_d. For (M3): if no trap exists in V_{d'} for J, separation in R^{d'+1} gives a
non-negative representation of P on V_{d'} whose mu has no mass in J. In both cases Lemma 3.2 (applied to
representations of P, same proof) bounds the masses locally uniformly, and test functions v = g_0 + eps g_{2j} of higher
degree control the tails, so a vague limit represents P on all of V; the tail step is not written out.
(Det) PROVED in the special case nu = 0 and sigma := mu - P a finite signed measure. Let sigma_e be its even extension.
Fubini gives (1/2pi) int Xi^2 B (sigmahat_e + 2 lambda) dt = 0 for every polynomial B. phi := sigmahat_e + 2 lambda is
bounded and continuous, so Psi(z) = int phi Xi^2 e^{izt} dt is analytic on |Im z| < pi/2 (Xi^2 = O(t^6 e^{-pi|t|/2})),
and all its Taylor coefficients at 0 vanish; so Psi = 0, phi Xi^2 = 0 a.e., phi = 0 (Xi has countably many real zeros,
phi is continuous). sigmahat_e = -2 lambda is the transform of -2 lambda delta_0, and sigma_e has no atom at 0: lambda = 0,
sigma = 0. []
The gap: P has infinite mass (P((0, U]) ~ 4 e^{U/2}), and a limit of KKT measures is only known to satisfy
int v dmu <= W(v) for v in the cone, which allows mu - P to grow like 1/g_0. V does not see a measure whose transform is
carried by the zeros of Xi (sum_rho a_rho cos(gamma_rho u) + b_rho u sin(gamma_rho u) integrates to 0 against every h in V).
A finite measure has no such transform; a growing one can (P minus the archimedean density is one, by the explicit
formula, but that density is not of the form lambda delta_0 + positive measure: it is singular and negative near u = 0).
That is the precise content of the brief's remark. A non-negative tempered f whose transform avoids (-14.13, 14.13)
is 0 (test against |phi_a|^2 with phihat_a supported in (-7, 7)), but mu - P may be negative at the prime powers, so
this does not settle (Det). (Det) is a uniqueness statement for positive solutions of the explicit formula on the
zero-killed class V; Johns II Theorem 2 concerns the full test class and does not imply it, as far as I can see.
The brief's diagnosis (lambda_d -> 0 and a uniform curvature bound) is right for the RATE, but the curvature bound is
not what is missing for convergence itself: (I) needs (M3) (existence of a zero near log n, which neither (M1) nor
(M2) provides) and (II) needs (M4).

## 6. Corrections to the brief

1. g_d denotes both the basis element and the minimiser; here the minimiser is G_d.
2. Phi is analytic only in |Im u| < pi/4, not entire (harmless).
3. T2: (iii) follows from (i); under (i) eta h(0) <= 1 always (w is negative definite), so (ii) reduces to
   eta W(h) > kt and the bound becomes (eta W(h) - kt + lambda_d (1 - eta h(0)))/sup_S w; lambda_d <= kt is not used.
4. T2 localises only if the top index of w is even; with eta b_{d'}(h) = b_{d'}(gt) the set S can contain [U, oo).
5. T3(a): the threshold is eta_0 = lambda'/(mu_k + lambda' h_k(0)) with h_k(0) < 0, not lambda'/mu_k.
6. T3(a): the trap of the newest point is wide and can contain two prime powers (d' = 6: e^u in [3.7, 5.2], 4 and 5;
   d' = 12: e^u in [7.67, 9.10], 8 and 9), so "the k-th zero stays near log n_k" holds for the older points only.
7. T4: the normal-family bound needs the explicit formula (pole term = W - A), contrary to "you do not need that
   identity"; the constants of the brief (ghat(i/2) = B(i/2)/4, B(i/2) = 2 x pole term) are right.
8. T4: the polynomial test class does determine finite signed perturbations of P; the obstruction is growth.
9. T1 at odd d: the b_d multiplier is unsigned (points at both infinities), 0 is in conv(Gamma).
10. T3(b) with gt = F_J: the stated facts of Props 5.6/5.7 give mass bounds but no localisation.

## 7. Files (this folder)

trap.py, trap_small.out: first scan, d' = 2, 4, 6, 8 (eta = 2 eta_0, KKT identity W(h_k) = mu_k + lambda' h_k(0)).
trap2.py, trap2_d8.out, trap2_d12_16.out, trap2_d24.out: traps, (H1) root check, grid margins, eps estimates, mass bounds
(usage: python trap2.py d' dps fac step u_extra; runs of 4 to 30 min at 60 to 100 digits).
moments.py, moments.out: -G_d''(0), G_d''''(0), B_d(i/2) for d = 2, ..., 24. asym.py, asym.out: check of Lemma 1.1.
