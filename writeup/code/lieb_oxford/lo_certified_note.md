# A certified Lieb–Oxford bound from the Lewin–Lieb–Seiringer reduction

Status labels: **[cited]** taken from Lewin–Lieb–Seiringer, *Improved Lieb–Oxford bound on the indirect and exchange
energies*, Lett. Math. Phys. 112 (2022), arXiv 2203.12473v3 (read in full, LaTeX source v07); equation numbers are
the paper's. **[proved]** proved in this note by hand or by exact symbolic computation (sympy, exact rationals).
**[CA]** computer-assisted: rigorous ball arithmetic (Arb via python-flint 0.9.0, 128-bit), every rounding enclosed.
**[float]** floating-point cross-check only, not part of any proof.

## Theorem A (uniform ball)

c_LO <= 1.604335.  More precisely, with mu = nu = B (uniform probability measure of the unit ball) and f = g of
Lemma 10, the ball-arithmetic enclosure gives J <= 1.6993342175 and c_LO <= 1.60433487882 (outward rounded: 1.604335).

## Theorem B (spheres)

c_LO <= 1.597856.  More precisely, with mu = sum_{j=1}^{48} mu_j sigma_{j/48}, nu = sigma_{2839/4096} (exact
rational weights listed in the script, MEASURE) and f = g of Lemma 10: J <= 1.37702535712, D(mu,mu) =
0.662496242827... (exact rational), c_LO <= 1.59785555399 (outward rounded: 1.597856).

## 1. The chain of statements

Notation as in the paper: D(n1,n2) = (1/2) iint dn1 dn2/|x-y|;
Phi_{mu nu}(a,b) = a^3 b^3 (1 - iint dmu(u) dnu(v)/|e1 + u/a - v/b|) (19);
Psi = Phi_{mu nu} + Phi_{nu mu} - Phi_{nu nu} (p. 8); F_{mu nu} = { f in C^0(R+,R+) : Psi(a,b) <= f(a)+f(b) } (24).

S1 **[cited]** (Section 2, unnumbered display on p. 9, derived from (18), (22), (23)): for radial probability measures
mu, nu with D(mu,mu) < oo and f in F_{mu nu},
  E_ind[rho] >= -( int f(|z|) |z|^-7 dz + D(mu,mu) ) int rho^{4/3}.
S2 **[proved here, also in the paper p. 9]** scaling: mu_t = t^3 mu(t.), nu_t likewise, f_t = t^-6 f(t.). Then
Psi_{mu_t nu_t}(a,b) = t^-6 Psi_{mu nu}(ta,tb) (substitute u -> u/t in (19)), so f_t is in F_{mu_t nu_t};
int f_t(|z|)|z|^-7 dz = t^-2 J_f and D(mu_t,mu_t) = t D(mu,mu). Minimising t^-2 J + t D over t > 0
(t = (2J/D)^{1/3}) gives c_LO <= (3/2)(2 J_f D(mu,mu)^2)^{1/3}, i.e. eq. (26) for one admissible triple.
S3 **[cited, Lemma 10, eq. (51); short proof repeated]** if mu, nu are supported in the closed unit ball and
mu({0}) = 0, then g(a) = max_{0<=b<=a} { Psi(a,b) - Psi(b,b)_+/2 } is in F_{mu nu}. Proof: Psi is continuous on
[0,oo)^2 and vanishes on the axes (Lemma 8), so g is continuous (maximum of a continuous function over [0,a]) and
g >= 0 (b = 0), g(b) >= Psi(b,b)_+/2 (take b' = b). Psi is symmetric: Phi_{mu nu}(b,a) = Phi_{nu mu}(a,b) by the
substitution (u,v) -> (-v,-u) and radial symmetry, and Phi_{nu nu} is symmetric. For b <= a,
Psi(a,b) <= g(a) + Psi(b,b)_+/2 <= g(a) + g(b); by symmetry for all a, b.
S4 **[classical]** D(B,B) = 3/5 (self-energy 3/(5R) of a uniform unit charge ball of radius R = 1); checked exactly
in the script.

Hence c_LO <= (3/2)(2 (9/25) J)^{1/3} with J = 4 pi int_0^oo g(r) r^-5 dr. Only J needs work.

Not used: the infimum over f (I(mu,nu)), Theorem 12, the truncated problem I_R, the discretisation (73)-(77),
the unproved tail rule (70). The paper's 1.58 is not used.

## 2. Psi for the ball: proof of (78) and (43) **[proved]**

Put R1 = 1/a, R2 = 1/b and W = iint dmu(u) dnu(v)/|e1 + u/a - v/b| = interaction of two unit charges, the first
(radius R1) centred at 0, the second (radius R2) centred at e1. Then Phi = a^3 b^3 (1 - W). Write
W = int_0^oo p(s) V1(s) ds, with V1 the (radial) potential of the first charge and p the density of s = |y| for y
distributed like the second charge. Archimedes: for y uniform on the sphere of radius R2 about e1, |y|^2 is
uniform on [(1-R2)^2,(1+R2)^2], so p(s) = s/(2 R2) on [|1-R2|, 1+R2]. For the uniform ball: p(s) = 3 s^2/R2^3 for
s <= R2 - 1 (whole sphere |y| = s inside) and p(s) = 3 s (R2^2 - (s-1)^2)/(4 R2^3) on [|1-R2|, 1+R2] (spherical cap
area pi s (R2^2-(s-1)^2) over the volume). Potentials: sphere 1/max(s,R1); ball (3R1^2 - s^2)/(2R1^3) for s <= R1,
1/s beyond (Newton). The integrals are polynomial; sympy evaluates them exactly in every configuration
(disjoint R1+R2 <= 1; partial overlap |R1-R2| < 1 < R1+R2 with R2 <= 1 and R2 >= 1; ball 1 inside ball 2;
ball 2 inside ball 1 with R2 <= 1 and R2 >= 1) and the differences with (78) simplify to 0. In each configuration
the inactive positive parts vanish: x1 = a+b-ab = ab(R1+R2-1), x2 = |a-b|-ab = ab(|R1-R2|-1). Same for (43)
(six configurations). Two further independent numerical routes agree with (78) at 16 points (Fourier route
W = (2/pi) int_0^oo F(kR1) F(kR2) sin(k)/k dk with F(z) = 3(sin z - z cos z)/z^3, to 1e-12 relative to a^3b^3)
and at 4 points (mixture of spheres with (43), scipy dblquad) **[float]**. Points next to the special lines
a+b-ab = 0, |a-b| = ab and a = b are included.

Consequences used below, for b <= a **[proved, exact]**:
 * x2 >= 0 (ball 1/a inside ball 1/b, i.e. b <= a/(1+a)): Psi = PA := a^3 b^3 (1 - 3b/2 + b^3/2) + 3 a b^6/10
   (T1 - T2 = PA identically; also directly: the average of (3R2^2 - |x-e1|^2)/(2R2^3) over a ball of radius R1
   about 0 is (3R2^2 - 1 - 3R1^2/5)/(2R2^3)).
 * x2 <= 0 <= x1: Psi = T1 = (a+b-ab)^4 (a^2b^2 - 5a^2 - 5b^2 + 4ab^2 + 4a^2b + 20ab)/160.
 * x1 <= 0 (b >= a/(a-1), a > 1): Psi = 0.
 * Psi(b,b)/2 = Delta(b) := b^6 (2-b)_+^4 (b^2 + 8b + 10)/320.
 * Psi is C^3 on (0,oo)^2: the pieces differ by x1^4 or x2^4 times a polynomial.

## 3. Structure of g for the ball

h(a,b) = Psi(a,b) - Delta(b) (Psi >= 0 here since mu = nu, by (21)). **[proved]** d/db h(a,b) = 0 on b = a
(symmetry of Psi), and h_bb(a,a) = -2 a^4 q_B(a,1) with q_B(a,1) = 9(a-2)^2 (a^4+4a^3-8a^2-8a+10)/80. So b = a is a
strict local maximum of h(a,.) for a < a_c and a local minimum for a_c < a < 1.624, where
a_c = 0.8890423748... is the smallest positive root of a^4+4a^3-8a^2-8a+10 (derived here; the analogous point for
the exchange problem in the paper is 3/8). For large a, Psi(a,b)/a^3 -> k(b) = b^3(1 - 3b/2 + b^3/2), whose maximum
on [0,1] is at b = (sqrt5-1)/2 with value kappa = (5 sqrt5 - 11)/4 = 0.0450849718... (eq. (53) made explicit;
k'(b) = 3b^2(b-1)(b^2+b-1)).

## 4. The integral J = 4 pi (S0 + S1 + S2 + S3)

**S0, a in [0, a0], a0 = 1/1024 [proved].** Psi(a,b) <= a^3 b^3 (W >= 0) and Psi(b,b) >= b^6 (1 - 3b/2)
(W <= b V_B(b e1) = b(3-b^2)/2 by (21)); a^3b^3 - b^6(1-3b/2)/2 has b-derivative 3b^2(a^3-b^3) + 21b^6/4 >= 0 on
[0,a], so g(a) <= a^6/2 + 3a^7/4 and S0 <= a0^2/4 + a0^3/4.

**S1, a in [a0, 7/8]: g(a) = Psi(a,a)/2 [CA].** With b = a u, D(a,u) := h(a,a) - h(a,au) =
Delta(a) - Psi(a,au) + Delta(au). Exactly: on u >= 1/(1+a) (piece T1), D = a^6 (1-u)^2 q_B(a,u); on
u <= 1/(1+a) (piece PA), D/a^6 = (1-u)^2 q_A(a,u) + r_A(a,u) (polynomial division, r_A linear in u). Ball
arithmetic on 1768 boxes of [1/1024, 7/8] x [0,1] (adaptive bisection, region-restricted u ranges, rounded
outward) shows q_B > 0 and (1-u)^2 q_A + r_A > 0. Hence h(a,b) <= h(a,a) for 0 <= b <= a, so g(a) = Psi(a,a)/2
and S1 = int_{a0}^{7/8} Delta(r) r^-5 dr is an exact rational (0.0761888789...).

**S2, a in [7/8, 100] [CA].** Geometric grid with ratio 1 + 2.5e-4 (18958 intervals, exact dyadic endpoints).
On I = [A0,A1] with midpoint am and half-width d, for every b in [0,A1] Taylor's theorem in a gives
H(a,b) <= H(am,b) + (a-am) H_a(am,b) + Haa+ (a-am)^2/2, where H is the piecewise formula (equal to h for b <= a)
and Haa+ >= max(0, sup_{I x [0,A1]} Psi_aa) is an enclosure computed on blocks of 64 intervals (factored form of
T1_aa, 4 x 1024 boxes, support cut b >= P/(P-1) where Psi = 0). G(s) = max_b [H(am,b) + s H_a(am,b)] is convex
in s (maximum of affine functions), so on I
   g(a) <= chord through (-d, G^(-d)), (d, G^(d)) + Haa+ d^2/2,
where G^(+-d) >= G(+-d) are certified by branch and bound in b: each box [lo,hi] gets the bound
F(c) + |F'(c)| r + max(sup_box F'', 0) r^2/2 (F is C^3; pieces selected by interval tests of x1, x2, hull of the
possible pieces, which agree to third order on the switching curves). The returned value is the exact maximum of
the bounds of all boxes left in the queue; relative gap to the best lower value <= 1e-10. The chord and the
remainder are integrated exactly against r^-5.

**S3, a >= R = 100 [proved].** t = 1/a. Region A (b <= 1/(1+t)): Psi/a^3 = k(b) + 0.3 t^2 b^6. Region B
(1/(1+t) <= b <= 1/(1-t)): the potential of ball 2 at distance <= 1+t is >= 1/(1+t), so Psi/a^3 <= t/((1-t)^3(1+t)).
Region C: Psi = 0. Since h <= Psi, k <= kappa, k decreasing on [0.618,1] and k(3/4) + 0.3/R^2 < kappa,
t/((1-t)^3(1+t)) < kappa at t = 1/R: g(a) <= kappa a^3 + 0.3 (3/4)^6 a, so S3 <= kappa/R + 0.1 (3/4)^6/R^3.

## 5. Numbers (script lo_ball_certify.py, run with rho = 2.5e-4)

 S0 <= 2.3865e-7, S1 = 0.076188878925 (exact rational), S2 <= 0.058588735243, S3 <= 0.000450867517,
 J <= 1.69933421745, c_LO <= 1.60433487882 -> 1.604335 (outward). Float cross-check (independent quadrature of a
 float g): J = 1.6993339650, c = 1.6043347994; the certificate exceeds it by 1.5e-7 relative. Paper's Table 2
 (discretised LP, same measures): 1.604311 ... 1.604377 (with f = g = the LP optimum there, per the paper's
 Fig. 2 caption). Sanity sampling [float]: Psi_aa <= Haa+ at 11880 random points; every 997th interval the
 certified G(+-d) exceeds a 4001-point float maximum. Output: 60 checks, 0 failed (out_ball_final.txt), run time
 115 s.

 The certified J exceeds the float J by 2.5e-7 (absolute). Rough split [heuristic estimates, not certified]:
 the S0 bound loses ~5e-9, the tail bound at R = 100 ~1.5e-7 (0.1(3/4)^6/R^3 against a true correction of about
 0.0056/R^3), the Taylor/chord remainder the rest; the float tail kappa/R is itself ~7e-8 too low. So the
 certificate is within ~1e-7 of the value of the method (mu = nu = B, f = g), which is 1.6043348.

## 6. Stage 2: mixtures of spheres (Theorem B)

**Measures.** mu = sum_j mu_j sigma_{j/K} (uniform probability measures on spheres of radius j/K), nu = sigma_s
(one sphere). Exact rational data in meas_K48.json (and inside lo_certify.py): K = 48, s = 2839/4096 =
0.693115234375, weights multiples of 2^-20 summing exactly to 1. They come from a float optimisation of the f = g
value (L-BFGS-B, then Powell, on a 210-point grid, warm start from K = 32; s2_opt.py) [float, heuristic: no claim
of optimality; two spheres for nu (K = 12) gave no gain]. Shape: the shell masses rise to a maximum 0.050 at
r = 9/16 and are then nearly flat (0.029-0.036 per shell, density ~ r^-2) up to r = 1; nu is one sphere of
radius 0.693 inside the support of mu.

**Psi [proved].** By bilinearity (the paper's (45)/(69), re-derived: Phi_{mu nu} = sum mu_j nu_k Phi_{sigma_j sigma_k}
since sum mu_j nu_k = 1), Psi = sum_{r,q} c_{rq} T_{rq} with c_{rq} = mu_r nu_q + nu_r mu_q - nu_r nu_q and
T_{rq}(a,b) = r^3 q^3 Phi_ss(a/r, b/q). From (43) (proved, six configurations), T_{rq} is
  I   a^3 b^3 (1 - b/q)                  if aq - br - ab >= 0   (sphere r/a about 0 inside sphere q/b about e1)
  II  a^3 b^3 (1 - a/r)                  if br - aq - ab >= 0   (the other way round)
  III a^2 b^2 (aq + br - ab)^2 / (4 r q)  if both are <= 0 <= aq + br - ab (overlap)
  0                                      if aq + br - ab <= 0  (disjoint).
Psi is C^1 with piecewise polynomial, bounded second derivatives (enough for every Taylor bound used). The
diagonal Psi(b,b) is piecewise polynomial with breakpoints |q - r|, q + r. Psi may be negative here (mu != nu), so
h = Psi - Psi(b,b)_+/2 = min(Psi, Psi - Psi(b,b)/2) and each box bound is the smaller of the two Taylor bounds.

**Certificate [CA] (cert_spheres.py).** Same scheme as for the ball, with: exact rational breakpoints in b at the
midpoint am of each a-interval (b = aq/(a+r), aq/(r-a), aq/(a-r) and the diagonal ones); on every sub-interval
between consecutive breakpoints Psi(am,.) and Psi_a(am,.) are single quartics (built by an event sweep, spot-checked
against direct exact evaluation); branch and bound inside sub-intervals; Haa+ on blocks of 32 a-intervals, with the
piece tests done at the box corners (x1, x2 are bilinear, so their extremes are at corners) and the factored form
b^2 (12 a^2 u^2 + 12 a u b r + 2 b^2 r^2)/(4rq), u = q - b, for piece III.
 S0 [proved]: for a <= 2^-12, Psi <= Phi_{mu nu} + Phi_{nu mu} <= 2 a^3 b^3 (Phi_{nu nu} >= 0 by (21)), so
 g <= 2a^6 and S0 <= a0^2.
 S3 [proved]: for a >= R = 600 and t = 1/a, piece I gives T/a^3 = L_q(b) := b^3 (1 - b/q)_+ exactly, piece 0
 gives T = 0 = L_q(b); only in the window b in [q/(1+rt), q/(1-rt)] (piece III) does T/a^3 differ from L_q, by at
 most r t b^3 max(1, b/q). Sum over pairs: Psi/a^3 <= Lambda(b) + t * (window terms), Lambda(b) =
 sum_k mu_k L_{rho_k}(b) = b^3 (1 - b V_mu(b e1)) (this is (53)). Windows are grouped by q and merged; on every
 cluster sup Lambda + C_cluster/R < kappa_mu (certified), so g(a) <= kappa_mu a^3 for a >= R and S3 <= kappa_mu/R,
 kappa_mu = sup Lambda = 0.0319135871783 (certified enclosure width 2e-15).

**Numbers (out_K48.txt; same in lo_certify_output.txt).** 48 atoms, 95 pairs with c != 0, 10101 a-intervals
(ratio 1 + 3e-3 below a = 1/4, 1 + 1e-3 above, up to R = 600), 372438 bisections, maximal relative
branch-and-bound gap 1e-10, 48 window clusters, none reaching kappa_mu = 0.031770351651336.
 S0 <= 5.96e-8, S2 <= 0.109527185984, S3 <= 5.2950586e-5, J <= 1.37702535712,
 **c_LO <= 1.59785555399 -> 1.597856** (outward). Float cross-check (independent numpy code built directly on (43),
 4800 a-points x 3001 b-points): J = 1.37702062, c = 1.59785372 (certified J is 4.7e-6 = 3.4e-6 relative higher, mostly the Haa
 remainder 4.2e-6 = 4 pi * 3.34e-7). 21 checks, 0 failed, 267 s.
 Other certified runs (same script): K = 32, s = 2843/4096: c_LO <= 1.597987 (out_K32.txt); K = 20: 1.598421
 (R = 200, older cruder tail bound, out_K20.txt); K = 12: 1.599768 (coarse steps). Float values of the f = g
 functional on the coarse optimisation grid: K = 12: 1.59923, 20: 1.59828, 32: 1.59788, 48: 1.59775 — the
 gain per doubling of K shrinks; with f = g and one sphere for nu this family seems to level off near 1.5975
 [heuristic extrapolation, not certified].

## 7. What is not proved here / what the theorems rest on

* Theorem 1 of the paper (the reduction (18) -> (22) -> unnumbered display p. 9 -> (26)) is cited, not re-proved.
  It is a refereed result (Lett. Math. Phys. 2022). Lemma 10 is cited and its short proof is repeated in S3 above.
  Nothing else from the paper is used: not the existence theory, not the numerics, not the unproved tail rule (70),
  not formulas (43)/(78) (both proved here).
* The symbolic steps rely on sympy (exact integration of polynomials, polynomial division, expansion); the
  computer-assisted steps rely on Arb/python-flint ball arithmetic and on the correctness of the scripts. Floating
  numbers enter only as exact dyadic endpoints or as heuristics for where to cut; every bound used is an arb upper
  endpoint.
* The measures of Theorem B are heuristic optima of the f = g functional on a coarse float grid; the certificate
  does not depend on how they were found. No optimality is claimed. The paper's 1.5765 uses the LP-optimal f (not g)
  and other measures; this note does not reach it and does not certify anything below 1.597856.
* Only f = g is used. For the ball the paper reports the discretised LP optimum equal to g (Fig. 2), so Theorem A is
  essentially the best this pair (B, B) gives; for sphere mixtures the LP f could be smaller than g (the paper
  says G^psi != G for its optimal measures).

## 8. Prior art

A web search (mathematical terms only) found a document titled "A certified explicit Lieb-Oxford constant below
1.64", dated 13 May 2026, on a private host (victoria.tail921ea4.ts.net) that did not resolve (DNS failure, two
URLs tried). Only search-engine snippets were seen, not the document. According to the snippets it works in the
LLS framework with ball smearing, certifies with rational interval arithmetic that 4 pi int g(r) r^-5 dr <= 1.81498
(< 1.8150) and concludes c_LO < 1.64 (with D = 3/5 this gives 1.5 (2 * 1.8150 * 0.36)^(1/3) = 1.63990), and states
that it does not claim optimality. Theorem A (1.604335) and Theorem B (1.597856) are below that value. Not found in
print otherwise: arXiv listing searches for Lieb-Oxford (done by the reader agent, files in ../lo/search) and two web
searches did not show another certified value below 1.64.

## 9. What a referee should attack first

1. The enclosure Haa+ of Psi_aa (an underestimate would break S2 silently): check the piece tests and the support
   cut b >= P/(P-1) (ball) / b >= max Pq/(P-r) (spheres). Random sampling passed (11880 and 3160 points), which is
   not a proof.
2. The convex-chord step: G(s) = max_b [H(am,b) + s H_a(am,b)] must be taken over a b-set that contains [0,a] for
   every a in I (here [0,A1]); the Taylor remainder must hold uniformly in b over that set.
3. The triangle certificate on [1/1024, 7/8] (ball): on a float grid the margin is smallest at (a,u) = (7/8, 1),
   where q_B(7/8,1) = 420633/20971520 = 0.0200573 > 0 (region A: D/(a^6(1-u)^2) >= 0.0889, at (7/8, 0)); also the
   boundary between the PA and T1 pieces at u = 1/(1+a).
4. The tails: the ball estimate Psi/a^3 <= t/((1-t)^3(1+t)) on the overlap window; for spheres the window bound
   |T/a^3 - L_q| <= r t b^3 max(1, b/q) and the clustering of windows.
5. Formula (78) and (43): re-derive one configuration by hand (e.g. partial overlap with R2 <= 1).
6. Theorem 1 as stated: c_LO is a sup over rho with integer mass; the bound E_ind >= -(J_f + D) int rho^{4/3}
   needs f in F_{mu nu} only on [0,oo)^2 including the axes (g(0) = 0 holds).

## 10. Files (all in /tmp/claude-1000/-home-bas-praime/9968315b-c905-4858-b310-2bb45fd31cf2/scratchpad/lo_cert/)

* lo_certify.py — the single self-contained verification script (Part A ball, Part B 48 spheres with the exact
  measure embedded); output lo_certify_output.txt ends with "81 checks, 0 failed" (60 + 21; 404 s). Generated by make_combined.py
  from the two working scripts below (identical code, parameters fixed).
* lo_ball_certify.py (+ out_ball_final.txt) — Theorem A; argument: a-step ratio (2.5e-4 used).
* cert_spheres.py (+ out_K48.txt, out_K32.txt, out_K20.txt) — Theorem B for any exact sphere-mixture JSON;
  arguments: measure.json rho_hi rho_lo R (1e-3 3e-3 600 used).
* meas_K48.json, meas_K32.json, meas_K20.json, meas_K12.json — exact rational measures; s2_rationalize.py makes them
  from the float optima s2_opt_K*_M1.json of s2_opt.py (heuristic optimiser, not part of any proof).
* explore.py, explore2.py, sym_test.py, s2_explore.py — first exploration (float), kept for the record.
