# Global uniqueness of feasible Beurling systems: proofs, a converse theorem, and the exact gap

Scratch folder: /tmp/claude-1000/prover_boot/ (checks: checks.py, check2.py, twist_check.py, lemmaG_check.py;
LP experiments: cont1.py, lp_ac*.py, lp_cert*.py; outputs *_out.txt).

## 0. Summary (labels)

| # | statement | label |
|---|---|---|
| R1 | Theorem 1: feasible <=> Xi_P'/Xi_P continues to Re s > 1/2 with Re >= 0; then Xi_P'/Xi_P(s) = int (s-1/2)/((s-1/2)^2+t^2) dmu_P(t) and int dmu_P/(1+t^2) < oo | PROVED (sec. 2) |
| R2 | Theorem 2: Lambda_P(1-s) = eps Lambda_P(s) (eps = +-1, growth condition) <=> F(tau_P + eps r delta_0) = eps (tau_P + eps r delta_0) | PROVED (sec. 3) |
| R3 | Lemma H (self-dual => primes) | PROVED (sec. 4), new short proof via Lemma G |
| R3' | your Lemma H' identity "<kappa', c sinc^2(pi c xi)> = 0 for 0 < c <= 1" | FALSE for c < 1: the value is r(1-c) for every system (for zeta: 1-c). The conclusion "kappa' >= 0 => primes" is PROVED (c = 1 suffices) |
| new | Lemma G (additive uncertainty): a positive measure tau on abs(x) >= 1 whose Fourier transform equals alpha delta_0 + beta on (-1,1) is alpha(Sha - delta_0) | PROVED (sec. 4) |
| new | odd functional equation Lambda_P(1-s) = -Lambda_P(s) is impossible | PROVED (Cor. 4.3) |
| R4 | for feasible P: P = primes <=> FE <=> mu_P pure point with integer masses | PROVED (sec. 5) |
| new | Theorem B: feasible + Xi_P meromorphic near the line + mu_P singular => primes | PROVED (sec. 5) |
| new | Prop. T (twisted summation formula) and Theorem E (converse theorem: (s-1)zeta_P entire of finite order and Lambda_P(1-s) = O(abs(s)^N abs(Lambda_P(s))) on a right half-plane => primes; no positivity of the Weil form used) | PROVED (sec. 6); formula checked numerically to 1e-16 on a non-self-dual example |
| new | Cor. F: feasible + Xi_P entire of finite order + finitely many zeros off the line => primes (includes "all zeros on the line") | PROVED (sec. 6) |
| R7 | corrected form: for feasible P with Xi_P entire of finite order, TFAE: P != primes; zeta_P has a nontrivial zero with Re < 1/2 (then infinitely many); Lambda_P(1-s)/Lambda_P(s) is not polynomially bounded on any Re s >= c' in (1,3); kappa' != 0 on (-1,1) | PROVED (sec. 7). Your direction "zeros left <=> faster growth on the left <=> kappa' below 1" holds in this form; the Carleman route is not needed |
| R5 | short-time form factor | HEURISTIC; needs multiplicity-free, separated systems (sec. 7) |
| R6 | finite-support relaxation numbers | not rechecked (your COMPUTED) |
| gap | (G1) feasible P without finite-order meromorphic continuation; (G2) Xi_P entire of finite order with infinitely many off-line zeros and Lambda_P(1-s)/Lambda_P(s) of positive exponential type; (G3) missing trivial zeros | OPEN (sec. 7) |

Integrality of the multiplicities is used nowhere. Every PROVED statement holds verbatim for real
multiplicities >= 0 and for the "Lambda-relaxation" (any positive von Mangoldt measure supported in u >= u0 > 0).
So, whatever carries global uniqueness, it is not integrality.

## 1. Setting and conventions

P: locally finite multiset of reals > 1. zeta_P(s) = prod (1-x^-s)^-1 = sum_{n>=0} a_n nu_n^-s, nu_0 = 1 < nu_1 <= ...,
a_0 = 1, a_n > 0. tau_P := sum_n a_n (delta_{nu_n} + delta_{-nu_n}): positive, even, no mass in (-1,1), mass 1 at +-1.
Gamma_R(s) = pi^(-s/2) Gamma(s/2), Lambda_P = Gamma_R zeta_P, Xi_P = (1/2) s (s-1) Lambda_P,
rho(s) := Gamma_R(s)/Gamma_R(1-s), r := residue of zeta_P at 1.

Two Fourier transforms:
* log-scale (Weil form): g^(t) = int g(u) e^{itu} du. D_P = A - C_P as in proofs.txt; feasible: D_P(f*f~) >= 0
  for f in C_c^oo. (K1) gives the spectral measure mu_P >= 0 with D_P(g) = int g^ dmu_P; mu_P is even, tempered.
* additive scale: (F phi)(xi) = int phi(x) e^{-2 pi i x xi} dx; F extends to tempered distributions by <F T, phi> = <T, F phi>.
  Triangle f_c(x) = (1 - abs(x)/c)_+, F f_c(xi) = c sinc^2(pi c xi), sinc(y) = sin(y)/y. Sha = sum_{m in Z} delta_m, F Sha = Sha.

Facts used from proofs.txt (refereed): Prop 5(i) psi_P(x) = O(x) for feasible P (uses C_P >= 0), so zeta_P,
zeta_P'/zeta_P converge absolutely on Re s > 1; Prop 5(ii); Lemma 1 (C_P determines P).
Also: zeta_P is absolutely convergent on Re s > 1 implies int dtau_P/(1+x^2) <= 2 zeta_P(2) < oo and tau_P tempered.

## 2. Theorem 1 (Herglotz form, R1)

Theorem 1. (a) If P is feasible, then int dmu_P/(1+t^2) < oo, and for Re s > 1
  Xi_P'/Xi_P(s) = H(s) := int (s-1/2)/((s-1/2)^2 + t^2) dmu_P(t)
  = (1/2) int [1/(s-1/2-it) + 1/(s-1/2+it)] dmu_P(t)                                  (2.1)
The right side is holomorphic on Re s > 1/2 and Re H(1/2+sigma+i tau) = pi (P_sigma * mu_P)(tau) >= 0,
P_sigma(x) = sigma/(pi (sigma^2+x^2)). So Xi_P continues to a zero-free holomorphic function on Re s > 1/2,
abs(Xi_P(sigma+it)) is non-decreasing in sigma > 1/2, and mu_P is the weak limit of (1/pi) Re H(1/2+sigma+it) dt.
(b) Conversely, if zeta_P'/zeta_P converges absolutely on Re s > 1 and Xi_P'/Xi_P has a holomorphic continuation H
to Re s > 1/2 with Re H >= 0, then P is feasible and its spectral measure is the Herglotz measure of H.

Lemma 2.1 (archimedean term). With a(t) = Re psi(1/4+it/2) - log pi and w > 0:
  (1/2 pi) int 2w/(w^2+t^2) a(t) dt = psi(1/4 + w/2) - log pi.
Proof. u(z) = Re psi(1/4 + z/2) - log pi is harmonic on Re z > -1/2, continuous on Re z >= 0, O(log(2+abs(z))).
Let v = u - (Poisson integral of u on iR); v is harmonic on Re z > 0, zero on iR, o(abs(z)); by reflection it is an
odd entire harmonic function of growth o(abs(z)), hence 0. At real z = w, psi(1/4+w/2) is real. []
(Checked: checks_out.txt (1), agreement to 1e-25 at w = 0.6, 1, 2.5.)

Proof of (a). Step 1. D_P is tempered (K1), so D_P(g) = int g^ dmu_P for all g in S (density of C_c^oo in S;
both sides continuous on S).
Step 2. Fix real s > 1, w = s - 1/2, g_s(u) = e^{-w abs(u)}, g_s^(t) = 2w/(w^2+t^2) > 0. Take k in C_c^oo even, >= 0,
int k = 1, and rho_eps := (k*k) scaled to support [-2eps, 2eps]; then rho_eps^ = k^(eps .)^2 is in [0,1] and -> 1.
g_{s,eps} := g_s * rho_eps is in S. With a cut-off chi(u/R), chi_R g_{s,eps} -> g_{s,eps} in S, and
D_P(chi_R g) = pole(chi_R g) + arch(chi_R g) - C_P(chi_R g). As R -> oo: pole and C_P terms converge by dominated
convergence (abs(g_{s,eps}(u)) <= C e^{-w abs(u)}, w > 1/2, and sum_{x,k} log x x^{-ks} < oo); arch is tempered.
So D_P(g_{s,eps}) = pole(g_{s,eps}) + arch(g_{s,eps}) - C_P(g_{s,eps}).
Step 3 (eps -> 0). pole(g_s * rho_eps) = c_eps pole(g_s) with c_eps = int rho_eps(v) cosh(v/2) dv -> 1 (cosh(./2) is an
eigenfunction of even convolutions). arch(g_{s,eps}) = (1/2pi) int g_s^ rho_eps^ a dt -> arch(g_s) (dominated:
g_s^ a = O(log t / t^2)). C_P(g_{s,eps}) -> C_P(g_s) (pointwise, dominated). Spectral side: rho_eps^ >= 0, so by Fatou
int g_s^ dmu_P <= lim D_P(g_{s,eps}) < oo; since g_s^ >= c/(1+t^2), int dmu_P/(1+t^2) < oo, and then dominated
convergence gives D_P(g_s) := lim D_P(g_{s,eps}) = int g_s^ dmu_P.
Step 4. pole(g_s) = 2 int_0^oo e^{-wu}(e^{u/2}+e^{-u/2}) du = 2/(s-1) + 2/s; C_P(g_s) = 2 sum log x x^{-ks} = -2 zeta_P'/zeta_P(s);
arch(g_s) = psi(s/2) - log pi (Lemma 2.1, 1/4 + w/2 = s/2). Sum: 2 Xi_P'/Xi_P(s). So Xi_P'/Xi_P = H on (1,oo); both sides
are holomorphic on Re s > 1 (zeta_P zero-free there), hence equal there. H converges locally uniformly on Re s > 1/2
because int dmu_P/(1+t^2) < oo. Re H: (1/2)[sigma/(sigma^2+(tau-t)^2) + sigma/(sigma^2+(tau+t)^2)] integrated against
the even mu_P equals pi (P_sigma * mu_P)(tau). Monotonicity: d/dsigma log abs(Xi_P) = Re H >= 0. []
(Check for zeta: checks_out.txt (4): xi'/xi(2) = 0.0690662 vs 400 zeros + smooth tail 0.0690649.)

Proof of (b). Herglotz representation in the half-plane (Nevanlinna; e.g. Garnett, Bounded Analytic Functions, I.3):
H(1/2+w) = b w + i beta + int [1/(w-it) - it/(1+t^2)] dnu(t), b >= 0, nu >= 0, int dnu/(1+t^2) < oo.
H(s-bar) = conj H(s) (zeta_P real on (1,oo)) makes nu even and beta = 0, and lets us symmetrise:
H(1/2+w) = b w + int w/(w^2+t^2) dnu. b = lim H(1/2+sigma)/sigma = 0, since Xi_P'/Xi_P(sigma) = O(log sigma).
For even real g in C_c^oo with support in [-a,a] put G(s) = int g(u) e^{(s-1/2)u} du (entire, rapidly decreasing on
vertical lines). For c > 1:
  (i) (1/2pi i) int_(c) 2 zeta_P'/zeta_P(s) G(s) ds = -2 sum_{x,k} log x x^{-k/2} g(k log x) = -C_P(g) (Fourier inversion,
      absolute convergence);
  (ii) (1/2pi i) int_(c) [2/s + 2/(s-1)] G(s) ds = int g(u) 2cosh(u/2) du (write 1/(s-1) = int_0^oo e^{-(s-1)u} du, Fubini);
  (iii) (1/2pi i) int_(c) [psi(s/2) - log pi] G(s) ds = arch(g) (shift to Re s = 1/2; Im psi(1/4+it/2) is odd, g^ even).
So D_P(g) = (1/2pi i) int_(c) 2 H(s) G(s) ds. Insert H = int (1/2)[1/(w-it') + 1/(w+it')] dnu(t'), w = c - 1/2 + it.
Exact identity: (1/2pi) int G(1/2+sigma+it) / (sigma + i(t-t')) dt = int_0^oo g(v) e^{it'v} dv (write the kernel as
int_0^oo e^{-(sigma+i(t-t'))v} dv and invert the Fourier transform; Fubini is justified since
int dnu/(1+t'^2) < oo and G decays). Summing the two halves: D_P(g) = int g^ dnu. Both sides only see the even part of g,
so this holds for all g in C_c^oo; then D_P(f*f~) = int abs(f^)^2 dnu >= 0. []

Corollary 1.2 (continuation across the line). Let P be feasible and let Xi_P extend meromorphically to an open set U
containing the line Re s = 1/2. Then Xi_P has no poles on the line, and
  mu_P = sum_gamma m_gamma delta_gamma + f(t) dt,  f(t) = (1/pi) Re Xi_P'/Xi_P(1/2+it) (continuous extension), f >= 0,
m_gamma = order of the zero of Xi_P at 1/2 + i gamma. (Near rho on the line H = m/(s-rho) + holomorphic and
Re m/(sigma + i(t-gamma)) -> pi m delta_gamma; a pole (m < 0) would give a negative atom.)

## 3. Theorem 2 (functional equation <=> self-duality, R2)

Theorem 2. Let zeta_P converge absolutely on Re s > 1, eps in {+1,-1}, sigma_eps := tau_P + eps r delta_0.
(a) If F sigma_eps = eps sigma_eps (as tempered distributions), then Lambda_P extends meromorphically to C,
Lambda_P(1-s) = eps Lambda_P(s), the only poles are simple at s = 1 (residue r) and s = 0 (residue -eps r),
Lambda_P(s) - r/(s-1) + eps r/s is entire and bounded in vertical strips, and zeta_P(-2k) = 0 for k >= 1.
(b) Conversely, if Lambda_P extends meromorphically to C with Lambda_P(1-s) = eps Lambda_P(s), its only poles are
simple poles at 0 and 1, and Lambda_P(sigma+it) -> 0 as abs(t) -> oo uniformly for sigma in compact sets
[this holds e.g. if s(s-1)Lambda_P(s) is entire of finite order: Phragmen-Lindelof applied to
s(s-1)Lambda_P(s) e^{-i pi s/4} in the upper half-strip, and symmetrically], then F sigma_eps = eps sigma_eps.

Proof. Theta_P(y) := sum_n a_n e^{-pi nu_n^2 y}; Lambda_P(s) = int_0^oo Theta_P(y) y^{s/2} dy/y for Re s > 1;
Theta_P(y) = O(e^{-pi y}) as y -> oo (nu_n >= 1). With G_y(x) = e^{-pi x^2 y}, F G_y = y^{-1/2} G_{1/y}, and
<sigma_eps, G_y> = 2 Theta_P(y) + eps r.
(a) Self-duality on G_y: 2Theta_P(y) + eps r = eps y^{-1/2} (2Theta_P(1/y) + eps r). Split the Mellin integral at y = 1:
Lambda_P(s) = int_1^oo Theta_P(y)(y^{s/2} + eps y^{(1-s)/2}) dy/y + r/(s-1) - eps r/s. Everything follows.
(b) Theta_P(y) = (1/2pi i) int_(c) Lambda_P(s) y^{-s/2} ds/2, c > 1. Shift to Re s = 1-c (the hypothesis kills the
horizontal segments; on the vertical lines Lambda_P decays exponentially by Gamma_R and the functional equation).
Residues: at s = 1, (r/2) y^{-1/2}; at s = 0, Res_0 Lambda_P = -eps r, giving -eps r/2. The remaining integral equals
eps y^{-1/2} Theta_P(1/y). So <sigma_eps - eps F sigma_eps, G_y> = 0 for all y > 0. T := sigma_eps - eps F sigma_eps is
even and tempered; y -> G_y is smooth into S, so T kills all derivatives in y at y = 1, i.e. x^{2k} e^{-pi x^2}, i.e. all
even Hermite functions h_k(x) = H_k(sqrt(2pi) x) e^{-pi x^2}; these span a dense subspace of the even part of S. T is even,
so T = 0. []

## 4. Lemma G and the additive lemmas (R3)

Lemma G. Let tau be a positive even measure with tau((-1,1)) = 0 and int dtau/(1+x^2) < oo. Suppose that on (-1,1)
F tau = alpha delta_0 + beta (alpha, beta real), i.e. <tau, F phi> = alpha phi(0) + beta int phi for even phi in
C_c^oo((-1,1)). Then supp tau is in Z\{0}, beta = -alpha and tau = alpha (Sha - delta_0).

Proof. Test function F2(x) := -sin(2 pi abs(x))/(2 pi) on [-1,1], 0 outside: continuous, even, F2(0) = 0, int F2 = 0, and
  (F F2)(xi) = G2(xi) := sin^2(pi xi) / (pi^2 (xi^2 - 1)),
which is >= 0 for abs(xi) >= 1 and vanishes there exactly on Z (G2(+-1) = 0), and G2(xi) = O(1/(1+xi^2)).
(Derivation: int_0^1 sin(2pi x) cos(2pi xi x) dx = sin^2(pi xi)/(pi(1-xi^2)); numerical check in lemmaG_check_out.txt.)
Approximate: F_{l,d}(x) := (F2 * psi_d)(x/l), psi_d an even bump of width d, l(1+d) < 1. Then F_{l,d} is in C_c^oo((-1,1)),
<tau, F F_{l,d}> = alpha F_{l,d}(0) + beta l int F2 = alpha F_{l,d}(0) -> alpha F2(0) = 0, while
F F_{l,d}(xi) = l G2(l xi) psi_d^(l xi) -> G2(xi) with a bound C/(1+xi^2); dominated convergence gives int G2 dtau = 0.
As G2 >= 0 on supp tau with zero set Z there, tau = sum_{m != 0} c_m delta_m, c_m >= 0. Then F tau = sum c_m e^{-2pi i m xi}
is 1-periodic; F tau - alpha Sha - beta is 1-periodic and vanishes on (-1,1), hence everywhere. So
tau = F^{-1}(alpha Sha + beta) = alpha Sha + beta delta_0; tau({0}) = 0 gives beta = -alpha. []

Remark. F2 = FG + f_1 with FG(x) = -[(1-abs(x)) + sin(2pi abs(x))/(2pi)] 1_{[-1,1]}, F FG = sinc^2(pi xi)/(xi^2-1) (>= 0 on
abs(xi) >= 1, value -1 at 0). Pairing with f_1 alone gives beta <= -alpha (your triangle), pairing with FG alone gives
beta >= -alpha; F2 does both at once. The LP that searched for such certificates (lp_cert2.py) returned exactly this FG
(value -1 with all g_k = 0).

Corollary 4.1 (Lemma H). If F(r delta_0 + tau_P) = r delta_0 + tau_P then zeta_P = zeta and P = primes.
Proof. On (-1,1): r + F tau_P = r delta_0, so Lemma G applies with alpha = r, beta = -r: tau_P = r(Sha - delta_0);
the mass at 1 is a_0 = 1 (nu_n > 1 for n >= 1), so r = 1 and zeta_P = zeta. Then -zeta_P'/zeta_P = -zeta'/zeta, the
generalized Dirichlet series sum Lambda_P(y) y^{-s} determines its coefficients, C_P = C_primes, and Lemma 1 gives
P = primes. [] (Your triangle-plus-periodicity proof is also correct; the mollification needed for the kink of f_1 is the
same as above, using Fatou with psi^ >= 0.)

Corollary 4.2 (corrected Lemma H'). For every Beurling system and every constant r, with
kappa' := F tau_P + r - r delta_0:  <kappa', c sinc^2(pi c xi)> = r(1-c) for 0 < c <= 1 (limit of mollified pairings;
<F tau_P, f_c> = 0 since f_c vanishes on abs(x) >= c). For zeta, kappa' = Sha - delta_0 and the value is 1 - c
(check2_out.txt; also sum_{n>=1} sin^2(n theta)/n^2 = theta(pi - theta)/2). So your "= 0 for 0 < c <= 1" holds only at c = 1.
If kappa' is a positive measure, c = 1 and Fatou give int sinc^2(pi xi) dkappa' = 0, so kappa' lives on Z\{0}; then on
(-1,1) F tau_P = r delta_0 - r, and Lemma G gives P = primes. []

Corollary 4.3 (no odd functional equation). If F(tau_P - r delta_0) = -(tau_P - r delta_0), then on (-1,1)
F tau_P = r delta_0 + r; Lemma G forces beta = -alpha, i.e. r = -r, r = 0 and tau_P = 0, contradicting a_0 = 1. []
(A second proof: pair with the +1 Hermite eigenfunction F = h_4 - 12.1 h_0, F(0) = -0.1, F > 0 on abs(x) >= 1,
checks_out.txt (3).)

## 5. R4: the equivalences and Theorem B

Theorem B. Let P be feasible. Assume Xi_P extends meromorphically to a neighbourhood of the line Re s = 1/2 and mu_P is
singular with respect to Lebesgue measure. Then P = primes.
Proof. By Cor. 1.2, f = 0 off the zero ordinates, i.e. Re H = 0 on the line. On a connected, conjugation- and
reflection-symmetric neighbourhood U of the line, K(s) := H(s) + H(1-s) is meromorphic, and on the line
1 - s = s-bar so K = 2 Re H = 0; hence K = 0 on U, d/ds log(Xi_P(s)/Xi_P(1-s)) = 0, and Xi_P(s) = c Xi_P(1-s) with c^2 = 1.
Define Xi_P(s) := c Xi_P(1-s) for Re s < 1/2: an entire function (no poles on the line by Cor. 1.2).
Growth: for 1/2 <= sigma <= c0, abs(Xi_P(sigma+it)) <= abs(Xi_P(c0+it)) (Thm 1(a), monotonicity), which decays like
abs(t)^A e^{-pi abs(t)/4}; by the functional equation the same holds for 1-c0 <= sigma <= 1/2. So Theorem 2(b) applies
(eps = c), and Cor. 4.1 (eps = +1) or Cor. 4.3 (eps = -1) finishes. []

Corollary 5.1. Let P be feasible. Then P = primes <=> Lambda_P(1-s) = Lambda_P(s) (with the growth of Thm 2(b)) <=>
mu_P is pure point with integer masses.
Proof. primes => (feasibility of the primes is RH) FE and mu_P = zeros. FE => primes: Thm 2 + Cor. 4.1.
Pure point with integer masses => primes: by Prop 5(iv) mu_P([T,T+1]) = O(log T), and masses >= 1, so the support is
locally finite; then H(1/2+w) = sum_gamma m_gamma w/(w^2+gamma^2) is meromorphic on C with integer residues m_gamma
at w = i gamma; Xi_P = exp(int H) is entire and mu_P is singular; Theorem B. []
Open inside R4: pure point with non-integer masses (then Xi_P has branch points on the line and the Mellin argument
fails) and singular continuous spectra (no continuation across the line).

## 6. The twisted summation formula and a converse theorem

Hypotheses.
(E1) (s-1) zeta_P(s) extends to an entire function of finite order.
(E2) For some c' in (1,3) and C, N: abs(Lambda_P(1-s)) <= C (1+abs(s))^N abs(Lambda_P(s)) for all Re s >= c'.
(Lambda_P(s) != 0 for Re s > 1 by the Euler product. (E2) forces Lambda_P to be holomorphic on Re s <= 1-c', so
zeta_P(-2k) = 0 for all k >= 1. Put B(s) := Lambda_P(1-s)/Lambda_P(s); (E2) says B is polynomially bounded on Re s >= c'.)

Proposition T. Under (E1), (E2), for every even F in C_c^oo((-1,1)):
  sum_{n>=0} a_n (F F)(nu_n) = (r/2) F(0) + zeta_P(0) (F F)(0),
i.e. F tau_P = r delta_0 + 2 zeta_P(0) on (-1,1).
Proof. M[phi](s) := int_0^oo phi(x) x^{s-1} dx. Tate's local functional equation (check with phi = e^{-pi x^2}):
M[F F](s) = rho(s) M[F](1-s), first on 0 < Re s < 1, then on Re s > 0 by continuation. For 1 < c < 3, Mellin inversion
and absolute convergence give
  S := sum a_n (F F)(nu_n) = (1/2pi i) int_(c) M[F](1-s) Lambda_P(s) / Gamma_R(1-s) ds.
On the strip 1-c' <= Re s <= c the integrand is meromorphic with poles only at s = 1 and s = 0: M[F](1-s) has poles at
s = 1, 3, 5, ... (only s = 1 in the strip), 1/Gamma_R(1-s) is entire with a simple zero at s = 1, Lambda_P has poles only at
0, 1 (and -2k, outside the strip). Residues: at s = 1, [r/(s-1)][-F(0)/(s-1)][-(s-1)/2] gives r F(0)/2; at s = 0,
Res Lambda_P = 2 zeta_P(0) and M[F](1) = (F F)(0)/2 give zeta_P(0)(F F)(0). Horizontal segments vanish: on Re s = c
zeta_P is bounded; on Re s = 1-c', zeta_P(1-s) = B(s) rho(s) zeta_P(s) is polynomially bounded by (E2); (E1) and
Phragmen-Lindelof give a polynomial bound inside; M[F](1-s) decays faster than any power of abs(t). After s -> 1-s the
remaining integral is
  I = (1/2pi i) int_(c') M[F](s) Lambda_P(1-s)/Gamma_R(s) ds = (1/2pi i) int_(c') M[F](s) B(s) zeta_P(s) ds.
If supp F is in [-a,a] with a < 1, integration by parts gives abs(M[F](sigma+it)) <= C_k a^sigma (1+abs(t))^{-k} for
sigma >= 1; with (E2) and abs(zeta_P) <= zeta_P(c') the integrand is holomorphic on Re s >= c' and the integral over
Re s = sigma is O(a^sigma) -> 0 as sigma -> oo. So I = 0. []
Check (twist_check_out.txt): zeta_P = zeta * (1 + eps/(s-rho1) + conj(eps)/(s - conj(rho1))), rho1 the first zeta zero, is
not self-dual (B = Phi(1-s)/Phi(s) != 1) and satisfies (E1), (E2); its Mellin measure is sum delta_n + D(x) dx. With
F = (1-x^2)^8: both sides agree to 5e-16 for three values of eps.

Theorem E (converse theorem). Let P be a Beurling system (real multiplicities >= 0 allowed) satisfying (E1) and (E2).
Then zeta_P = zeta and P = primes.
Proof. Prop. T and Lemma G (alpha = r, beta = 2 zeta_P(0); tau_P >= 0 lives on abs(x) >= 1; int dtau_P/(1+x^2) < oo)
give tau_P = r(Sha - delta_0) and zeta_P(0) = -r/2; a_0 = 1 gives r = 1, zeta_P = zeta; then P = primes as in Cor. 4.1. []
Remarks. (1) No positivity of the Weil form is used: this is a Hamburger-type theorem in which the functional equation is
replaced by the one-sided growth bound (E2), at the price of positive coefficients. (2) In the example above Lemma G
predicts that D must change sign; it does (min D = -0.09; lemmaG_check_out.txt also confirms int G2 D dx = 0 to 1e-14).
(3) Sharpness of the method: if only abs(B) <= C(1+abs(s))^N e^{h sigma} with h > 0, Prop. T holds on (-e^{-h}, e^{-h})
only, and the scaled certificate is negative on (1, e^h), where tau_P lives. (4) Lemma G alone gives an additive
characterisation without any continuation: P = primes <=> F tau_P is in span(delta_0, 1) on (-1,1) <=> kappa' = 0 on (-1,1).
(5) Prior art: no literature search was done for Theorem E or Lemma G; Lemma G is close to the equality case of the
one-dimensional Cohn-Elkies bound and may be known.

Corollary F. Let P be feasible, let Xi_P extend to an entire function of finite order (equivalently (E1) together with
zeta_P(-2k) = 0 for k >= 1), and let Xi_P have only finitely many zeros off the line Re s = 1/2. Then P = primes.
Proof. Feasibility puts all zeros in Re s <= 1/2 (Thm 1). Let rho_1..rho_m (a conjugation-closed multiset) be the
zeros with Re < 1/2. B(s) = Xi_P(1-s)/Xi_P(s) (= Lambda_P(1-s)/Lambda_P(s)): the on-line zeros cancel (1 - rho = conj rho
on the line, equal multiplicities), so B = R e^Q with R(s) = prod_j (s - 1 + rho_j)/(s - rho_j) and Q a polynomial
(Hadamard). abs(R) = 1 on the line (pair rho with conj rho) and abs(R) <= 1 on Re s >= 1/2 (each factor compares the
distances to a point and to its mirror image). Hence Re Q(1/2+it) = 0; with B(s-bar) = conj B(s), Q(1/2+w) =
i pi k + odd real polynomial in w. On the line B'/B = -2 Re H, so by Cor. 1.2
  f(t) = (1/2pi) [ 2 sum_j P_j(t) - p(t) ],  P_j(t) = (1/2 - beta_j)/((1/2-beta_j)^2 + (t-gamma_j)^2),  p(t) = Re Q'(1/2+it),
p an even real polynomial. f >= 0 and int f/(1+t^2) <= int dmu_P/(1+t^2) < oo (Thm 1) force p = const = p_0, and
p_0 <= 2 sum P_j(t) -> 0 gives p_0 <= 0. So B = +-e^{p_0 (s-1/2)} R(s), abs(B) <= 1 on Re s >= 1/2: (E2) with N = 0.
Theorem E. [] (With m = 0 this is "feasible + Xi_P entire of finite order + all zeros on the line => primes"; the
exponential twist e^{-b(s-1/2)}, b > 0, that positivity alone would allow is killed by Lemma G.)

## 7. The exact gap, and remarks on R5-R7

Corrected R7 (PROVED). For feasible P with Xi_P entire of finite order the following are equivalent:
(i) P != primes; (ii) zeta_P has a zero with real part < 1/2 other than -2, -4, ... (and then infinitely many);
(iii) for no c' in (1,3) is Lambda_P(1-s)/Lambda_P(s) polynomially bounded on Re s >= c';
(iv) kappa' is not 0 on (-1,1).
(i)<=>(iii): Theorem E (and B = 1 for zeta). (i)=>(ii) infinitely many: Cor. F. (ii)=>(i): the primes are feasible only
under RH. (i)<=>(iv): Lemma G (no hypotheses needed). Note the direction of growth: a non-prime feasible system of
this kind must be LARGER on the left (not dominated by the right side up to a polynomial). In particular, with finitely
many off-line zeros positivity gives abs(B) <= 1 (smaller on the left), and that case is impossible.

What remains open (the exact gap for "feasible => primes"):
(G1) feasible P whose zeta_P has no meromorphic continuation of finite order to C (e.g. a natural boundary on
Re s = 1/2; mu_P with a singular continuous part or non-integer atoms). Nothing proved here beyond Thm B and Cor. 5.1.
(G2) Xi_P entire of finite order with infinitely many off-line zeros and B of positive exponential type.
HEURISTIC (needs a Blaschke-type factorisation of B, not written out): B = +-B_0 e^{h(s-1/2)} with h > 0 and
f = (1/2pi)(2 sum_rho P_rho - h) >= 0, so the off-line zeros must keep sum_rho P_rho(t) >= h/2 at every t: a positive
density (>= h/2pi per unit height) of zeros of zeta_P at bounded distance left of the line, and by Prop. T (on the
smaller interval) the Beurling integers carry Fourier mass in e^{-h} <= abs(xi) < 1.
(G3) (s-1)zeta_P entire of finite order but zeta_P(-2k) != 0 for some k (Xi_P has poles off the line): Prop. T acquires
polynomial terms (residues at -2k give moments of F) and Lemma G does not apply as stated.
Equivalent additive form of the whole question: does feasibility force F tau_P to lie in span(delta_0, 1) on (-1,1)?

Route (a) of the brief (feasibility + integrality => kappa' >= 0): not achieved; nothing in the proofs uses
integrality, and Lemma G shows that the weaker "kappa' = 0 on (-1,1)" is already equivalent to P = primes.
Route (b): done in the forms Theorem E / Cor. F; the theta function enters through Theorem 2 and Prop. T.
Route (c) (counterexample): none found. Real-weight and continuous relaxations: every PROVED statement holds there,
so a relaxed counterexample would have to sit in (G1)-(G3). LP experiment (lp_ac.py, lp_ac2_out.txt; COMPUTED, float):
absolutely continuous von Mangoldt densities e^u(1+d(u)) on [log x0, log U] (d piecewise constant, >= -1), gap below x0,
exact e^u beyond U. The optimum of min over t in [0,60] of Re H(1/2+it) is negative in every run (-1.1e-3 to -3.7e-4 for
x0 = 1.1..1.7, U = 10, 30; -0.17 for x0 = 2.5), slowly rising towards 0 with U, and the optimiser concentrates its mass
in narrow peaks (largest near 2.9, 4.6, 7.7, 12-26). Since only t <= 60 is constrained (the full constraint fails beyond
t = 60 at the optimum), these values are upper bounds: no feasible system exists in those finite families. This is
consistent with sup = 0 attained only by atomic systems, and proves nothing about it.

R5 (HEURISTIC, not used). The mean-value step is Montgomery-Vaughan for sum_x m_x b(x) x^{-it}. Its diagonal is
sum m_x^2 abs(b(x))^2, so with multiplicities the ramp is weighted by m_x^2, not m_x; and the off-diagonal error is
O(1/(T delta)) with delta the minimal log-gap, which for Beurling systems can be arbitrarily small. The statement holds
as written for multiplicity-free systems with abs(x - y) >= c > 0 (then delta >= c e^{-L} and T >= e^{(1+eps)L} suffices).
The constant L^2/(6 pi^2) was not re-derived.

R6: your finite-support numbers were not recomputed. They concern finite support; the theorems here concern L = oo.

## 8. Numerical checks (all small, foreground)

checks_out.txt: (1) Lemma 2.1 to 1e-25; (3) Hermite +1 eigenfunction for Cor. 4.3; (4) Herglotz form for zeta at s = 2, 1.2.
check2_out.txt: Cor. 4.2 value 1 - c for zeta (c = 0.3, 0.77, 1).
twist_check_out.txt: Prop. T on a non-self-dual example, agreement 5e-16.
lemmaG_check_out.txt: F F2 = G2 at five points; int G2 D = 0 to 1e-14 for the example; D changes sign.
lp_ac2_out.txt, cont1.py: continuous relaxation experiments (Sec. 7).


## 9. Corrections after the referee report ([referee_report.md](referee_report.md); no failing step found)

1. Theorem E: add the hypothesis that zeta_P converges in some half-plane (Landau's theorem with a_n >= 0 and finite order then gives absolute convergence on Re s > 1).
2. Theorem 1, part (b): Fubini is applied to the symmetrised kernel w/(w^2 + t'^2), not to each half-kernel.
3. Remark after Lemma G: the two inequalities are swapped (triangle: alpha + beta >= 0; F2: -alpha - beta >= 0). Not used in the proof.
4. R7 (ii): zeros of Xi_P, counted with multiplicity.
5. The open cases: G3 reads "Xi_P meromorphic of finite order with poles" (missing trivial zeros, or poles of zeta_P in Re s < 1/2 away from -2k, shape zeta(s)/(1 - 2^(-s))). G2 reads "Xi_P entire of finite order and Lambda_P(1-s)/Lambda_P(s) not polynomially bounded on any half-plane Re s >= c'"; that this means positive exponential type and a positive density of left zeros is heuristic.
6. The LP experiment used a cap d <= 50 on the density and a finite t-grid.
