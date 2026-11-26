# The gap G1-G3 of global.md: second attempt (26 Nov 2026)

Notation as in global.md. B(s) := Lambda_P(1-s)/Lambda_P(s). F2(x) = -sin(2pi abs x)/(2pi) on [-1,1] (transform
G2 = sin^2(pi xi)/(pi^2(xi^2-1))), FG = F2 - f_1 (transform sinc^2(pi xi)/(xi^2-1), FG <= 0, FG(0) = int FG = -1).
m_j(F) := int x^j F(x) dx. For k >= 1: m_2k(F2) > 0 (int_0^1 x^2k sin 2pi x dx < 0) and m_2k(FG) < 0.

## Outcome

G1 and G2 are not closed. G3 is closed only in narrow subcases. What moved: one obstruction (PROVED), one sign
theorem for G3 (PROVED), and some LP numbers (COMPUTED, uncertified).

## 1. G2: additive arguments cannot close it (PROVED)

Example: zeta_P(s) = zeta(s)(1 + lambda L^-s) with L > 1 and 0 < lambda < sqrt(L). It has every property used by
Lemma G, Prop. T and their extensions:
- positive coefficients; nu_0 = 1, a_0 = 1, all other frequencies (n and Ln) >= 1, no mass in (-1,1);
- (s-1) zeta_P entire of order 1, and zeta_P(-2k) = 0, so Xi_P is entire;
- infinitely many off-line zeros (a vertical row at Re s = log lambda/log L < 1/2);
- B = (1 + lambda L^(s-1))/(1 + lambda L^-s) has exponential type h = log L > 0 on the right;
- F tau_P = (1 + lambda/L) delta_0 - (1 + lambda) = r delta_0 + 2 zeta_P(0) on (-1/L, 1/L), exactly Prop. T on
  (-e^-h, e^-h) (direct Poisson computation).

It is excluded only by two things.
- Weil positivity fails: on the line, Re Xi_P'/Xi_P = -log L Re[z/(1+z)] with z = lambda L^(-1/2-it) and abs(z) < 1.
  This changes sign, and its mean is 0.
- It has no positive Euler product: log(1 + lambda L^-s) has alternating coefficients.

So any proof of G2 must use Weil positivity together with the Euler product. Positivity of tau_P plus a spectral
identity on a shrunken interval is not enough. The same holds at a = e^-h for every a < 1: a positive tau with an atom
at 1 can have F tau in span(delta_0, 1) on (-a, a) without being the comb. Example: (Sha - delta_0) +
lambda(Sha_L - delta_0) with 1 < L < 1/a. Lemma G is sharp at a = 1.

On your idea (1b), the claim that the kernel is "delta_1 minus positive" (PROVED by computation):
- A single off-line pair (rho, conj rho) has Blaschke kernel delta_1 - 2 delta y^(-beta-1)[2Re((1 + i delta/gamma)
  y^(-i gamma))] dy on (0,1). This oscillates in sign.
- "delta_1 minus a positive part" holds for real off-line zeros.
- It also holds for one conjugation-symmetric vertical row. There Bl = c_0 - sum_{m>=1} c_m lambda^(-m(s-1/2)) with
  c_m > 0, and the kernel is atomic.

Positivity side (PROVED, Xi_P entire): 2 pi f(t) = -B'/B(1/2+it) = -(d/dt) arg B(1/2+it). Feasibility therefore means
arg B is non-increasing on the line, and the mean of f equals lim (N(T) - N_on(T))/T.

In every example above (Dirichlet-polynomial twists, and removing a prime: rows of density exactly h/2pi), arg B is
bounded. Then f is almost periodic with mean 0, so f = 0 and Thm B applies. Conversely (HEURISTIC, needs the
Littlewood count N_P = N + hT/2pi + o(T)), a G2 system needs off-line zeros of density strictly above h/2pi, or a
compensation that is not almost periodic.

## 2. G3: missing trivial zeros (PROVED in the stated class)

Class: P feasible; (s-1) zeta_P entire of finite order; finitely many off-line zeros; zeta_P(-2k) != 0 exactly for k in
a finite set S.

Structure (proof as in Cor. F):
- B = eps e^(-b(s-1/2)) R(s) R_S(s), with R_S = prod_S (s+2k)/(s-1-2k) and R the finite Blaschke product of the
  off-line zeros.
- 2 pi f = b + 2 sum_j P_j(t) - 2 sum_S P_(1/2+2k)(t), so b >= 0.
- Pushing the contour of Prop. T past the poles of B at 1+2k gives, for even F in C_c^oo((-e^b, e^b)):

    <tau_P, F^> = r F(0) + 2 zeta_P(0) F^(0) - sum_{k in S} rho_k zeta_P(1+2k) m_2k(F),
    rho_k = Res_{1+2k} B = -2 pi^k (-1)^k zeta_P(-2k)/(k! Lambda_P(1+2k)).

Theorem G3. In this class:

    sum_S (-1)^k zeta_P(-2k) zeta_P(1+2k) m_2k(F2) pi^k/(k! Lambda_P(1+2k)) >= 0   (pairing with F2).

For S = {k}:
- (-1)^k zeta_P(-2k) > 0;
- zeta_P(0) < -r/2, i.e. B(1) > 1 (pairing with FG).

If moreover there are no off-line zeros:
- b >= 4/(4k+1) (f(0) >= 0);
- B(1) = e^(-b/2)(2k+1)/(2k) > 1 forces b < 2 log(1 + 1/2k);
- this leaves a window [0.8, 0.811) for k = 1 and [0.444, 0.446) for k = 2. It is nonempty for every k
  (2/(4k+1) < log(1 + 1/2k)).

LP (g3_cert2.py; COMPUTED, float, triangle basis, F^ >= 0 enforced on [1,40]): looks for F supported in [-A, A] with
F^ >= 0 on abs(xi) >= 1, F^(0) = m_2k(F) = 0 and F(0) < 0. Such an F kills the case whatever the unknown constants.
- Found for S = {1} at A = 1.5 (F(0) = -0.97) and S = {2} at A = 1.5 (F(0) = -3.2).
- Not found at A <= 1.2.
- Since e^b >= 2.23 (k = 1) and >= 1.56 (k = 2), this would close the cases S = {1} and S = {2} without off-line zeros.
- Uncertified: the solutions violate F^ >= 0 by 1e-4 to 1e-3 between grid points and in the tail on [1, 400].
- For large k the window pushes A towards 1 + 1/k, where no certificate was found.

## 3. G1: nothing moved

Positivity gives a Herglotz function. Herglotz functions can have singular continuous boundary measures, hence no
continuation. Your split idea needs F(mu_c) = A - C_P - F(mu_pp) to vanish near 0. Nothing ties the pure point part of
mu_P to A near 0 when the spectrum is not discrete. I found no way to get continuation from the comb plus positivity.

## 4. Sharpest statement of what is missing

Close G2 in two steps, both open:
- (a) Show that a feasible Euler product cannot have arg B(1/2+it) decreasing at a positive average rate, i.e. f with
  positive mean, i.e. N_on(T) <= N(T) - eps T.
- (b) Show that f >= 0 with mean 0 forces f = 0. This is true when f is almost periodic, and unproved otherwise.

After that only infinite S (G3) and G1 would remain.

Files: gap.md, g3_cert2.py (output in the text above), g3_cert.py (first, unreliable run: non-nested grids).
