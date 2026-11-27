# Flowed characteristic polynomial of a Haar unitary: cluster expansion, exact variance, scaling limit

Status: Claims 1, 2 and 3 are **proved** as stated below, after the sign, the symmetry factor and the measure
d-sigma are fixed (Sections 1, 3 and 4). The circle property for t >= 0 is **proved** twice: once by citing
Lee-Yang, and once by a self-contained argument through the dynamics of the zeros (Section 2). Every number
below comes from `check.py`, which prints "41 checks, 0 failed" in 24 s.

Labels: [proved], [computed], [heuristic], [classical, cited]. "Derived here" means that I wrote the
argument in this note. It does not mean "new": no literature search for priority was made (Section 8).

---

## 0. Normalisation

- U is Haar on U(N) with eigenvalues z_1..z_N.
- E(x) = prod_i (1 + x z_i) = sum_j e_j x^j.
- E_t(x) = sum_{j=0}^N e^{-t j(N-j)} e_j x^j.
- Q := e^t.
- The two factors of the multiplier are e^{-t j (N-j)} = e^{-tNj} * Q^{j^2}. Hence

      E_t(x) = F_t(e^{-tN} x),   F_t(y) := sum_{j=0}^N Q^{j^2} e_j y^j.                      (0.1)

- E_t(0) = 1, and the x^N coefficient is e_N = det U, which is not 0. So E_t(x) = prod_{i=1}^N (1 + zeta_i(t) x)
  with zeta_i(t) = -1/(root_i). Then

      log E_t(x) = sum_{k>=1} (-1)^{k-1} p_k(t) x^k / k,    p_k(t) = sum_i zeta_i(t)^k.

- Define q_k by log F_t(y) = sum_k (-1)^{k-1} q_k y^k / k. By (0.1), **p_k(t) = e^{-tNk} q_k** [proved]. So
  q_k is the k-th power sum of the negated reciprocal roots of F_t. Those numbers have modulus e^{tN} when
  t >= 0 (Section 2): this is the "radius factor".
- For a partition lambda: z_lambda = prod_i i^{m_i} m_i!, r = l(lambda) is the number of parts,
  p_lambda = prod_s p_{i_s}, and eps_lambda = (-1)^{|lambda| - r}.

## 1. Claim 1: the cluster (Mayer) expansion [proved]

**Theorem 1.** Take any N, any complex t and any complex numbers z_i (Haar measure and |z_i| = 1 are not
used). For every k >= 1, including k > N:

    p_k(t) = e^{-tNk} sum_{lambda |- k} c_lambda(t) p_lambda,

    c_lambda(t) = (-1)^{r-1} (k / z_lambda) prod_s e^{t i_s^2} Phi_lambda(t),

    Phi_lambda(t) = sum_{G connected graph on {1..r}} prod_{(s,s') in E(G)} (e^{2 t i_s i_s'} - 1).

Equivalently, the coefficient of prod_s (p_{i_s}/i_s) is (-1)^{r-1} (k / prod_i m_i!) prod_s e^{t i_s^2} Phi_lambda.

- The **sign** is (-1)^{l(lambda)-1}.
- The **symmetry factor** is 1/prod_i m_i(lambda)!.
- Repeated parts are distinct labelled vertices with equal weights. Phi depends only on the multiset of
  parts, because relabelling the vertices permutes the connected graphs.
- The c_lambda are elements of the free algebra Lambda = C[p_1, p_2, ...]. There they are unique and do not
  depend on N. In N variables with k > N the p_lambda are linearly dependent, so only the identity holds;
  the coefficients are no longer unique.

**Proof.**

1. *Universal series.* Let R = C[p_1, p_2, ...] with the p_i independent indeterminates. Put
   x_i = (-1)^{i-1} p_i / i and e_j^u := [y^j] exp(sum_i x_i y^i) = sum_{lambda |- j} eps_lambda p_lambda / z_lambda.
   (Expand the exponential: the term prod_i x_i^{m_i}/m_i! is eps_lambda p_lambda / z_lambda.)
   Let G_t(y) = sum_{j>=0} Q^{j^2} e_j^u y^j, an element of 1 + yR[[y]].
   The evaluation map ev: p_i -> sum_l z_l^i is a ring homomorphism. Applied coefficientwise in y, it
   commutes with exp and log, because the y^k coefficient of exp or log is a universal polynomial in
   finitely many coefficients. The identity prod_l (1 + z_l y) = exp(sum_i (-1)^{i-1} p_i(z) y^i / i) in C[[y]]
   gives ev(e_j^u) = e_j(z), which is 0 for j > N. So ev(G_t) = F_t and log F_t = ev(log G_t).

2. *Sequences.* A multiplicity vector m with sum i m_i = j is reached by r!/prod m_i! sequences
   (i_1..i_r) of positive integers, r = sum m_i. Hence

       G_t = sum_{r>=0} (1/r!) sum_{(i_1..i_r)} Q^{(i_1+..+i_r)^2} prod_s x_{i_s} y^{i_s}.

   Write Q^{(sum i_s)^2} = prod_s Q^{i_s^2} prod_{s<s'} (1 + f_{ss'}), where f_{ss'} = Q^{2 i_s i_s'} - 1. Expand
   the product over pairs as a sum over all graphs G on [r]. Put w(i) = Q^{i^2} x_i y^i. Then

       G_t = sum_r (1/r!) g_r,   g_r = sum_{i in Z_{>0}^r} prod_s w(i_s) sum_{G on [r]} prod_{e in G} f_e.

   Each y-coefficient is a finite sum, because r <= degree.

3. *Exponential formula.* A graph on [r] is the same as a set partition pi of [r] (its components) together
   with a connected graph on each block. The vertex weights and the edge weights f_{ss'} (which depend only on
   i_s and i_s') factor over blocks. So g_r = sum_{pi in Pi_r} prod_{B in pi} c_{|B|}, where

       c_n = sum_{i in Z_{>0}^n} prod_s w(i_s) Phi(i),   Phi(i) = sum_{G connected on [n]} prod_e f_e.

   Then sum_r g_r / r! = exp(sum_n c_n / n!). The proof: the number of set partitions of [r] into blocks of
   given ordered sizes (r_1..r_m) is (1/m!) r!/prod r_l!, so

       sum_r g_r/r! = sum_m (1/m!) sum_{r_1..r_m} prod_l c_{r_l}/r_l! = exp(C),   C = sum_n c_n/n!.

   This is the classical exponential formula (Stanley, EC2, Cor. 5.1.6). The steps are legitimate in R[[y]]
   because c_n has y-order >= n. Hence log G_t = C.

4. *Coefficient.* In [y^k] C, group the sequences by their multiset lambda: there are r!/prod m_i! of them,
   and Phi is symmetric. This gives

       [y^k] C = sum_{lambda |- k} (1/prod m_i!) prod_s (Q^{i_s^2} (-1)^{i_s-1} p_{i_s}/i_s) Phi_lambda.

   Multiply by (-1)^{k-1} k. The sign is (-1)^{k-1} (-1)^{k-r} = (-1)^{r-1}, and (prod_s i_s) prod m_i! = z_lambda.
   Then apply step 1 and (0.1). QED.

**Second routes, checks and special cases** [proved]:

- t = 0: f = 0, so only the one-vertex graph survives and c_lambda = delta_{lambda,(k)}.
- *c_(k) = e^{tk^2} exactly, by a route without graphs.* Reduce modulo the ideal generated by p_i (i != k)
  and p_k^2. Then e_j^u = 0 except e_0 = 1 and e_k = (-1)^{k-1} p_k / k. So log G_t = (-1)^{k-1} Q^{k^2} p_k y^k/k
  modulo the ideal, and q_k = Q^{k^2} p_k modulo the ideal.
- *First order, by a route without graphs.* d/dt G_t at t = 0 is sum_j j^2 e_j y^j = D^2 E with D = y d/dy.
  Also D^2E/E = D(L) + L^2, where L = D log E = sum (-1)^{k-1} p_k y^k. So

      q_k = (1 + t k^2) p_k - t k sum_{a=1}^{k-1} p_a p_{k-a} + O(t^2).

  This agrees with Theorem 1: c_(a,b) = -k e^{t(a^2+b^2)} (e^{2tab} - 1)/z_(a,b), and c_lambda = O(t^{r-1}).
- k = 2: q_2 = Q^4 p_2 - (Q^4 - Q^2) p_1^2.
- check.py [computed]:
  - Route A (log recursion in the p-basis, exact polynomials in Q) equals route B (literal enumeration of
    connected graphs) for all 63 partitions of k <= 8 with r <= 6.
  - For the 3 partitions with r = 7, 8, route A equals a set-partition recursion for the connected sum.
  - The counts of connected graphs on 1..6 vertices are 1, 1, 4, 38, 728, 26704.
  - At N = 6, for a random U and k = 1..10 (k > N included) at t = 0.03, -0.02, 0.2, the formula matches
    p_k(t) computed from E_t by Newton's identities to relative error 3.7e-12.

## 2. Where the roots of E_t lie

**Theorem 2.**
(a) For every t >= 0 and every U, all roots of E_t lie on |x| = 1. So |zeta_i(t)| = 1 and |p_k(t)| <= N.
The numbers -1/y_i, where y_i are the roots of F_t, have modulus e^{tN}.
(b) For every t < 0 and N >= 2 the statement fails for some U, for example U = I.

*Proof of (a), route 1* [classical, cited]. E_t(x) = sum_{S subset of [N]} prod_{i in S, j not in S} q * prod_{i in S} w_i,
with q = e^{-t} in (0,1] and w_i = x z_i. This is the Lee-Yang polynomial of the Curie-Weiss ferromagnet. By
the Lee-Yang circle theorem in its multivariate form, with real symmetric couplings in [-1,1], it does not
vanish when all |w_i| < 1. Sources: Lee-Yang, Phys. Rev. 87 (1952) 410, Appendix; Asano 1970; Ruelle 1971.
All three are cited from memory and were not re-read. The complementation S -> [N] \ S gives
P(w) = prod w_i P(1/w), so P also does not vanish when all |w_i| > 1. Here |w_i| = |x|, so every root has |x| = 1.

*Proof of (a), route 2* [proved; derived here].
- Assume first that the eigenvalues are distinct. With e_N = e^{i phi}, the identity e_{N-j} = e_N conj(e_j)
  holds on the torus, and the multiplier is symmetric under j <-> N - j. Hence
  f_t(theta) := e^{-i phi/2} e^{-iN theta/2} E_t(e^{i theta}) is real for real theta.
- Write m = j - N/2. Then -j(N-j) = m^2 - N^2/4, so d/dt f = -f'' - (N^2/4) f.
- At a simple zero theta_i: d theta_i/dt = -(d f/dt)/f' = f''/f'. Since f = C prod_j sin((theta - theta_j)/2),

      d theta_i / dt = sum_{j != i} cot((theta_i - theta_j)/2) = d W / d theta_i,
      W = sum_{i<j} log |e^{i theta_i} - e^{i theta_j}|^2.

- Let S be the set of t >= 0 at which E_t has N distinct roots on the circle. Then 0 is in S.
  - *S is open.* f_t is real and changes sign at each simple zero, so for nearby t each of N disjoint small
    intervals still holds a zero (intermediate value theorem). The degree is N, so there are no others.
  - *Let T = sup{T': [0,T'] in S}, and suppose T < infinity.* On [0,T) the zeros satisfy the ODE, so
    dW/dt = sum (dtheta_i/dt)^2 >= 0 and W >= W(0) > -infinity. The leading coefficient e_N does not depend
    on t, so the roots at T are limits of unimodular roots and hence unimodular. A double root at T would
    force W -> -infinity, which is impossible. So T is in S, and by openness S contains [0, T + eps].
    This contradicts the definition of T.
  - Hence S = [0, infinity).
- For U with repeated eigenvalues, the roots depend continuously on U. QED.

Byproduct: the zeros perform gradient ascent of log|Vandermonde|^2. They repel each other and tend to equal
spacing; as t -> infinity, E_t -> 1 + e_N x^N.

*Proof of (b).* For U = I: sum_i zeta_i = e_1(E_t) = N e^{|t|(N-1)} > N. If all |zeta_i| = 1, the sum would
have modulus <= N.

check.py [computed]:
- For 2000 Haar U at N = 8: max ||x| - 1| = 1.7e-14, 5.1e-15 and 3.9e-15 at t = 0.001, 0.03 and 0.3.
- At t = -0.03: U = I has a root off the circle by 1.816, and 1548/2000 Haar U have a root off the circle.
  Openness also holds for t slightly negative, so each U with distinct eigenvalues keeps its roots on the
  circle for t > -eps(U). For t < 0 the statement therefore fails for "some U", not for "all U".
- The ODE checked by finite differences: relative error 1.3e-10.

## 3. Claim 2: exact variance [proved]

**Theorem 3.** For real t and 1 <= k <= N:

    E |p_k(t)|^2 = e^{-2tNk} H_k(t),
    H_k(t) := sum_{lambda |- k} c_lambda(t)^2 z_lambda = k^2 sum_{lambda |- k} e^{2t sum i_s^2} Phi_lambda(t)^2 / z_lambda.

For complex t, replace c^2 by |c|^2 and e^{-2tNk} by e^{-2 Re(t) N k}. For k <= N the right-hand side depends
on N only through e^{-2tNk}.

*Proof.* By Theorem 1, E|q_k|^2 = sum_{lambda,mu} c_lambda conj(c_mu) E[p_lambda conj(p_mu)]. Frobenius gives
p_lambda = sum_{nu |- k} chi^nu(lambda) s_nu, and s_nu(z_1..z_N) = 0 exactly when l(nu) > N (Macdonald I.7,
cited). For l(nu) <= N the s_nu are the characters of distinct irreducible representations of U(N), so
E s_nu conj(s_nu') = delta by Schur orthogonality. When k <= N every nu |- k has l(nu) <= N, so
E p_lambda conj(p_mu) = sum_nu chi^nu(lambda) chi^nu(mu) = delta_{lambda mu} z_lambda (column orthogonality of S_k). QED.

This is Diaconis-Evans, Trans. AMS 353 (2001) 2615, **Theorem 2.1(a)**. I read it in the Berkeley technical
report version (No. 577), text extracted from the PDF:

> for n >= (sum_j j a_j) v (sum_j j b_j), E[prod_j (Tr M^j)^{a_j} conj(Tr M^j)^{b_j}] = delta_ab prod_j j^{a_j} a_j!

The paper adds: "Part (a) is Theorem 2 in [DS94] with a slightly different (and incorrect) condition on a, b."
I did not read Diaconis-Shahshahani, J. Appl. Probab. 31A (1994) 49, myself. Here sum j a_j = |lambda| = k <= N.

**Extension to all k** [proved, same proof]:

    E |p_k(t)|^2 = e^{-2tNk} sum_{nu |- k, l(nu) <= N} ( sum_lambda c_lambda chi^nu(lambda) )^2.

For k > N the truncation l(nu) <= N enters, and the full Hall norm H_k is strictly larger than the expectation.

check.py [computed]:
- A route independent of c_lambda: q_k from Newton's identities in N variables, then exact Weyl integration
  (constant term of q(z) q(1/z) prod_{i != j}(1 - z_i/z_j)/N!), in exact rationals at Q = 11/10 and Q = 4/5.
  For N = 2, 3, 4 it equals H_k in all 18 cases with k <= N. It equals the truncated sum in all 30 cases.
  For every k > N it differs from H_k: for example N = 4, k = 6, Q = 11/10 gives 428.570651 against
  14529.767197.
- Small-t expansion (exact Taylor coefficients for k <= 8, and proved for all k by the two graph-free routes
  of Section 1):

      H_k = k e^{2tk^2} + t^2 k^2 (k^3 - k)/3 + O(t^3).

  Derivation: H_k = k c_(k)^2 + sum_{lambda != (k)} c_lambda^2 z_lambda, with c_(k) = e^{tk^2}. For lambda != (k),
  c_lambda(0) = 0, so the second sum is t^2 k^2 || sum_a p_a p_{k-a} ||^2 + O(t^3). In the Hall norm,
  <p_a p_{k-a}, p_b p_{k-b}> is nonzero only for b in {a, k-a}, which gives
  || sum_a p_a p_{k-a} ||^2 = 2 sum_{a=1}^{k-1} a(k-a) = (k^3 - k)/3.

- A consequence that is not an identity, from Theorem 2 with N = k: H_k(t) <= k^2 e^{2tk^2} for t >= 0. The
  ratio is 0.897 at (t, k) = (0.3, 6) and 0.997 at (0.3, 12), so the bound is nearly attained. Checked exactly
  at Q = 11/10, 3/2 and 5 for k <= 6.
- Monte Carlo: Section 7.

## 4. Claim 3: the scaling limit [proved]

**Theorem 4.**
- Hypotheses: tau is real (either sign), k = k_N with k_N <= N and k_N/N -> u in (0, 1], and t = tau/N^2.
- The measure: d sigma = da_1 ... da_{r-1}, Lebesgue measure on the projection of
  Delta_r(u) = {a in R_{>0}^r : sum a = u} to the first r-1 coordinates. It is **not** the surface measure,
  which is sqrt(r) times larger. For r = 1, the integral is evaluation at a_1 = u.
- The statement: the series K_tau(u) converges absolutely and E|p_k(t)|^2 / N -> K_tau(u), where

      K_tau(u) = u^2 e^{-2 tau u} sum_{r>=1} (1/r!) I_r,   I_r = int_{Delta_r(u)} h_r d sigma,
      h_r(a) = prod_s (e^{2 tau a_s^2}/a_s) Phi_r(a)^2,   Phi_r(a) = sum_{G conn on [r]} prod_{edges} (e^{2 tau a_s a_s'} - 1).

- K is entire in tau, and K_0(u) = u.

*Proof.*

4.1 *Sum over compositions.* Use 1/z_lambda = (prod_s 1/i_s)(1/prod m_i!), the fact that there are
r!/prod m_i! orderings, and the symmetry of Phi. Theorem 3 becomes

    E|p_k(t)|^2 = e^{-2tNk} k^2 sum_{r=1}^{k} (1/r!) sum_{(i_1..i_r) composition of k} prod_s (e^{2t i_s^2}/i_s) Phi(i; t)^2.

4.2 *Exact rescaling.* Put a = i/N and u_N = k/N. Then t i_s^2 = tau a_s^2, 2t i_s i_s' = 2 tau a_s a_s',
e^{-2tNk} = e^{-2 tau u_N}, and prod 1/i_s = N^{-r} prod 1/a_s. So, exactly,

    E|p_k|^2 / N = u_N^2 e^{-2 tau u_N} sum_{r=1}^{k} (1/r!) S_r(N),   S_r(N) = N^{1-r} sum_{compositions} h_r(i/N).

4.3 *h_r is entire for r >= 2.* Write e^{2 tau x} - 1 = x g(x) with g entire. A connected graph on r >= 2
vertices has every degree >= 1, so Phi_r = (prod_s a_s) Psi_r with Psi_r = sum_G prod_s a_s^{deg(s)-1} prod_e g(.)
entire. Hence h_r = prod_s (a_s e^{2 tau a_s^2}) Psi_r^2 is entire. The apparent singularity prod 1/a_s cancels.

4.4 *Fixed r.*
- For r = 1, S_1 = h_1(u_N) -> e^{2 tau u^2}/u.
- For r >= 2, attach to the point b = (i_1..i_{r-1})/N the cube prod_s ((i_s - 1)/N, i_s/N]. These cubes are
  disjoint. They lie in D_v = {b > 0, |b| <= v} with v = u_N, and they cover D_{v - (r-1)/N}. This follows from
  i_s = ceil(N b_s), because |i| < N|b| + r - 1 <= k.
- So S_r(N) is a Riemann sum of H(b) = h_r(b, v - |b|) over a region whose volume differs from that of D_v by
  O(r/N), with mesh sqrt(r)/N.
- By the uniform continuity of h_r on [0,1]^r, S_r(N) - int_{D_{u_N}} H -> 0, and int_{D_{u_N}} -> int_{D_u} = I_r.

4.5 *Domination uniform in N (the interchange of lim_N with sum_r).*
- Take tau >= 0, so all f >= 0.
- (i) Every connected graph contains a spanning tree T, so Phi <= sum_T prod_{e in T} f_e prod_{e not in T} (1 + f_e).
- (ii) prod_{all pairs} (1 + f) = e^{t(k^2 - sum i^2)} <= e^{tk^2}.
- (iii) f <= 2t i i' e^{2t i i'}.
- (iv) The weighted Cayley formula sum_T prod_s x_s^{deg_T(s)} = prod x_s (sum x_s)^{r-2} (Pruefer code; classical).
- Together:

      Phi(i; t) <= e^{2tk^2} (2t)^{r-1} k^{r-2} prod_s i_s.

- The compositions satisfy sum prod_s i_s = C(k+r-1, 2r-1) <= (2k)^{2r-1}/(2r-1)!. This is [x^k] of
  (x/(1-x)^2)^r.
- Collect the powers of N. The r-th term of 4.2 is at most

      M_r = e^{6 tau} 2^{2r-1} (2 tau)^{2r-2} / (r! (2r-1)!),

  independently of N, and sum_r M_r < infinity.
- For tau < 0, use |e^w - 1| <= e^{|w|} - 1, e^{2t i^2} <= 1 and e^{-2tNk} = e^{2|tau| u_N}. The same bound
  holds with |tau| and an extra factor e^{2|tau|}.
- Tannery's theorem (dominated convergence for series) then gives lim sum_r = sum_r lim. The same bound shows
  that the series for K converges absolutely.
- For complex tau, the analytic function G_N(tau) := e^{-2tNk} sum c_lambda^2 z_lambda / N (no conjugate) obeys
  the same domination with |tau|. So G_N -> K locally uniformly, and K is entire (Weierstrass). QED.

4.6 *Properties of K* [proved].
- For tau >= 0 the terms T_r := u^2 e^{-2 tau u} I_r / r! satisfy, with B_r = (2 tau)^{2r-2} u^{4r-3} / (r! (2r-1)!),

      e^{-2 tau u} B_r <= T_r <= e^{6 tau u^2} B_r.

  This uses the tree bounds and int_{Delta_r} prod a_s d sigma = u^{2r-1}/(2r-1)!.
- Small-tau expansion, by two routes:

      K_tau(u) = u - 2 tau u^2 (1-u) + tau^2 [ 2u^3 (1-u)^2 + u^5/3 ] + O(tau^3).

  - Route (i), from the integral: r = 1 gives u e^{-2 tau u(1-u)}. r = 2 gives
    (u^2/2) int_0^u 4 tau^2 a (u-a) da = tau^2 u^5/3. r >= 3 is O(tau^4).
  - Route (ii), from the finite-N expansion of H_k in Section 3: rescale, then pass to the limit coefficientwise,
    which is allowed by local uniform convergence in tau. The 1/N^2 term -u^3/(3N^2) disappears.
- Speed of convergence [computed; heuristic reason]. The observed error is O(1/N^2): the fitted order is
  1.98-2.00 in all 9 cases. Reason: h_r vanishes on the faces of the simplex, so Euler-Maclaurin has no 1/N
  term. The step 4.4 proof gives only o(1).

check.py [computed]:
- Simplex quadrature with r <= 7; the r = 7 term is below 4.2e-8 of K.
- The finite-N values are exact via Theorem 3, using the p-basis log recursion. The double-precision H_k
  agrees with a 60-digit computation to 1e-12.

| tau | u | K_tau(u) | largest-N value | rel. error | Richardson in 1/N^2, rel. error |
|---|---|---|---|---|---|
| 0.5 | 0.50 | 0.39133216 | 0.39132863 (N=48) | 9.0e-6 | 2.0e-9 |
| 0.5 | 1.00 | 1.06959634 | 1.06944168 (N=24) | 1.4e-4 | 3.3e-7 |
| 1.0 | 0.75 | 0.55997191 | 0.55986944 (N=32) | 1.8e-4 | 4.9e-7 |
| 1.0 | 1.00 | 1.24033191 | 1.23961527 (N=24) | 5.8e-4 | 3.6e-6 |
| 2.0 | 0.50 | 0.19674125 | 0.19671280 (N=48) | 1.4e-4 | 1.9e-7 |
| 2.0 | 1.00 | 1.76699977 | 1.76292920 (N=24) | 2.3e-3 | 4.4e-5 |

Three more cases are in the script output.
- A failed prediction, recorded: I first required a raw error below 0.2% at the largest N. The case tau = 2,
  u = 1 gives 0.23%.
- A second recorded failure: the first run of the domination check "failed" with term/M_r = 1.9e16. This was
  double-precision rounding noise in the r-th terms at large r (around 1e-33, against true values down to
  1e-69). In 60 digits the maximum of term/M_r is 3.6e-2.

## 5. Extending the parameters: what enters

- **k > N:** the Schur truncation l(nu) <= N enters (Section 3, proved; checked exactly and by Monte Carlo).
  Claim 2 as stated is false there.
- **u > 1 in the limit:** not done. At tau = 0 the limit is 1 (Diaconis-Evans 2.1(b): E|Tr M^k|^2 = min(k, N)).
  For tau != 0 the truncated sum would need an asymptotic analysis of the characters with l(nu) <= N. **Open.**
- **t < 0:** Claims 1 to 3 are unchanged; Theorem 4 holds with |tau| in the bounds. The roots leave the circle
  (Section 2(b)): in the zero dynamics, time reversal makes the zeros attract each other and collide.
- **Complex t:** Claim 1 is unchanged. Claim 2 holds with |c|^2 and e^{-2 Re(t) N k}. This is not tested
  numerically.

## 6. Hypotheses: where each is used

| hypothesis | used in |
|---|---|
| Haar measure | Theorem 3 only (Schur orthogonality); hence Theorem 4 |
| unimodular z_i | Theorem 2 (self-inversive symmetry, Lee-Yang with \|w_i\| = \|x\|); Theorem 3 (conj p = p(1/z)) |
| k <= N | Theorem 3 (every nu \|- k has l(nu) <= N); Theorem 4 (k_N <= N) |
| t real | Theorem 3 (c real); otherwise the |c|^2 form |
| t >= 0 | Theorem 2(a) only; the bound H_k <= k^2 e^{2tk^2} |
| distinct eigenvalues | route 2 of Theorem 2 only (removed by continuity) |
| t = tau/N^2, k/N -> u | Theorem 4 (4.2 is an exact identity under this scaling) |

## 7. Monte Carlo (check.py, 200,000 Haar matrices per N, scipy.stats.unitary_group, fixed seed)

The formula is e^{-2tNk} H_k for k <= N, and the truncated form for k > N. "r=1" is k e^{-2tk(N-k)}, the
single-vertex term alone. The z-score is (MC - formula)/SE. The cells in one (N, t) block share samples.
Since p_1(t) = e^{-t(N-1)} p_1, the k = 1 rows repeat one statistic.

| N | t | k | MC | SE | formula | z | z vs r=1 |
|---|---|---|---|---|---|---|---|
| 6 | 0.01 | 4 | 3.43278 | 0.00700 | 3.43409 | -0.19 | 3.5 |
| 6 | 0.01 | 6 | 6.21458 | 0.01260 | 6.21939 | -0.38 | 17.0 |
| 6 | 0.03 | 3 | 1.78960 | 0.00401 | 1.78186 | 1.93 | 10.3 |
| 6 | 0.03 | 5 | 4.20991 | 0.00816 | 4.21002 | -0.01 | 62.0 |
| 6 | 0.03 | 6 | 7.55790 | 0.01443 | 7.56663 | -0.60 | 108.0 |
| 6 | 0.03 | 7 (k>N) | 7.38361 | 0.01415 | 7.38246 | 0.08 | full Hall norm 15.362: z = -564 |
| 6 | 0.03 | 8 (k>N) | 7.11821 | 0.01382 | 7.11193 | 0.45 | full Hall norm 35.232: z = -2035 |
| 6 | -0.01 | 8 (k>N) | 6.49068 | 0.01357 | 6.47277 | 1.32 | full Hall norm 6.860: z = -27 |
| 8 | 0.01 | 6 | 4.89681 | 0.01004 | 4.89235 | 0.44 | 17.6 |
| 8 | 0.03 | 4 | 1.61972 | 0.00361 | 1.62277 | -0.85 | 24.4 |
| 8 | 0.03 | 6 | 3.68657 | 0.00719 | 3.68307 | 0.49 | 106.5 |
| 8 | -0.01 | 6 | 8.01380 | 0.01745 | 8.00297 | 0.62 | 22.1 |

- All 42 cells (N = 6: k = 1..8; N = 8: k = 1..6; t = 0.01, 0.03 and -0.01, which the claim did not name)
  agree within 3 SE. The largest |z| is 1.93.
- The single-vertex term alone is rejected (|z| > 5) at t = 0.03 for all 4 <= k <= N, so the cluster
  terms are being seen.
- For k > N the full Hall norm is rejected with |z| >= 24.
- The full table is in the output of check.py.

## 8. Classical versus derived here; sources

- **Classical:**
  - exponential formula (Stanley EC2 5.1.6);
  - Frobenius formula, Schur orthogonality, column orthogonality (Macdonald, Ch. I);
  - moments of traces (Diaconis-Shahshahani 1994; Diaconis-Evans 2001 Thm 2.1, **read**);
  - Lee-Yang circle theorem (cited from memory, **not re-read**);
  - Cayley/Pruefer weighted tree count;
  - the tree bound on connected sums with f >= 0.
- **Derived here:**
  - the explicit cluster formula for p_k(t), with its sign and symmetry factor;
  - the zero dynamics d theta/dt = grad log|Delta|^2 and the second proof of the circle property built on it.
    The real-line analogue for the heat flow (de Bruijn; Csordas-Smith-Varga 1994) is known to me only from
    memory; I did not search for the circle version;
  - the truncated-Schur formula applied to the flow;
  - the uniform domination and the limit K_tau(u), with its two-sided term bounds and its small-tau expansion.
- **Not searched:** no priority search was made. One web search, for the Diaconis-Evans statement, plus one
  fetch of the technical report PDF, whose text was extracted locally.
