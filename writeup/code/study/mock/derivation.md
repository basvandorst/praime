# Cumulants of the low-zero linear statistic over all characters mod q: theorem, conjecture, finite-q structure

Scope: Hughes–Rudnick, "Linear statistics of low-lying zeros of L-functions" (arXiv math/0208230, Q. J. Math. 54
(2003)); companion "Mock-Gaussian behaviour for linear statistics of classical compact groups" (math/0206289,
J. Phys. A 36 (2003)). Both read from the arXiv TeX sources. All numbers below come from the scripts in this folder
(file names in brackets). Tags: [proved] (in the paper or here), [computed], [heuristic], [conjecture].

## 1. What the paper states

Family. q an odd prime; all q-2 non-principal characters mod q, even and odd together (no split by parity).
Zeros 1/2 + i gamma (GRH not assumed), scaled x = gamma log q/(2 pi). Test functions: f real, even,
f-hat compactly supported, f(x) << (1+|x|)^(-1-delta); W_f(chi) = sum_j f(x_{chi,j}). Fourier convention
f-hat(u) = int f(x) e^(-2 pi i x u) dx.

Explicit formula (their eq. 2.1, our check in sect. 2): W_f(chi) = Wbar(chi) + W^osc(chi),
W^osc(chi) = -(1/log q) sum_n Lambda(n) n^(-1/2) f-hat(log n/log q) (chi(n) + conj chi(n)).
Wbar depends on chi only through the parity a(chi) and equals int f + O(1/log q).
(Their printed gamma term "Gamma'/Gamma(1/2 + a + ir) + Gamma'/Gamma(1/2 + a - ir) - (1/2) log pi" is a misprint; the
correct term, checked against zeros to 1.5e-8, is log(q/pi) + Re psi(1/4 + a/2 + ir/2). Asymptotics unaffected.)

* Thm 3.1 [proved]: supp f-hat in [-2, 2] => <W_f>_q = int f + O(1/log q).
* Thm 4.1 [proved]: supp f-hat in [-1, 1] => variance -> sigma_f^2 = int_{-1}^{1} |u| f-hat(u)^2 du.
* Thm 5.1 [proved]: supp f-hat in [-alpha, alpha], m < 2/alpha => lim <(W^osc)^m>_q = Gaussian moment
  (0 for odd m, (m-1)!! sigma_f^m for even m). The statement is about the prime sum W^osc; for W_f itself the
  parity-dependent Wbar adds O(1/log q). The error terms are O(q^(m alpha/2 - 1 + eps)) (congruences) and O(1/log q)
  (non-matching equalities such as p * p = p^2).
* RMT side [proved]: Z_f(U) = sum_n (1/N) f-hat(n/N) Tr U^n on U(N): E Z_f = f-hat(0) for all supports; variance ->
  int min(|u|, 1) f-hat^2 for all supports; supp f-hat in [-2/m, 2/m] => first m centred moments Gaussian
  (Diaconis–Shahshahani); all moments for any support via Dyson's correlation functions (their Thm 7.4); U(N) mock-Gaussian = their Thm 6.5.
* Conjecture (abstract and sect. 6): the moments of W_f follow the U(N) moments beyond the proven range, so the
  higher moments are not Gaussian there and W_f is not asymptotically normal. Remark after Thm 4.1: the variance
  int min(1,|u|) f-hat^2 should hold for all admissible f. No explicit formula is written for L-functions beyond
  "the U(N) moments".
* Companion paper (math/0206289, Theorem 7, credited to Soshnikov): C_l^{U(N)}(g) = sum mu_l(k) prod g_{k_j}; C_2 = sum min(|n|,N)
  g_n g_{-n}; for l >= 3, mu_l(k) = 0 whenever sum |k_j| <= 2N.

### 1.1 The U(N), N -> infinity cumulants (derived here; [proved], checked two ways)

log E exp(t Z) = log det(I + (e^{tf} - 1) K_sine), K_sine-hat = 1_[-1/2, 1/2]. Expanding the log and the exponential:

  kappa_m^U(f) = int_{u_1+...+u_m = 0} prod_j f-hat(u_j) Q_m(u) du_1..du_{m-1},
  Q_m(u) = sum over ordered set partitions (B_1..B_n) of {1..m} of (-1)^(n+1)/n * (1 - R)_+,
  R = max - min of the partial sums S_0 = 0, S_i = w_1 + .. + w_i (i < n), w_i = sum_{l in B_i} u_l.

Check: Q_2(u, -u) = 1 - (1 - |u|)_+ = min(|u|, 1). Independent check: exact finite-N cumulants from Heine's identity
E e^{tZ} = det T_N(e^{tF}) (traces of Toeplitz products; [rmt.py], `toeplitz_cumulants`), N = 100..400: agree with the
sine-kernel quadrature to 1e-4 relative for kappa_3 (s = 1.5: 0.025154 vs 0.025152) and 1.2e-3 for kappa_4 (s = 1.5:
-0.105637 on the 3-D grid h = s/40 vs -0.105511 at N = 400; the gap is the grid error).

Consequence (from the Soshnikov support condition): kappa_m^U = 0 whenever floor(m/2) * s <= 1 (the positive frequencies
of a balanced m-tuple with |u_j| <= s sum to at most floor(m/2) s). So the U(N) threshold is s = 2/m for even m but
s = 2/(m-1) for odd m: **kappa_3^U = 0 up to s = 1**, beyond Hughes–Rudnick's proven 2/3. Numerically kappa_3^U = 0
to 1e-16 for s <= 1.

Test function: f-hat(u) = (1 - |u|/s)_+ (f = Fejer kernel s (sin(pi s x)/(pi s x))^2 >= 0, decay x^-2, admissible; the
standard choice; linear in u, so the binned prime sum needs one accumulator per s). Any admissible f-hat is continuous,
so f-hat -> 0 at the edge and the RMT excess starts at high order beyond each threshold (table below): the thresholds
themselves are invisible; readable cells lie well beyond them.

kappa_m^U for this f-hat [computed, N = 300 Toeplitz = sine kernel]:

| s | kappa_2^U | diag int|u|f^2 = s^2/6 | U_2 = -int(|u|-1)_+ f^2 | kappa_3^U | kappa_4^U | k4/k2^2 |
|---|---|---|---|---|---|---|
| 0.6 | 0.06000 | 0.06000 | 0 | 0 | -5.9e-8 | -0.0000 |
| 0.8 | 0.10667 | 0.10667 | 0 | 0 | -1.22e-4 | -0.011 |
| 1.0 | 0.16667 | 0.16667 | 0 | 0 | -2.976e-3 | -0.107 |
| 1.2 | 0.23981 | 0.24000 | -1.85e-4 | +1.567e-3 | -2.029e-2 | -0.353 |
| 1.4 | 0.32449 | 0.32667 | -2.177e-3 | +1.294e-2 | -6.783e-2 | -0.644 |
| 1.6 | 0.41823 | 0.42667 | -8.44e-3 | +4.296e-2 | -0.1526 | -0.873 |
| 2.0 | 0.62500 | 0.66667 | -4.167e-2 | +0.1833 | -0.4204 | -1.076 |

Finite N: the exact U(N) cumulants at N = 9, 12, 14 (N ~ log q) differ from N = infinity only through the sampling of
f-hat at n/N, by 1-10 %: kappa_4(s = 1) = -0.00268, -0.00281, -0.00285 vs -0.00298; kappa_4(1.4) = -0.0673 vs -0.0678
[reg/extrap_check_out.txt]. Nothing in the arithmetic corresponds to that discreteness (log p is continuous), and the
finite-q effects found below are 10-100 times larger, so no finite-N matching is used; the finite-q structure is
handled exactly instead (sect. 3-4).

## 2. The statistic from primes, all characters at once [proved; validated]

L = log q. g = least primitive root, l(n) = discrete log of n mod q, chi_k(g) = e(k/(q-1)), k = 1..q-2 (k = 0 is chi_0);
chi_k(-1) = (-1)^k, so even characters = even k ((q-3)/2 of them), odd = odd k ((q-1)/2).

  P_s(chi) := W_f(chi) - Wbar_{a(chi)} = -(2/L) sum_{n < q^s, q not| n} Lambda(n) n^(-1/2) f-hat_s(log n/L) Re chi(n)
  Wbar_a = (1/L) int f(x) [log(q/pi) + Re psi(1/4 + a/2 + i pi x/L)] dx     (a = 0 even, 1 odd).

n = 0 mod q: only n = q^k (Lambda = log q, chi = 0) -- dropped; prime powers p^k included with weight log p p^(-k/2)
f-hat(k log p/L). Binning: B_s[l] = sum_{n: l(n) = l} Lambda(n) n^(-1/2) f-hat_s(log n/L); then
S_s(k) = sum_l B_s[l] e(k l/(q-1)) is one FFT of length q-1 (numpy ifft times (q-1)), P_s(chi_k) = -(2/L) Re S_s(k).
P(chi) = P(conj chi): complex characters come in identical pairs, the quadratic character once; the family moments
weight each character once (as in the paper), the sampling model uses n_eff = (q-1)/2 pairs.

Validation:
* [validate_fft.py -> pilot/validate_fft_out.txt] q = 10007, 200 random characters, s = 0.5, 1, 1.5 (+25 at s = 2):
  direct sums (plain full sieve, count = PARI primepi; PARI znlog) agree with the FFT to <= 2.1e-14; family size
  q-2; conjugate symmetry; parity via chareval.
* [validate_zeros.py -> pilot/validate_zeros_out.txt] all 99 (q = 101) and 149 (q = 151) characters, f-hat = T*T
  (cubic B-spline shape, support s = 1, 2; f ~ x^-4), PARI lfunzeros on |gamma| <= 150: zero sum = Wbar_a + P to
  <= 1.5e-8 (tail bound 5e-9). The parity term is large at small q: at q = 101, s = 1: Wbar_even = 0.0831,
  Wbar_odd = 0.1887 while f-hat(0) = 1/3; mixing parities would add (Wbar_o - Wbar_e)^2/4 to kappa_2 etc.; using P
  (zeros minus the exact parity term) removes it.
* Speed: the segmented odd-only sieve does 1e7 numbers in ~0.013 s; q = 10007 with s up to 2 (all n < 1.0e8),
  19 supports: 0.6 s total [mock_sum.py].

## 3. Exact finite-q moments

Orthogonality on the family (A, B products of prime powers coprime to q):
  <chi(A) conj chi(B)>_F = [(q-1) 1{A = B mod q} - 1]/(q-2).
With P = -(1/L)(S + conj S), S = sum a_n chi(n), a_n = Lambda(n) n^(-1/2) f-hat(log n/L):
  <P^m> = (-1/L)^m sum_{S subset {1..m}} sum_{n_1..n_m} prod a_{n_j} [(q-1) 1{A_S = B_S mod q} - 1]/(q-2),
  A_S = prod_{j in S} n_j, B_S = prod_{j not in S} n_j.
Every tuple falls in exactly one class [proved, elementary]:
 (E) A_S = B_S as integers: weight exactly 1. By unique factorisation these are exactly the moments of the
     independent-phase (Steinhaus) model chi(p) = z_p iid uniform, chi(p^k) = z_p^k. Its cumulants are additive over
     primes: kappa_m^St = sum_{p != q} kappa_m(X_p), X_p = -(2/L) sum_k a_{p^k} cos(k theta), computed exactly per prime
     (theta-grid with M > 4K points). For p > q^(s/2): kappa_2 = 2a_p^2/L^2, kappa_3 = 0, kappa_4 = -6 a_p^4/L^4;
     kappa_3 comes only from p, p, p^2 (leading -(6/L^3) sum_p a_p^2 a_{p^2}). [mock_sum.py, mock_test.py:steinhaus]
 (N) 0 < |A_S - B_S| < q: never congruent, weight exactly -1/(q-2) ("no-wrap" deficit).
 (F) |A_S - B_S| >= q: weight [(q-1) 1{cong} - 1]/(q-2); zero on average if A_S mod q were uniform on units.
Model R := (E) + (N) exactly, (F) replaced by its uniform average 0. Then
  mu_m^R = mu_m^St + Def_m,  Def_m = (-1/L)^m (-1/(q-2)) sum_S N_S,  N_S = sum over (N)-tuples of prod a.
Types and multiplicities (ordered tuples): m = 1: (1|0) x2 [2 <= n <= q]; m = 2: (2|0) x2 [n1 n2 <= q], (1|1) x2;
m = 3: (3|0) x2, (2|1) x6; m = 4: (4|0) x2, (3|1) x8, (2|2) x6. All sums exact (window sums on sorted lists) except the
(2|2) sum for pair products A >= 50q (20q for q ~ 1e6), where the window sum is 2q times the local density of pair
products (log-histogram convolution); changing the cut from 10q to 200q moves R_4 by < 5e-7 [checked at q = 300007].
Cumulants: R_m := kappa_m(mu^R) - kappa_m^St (mu_1^R = Def_1 included). Structured part S_m := E_m - R_m, where
E_m = kappa_m - kappa_m^St is the measured excess over the diagonal. S_m is exactly the contribution of the
wrap-around tuples (F): the excess of congruent wrap-around tuples over their uniform share 1/(q-1).

Where each regime sits [proved for the bookkeeping; the limits are heuristic/conjectural as labelled]:
* If every balanced product is < q (m s/2 < 1, the paper's range) the balanced types have no (F) tuples; only
  unbalanced types (e.g. n1 n2 n3 = n4 mod q, n1..n4 = 1) wrap, and they are equidistributed (HR's Lemmas 5.2-5.3 bound
  them by q^(m s/2 - 1 + eps)). Then kappa_m = St_m + R_m + (fluctuation). R_m is NOT small at accessible q: it
  carries the principal-character subtraction -(sum a)^2/(q-1) that nothing cancels when products are below q
  (m = 2: -(2/L^2)(sum_{n<q} a_n)^2/(q-1) + ...). It decays like q^(m s/2 - 1) below and like 1/log q above the range.
* m = 2, s > 1 (F) = pairs n = m mod q, n != m: 2 sum_{h>=1} sum_n a_n a_{n+hq}. Hardy–Littlewood gives the average
  S(hq) ~ S(h); with sum_{h<=H} (S(h) - 1) = -(1/2) log H + O(1) (Montgomery–Soundararajan; Goldston–Montgomery) the
  structured part -> -int (|u|-1)_+ f-hat^2 = U_2 [heuristic derivation, reproduces Montgomery–Hooley]: this is the
  variance of Lambda in progressions to one prime modulus, Hooley's conjecture / Friedlander–Goldston (asymptotic
  for a single q under GRH + a uniform Hardy–Littlewood hypothesis), proved in F_q[t] by Keating–Rudnick. The first
  lower-order terms (constant in the singular-series sum, Euler–Maclaurin boundary of the h-sum, O(1) part of the
  pair integral, the n m = 1 mod q pairs, prime powers) are all O(1/log q); I did not complete that bookkeeping.
* m = 3, (F) of type (2|1): n1 n2 = n3 + j q -- the covariance over residues of the Lambda*Lambda sum (products up to
  q^(2s)) and the Lambda sum (up to q^s): a shifted-convolution (Hardy–Littlewood type) correlation of Lambda*Lambda
  with Lambda along the progression. RMT says its structured part -> 0 for s <= 1 and -> kappa_3^U > 0 beyond.
* m = 4, (F) of type (2|2): n1 n2 = n3 n4 mod q, |n1 n2 - n3 n4| >= q -- the variance of Lambda*Lambda in progressions
  mod q beyond the diagonal range, the analogue for Lambda*Lambda of the Keating–Rodgers–Roditty-Gershon–Rudnick
  matrix-integral conjecture for d_k in progressions (proved in F_q[t] as q -> infinity; Rodgers–Soundararajan proved
  an average-over-moduli version for d_k in a range). Plus (3|1) (Lambda*Lambda*Lambda against Lambda) and (4|0).
* Deficit asymptotics [heuristic, checked numerically]: for A, B >> q the (2|2) deficit is ~ -(12/L) int_{v>1}
  (f-hat*f-hat)(v)^2 dv, the (3|1) deficit ~ -(16/L) int_1^s f-hat (f-hat*f-hat*f-hat), plus O(1/L^2) from products < q.
  Computed (prediction side only, no characters): R_4(s = 1) = -0.0185, -0.0140, -0.0106 at q = 10007, 100003,
  1000003; R_3(s = 1) = 0.0114, 0.0064, 0.0038; R_2(1) = -0.0047, -0.0019, -0.0009 [explore_def_scaling.py ->
  pilot/deficit_scaling_out.txt; same values in reg/predictions_frozen.json].
  So the deficit -> 0, but only like ~1/log q; at s = 1.4, q ~ 1e4-1e5 it is as large as kappa_4^U itself.

## 4. The pilot, q = 10007 (development modulus; counts as used) [computed]

mock_sum.py: 0.6 s for 19 supports s = 0.2..2.0 (all n < 1.0e8). Reproduces the scout: at s = 0.4,
kappa_3/kappa_2^1.5 = -0.236, kappa_4/kappa_2^2 = -0.20. kappa_3 = -0.000723 is the Steinhaus value (-0.000750) plus
the deficit (+0.000027); kappa_4 = -0.000086 = Steinhaus (-0.000071: small primes, the bounded single-cosine kurtosis
-3/8 per prime) + deficit (-0.000015).

Readings (pilot/pilot_read_out.txt; se = delta method from moments to order 8, n_eff = 5003; a bootstrap over conjugate
pairs agrees within ~10 %, pilot/pilot_se_out.txt):

| s | E_2 | R_2 | S_2 | U_2 | se_2 | E_3 | R_3 | S_3 | U_3 | se_3 | E_4 | R_4 | S_4 | U_4 | se_4 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.4 | -0.000045 | -0.000045 | +0.000000 | 0 | 0.0004 | +0.000028 | +0.000027 | +0.000000 | 0 | 0.00009 | -0.000016 | -0.000015 | -0.000000 | 0 | 0.00002 |
| 0.6 | -0.00031 | -0.00030 | -0.00001 | 0 | 0.0009 | +0.00043 | +0.00044 | -0.00001 | 0 | 0.0003 | -0.00048 | -0.00053 | +0.00005 | -0.0000 | 0.00013 |
| 0.8 | -0.00098 | -0.00126 | +0.00028 | 0 | 0.0016 | +0.0032 | +0.0029 | +0.0003 | 0 | 0.0008 | -0.0044 | -0.0047 | +0.0003 | -0.0001 | 0.0004 |
| 1.0 | -0.0032 | -0.0047 | +0.0014 | 0 | 0.0022 | +0.0132 | +0.0114 | +0.0018 | 0 | 0.0016 | -0.0184 | -0.0185 | +0.0001 | -0.0030 | 0.0011 |
| 1.2 | -0.0120 | -0.0130 | +0.0010 | -0.0002 | 0.0030 | +0.0369 | +0.0319 | +0.0051 | +0.0016 | 0.0027 | -0.0437 | -0.0443 | +0.0005 | -0.0203 | 0.0026 |
| 1.4 | -0.0302 | -0.0250 | -0.0052 | -0.0022 | 0.0039 | +0.0788 | +0.0648 | +0.0140 | +0.0129 | 0.0040 | -0.0764 | -0.0827 | +0.0062 | -0.0678 | 0.0050 |
| 1.6 | -0.0587 | -0.0398 | -0.0188 | -0.0084 | 0.0051 | +0.1415 | +0.1097 | +0.0317 | +0.0430 | 0.0055 | -0.1068 | -0.1323 | +0.0256 | -0.1526 | 0.0092 |
| 1.8 | -0.0993 | -0.0570 | -0.0423 | -0.0211 | 0.0067 | +0.2255 | +0.1664 | +0.0591 | +0.0982 | 0.0071 | -0.1251 | -0.1917 | +0.0666 | -0.2728 | 0.0150 |
| 2.0 | -0.1517 | -0.0761 | -0.0755 | -0.0417 | 0.0085 | +0.3331 | +0.2343 | +0.0988 | +0.1833 | 0.0089 | -0.1196 | (n.c.) | | -0.4204 | 0.0222 |

Findings:
1. In the proven range (s < 2/m) the measured cumulants equal Steinhaus + model R to < 0.05 se (e.g. s = 0.5:
   E_4 = -0.000113, R_4 = -0.000111). The exact finite-q bookkeeping works.
2. Below the U(N) thresholds the measured excess over the diagonal is large and entirely the no-wrap deficit:
   E_4(s = 1.0) = -0.0184 = R_4 (-0.0185) while kappa_4^U = -0.0030; E_3(1.0) = +0.0132 with R_3 = +0.0114.
   Against the naive reading "excess over the diagonal = random-matrix excess" the pilot is off by factors 3-30.
3. Beyond the thresholds the structured (wrap-around) part, which must carry the RMT term asymptotically, is
   * m = 3: S_3 = U_3 within 1 se at s = 1.3-1.5 (S/U = 1.6, 1.08, 0.87), then falls below (s = 2: 0.54 U_3);
   * m = 2: S_2 has the RMT sign and 1.8-2.4 x its size for s = 1.4-2.0 (lower-order terms of the same sign);
   * m = 4: S_4 ~ 0 or positive (+0.006 at s = 1.4) where U_4 = -0.068: the RMT fourth cumulant is absent from the
     structured part at this q (15 se), while the total E_4 = R_4 + S_4 happens to be close to U_4 near s = 1.4-1.5.
   No single decomposition (total, structured, or an effective-N rescaling) matches RMT for m = 2, 3, 4 at once.
   The finite-q terms are O(1) relative to the RMT term at q ~ 1e4 and decay (where computable) like ~1/log q.

## 5. Test design, power, and what can be read

Statistical power is not the problem. Expected pooled se (4 fresh moduli per group; scaled from the pilot,
reg/power_out.txt) and the z of the N = infinity RMT term:
* kappa_4: z(U_4) >= 5 for s >= 0.9 at q ~ 3e4 (z = 9 at s = 1.0, 46 at 1.4), s >= 0.8 at q ~ 1e6;
* kappa_3: s >= 1.3 at q ~ 3e4 (5.6), s >= 1.2 at q ~ 3e5 (6.1; 34 at 1.4);
* kappa_2: s >= 1.4 at q ~ 3e5 (6.0).
But the pilot deviations D_m = S_m - U_m are of the same size as U_m (D_4(1.4) = +0.074 +- 0.005), and their decay
with q is not predicted by any available theory: the 1/log q terms of the 3- and 4-level densities for this family
(they would follow from the ratios conjecture of Conrey–Farmer–Zirnbauer for the family mod q, carried to 3 and 4
levels with all arithmetic factors) are not in print as far as found. Hence:

* No definite test of the beyond-range RMT prediction is possible at q <= 1e6 with the N = infinity values (they would
  be "falsified" by finite-q terms we can already see) -- recommendation: do not register it as a pass/fail test.
* Definite readings that remain (round_draft.txt):
  T0 the exact finite-q identity in the proven range (S_m = 0 within sampling);
  T1 the conjectured Gaussian extension beyond HR's proof: S_3 = 0 on s in (2/3, 1] (U(N) says kappa_3 = 0 up to 1),
     S_4 = 0 on (1/2, 0.7] (kappa_4^U < 1e-5 there) -- with the caveat that a nonzero S could be a lower-order term;
     the pilot has S_3 = +0.0003, +0.0010, +0.0018 at s = 0.8, 0.9, 1.0 (0.4-1.2 se), the fresh set resolves 1e-4;
  T3 the direction of convergence: D_m(q) = S_m - U_m against 1/log q across fresh moduli spanning log q = 10.3-13.8;
     the conjecture needs |D| -> 0. Linear extrapolation in 1/log q is checked on the deficit R_m, which is known to
     vanish: for m = 4 it extrapolates to +0.006, +0.004, -0.0006, -0.007, -0.016 (s = 1.0, 1.1, 1.2, 1.3, 1.4)
     instead of 0 [reg/extrap_check_out.txt] -- registered as the allowance; for m = 2 and 3 the same check gives
     |alpha_R| = 0.005-0.024, larger than U_2, U_3 at s <= 1.4, so T3 is restricted to m = 4, s = 1.2-1.4
     (expected separation of "converges like 1/log q" from "pilot deviation persists": 4.2, 4.9, 4.3 combined errors).
  T1 power [reg/t1_power_out.txt]: pooled se over the 16 fresh moduli 2e-5..1e-4; a pilot-calibrated 1/log q
     lower-order term (S_pilot * 9.21/log q) would show at 6-18 se (m = 3, s >= 0.8; m = 4, s >= 0.6; m = 2, s = 1).
  T2 everything else beyond the thresholds: reported, no reading.

## 6. Prior art (web search of the named papers and related terms; not exhaustive)

* Hughes–Miller (Duke 2007, math/0507450): orthogonal family (newforms of prime level, split by sign), n-th centred
  moments Gaussian for supp in (-1/n, 1/n), and computed beyond (to 1/(n-1)) where they agree with RMT and odd moments
  become non-zero -- analytic, no numerics.
* Entin–Roditty-Gershon–Rudnick (GAFA 2013, 1208.5962): function fields, large-q limit, Katz–Sarnak equidistribution
  identifies Gao's combinatorial n-level factor with RMT (symplectic family) -- theorem, no number-field numerics.
* Chandee–Lee / Chandee–Klinger-Logan–Li (1706.02848): n-level density of primitive Dirichlet L-functions averaged
  over moduli, total support < 2 (large sieve) -- theorem, averaged over q.
* Fiorilli–Miller (ANT 2015, 1111.3896): 1-level density over all chi mod q, q in (Q/2, Q], lower-order terms, a term
  not predicted by the ratios conjecture, support up to (-4, 4) under hypotheses -- 1-level only, averaged over q.
* Drappeau–Pratt–Radziwiłł (2002.11968): unconditional 1-level density beyond support 2, averaged over q.
* Variance of primes in progressions (the m = 2 regime): Hooley, Friedlander–Goldston, Fiorilli (1301.5663,
  distribution of the variance), Fiorilli–Martin (2008.05837, disproof of Hooley's upper bound for q ~ loglog x);
  Keating–Rudnick (function fields). Divisor sums in progressions (the m = 4 (2|2) regime): KRRR, Rodgers–Soundararajan
  (1610.06900).
Verdicts:
* Numerical test of the moments/cumulants of W_f for chi mod a prime q beyond the mock-Gaussian range: NOT FOUND.
* Finite-q diagonal formulas: the split into equality (J_eq, J_diag, J_non) and congruence (J_cong) terms is Hughes–
  Rudnick's own proof -- KNOWN IN SUBSTANCE; the exact per-prime Steinhaus cumulants and the exact no-wrap deficit
  (model R) as an additive finite-q correction: NOT FOUND as such (for m = 2 the same boundary terms are part of the
  Montgomery–Hooley / Friedlander–Goldston variance computation -- KNOWN IN SUBSTANCE there).
* kappa_3^U = 0 up to support 1 (odd-m threshold 2/(m-1)): follows directly from the Soshnikov/Hughes–Rudnick support
  condition; not seen stated for the L-function family in the sources read.

## 7. Doubts

1. The family average is the whole population; "se" is a model of the arithmetic q-to-q fluctuation (characters as if
   independent draws, n_eff = (q-1)/2). Delta method = bootstrap within ~10 %; whether it matches the scatter between
   moduli is untested (no other modulus looked at); the registration reads chi^2/dof across moduli as a check.
2. That the RMT term lives in the wrap-around part and that the deficit -> 0 rests on my O(1/log q) estimate plus the
   computed decay at q = 1e4, 1e5, 1e6; the decay is slow (~1/L) and its exact form unknown.
3. m = 4 at the pilot: S_4 ~ 0 against U_4 = -0.068 at s = 1.4 is either a large positive lower-order term or a gap in
   the identification; one modulus cannot tell.
4. The m = 2 lower-order bookkeeping (singular-series constant, Euler–Maclaurin, n m = 1 pairs) is incomplete, so even
   m = 2 has no frozen finite-q amplitude.
5. Only the Fejer f-hat was used; its edge zero makes near-threshold cells uninformative. A band-pass f-hat
   concentrated near |u| = u0 would sharpen the 2/m thresholds (kappa_4 at u0 ~ 1/2) but not remove the finite-q terms.
6. Supports are capped by the (2|2), (3|1) deficit cost, not by the sieve: q^s up to ~1e11 is feasible for P itself.
