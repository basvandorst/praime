# Murmuration density of the symmetric-square family at the conductor scale p ~ N^2 (N prime)

Folder: `writeup/code/study/lowzero/universal/sym2murmur/ (jobs run in a scratch folder) `.
Labels: **proved** (exact, with a proof here or cited), **derived** (exact algebra done here, checked numerically),
**computed**, **heuristic**, **conjecture**.

## 0. Result

Family: f runs over the newforms of S_k(Gamma_0(N)), N prime, k even with S_k(SL_2(Z)) = 0 (k = 2, 4, 6, 8, 10, 14; tested k = 2, 4),
natural weights (each form weight 1). Statistic: the coefficient of L(s, Sym^2 f) at a prime p != N, lambda_f(p^2) = lambda_f(p)^2 - 1
(root number of Sym^2 f is +1; conductor N^2). Variable y = p/N^2, normalization M(y) = lim_{N -> oo} E_nat[ sqrt(p) lambda_f(p^2) ].

    M_k(y) = 24 [k=2] y^{3/2}  -  6 sqrt(y)/(k-1)  -  (12/(pi (k-1) sqrt y)) * sum_{1 <= m < 4y} nu(m) sqrt(m(4y - m)) U_{k-2}(1 - m/(2y))

U_{k-2} the Chebyshev polynomial of the second kind (U_0 = 1, U_2(x) = 4x^2 - 1); nu multiplicative-type, nu(m) = prod_l nu_l(v_l(m)) over
all primes l:

    nu_l(0) = 1 - u_l,  u_l = 1/((l-1)(l^2-1));   nu_l(1) = 1;
    nu_l(e) = 1 + 1/(l^2-1) + (l/(l+1)) sum_{j=1}^{floor(e/2)} l^{-j}  - [e odd] l^{-(e-1)/2}/(l^2-1)   (e >= 1)

(nu_l(2) = 1 + l/(l^2-1), nu_l(3) = 1 + 1/l, nu_l(oo) = l/(l-1)). So nu(m) = C0 * prod_{l | m} nu_l(v_l(m))/(1-u_l),
C0 = prod_l (1 - u_l) = 0.615133 (the same Euler product as the gate's M_prime), and the mean of nu over m is exactly 1
(local check: (1-1/l)(1-u_l) + (1/l-1/l^2) + sum_{e>=2} ... = 1; computed mean over m <= 1e5: 0.99997).
nu(1..16) = 0.6151 0.9227 0.6561 1.5378 0.6216 0.9842 0.6173 1.3840 0.9022 0.9324 0.6156 1.6404 0.6154 0.9259 0.6630 1.6916.

Finite N: E_nat[sqrt(p) lambda_f(p^2)] = rho_N (1/sqrt(p) + M_k(y)) + o(1), rho_N = (k-1)N/(12 dim S_k(N)) (exact factor of every
trace-formula term; the 1/sqrt(p) is the Plancherel term 1/p).

Shape: for y < 1/4 the sum is empty and M_k(y) = 24[k=2]y^{3/2} - 6 sqrt(y)/(k-1) exactly; at each y = m/4 a new term
sqrt(m(4y-m)) switches on (square-root corner); M_k stays bounded as y -> oo (the 6 sqrt(y)/(k-1) cancels the growth of the
class-number sum: computed, sect. 4.4). Values (k=2): M(0.1) = -1.138, M(0.5) = 0.920, M(1) = 2.540, M(2) = 2.768, M(3) = 3.523;
(k=4): -0.632, -0.307, 0.350, 2.382, 2.338.

Equivalent GL(3) form (derived; formal): M_k(y) = c_0 + sum_rho A_3(rho) K_3(rho y), K_3 = the inverse Mellin transform on Re w = 0 of
the gamma ratio of Q^{it} X_{Sym^2 f}(1/2+it) (oscillates like (rho y)^{1/6} e(3 (rho y)^{1/3})), main line
A_3(1) = prod_l (1 - l^-2 + l^-5) = 0.636446 for natural weights and A_3(1) = 1 for harmonic weights: exactly the Plancherel and
Sato-Tate means of the one-point-law constant c_f = prod G_p(1/p) of fact 2057. The (P/N^2)^{1/3} regularity belongs to y >> 1;
at y ~ 1 the corners at y = m/4 dominate.

Harmonic weights (derived from Petersson, heuristic averaging, NOT tested):
M^h_k(y) = sum_{s squarefree} (1/(sqrt(s) phi(s))) sum_{q >= 1, (q,s)=1} K21(y/(s q^2)),  K21(z) = 4 pi (-1)^{k/2} sqrt(z) cos(4 pi z) J_{k-1}(4 pi z);
GL(3) form: lines rho = j^2/(s d^2), d | s, amplitude mu(d)/(j sqrt(s) phi(s)); main line 1. M^h ~ y^{k-1/2} as y -> 0.

Honest size: a closed-form, data-checked (discovery only) density for one degree-3 family at its conductor scale, with every local
factor, plus an exact kernel identity that ties it to the GL(3) one-point law; the key evaluation step (sect. 2.4) is heuristic;
not found in print; holdout not yet run.

## 1. Setup

lambda_f(n) = a_f(n)/n^{(k-1)/2}; p != N prime: lambda_f(p^2) = lambda_f(p)^2 - 1 = lambda_{Sym^2 f}(p). Gamma factor of Sym^2 f:
Gamma_R(s+1) Gamma_C(s+k-1) = Gamma_R(s+1)Gamma_R(s+k-1)Gamma_R(s+k), mu = (1, k-1, k); conductor N^2 (f is Steinberg at N),
root number +1. The normalization sqrt(p)*lambda is the gate's dictionary (3.1) with Q = N^2: E[sqrt(p) lambda_L(p)] = M(p/Q) makes the
prime side of the zero statistic O(1) (sect. 5).

## 2. Route C: Eichler-Selberg trace formula (natural weights)

2.1 (proved; Cohen's form of the trace formula, N odd prime, gcd(n,N)=1, k as in sect. 0). Tr T_n | S_k(Gamma_0(N)) = A1 + A2 + A3 + A4,
A1 = n^{k/2-1}(k-1)(N+1)/12 [n square], A3 = -sum_{d|n} min(d,n/d)^{k-1} (the tau-sum over {1,N} gives 2), A4 = sigma(n)[k=2],
A2 = -(1/2) sum_{t^2<4n} U_{k-2}(t,n) [ H(4n-t^2)(1 + chi_N(t^2-4n)) + N H((4n-t^2)/N^2) [N^2 | 4n-t^2] ],
since mu(t,f,n) = 1 + chi_N(D/f^2) for N not dividing f and N + 1 for N | f (then N^2 | D and chi_N(D) = 0).
Implemented in `code/trace.py`; checked against PARI trace(mfheckemat) for N in {11,23,37,43,101}, k in {2,4}: 1110 checks, 0 mismatches
(`code/validate_trace.py`).

2.2 (proved) n = p^2. Subtracting the level-1 formula (Tr = 0 at level 1):
Tr T_{p^2} = p^{k-2}(k-1)N/12 - (2 + p^{k-1})/2 + A2^chi + A2^{N^2},
A2^chi = -(1/2) sum_{|t|<2p} U_{k-2}(t,p^2) H(4p^2-t^2) chi_N(t^2-4p^2),
A2^{N^2} = -N sum_{1 <= m < 4y} U_{k-2}(2p - mN^2, p^2) H(m(4p - mN^2))   (t = +-(2p - mN^2); the two signs coincide).
Normalized by p^{k-1} dim and multiplied by sqrt(p): Plancherel rho_N/sqrt(p); hyperbolic -6 rho_N sqrt(y)/(k-1) (1 + 2p^{1-k}).

2.3 (heuristic, Zubrilina-type class-number averaging; checked) H(D) = (sqrt(D)/pi) L(1, -D) with Zagier's L(1,-D) (local factors
(l/(l+1)) sum_e N_{l^e}(-D) l^{-e}). With p prime in a window and N prime, the average of pi H(D)/sqrt(D), D = m(4p - mN^2), is nu(m):
computed directly from class numbers over N in [101,199], p <= 3N^2 (`code/check_nu.py`):
m = 1..11: 0.6149(6) 0.9229(10) 0.6559(6) 1.5368(20) 0.6220(7) 0.9838(11) 0.6171(9) 1.3848(25) 0.9025(15) 0.9324(21) 0.6177(22)
against nu(m) = 0.6151 0.9227 0.6561 1.5378 0.6216 0.9842 0.6173 1.3840 0.9022 0.9324 0.6156 (agreement within errors; nu itself was
derived on route P, sect. 3, so this is the route-B/C cross-check asked for). With U_{k-2}(2p - mN^2, p^2) = p^{k-2} U_{k-2}(1 - m/(2y)) and
sqrt(D) = N sqrt(m(4y-m)), A2^{N^2} gives the m-sum of sect. 0.

2.4 (heuristic: the one unproved evaluation) A2^chi contributes only its smooth main term -(1/(2 pi)) S int U sqrt(4p^2-t^2) dt, where S is
the average over t of chi_N(D_t) L(1,-D_t)'s local factors: at N it is avg chi_N + (1/N) avg chi_N^2 + ... = -1/N + (N-2)/N^2 + ... = -2/N^2
(n = p^2 is a square mod N, so D_t == 0 mod N for two classes of t), elsewhere 1; int U sqrt(4p^2-t^2) dt = 2 pi p^2 [k=2] since
int_0^pi sin((k-1)theta) sin(theta) d theta = (pi/2)[k=2]. Normalized: 2p^2/N^2 / (p (N/12)) * sqrt(p) = 24 y^{3/2} for k = 2, 0 for k >= 4
(the same coefficient 24 = (6/pi^2)(4 pi^2) comes out of route P, 3.5).
In particular A2^chi has no lower-order term of size p^{k-1}, so the hyperbolic -(1/2)p^{k-1} of A3 survives. Evidence: (a) for y < 1/4,
where A2^{N^2} = 0, the data follow 24[k=2]y^{3/2} - 6 sqrt(y)/(k-1) (k=4 bins y in (0,0.25]: -0.21(7) -0.57(15) -0.72(21) -1.10(26) -0.83(29)
against -0.32 -0.55 -0.71 -0.84 -0.95); (b) route P without this term grows like +6 sqrt(y)/(k-1) (sect. 4.4), so only with it is M bounded.

## 3. Route P: Petersson formula (harmonic weights), then natural weights

3.1 (proved: classical) Sum^h lambda_f(m^2) = delta_{m,1} + 2 pi i^{-k} sum_{c == 0 (N)} S(m^2,1;c)/c J_{k-1}(4 pi m/c).
3.2 (proved: Salie evaluation, IK Lemma 12.3 type) for l odd, (a,l)=1: S(a^2,1;l^{2j}) = l^j (e(2a/l^{2j}) + e(-2a/l^{2j})); hence
S(p^2,1;N^2) = 2N cos(4 pi p/N^2): the modulus c = N^2 resonates with J_{k-1}(4 pi p/N^2) at p ~ N^2.
3.3 (heuristic: average over p in windows; Kloosterman sums with N not dividing c/N are dropped, size N^{-3/2}) For c = N^2 u, (u,N)=1,
twisted multiplicativity + reciprocity ubar/N^2 + (N^2)bar/u == 1/(uN^2) give
E_p[ S(p^2 n^2,1;N^2u) J(4 pi p n/(N^2u)) ] = 2N a(u,n) cos(4 pi y n/u) J_{k-1}(4 pi y n/u), a(u,n) = prod alpha(l, v_l(u), v_l(n)),
alpha(l,e,v) = phi(l^e)^{-1} sum_b e(-2 l^v b/l^e) S(l^{2v} b^2, 1; l^e) (computed, `code/local_a.py`, `code/local_a2.py`, l^e <= 729):
alpha(l,1,0) = 1/(l-1), alpha(l,1,v>=1) = -1, alpha(l,2j,0) = l^j, alpha(l,2j,v>=1) = 0, alpha(l,odd>=3,.) = 0 (l = 2 included).
So E^h[sqrt(p) lambda(p^2)] -> M^h of sect. 0 (n = 1; u = s q^2).
3.4 (heuristic: E_nat[X] ~ (6/pi^2) Sum^h X L(1,Sym^2 f), L(1,Sym^2f) = zeta(2) sum lambda(r^2)/r; Hecke lambda(p^2)lambda(r^2) =
lambda((pr)^2)) natural lines rho = r/(s q^2), gcd(r,sq) = gcd(s,q) = 1, s squarefree, q >= 1, amplitude
A(rho) = (6/pi^2) r^{-3/2} w(s) v(q), w(s) = prod_{l|s} l^{3/2}/((l-1)(l^2-1)), v(q) = prod_{l|q} l^2/(l^2-1); M^P = sum A K21.
3.5 (derived, exact) Poisson in r: sum_{r>=1} r^{-3/2} K21(rz) = 4 pi^2 z^{3/2}[k=2] - (4 pi/(k-1)) sqrt(z) sum_{1<=j<4z} sin((k-1)theta_j),
cos theta_j = 1 - j/(2z) (the r-summand is entire of exponential type, so only frequencies |j| <= 4z survive; Weber-Schafheitlin
int_0^oo J_nu(u) cos(bu) du/u = cos(nu arcsin b)/nu, = 0 for b > 1 and nu odd). Regrouping by m = j s q^2/d gives exactly
M^P = 24[k=2] y^{3/2} - (12/(pi(k-1) sqrt y)) sum_m nu(m) sqrt(m(4y-m)) U_{k-2}(1 - m/(2y)) with the nu of sect. 0 (closed form vs direct line
sum agree to the truncation error, `code/model.py` main). Route P therefore reproduces route C's class-number terms and their local
factors exactly, but not the hyperbolic term.

## 4. Route B: the one-point law in degree 3, and the kernel identity

4.1 Kernel (derived; Meijer-G form checked against GL(2)): K_d(z) = (1/2pi) int G(t) z^{it} dt, G(t) = prod_j Gamma_R(1/2-it+mu_j)/Gamma_R(1/2+it+mu_j),
= 2 G^{d,0}_{0,2d}((pi^d z)^2 | a; 1-a), a_j = (1/2+mu_j)/2 (`code/kernels.py`: reproduces 2 pi sqrt(z) J_{k-1}(4 pi sqrt z) for GL(2) to 6 digits).
A family density sum_rho A(rho) K_d(rho y) is exactly the one whose zero statistic E[sum_rho h X^{-1} r^{1/2-rho}] has prime side
-(1/2pi) sum A(rho) Hhat(ln(r/rho)) (gate dictionary 3.1).
4.2 (proved here: gamma-function identity; checked at 48 random complex w, k = 2..8, to 5e-31, `code/constants.py`; on the imaginary axis to
1e-16, `code/duality_check.py`)  chi(-2w) Ktilde21(w) = G_3(t = i w),  chi(s) = zeta(s)/zeta(1-s). Proof: Ktilde21 by Weber-Schafheitlin
with J_{-1/2}; Gamma(-w) cancels; duplication for Gamma(1/2+w), reflection for Gamma((3/2-k-w)/2) and Gamma((1/2+w)/2); the sine ratio is
(-1)^{k/2}, which cancels the (-1)^{k/2} of K21.
4.3 (derived, formal) Hence sum_{q} K21(y/(s q^2)) = (pole terms) + sum_j j^{-1} K_3(j^2 y/s) (Mellin: zeta(-2w) = chi(-2w) zeta(1+2w));
the q-lines of route P are the Poisson duals of GL(3) lines. Main GL(3) line: harmonic: only s = 1, j = 1: A_3(1) = 1. Natural: local
enumeration (r exponent a, s, Moebius e | rs, d from v(q) = sum_{d|q} h(d), h(l) = 1/(l^2-1), j exponent b; a + 2b - [l|s] - 2[l|e] - 2[l|d] = 0)
gives the local factor (1-l^-2)(1 + 1/(l(l^2-1)) - l^-3) = 1 - l^-2 + l^-5, product 0.636446: the gate's lead (Plancherel mean of
G_p(1/p)) is VERIFIED; harmonic constant 1 = Sato-Tate mean of G_p(1/p). So the per-form one-point law averaged in the family's local
model and the Petersson/trace computation give the same main GL(3) line. Convergence of the j-series is only conditional (K_3 grows
like z^{1/6}), so the GL(3) form is formal.
4.4 (computed) The Mellin transform of M^P has a double pole at w = -1/2 (zeta(3/2+w) from r, zeta(-2w) from q) times the simple zero
Ktilde21(-1/2) = 0 (since int_0^oo J_nu(x) cos x dx/x = cos(nu pi/2)/nu = 0): a sqrt(y) term. Numerically M^P(y)/sqrt(y) at y = 300.137:
5.969 (k=2), 1.980 (k=4), 1.162 (k=6) against 6/(k-1) = 6, 2, 1.2 (`code/largey.py`). The trace formula's hyperbolic term is exactly
minus this residue: the natural density M = M^P - 6 sqrt(y)/(k-1) has no pole at -1/2 and stays bounded (M(300.137) = -0.54, -0.34, -0.65).
The Petersson route loses this term in step 3.4 (interchange of the r-sum with the family average, uniform only for r << N); the trace
formula keeps it. For harmonic weights there is no double pole, so no such term is expected (untested).

## 5. Which step is the conjecture; regimes

- Proved: 2.1, 2.2 (exact decomposition), 3.1-3.2, 4.2 (identity), 3.5 (Poisson algebra).
- Heuristic: 2.3 (class-number averages; a theorem in this style needs averaging over N in [X, X+Y], as Zubrilina and Kundu-Mueller do; prime
  N only plus p-windows would need equidistribution of p mod d^2 jointly with N), and 2.4 (A2^chi has no p^{k-1} term: the main open step; via
  Petersson it is the statement that the Kloosterman sums S(p^2,1;Nm), N not dividing m, average to o(1/N) over p, plus the correct
  harmonic-to-natural conversion with the dual part of the AFE). The per-form one-point law is NOT used for the family formula; it reappears
  only as the identification of the main GL(3) line (4.3), which is a consequence, not an input.
- Regimes: the formula is for p ~ y N^2, y fixed, N -> oo. y < 1/4 (N << p < N^2/4): M = 24[k=2] y^{3/2} - 6 sqrt(y)/(k-1), i.e.
  E_nat[lambda_f(p^2)] = 1/p - 1/(2 dim) + O(...) with a p-independent bias -1/(2 dim) (the hyperbolic term) all the way down to p ~ N^{1+delta};
  the natural data at y in (0, 0.05] agree (k=2: -0.76(11) against -0.79). Kundu-Mueller (2507.00738) live at P^2 ~ N (P ~ sqrt N) and their
  statistic is root-number weighted, eps_f lambda_f(P^2): it is the GL(2) murmuration at the square index n = P^2, a different statistic at
  a different scale; it neither covers p ~ N nor p ~ N^2 (correcting the brief). In between (sqrt N << p << N^2/4) the unweighted average has no
  oscillating term at all in this model. y >> 1: the GL(3) lines take over (y^{1/6} e(3 y^{1/3})), the per-form one-point-law regime.

## 6. Discovery data (computed; looked at freely)

Data: own Eichler-Selberg code + Gram identity s(p) = sum_f lambda_f(p)^2 = w^T G^{-1} w, G_{mm'} = sum_f lambda_f(m)lambda_f(m') from traces
of T_{mm'/d^2}, m <= dim (`code/gen_data.py`); validated against exact PARI eigenform coefficients (all p <= 3N^2) at (k,N) = (2,101),
(4,37), (2,53): max relative error 1.1e-15 (`code/validate_data.py`). Class numbers: own sieve over reduced forms (`code/hurwitz.py`,
0 mismatches against PARI for m <= 1e5). Discovery range: k = 2, all 21 primes N in [101,199]; k = 4, all 14 primes N in [101,167];
every prime p <= 3N^2, p != N (147,334 and 76,756 pairs). (Also computed, not used: k=2 N = 53, 293; k=4 N = 37.)

Fits, 60 bins of width 0.05 on (0,3], dim-weighted bin means, empirical bin errors (`code/compare.py`, `code/fit2.py`, `test_dry.out`):

| k | amplitude of M (pred. 1) | chi2(M)/60 | chi2 null | no hyperbolic | smooth only | nu = 1 | harmonic M^h | 0.6364 K3 | 0.6364 K2 | Zubrilina shape |
|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 0.983 +- 0.075 | 56.0 | 228 | 2554 | 89565 | 1304 | 143 | 230 | 333 | 649 |
| 4 | 1.026 +- 0.062 | 55.6 | 327 | 617 | 411 | 165 | 145 | 372 | 438 | 808 |

Two-component fit (k=4): smooth part a = 1.026 +- 0.063, class-number part b = 1.025 +- 0.078 (corr 0.73); for k = 2 the parts are
degenerate (corr ~1). Noise: residual variance / (p/dim) = 0.969 (k=2), 0.953 (k=4) (independent Sato-Tate forms: 1); empirical/predicted
bin error (c = 0.96) median 1.007 and 1.001; correlation of k=2 and k=4 averages at the same (N,p): -0.022. No dependence on N mod 4
visible (bin means for N = 1 and 3 mod 4 agree within errors).

## 7. Not done / not checked

Squarefree levels; weights 6-14 (int64 limit in trace.py); the harmonic formula against data; the GL(3) form summed numerically; the
constant c_0; a proof of 2.3-2.4; the elliptic-curve and genus-2 families (task G: not done).

## 8. Prior art (verdict for the formula of sect. 0: NOT FOUND in the sources read; open while the sources marked unknown are unread)

| source | status | what it has |
|---|---|---|
| Kundu-Mueller, arXiv:2507.00738 | relevant sections read (abstract, intro incl. Outlook, Sect. 2 trace set-up, start of Sect. 3) | root-number-weighted eps_f P lambda_f(P^2), y = P^2/N, squarefree N: Zubrilina's shape; Sym^2 (root number +1, conductor N^2) posed as "more delicate"; Sutherland's (P/N)^{1/3} remark |
| Lowry-Duda, arXiv:2506.01640 | Sect. 5.3 read (gate: read in full) | Sym^2 murmurations proposed via GL(2) trace formulas on lambda_f(p^2); no formula; "initial heuristic estimates suggest this would be very challenging" (for lambda_f(p^m)) |
| Zubrilina, arXiv:2310.07681 | Thm 1 and its footnotes read; full source grepped | squarefree levels only; footnote: densities "slightly different when one considers all levels". The prime-level density M_prime of the gate is NOT stated there (settled) |
| Cowan, arXiv:2408.12723 | grepped in full (gate read Sect. 1-3, 5-6) | Sym^2 appears only through L(1, sym^2 E) as a normalizer |
| Sawin-Sutherland, arXiv:2504.12295 | grepped in full (gate read Sect. 1-2) | nothing on degree 3 |
| Hua-Mok, arXiv:2609.00777 (Sep 2026) | abstract only | root-number-weighted prime powers at p^k/X (squarefree, GL(2)); extends Kundu-Mueller; not Sym^2, not p ~ N^2 |
| arXiv API: all 101 entries matching "murmurations" (to 26 Nov 2026) | titles and abstracts scanned (`dl/arxiv_murm.xml`) | no abstract on Sym^2, GL(3), Rankin-Selberg or degree >= 3 densities; the only "cubic" hits are cubic characters over function fields (2608.01337) and Cardi's "cubic pointwise range" K up to X^3 (2607.08982), both GL(1)/GL(2) |
| Tomczak 2606.08353; Burrin et al. 2603.25564; Gao-Zhao 2607.20853; Cowan 2504.09944; Wachs 2603.22807, 2603.04604 | abstracts only | GL(1)/GL(2) |
| Sutherland, letter to Rubinstein-Sarnak (math.mit.edu/~drew/RubinsteinSarnakLetter.pdf); IHP2026 slides | could not open (image PDFs, no text tool in this environment) | unknown |
| ILS 2000 (symmetric-square sections); Ricotta-Royer arXiv:math/0703760 and 0806.2908 (Sym^r one-level densities, level aspect, prime level, harmonic weights); Dueñez-Miller; Güloğlu | not read (abstract of math/0703760 only) | unknown: these evaluate exactly the Petersson sums Sum^h lambda_f(p^2) and could contain the c = N^2 Salie main term in some support range |

Ingredients that are classical: Eichler-Selberg; Salie's evaluation of S(a^2,1;l^{2j}); Petersson; Weber-Schafheitlin. The kernel identity
chi(-2w) Ktilde21(w) = G_3 is a gamma-factor identity of the kind used for Voronoi formulas of Sym^2 lifts; not searched for.
