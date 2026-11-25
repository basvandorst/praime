# Prior-art gate: murmurations vs the family theorem (fact 2054) and the one-point law (fact 2055)

Folder: this folder
(the downloaded sources are not stored; reading log `reading_log.txt`). Labels: **proved**, **proved (cited)**,
**derived** (algebra done here, exact), **computed**, **heuristic**, **conjecture**.

## 1. Verdicts

- **C1 (murmuration densities = root number x coefficient correlation carried by the Petersson term c = N and the kernel
  J_{k-1}(4 pi sqrt(p/N))): KNOWN.** Miller–Montague, arXiv:0911.1830, Lemma 3.5 (label `lem:V1`) and its proof
  (equations labelled `ln:JPH`, `ln:INT`); ILS 2000 Lemmas 6.5–6.6 (cited there); restated as a murmuration theorem in
  Cowan, arXiv:2408.12723, Theorem 1.4; level-1 weight-aspect version in Sarnak's letter (Aug 2023), eqs. (15)–(17), and
  Lowry-Duda arXiv:2506.01640 eq. (ils_murmuration). Natural-weight (Eichler–Selberg) version: Zubrilina arXiv:2310.07681,
  Theorem 1 and its Bessel form Theorem 3 (numbering of the v2 source; labels mainthm, tm3).
- **C2 (murmurations <-> zero sums twisted by the functional-equation factor; one-level density with root number):
  KNOWN IN SUBSTANCE.** The murmuration = transition-range term of the root-number-split one-level density (Sarnak letter
  eqs. (7'), (19): "the third term ... is the murmuration term"); its source is the functional-equation factor X of the
  approximate functional equation / ratios recipe (Cowan 2408.12723, Remark 3.3 and Sec. 6; Miller–Montague Lemma 3.5 gives
  the root-number term as 2 lim int phi(x log R) chi(4 pi i x) X_L(1/2 + 2 pi i x) dx with chi(u) = sum_b mu^2(b)/(phi(b) b^u)).
  Fact 2054 is this formula evaluated (formally: MM need compact Fourier support) at the test function h * X^{-1}: the X
  cancels, the b = 1 term is the main term, the b >= 2 terms are the lines r = 1/b^2 of Theorem D; the proof route (c = Nb,
  mu^2(b)/phi(b), GR 6.561.14) is the same as the repo's Lemmas 2–4. **Extra in the repo, not found:** the statistic
  sum_rho h(rho) X(rho)^{-1} at height T -> infinity (not low-lying, no compact Fourier support), the exact main term
  -(1/2pi) int H, the unconditional range (Thm C), Hecke twists lambda_f(l) with Ramanujan sums (Thm D for l > 1), and the
  per-form reading (Prop D').
- **C3 (single L-function statement of murmuration type): KNOWN IN SUBSTANCE on the prime side, NOT FOUND on the zero side
  and NOT FOUND as "murmuration of one form".** ILS 2000 Appendix C, (C.26)–(C.33) (level 1, prime-side resonance with
  constant Z, conditional on Hypotheses B, C; plus (C.8') and the closing remark on zeros of L(s, Delta), Rubinstein's
  numerics). Sawin–Sutherland arXiv:2504.12295, proof of Prop. 1.6 (label `primes`, lines 339–348 of the source): an exact identity for
  ONE elliptic curve (any conductor, root number eps(E)): Ramanujan-sum sieve + Voronoi with additive twists, then a
  Cramér heuristic; reading off its dual terms n = q^2 gives the repo's constant c_E = sum_q mu(q) lambda_E(q^2)/q
  (this read-off is not written there). The murmuration papers say the opposite of a one-form law (Martin 2505.01093,
  Sec. 2 remark (1): "something one only sees averaging over large families"; BBLL 2310.07746 remark 2).
- **C4 (predicted densities for degree >= 3 / non-GL(2) families; the constants B_L(1)/L(1, Ad L), prod G_p(1/p)):
  NOT FOUND.** Sym^2: proposed as a problem only (Lowry-Duda 2506.01640 Sec. 5.3; Kundu–Müller 2507.00738 intro, who also
  report Sutherland's observation that degree 3 looks most regular in (P/N)^{1/3}); genus 2 and 3: data only (Sutherland's
  pages and slides); Sawin–Sutherland Sec. 1.3 suggests their method for geometric families and hypergeometric motives,
  no formula. No source contains B_L(1)/L(1, Ad L) or prod_p G_p(1/p).


## 2. Sources (status)

| source | status | nearest formula |
|---|---|---|
| Sarnak, Letter to Sutherland and Zubrilina on murmurations and root numbers (Aug 2023), IAS pdf "Nina and Drew letter_0.pdf" (landing page publications.ias.edu/sarnak/paper/2726 refused: Cloudflare) | read in full (26 handwritten pages, as images) | (7') transition of root-number-split 1-level density; (12) Zubrilina density; (13) quadratic characters (with Rubinstein); (14)–(17) arithmetic (harmonic) weights, level 1: M^arith(x) = 4 pi W sum_c mu^2(c)/(c^2 phi(c)) Phi(16 pi^2 x/c^2), Z^arith = point masses; (19) "the third term is the murmuration term" |
| Zubrilina, Murmurations, arXiv:2310.07681 (Invent. Math. 2025) | read: intro, Thm 1.1–1.5, Sec. 2, Sec. 3.4, Sec. 4–6 | Theorem 1 (density, natural weights, squarefree N); Theorem 3: M_k(y) = alpha sqrt(y) sum_{d,s} Q(d) J_{k-1}(4 pi s sqrt(y)/d)/s, Q(d) = mu^2(d) prod_{p|d} p^2/(p^4-2p^2-p+1); line 373: harmonic weights grow like y at 0 |
| "On Murmurations", arXiv:2609.29631 | intro read | corrects Zubrilina's proof (U_{k-2} inside the r-sum; Lemma 7.3); the density of Thm 1.1 stands |
| Miller–Montague, An orthogonal test of the L-functions Ratios Conjecture II, arXiv:0911.1830 | read: Sec. 2.3, Sec. 3 (Lemmas 3.3–3.5 with proofs), App. A | Lemma 3.5 (lem:V1): root-number-twisted harmonic prime sum = 2 lim int phi(x log R) chi(eps + 4 pi i x) X_L(1/2 + 2 pi i x) dx + O(N^{sigma/2-1+eps}), chi(u) = prod_p (1 + 1/((p-1)p^u)); proof via c = Nb, mu^2(b)/phi(b), GR 6.561.14 |
| Cowan, Murmurations and explicit formulas, arXiv:2306.10425 | read in full | Heuristic 2.2: averaged explicit formula, zero density rho_F(gamma)/(1/2 + i gamma); unweighted zeros |
| Cowan, Murmurations and ratios conjectures, arXiv:2408.12723 | read: Sec. 1, 2 (Thm 2.1), 3 (Thm 3.1, Remarks 3.2–3.3), 5.3, 5.4, 6 | Thm 1.4 (newspaces, GRH, from Miller–Montague); Remark 3.3: the factor (d/pi)^{-r} of the functional-equation term yields murmurations; Thm 2.1: conversion to Sarnak's density |
| Bober–Booker–Lee–Lowry-Duda, weight aspect, arXiv:2310.07746 | read: Sec. 1 (Thm 1.1, Remarks 1–9), Sec. 2 | nu(E) with point masses at squares; Remark 4: murmurations visible in ILS (Petersson); Remark 9: integer murmurations |
| Booker–Lee–Lowry-Duda–Seymour-Howell–Zubrilina, Maass, arXiv:2409.00765 | read: Sec. 1, start of Sec. 2 | same nu(E) as BBLL |
| Lee–Oliver–Pozdnyakov, Dirichlet characters, arXiv:2307.00256 | read: abstract, Thms 1.1–1.2, Cor 1.3, Sec. 5 theorems | Gauss-sum-normalized characters: density cos(2 pi y) (degree-1 resonance) |
| Wang, Hecke L-functions of imaginary quadratic fields, arXiv:2503.17967 | read: abstract, Sec. 1.1–1.4 (Thm 1) | almost periodic density, Sp symmetry |
| Sawin–Sutherland, elliptic curves by height, arXiv:2504.12295 | read: Sec. 1, Sec. 2 (Lemma 2.1, proofs of Props 1.5, 1.6) | per-curve exact Voronoi + Ramanujan-sum identity (proof of Prop 1.6); Conj 1.1; prime-conductor prediction (Sec. 1.5) |
| Lowry-Duda, On murmurations and trace formulas, arXiv:2506.01640 | read in full | eq. (ils_murmuration); Sec. 5.3 Sym^2 murmurations proposed |
| Kundu–Müller, p-power coefficients, arXiv:2507.00738 | intro read | lambda_f(P^2) density = same shape as Zubrilina; Sym^2 "more delicate"; degree 3 regular in (P/N)^{1/3} (Sutherland, unpublished) |
| Martin, Variations on murmurations, arXiv:2505.01093 | read: abstract, Sec. 1–2 | remark: correlation only visible in large families |
| Kuan, Murmurations using Petersson, arXiv:2507.11418 | intro + Thm 1.1 read | level 1, weight aspect, harmonic |
| ILS 2000, Publ. IHES 91, Appendix C | relevant section read (OCR text of the published version saved by an earlier session) | (C.26)–(C.33): prime-side resonance with constant, level 1, under Hypotheses B, C; final remark on zeros (Rubinstein) |
| Martin, Root number bias, arXiv:2207.08121; Distribution of local signs ..., arXiv:2409.02338 | abstracts (+ grep of source) | family statements only |
| He–Lee–Oliver–Pozdnyakov arXiv:2204.10140; depth aspect arXiv:2603.25564; Carmichael arXiv:2502.03437 | grep / abstract only | no single-form or zero statements found by grep |
| Sutherland: murm_g2c.html, murm_nsusu.html, slides MurmurationsSimons/Harvard/NTW.pdf | text layer only (figures not inspected) | genus 2/3 data, no formula |

## 3. Mathematics: the transform (derived / proved / computed as labelled)

**3.1 Dictionary (derived; the integral is GR 6.561.14 = repo Lemma 4).** Let M(y) be a family murmuration density in the
normalization E_F[eps_f a_f(p)] ~ M(p/N), a_f(p) = sqrt(p) lambda_f(p), weight k. The prime part of the family average of
Z_f(H; r) = sum_rho h(rho) X_f(rho)^{-1} r^{1/2-rho} is (explicit formula, fact 2052, PNT smoothing; heuristic interchange)

    E_F[eps_f Z_f(H; r)] ~ -(1/2pi) int_0^inf y^{-1} M(y) I_{Nyr} dy = -(1/2pi) int H(t) e^{2 i theta(t)} (Nr)^{-it} Mtilde(-it) dt,

Mtilde(w) = int_0^inf y^{w-1} M(y) dy. For one Bessel term M = c * sqrt(y) J_{k-1}(4 pi sqrt(rho y)):
Mtilde(-it) = c (2 pi sqrt(rho))^{-1+2it} Gamma(k/2 - it)/Gamma(k/2 + it), and with e^{2i theta(t)} = (N/4pi^2)^{it} Gamma(k/2+it)/Gamma(k/2-it)
all gamma factors and N cancel:

    Bessel term  c * sqrt(y) J_{k-1}(4 pi sqrt(rho y))   <->   line r = rho with amplitude A(rho) = c/(2 pi sqrt(rho)),
    E_F[eps_f Z_f(H; r)] = -(1/2pi) sum_rho A(rho) Hhat(ln(r/rho)),   i.e.  M(y) = 2 pi sum_rho A(rho) sqrt(rho y) J_{k-1}(4 pi sqrt(rho y)).

(Degree d: the Bessel kernel is replaced by the inverse Mellin transform of the gamma ratio of X, the GL(d) Voronoi kernel,
oscillating like e(d (rho y)^{1/d}); same cancellation.)

**3.2 Harmonic weights, prime level (proved modulo the standard averaging of Kloosterman sums over primes).** Repo Lemma 2 +
Lemma 3(i),(ii): sqrt(N) Sum^h eps_f lambda_f(p) = 2 pi sum_{(m,N)=1} S(p, Nbar; m) m^{-1} J_{k-1}(4 pi sqrt(p/N)/m); averaging p
in short windows, S(p, Nbar; m) -> mu^2(m)/phi(m). Hence the harmonic density

    M^h(y) = 2 pi sqrt(y) sum_m mu^2(m)/(m phi(m)) J_{k-1}(4 pi sqrt(y)/m),   lines rho = 1/m^2, A = mu^2(m)/phi(m),

which is exactly Theorem D at l = 1 (uniform form) and exactly Miller–Montague's chi(u) = sum_b mu^2(b)/(phi(b) b^u) times
X_L: the m = 1 term gives the main term -(1/2pi) int H of fact 2054. With a Hecke twist l: M^h_l(y) = E^h[eps lambda(l) a(p)]
= 2 pi sqrt(y) sum_m mu(m) c_m(l)/(m phi(m)) J_{k-1}(4 pi sqrt(l y)/m) <-> Theorem D's lines rho = l/m^2 with A = mu(m)c_m(l)/(phi(m) sqrt(l)).
So fact 2054 (with its last formula) is the Mellin/Bessel transform of the (Hecke-twisted) harmonic murmuration density,
and conversely the line amplitudes plus the gamma factor determine the density (Mellin inversion on Re w = 0).

**3.3 Zubrilina's density (natural weights, squarefree N) <-> the one-point law with lines (derived exactly; the averaging
step is heuristic).** Zubrilina, Theorem 3 (quoted): "For any even k > 0, the weight k murmuration density function is given by
M_k(y) = alpha sqrt(y) sum_{d,s in N} Q(d) J_{k-1}(4 pi s sqrt(y)/d)/s, where Q(d) := mu^2(d) prod_{p|d} p^2/(p^4 - 2p^2 - p + 1)",
alpha = 2 pi prod_p (1 - p - 2p^2 + p^4)/(p^4 - 2p^2 + p) (Theorem 1). Apply 3.1 to it: the term alpha Q(d) s^{-1} sqrt(y) J_{k-1}(4 pi s sqrt(y)/d) is the line rho = (s/d)^2
with amplitude (alpha/2pi) Q(d) d/s^2; summing pairs with the same ratio, A(rho) = prod_p L_p(v_p(rho)) with
L_p(0) = (p-1)(p+1)^2/(p(p^2+p-1)), L_p(-2) = p^2/((p-1)(p^2+p-1)), L_p(2e) = p^{-2e} L_p(0), L_p = 0 otherwise.
The conjecture of fact 2055 with lines gives per form the local factor S_p(v)/(1-1/p) (Prop D'). Averaging with natural
weights over squarefree N: Plancherel measure at p not dividing N (E[lambda(p^{2j})] = p^{-j}), lambda(p)^2 = 1/p with random sign
at p || N, and P(p | N) = (p-1)/(p^2+p-1) (squarefree N weighted by phi(N)). **Result (sympy, exact, v = -6..6):
the averaged local factors equal L_p(v) for every v.** In particular the r = 1 line: E_nat[c_f] = prod_p L_p(0) = gamma/(2 pi^2)
= 0.862991 (Zubrilina's own identity alpha sum Q(d)/d = gamma/pi, checked). So Zubrilina's density is the natural-weight
family average of the per-form one-point law with lines, Bessel coefficient by Bessel coefficient.

**Answer to the question asked.** Transforming Zubrilina's density against the kernel of fact 2054 does NOT give
-(1/2pi) int H: it gives -(1/2pi)(gamma/(2 pi^2)) int H = -(0.862991/2pi) int H plus the lines rho = (s/d)^2 != 1 (only d = s
resonates at r = 1). The coefficient 1 of fact 2054 belongs to harmonic weights at prime level (3.2, the m = 1 Bessel term);
Zubrilina's family has natural weights and squarefree levels, and its coefficient 0.862991 is exactly the value the one-point
law predicts for that weighting (E_nat[c_f] over squarefree levels). So hypothesis (a) holds in the precise form "fact 2054
(with Hecke twists and lines) = Mellin/Bessel transform of the harmonic murmuration density; natural-weight densities are its
t^2-twisted averages", and hypothesis (b) holds as a derived (heuristic-averaging) statement: the per-form one-point law with
lines, averaged in the family's local model, reproduces Zubrilina's density exactly.

**3.4 Prime level, natural weights (derived; a theorem-level route also exists).** The same averaging at prime N gives
A((s/d)^2) = (6/pi^2) s^{-2} prod_{p|d} p^2/((p-1)^2(p+1)) ((s,d) = 1, d squarefree), and by Poisson in s
(the s-sum of J_1(cs)/s has f-hat supported on |xi| <= 1) the closed form, k = 2:

    M_prime(y) = 2 pi sqrt(y) - 12 y + (12/pi) sum_{1 <= r <= 2 sqrt y} nu'(r) sqrt(4y - r^2),
    nu'(r) = sum_{delta | r, delta sqfree} c(delta),  c(delta) = (C0/delta) prod_{p|delta} [p^2/((p-1)^2(p+1))]/(1 - u_p),
    u_p = 1/((p-1)^2(p+1)),  C0 = prod_p (1 - u_p) = 0.615133.

(Zubrilina's k = 2 formula has the same shape with beta, gamma, alpha nu(r) in place of 2 pi, 12, (12/pi) nu'(r).)
The same amplitudes follow WITHOUT the conjecture from Theorem D with l = t^2 and the harmonic -> natural conversion
(natural weight = harmonic x sum_t lambda_f(t^2)/t): A_nat(rho) = sum_{t/m = sqrt(rho)} mu(m) c_m(t^2)/(phi(m) t^2);
computed equal to the formula above for 12 pairs (s, d) to 1e-6 (hecke_to_natural_check_out.txt). This is the reverse
direction: fact 2054's Hecke-twisted lines determine the natural-weight density.

**3.5 Transforms (computed, transform_check_out.txt).** -(1/2pi) int y^{-1} M(y) G(y) dy divided by -(1/2pi) int H, Gaussian H:

| t0, sigma | harmonic M^h | predicted lines | Zubrilina M_2 | predicted lines | prime natural M_prime | predicted lines |
|---|---|---|---|---|---|---|
| 10, 2 | 1.00613-0.02073i | 1.00577-0.02062i | 0.84170-0.00906i | 0.84141-0.00897i | 0.58204-0.00777i | 0.58177-0.00769i |
| 14, 2 | 1.01948-0.01089i | 1.01819-0.01133i | 0.89478+0.00987i | 0.89383+0.00957i | 0.63932+0.00770i | 0.63835+0.00727i |
| 18, 2.5 | 1.00253+0.00047i | 1.00243+0.00044i | 0.86364-0.01647i | 0.86359-0.01650i | 0.61201-0.01632i | 0.61195-0.01634i |

Main terms (r = 1 line only): 1, 0.862991, 0.607927. Agreement transform vs line sum: <= 1.3e-3 (quadrature). The
Bessel–Mellin integral itself (mpmath quadosc): 7.6e-10 to 9.8e-6 (nu = 1), <= 1e-11 (nu = 3).

## 4. Numerical check on real data (computed)

Data: the class-number formula T(p,N) = H(4pN)/2 + sum_{r>=1, r^2 N < 4p} H(4pN - r^2 N^2) - (p+1) for
sum_{f in S_2^new(N)} eps_f a_f(p) (k = 2, N squarefree, p not dividing N) was first checked EXACTLY against PARI
(-trace(T_p W_N) from mfheckemat and mfatkininit) at N = 389, 390, 431, 101: 1014 pairs (N, p), 0 mismatches, 931 of them with
r-terms (trace_validate_out.txt). Then natural averages E[eps_f a_f(p)] = sum T/sum dim over (i) all 103 prime levels
1500 <= N <= 2297 (total dim 16276) and (ii) all 120 squarefree levels 1500 <= N <= 1699 (total dim 11288), primes p <= 12N
(summarize_data_out.txt; 48 bins of width 0.25 in y = p/N):

| family | RMS(data - M_prime) | chi2 vs M_prime | RMS(data - M_zub) | chi2 vs M_zub |
|---|---|---|---|---|
| prime levels | 0.069 | 16.1 | 0.419 | 907.6 |
| squarefree levels | 0.408 | 605.7 | 0.065 | 10.1 |

(chi2 over 48 bins with bin standard errors; the errors are conservative, they include the slope of the density inside a bin.)
So: Zubrilina's Theorem 1 (k = 2) fits the squarefree family and not the prime family; the prime-level density M_prime derived
in 3.4 from the one-point law with lines (equivalently from Theorem D with l = t^2) fits the prime family and not the squarefree one.

Transform of the data = natural-weight family average of the prime part of sum_gamma H(gamma) e^{2 i theta_f(gamma)} (eps_f-weighted),
divided by -(1/2pi) int H (Gaussian H):

| window | family | data | model at the same primes | model integral (3.5) | main term |
|---|---|---|---|---|---|
| t0 = 10, sigma = 2 | prime | 0.5798-0.0079i | 0.5709-0.0150i (M_prime) | 0.5820-0.0078i | 0.6079 |
| t0 = 10, sigma = 2 | squarefree | 0.8327-0.0288i | 0.8354-0.0239i (M_zub) | 0.8417-0.0091i | 0.8630 |
| t0 = 14, sigma = 2 | prime | 0.7014+0.0572i | 0.6961+0.0550i (M_prime) | 0.6393+0.0077i | 0.6079 |
| t0 = 14, sigma = 2 | squarefree | 0.9487+0.0510i | 0.9488+0.0515i (M_zub) | 0.8948+0.0099i | 0.8630 |

1-sigma noise of the data transform (from the bin errors): 0.019–0.025. Data vs model at the same primes: <= 0.011
in every case (within noise); the wrong model is off by about 0.25 (M_zub at prime levels gives 0.826 and 0.960, M_prime at
squarefree levels 0.579 and 0.686: model_at_samples_*_out.txt). The common shift of data and sampled model away from the
integral at t0 = 14 (about +0.06 + 0.05i) is a prime-sampling effect present in the model alone (the PNT remainder, term E_1 of
Theorem A; that it appears at t0 = 14 and not at t0 = 10 suggests the zeta zero 14.13, not checked). The data are at
T <= 20 with N about 1500–2300, i.e. T of order N^{1/3}: the regime of Theorem B, not the limit.

## 5. Leads: testable predictions not used by the repo (heuristic unless stated)

1. **Weighting changes the family constant and creates new lines (derived; prime-level part follows from Theorem D +
   harmonic->natural conversion, see 3.4).** Natural weights at prime level: E_nat[eps_f Z_f(H; r)] = -(1/2pi) sum_rho A(rho) Hhat(ln(r/rho)),
   A((s/d)^2) = (6/pi^2) s^{-2} prod_{p|d} p^2/((p-1)^2(p+1)): A(1) = 0.607927, A(1/4) = 0.810569, A(4) = 0.151982,
   A(1/9) = 0.341959, A(1/36) = 0.455946. Squarefree levels (Zubrilina's family): A(1) = 0.862991, A(1/4) = 0.767103,
   A(4) = 0.215748, A(1/9) = 0.364074. Harmonic weights (Theorem D): A(1) = 1, A(1/m^2) = mu^2(m)/phi(m), A(4) = 0.
   The lines r = s^2 (s >= 2) exist only with natural weights: a sharp zero-side test (e.g. N = 10007 data of the repo's
   referee, weighted by 1 instead of omega_f).
2. **Elliptic curves of prime conductor (Sawin–Sutherland local model, Sec. 1.5 of 2504.12295).** The one-point law with
   lines averaged in their model gives A(1) = E[c_E] = zeta(3)/zeta(6) = 1.181565, A(1/4) = 1.050280, A(1/9) = 0.569683,
   A(1/25) = 0.293047, A(1/36) = 0.506385, A(1/49) = 0.196355 (lead_constants_out.txt); these are exactly the Bessel
   coefficients of THEIR predicted density (m = q term = zeta(3)/zeta(6), checked to 10 digits; A(1/4) from their double
   sum: 1.050280). So the repo's 138 curves (13 lines each) can be compared with Sawin–Sutherland's density: the mean of
   c_E over prime-conductor curves should tend to 1.1816, the r = 1/4 line mean to 1.0503 (harmonic GL(2): 1).
3. **Sym^2 families (degree 3).** Constant c = prod_p G_p(1/p) of fact 2057: Sato–Tate mean of G_p(1/p) = 1 exactly
   (traces: E[|lambda(p)|^2 - 1] = 0, E[lambda(p^3)] = 0 since Sym^3(Sym^2) = Sym^6 + Sym^2), so a harmonic Sym^2 family has
   constant 1; Plancherel mean = 1 - p^{-2} + p^{-5}, natural Sym^2 family constant prod_p(1 - p^{-2} + p^{-5}) = 0.636446
   (good primes only). The density is then sum_rho A(rho) K_3(rho y), y = p/N^2, K_3 the inverse Mellin transform (Re w = 0)
   of the gamma ratio of X_{Sym^2 f}, oscillating like e(3 (rho y)^{1/3}) (fact 2056), i.e. regular in (p/N^2)^{1/3}, which
   matches Sutherland's unpublished observation quoted by Kundu–Müller. The harmonic version is checkable by Petersson
   (sum^h lambda_f(p^2) at p ~ N^2 is a Kloosterman series) without any zeros.
4. **General recipe (degree d, any family with a local model):** M_F(y) = sum_rho A_F(rho) K_d(rho y),
   K_d(y) = (1/2pi) int y^{it} prod_j Gamma_R(1/2 - it + conj(mu_j))/Gamma_R(1/2 + it + mu_j) dt (inverse Mellin on Re w = 0 of
   Q^{it} X(1/2+it)), A_F(rho) = E_F[c_L rho_L(rho)/sqrt(ab)] with the per-form local factors of facts 2055/2057 averaged in
   the family's local model (Sato–Tate, Plancherel, point counts). Degree 1 reproduces the shape cos(2 pi y) of
   Lee–Oliver–Pozdnyakov (Gauss-sum normalized characters) up to normalization (not checked in detail).
5. **Arithmetic -> natural weights.** Lowry-Duda (2506.01640 Sec. 3) writes that "it's not apparent how to translate an
   arithmetically weighted murmuration to a non-arithmetically weighted normalization". At the level of lines this is
   explicit: A_nat(rho) = sum_{t/m = sqrt(rho)} mu(m) c_m(t^2)/(phi(m) t^2) (3.4, checked); for level 1 this would turn
   Sarnak's (16)–(17) into a natural-weight weight-aspect density to compare with BBLL's nu(E) (not done).
6. **Local root numbers (Martin 2409.02338 conjectures murmurations for combinations of local signs).** At p || N the per-form
   line factor of the conjecture is S_p(v)/(1-1/p) with lambda_f(p) = -w_p(f) p^{-1/2} (k = 2, w_p the Atkin–Lehner sign):
   v = 0: 1 + 1/p; v = +1: -w_p (1 + 1/p)/p; v = -1: w_p/(p-1); v = +2: p^{-2}(1 + 1/p); v = -2: 0 (derived, not checked). So the
   family average of eps_f w_q(f) a_f(p) over levels divisible by a fixed prime q should have Bessel components only at
   rho = q^{+-1} s^2/d^2, with these local factors: a concrete formula for the local-sign murmurations Martin conjectures.
7. **Saturation exponent.** Sarnak's letter (9)–(10): averaging over p in windows of length N^gamma is needed, with the guess
   gamma + delta > 1 (delta = conductor dimension of the family). A single L-function has delta = 0, so gamma > 1: the
   signal appears only after averaging over p across many multiples of N, i.e. only in the large-y (T -> infinity) regime,
   which is exactly where the one-point law lives (prime sum up to x = N T^2/4 pi^2, signal ~ c x^{3/4} vs noise ~ x^{1/2}).
   This is the precise sense of "murmuration of one form" and is consistent with Martin's remark that fixed-p correlations are
   invisible for one form.

## 6. What I could not open or check

- Not opened: the published versions of Zubrilina (Invent. Math. 2025), BBLL (ANT 2026, msp.org/ant/2026/20-2/p03, PDF link
  seen in a search, not fetched), Cowan's companion paper on quadratic twists ([quadratictwists] in 2408.12723), Conrey–Snaith and
  Huynh–Miller–Morrison (ratios conjectures for the families), ILS 2000 main text Sec. 6–7 (only through Miller–Montague's
  quotation and the earlier session's notes), Fiorilli–Miller, Hughes–Miller, Ricotta–Royer (lower-order terms of 1-level
  densities with root-number splitting: these could contain the X-weighted statistic in another guise; **unknown: not read**),
  Sutherland's letter to Rubinstein–Sarnak (math.mit.edu/~drew/RubinsteinSarnakLetter.pdf), Lowry-Duda's Maass-murmurations
  blog post cited by Sarnak, Sutherland's slide figures (only the text layer was read), the depth-aspect paper 2603.25564 beyond
  its introduction, HLOP 2204.10140 beyond grep, Martin's papers beyond abstracts.
- Refused: publications.ias.edu/sarnak/paper/2726 (Cloudflare challenge); the PDF of the letter itself downloaded fine.
- C2's "not found" part (the zero sum weighted by X(rho)^{-1} at height T, averaged over a family) is therefore conditional on
  the unread lower-order-term papers listed above.
- Not checked: the degree-d kernel K_d and the Sym^2 / local-sign leads (5.3, 5.6) numerically; the harmonic-weight density
  M^h against eigenform data (it is a consequence of the Petersson formula, which the repo already checked); squarefree-level
  version of fact 2054 (the repo proves N prime only); the identification of the t0 = 14 prime-sampling shift with the
  zeta zero 14.13; uniformity in l for Theorem D with l = t^2 (repo GAP 4), used in 3.4.
- Heuristic steps: interchange of the y-integral with the family/prime limits (3.1); Plancherel/local-model independence
  across primes and the interchange with the infinite product (3.3, 3.4, leads 2–3); Kloosterman averaging over primes in short
  windows (3.2, standard under GRH as in ILS/Miller–Montague).

## 7. Files

All in this folder:
`murmur_gate.md` (this note), `reading_log.txt`, `densities.py` (M_harm, M_zub, M_prime, kernel G), `transform_check.py` + `_out.txt`,
`local_factors_check.py` + `_out.txt`, `hecke_to_natural_check.py` + `_out.txt`, `mprime_crosscheck.py`, `trace_validate.py` + `_out.txt`,
`trace_density_data.py`, `data_prime_1.npz`, `data_squarefree_1.npz`, `summarize_data.py` + `summarize_data_out.txt`,
`model_at_samples.py` + `model_at_samples_prime_out.txt`, `model_at_samples_squarefree_out.txt`, `lead_constants.py` + `_out.txt`,
`pdf2txt.py`, `pdfimgs.py` (PDF text / page-image extraction), `dl/` (downloaded sources, one folder each; Sarnak letter PDF in `dl/sarnak/` (page images deleted after reading; regenerate with pdfimgs.py)).
