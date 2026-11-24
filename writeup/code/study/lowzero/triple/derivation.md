# The three-point function of Dirichlet L-zeros, with missing orbits: the third cumulant of the window count

Folder: `writeup/code/study/lowzero/triple/` (all files named
below are there). Nothing in the repo was changed; the repo's frozen `bklib.py` (sha256 9bfc1a3b…) and its cached grid
`amp/grid.npz` were copied into `lib/`.

| script | output | what it does |
|---|---|---|
| `check_B1.py` | `check_B1_out.txt` | the arithmetic factor B1 of the triple correlation recomputed from the ratios recipe's local factor (θ-integral) |
| `lib/gtab.py` | `lib/gtab.npz` | table of ζ′/ζ(1 + iε) + 1/(iε) (mpmath, 25 digits) |
| `lib/tplib.py` | | the three-point function T3 for χ mod k, the quadrature for κ[3](ℓ), the closed form of the diagonal part |
| `check_forms.py` | `check_forms_out.txt` | GUE limit; R3 = 0 at a coincidence; ratios form against Bogomolny–Keating 2013's form |
| `kf3_test.py` | `kf3_test_out*.txt` | quadrature convergence, prime cutoff, closed form against quadrature |
| `kappa3_tables.py` | `kappa3_tables_out.txt`, `kappa3_tables.json` | κ3(ℓ) and Δκ3(ℓ; p) at L = 8.1, 9.4, 10.2, 10.7 |
| `cue.py` | `cue_out_11.txt`, `cue_11.npz` | 150 Haar U(800) matrices: the GUE value by simulation; estimator bias |
| `noise.py` | `noise_out.txt`, `noise_clean_out.txt`, `noise_*resid.npz` | NOISE ONLY on the used zero sets (residuals of the indicator regression) |
| `gapcheck.py` | `gapcheck_out.txt` | the three sets whose κ3 residual is an outlier: missing pairs |
| `predict3.py` | `pred3_*.json`, `predictions_kappa3.json`, `predict3_out.txt`, `b1grid.npz` (cache) | frozen per-L-function predictions (full law and alternatives) |
| `interp_check.py` | `interp_check_out.txt` | the cached/interpolated path of predict3 against direct evaluation |
| `kf3_grid_check.py` | `kf3_grid_check_out.txt` | quadrature convergence of the arithmetic part |
| `power3.py` | `power3_out.txt`, `power3.json` | expected z of every reading |
| `registry.py` | `registry_out.txt` | expected value of every reading under every hypothesis (numbers of the draft) |
| | `round_draft.txt` | the pre-registration draft |

Labels: **proved** (algebra/calculus from stated formulas), **conjecture** (the ratios conjecture / BK heuristic itself),
**computed**, **checked** (an independent numerical comparison), **heuristic**, **literature**.

## 0. Answers in short

1. **R3 with missing orbits** (§3; conjecture, two routes that agree to 10⁻⁷). In unfolded units the connected
   three-point function of the zeros of L(s, χ), χ primitive mod k, is T3 = (1/L³) Σ of six terms
   J*_k(a; b1, b2) = Q_k(a+b1, a+b2) − P1_k(a+b1) P3_k(a, b1, b2) − P1_k(a+b2) P3_k(a, b2, b1), with
   P1_k(x) = e^{−Lx}(k/φ(k))² L(1+x, χ₀)L(1−x, χ₀)A_k(x), P3_k = B1_k + (L′/L)(1+a+c, χ₀) − (L′/L)(1+c−b, χ₀),
   L = ln(kt/2π): Conrey–Snaith's ζ formula with every prime p | k deleted from every sum and product. Per bad prime:
   the triple-orbit term q_p is removed from Q; P1 is multiplied by BK's two-point factor conj F_p(ε); P3 changes by
   −b1_p + ln p/(p^{1+a+c} − 1) − ln p/(p^{1+c−b} − 1). Not additive over bad primes (P1·P3 is a product).
2. **κ3 and Δκ3** (§2, §4, §5). κ3(ℓ) = κ[3](ℓ) + 3Σ²(ℓ) − 2ℓ, κ[3] = ∫_{[0,ℓ]³} T3. Per bad prime,
   Δκ3 = Δκ3_diag + Δκ3_off with the closed form Δκ3_diag = (6/π³) Σ_{m,n≥1} sin mθ sin nθ sin(m+n)θ/(mn(m+n)p^{m+n}),
   θ = πℓ ln p/L (the triple-orbit term; equal to −E[X_p³] of the random-phase model, i.e. Berry's diagonal one order
   up), and Δκ3_off = 3ΔΣ² + Δκ[3]_P, in which the two-point drop (3ΔΣ² down to −0.33 for p = 2) is cancelled by the
   three-point part to −0.007 … +0.006. For p = 2 (L = 9.4–10.7) Δκ3 rises to +0.016 near ℓ = 4–5.5 with a period-1
   sawtooth beyond; p = 3: ±0.008; p = 5: ±0.004. The level for k = 1 is 0.006–0.025 below GUE for ℓ = 0.5–6.
3. **GUE κ3**: +0.0453, +0.0112, +0.0119, +0.0068, …, +0.0019 at ℓ = 0.5, 1, 1.5, 2, …, 8; 150 CUE(800) matrices agree
   at all 16 lengths (|z| ≤ 1.4).
4. **Power** (§7–8; noise only from the used sets). Per-L-function sd of κ3 is 0.0015–0.005, far below CUE's. On a fresh
   set of 186 L-functions: level vs GUE ≈ 200 se (combined); p = 2 coefficients vs GUE 68 se, vs "triple-orbit term not
   removed" 84 se; p = 3: 25/38; p = 5: 10/12; the triple-orbit amplitude γ_Q to ±0.013. On the 546 pooled sets
   roughly √3 more. Recommendation: **run** (test Fn, `round_draft.txt`; predictions `predictions_kappa3.json`).
5. **By-product** (§6): κ3 exposes three used zero sets with a missing pair near the top of their band (C: D = −195;
   Dn: D = −939; complex set: q = 173), missed by round 136's band-mean completeness check.
6. **In print** (§9): the ζ triple correlation with lower-order terms was compared with 100,000 ζ zeros by Conrey–Snaith
   (2007); nothing found for Dirichlet L-functions or for the conductor's primes, nor for the third cumulant of the
   count as a test statistic.
7. **Doubts** (§10): the conjecture's finite-height accuracy at the 10⁻³ level of the κ3 level; the estimator's edge bias
   (< 5·10⁻⁴, unresolved); completeness (κ3 is very sensitive to missing zeros).

## 1. Sources read

- J. B. Conrey, N. C. Snaith, *Correlations of eigenvalues and Riemann zeros*, Commun. Number Theory Phys. 2 (2008),
  arXiv:0803.2795. TeX read: the ratios conjecture (Conj. 2), Theorem 3 (averages of log-derivatives), Theorem 5 (n-point
  correlation), §5 examples: the triple correlation of ζ, eqs. (Jzabb), (Qprime), (Pprime), and R_{ζ,t,3}.
- J. B. Conrey, N. C. Snaith, *Triple correlation of the Riemann zeros*, arXiv:math/0610495 (J. Théor. Nombres Bordeaux
  20 (2008)). TeX read: Theorem 1 and the numerical section. Its P(x, y) has a misprint, 1/(1 − p^(1−x+y)) where the 2008
  paper has (1 − p^(−1+x−y)); `check_B1.py` settles it in favour of the 2008 form (below).
- E. Bogomolny, J. P. Keating, *A method for calculating spectral statistics based on random-matrix universality with an
  application to the three-point correlations of the Riemann zeros*, J. Phys. A 46 (2013) 305203, arXiv:1307.6012. TeX read
  in full: kernel K(x, y) = sin(π(N̄(x) − N̄(y)))/(π(x − y)) with the short primes in N̄, R_n = ⟨⟨det K⟩⟩ over independent
  prime phases; their R3^diag (eq. rdiag) and R3^osc (eq. rosc). The TeX contains no figure and no numerics.
- Bogomolny–Keating, Nonlinearity 8 (1995) 1115 and 9 (1996) 911, and PRL 77 (1996) 1472: not obtained (not on arXiv);
  what they contain is taken from the introduction of Conrey–Snaith 2007 (BK96 "outlined" the three-point lower-order terms
  in semiclassical language; the full expression was a private communication, "bk06", which agrees with CS).
- Abstracts only: Fazzari–Gerspach, arXiv:2412.20099 (third moment of log ζ, assuming a triple-correlation conjecture; no
  numerics); Forrester–Mays, Proc. R. Soc. A 471 (2015), arXiv:1506.06531 (finite-size corrections, two-point function and
  spacings of Odlyzko's zeros; no three-point statistic).

## 2. Definitions; the cumulant identity (proved)

Zeros 1/2 + iγ of L(s, χ), χ primitive mod k; unfold by u = (γ/2π) ln(kγ/(2πe)) (density 1, as in r134.py); L := ln(kE/2π)
at height E; raw separations ε = 2πr/L for unfolded separations r. n(x) = #{j : x ≤ u_j < x + ℓ}.

With unit density, write R2 = 1 + T2, R3 = 1 + T2(12) + T2(13) + T2(23) + T3 (T_k the connected, "cluster" or Ursell
functions; T2 = R(r) − 1 of `bk/derivation.md`). The factorial cumulants of n are κ[k](ℓ) = ∫_{[0,ℓ]^k} T_k (standard;
e.g. Daley–Vere-Jones), and the ordinary cumulants follow from κ1 = κ[1], κ2 = κ[2] + κ[1], κ3 = κ[3] + 3κ[2] + κ[1]:

    Σ²(ℓ) = ℓ + ∫∫_{[0,ℓ]²} T2,      κ3(ℓ) = E[(n − ℓ)³] = κ[3](ℓ) + 3Σ²(ℓ) − 2ℓ,      κ[3](ℓ) = ∫∫∫_{[0,ℓ]³} T3.

T3 is symmetric, translation invariant and (for these L-functions) even, so with r1 = 0, r2 = t, r3 = s:

    κ[3](ℓ) = 6 ∫_0^ℓ (ℓ − s) F(s) ds,    F(s) = ∫_0^s T3(0, t, s) dt = 2 ∫_0^{s/2} T3(0, t, s) dt.

Check (proved): for ℓ → 0, T3 → 2 (R3 → 0, T2 → −1), so κ[3] ≈ 2ℓ³ and κ3 ≈ ℓ − 3ℓ² + 2ℓ³, the Bernoulli value.
κ3 is the three-point statistic measured; κ[3] is its purely three-point part. κ3 mixes the two- and three-point
functions; that is a feature, not a defect (§5.3): the bad-prime changes of the two parts cancel to a few per cent.

## 3. The three-point function for primitive χ mod k

### 3.1 Route (b): the ratios recipe (Conrey–Snaith 2008, adapted; conjecture → proved consequences)

For ζ, CS 2008 (§5, "Triple correlation for ζ") give, with ℓ_t = ln(t/2π),

    R_{ζ,t,3}(u, v, w) = ℓ_t³ + ℓ_t Σ_{6 ordered pairs} J*(iu; −iv) + Σ_{6} J*(a; b1, b2),

the six three-point terms being J*(−iw; iu, iv), J*(−iv; iw, iu), J*(−iu; iw, iv) and the three with all signs flipped, and

    J*(a; b1, b2) = Q(a+b1, a+b2) − P1(a+b1) P3(a, b1, b2) − P1(a+b2) P3(a, b2, b1),
    P1(x) = e^{−ℓ_t x} ζ(1+x) ζ(1−x) A(x),    P3(a, b, c) = B1(a+b, a+c) + (ζ′/ζ)(1+a+c) − (ζ′/ζ)(1+c−b),
    A(x) = Π_p (1 − p^{−1−x})(1 − 2/p + p^{−1−x})/(1 − 1/p)²,
    Q(x, y) = −Σ_p ln³p / ((p^{1+x} − 1)(p^{1+y} − 1)),
    B1(x, y) = Σ_p (1 − p^{−x})(1 − p^{−x} − p^{−y} + p^{−1−y}) ln p / ((1 − p^{−1+x−y})(1 − p^{−1−y})(1 − 2/p + p^{−1−x}) p^{2−x+y}).

R/(2π)³ is the density of distinct triples in raw height, so in unfolded units R3/L³; subtracting the pair terms
(which are the two-point function of `bk/derivation.md`) leaves the connected function

    T3(r1, r2, r3) = (1/L³) Σ_{6} J*(…),   evaluated at raw positions 2πr_j/L.               (proved from the conjecture)

**Adaptation to χ mod k** (the step that is ours; the same three rules as for two points in `amp/derivation.md` §1.1):
(i) the diagonal Dirichlet series carries |χ(n)|² = 1_{(n,k)=1}, so the local factor of every Euler product at p | k is 1:
every prime sum and product runs over p ∤ k; (ii) the global ζ-factors become L(·, χ₀), and each pole/zero pair
evaluated in a swap contributes the derivative of 1/L(1+z, χ₀) at 0, which is k/φ(k) instead of 1, so
ζ(1+x)ζ(1−x) → (k/φ(k))² L(1+x, χ₀)L(1−x, χ₀) and ζ′/ζ → L′/L(·, χ₀); (iii) X_χ(s+α)X_χ̄(1−s+β) ≈ (kt/2π)^{−α−β} and
−X′_χ/X_χ ≈ ln(kt/2π), so ℓ_t → L. Root numbers cancel (ε_χ ε_χ̄ = 1). Every primitive χ, real or complex, any exponent of
p in k, gives the same function (only |χ|² enters). So:

    T3_χ = (1/L³) Σ_6 J*_k,   J*_k(a; b1, b2) = Q_k(a+b1, a+b2) − P1_k(a+b1) P3_k(a, b1, b2) − P1_k(a+b2) P3_k(a, b2, b1),
    P1_k(x) = e^{−Lx} (k/φ(k))² L(1+x, χ₀) L(1−x, χ₀) A_k(x),
    P3_k(a, b, c) = B1_k(a+b, a+c) + (L′/L)(1+a+c, χ₀) − (L′/L)(1+c−b, χ₀),
    Q_k, B1_k, A_k: the sums and products over p ∤ k only.

### 3.2 Route (a): Bogomolny–Keating 2013 (random-matrix universality with short orbits; heuristic)

BK take the sine kernel with the smooth count N̄(E) replaced by the count including all orbits (prime powers) up to a
cut-off p*, R_n = ⟨⟨det K(E + e_i, E + e_j)⟩⟩, the average being an independent uniform average over the phases E ln p.
For L(s, χ) the oscillating density is −(1/2π) Σ_p Σ_m ln p p^{−m/2} (χ(p)^m e^{−imE ln p} + c.c.). For p ∤ k, χ(p) is a
phase that the uniform average over E ln p absorbs; for p | k the orbit is absent. So, *in this route the rule is
literally "the orbits p | k are missing"*: delete p | k from every average, and the divergent product
Π_{p<p*}(1 − 1/p)/(1 − p^{−1−s}) → s ζ(1+s) acquires Π_{p|k}(1 − p^{−1−s})/(1 − 1/p), i.e. ζ → (k/φ(k)) L(·, χ₀); the mean
density becomes L/2π. Their result, written in our notation, is

    T3 = (2π/L)³ [R3^diag + R3^osc],
    R3^diag = −(1/(2π)³) Σ_{p∤k} ln³p [1/((p^{1−ie12}−1)(p^{1−ie13}−1)) + 1/((p^{1−ie21}−1)(p^{1−ie23}−1))
                                    + 1/((p^{1−ie32}−1)(p^{1−ie31}−1))] + c.c.,
    R3^osc = −(e^{iL e12}/(2πi)³) (k/φ)²|L(1+ie12, χ₀)|² Π_{p∤k}(1 − (1−p^{ie12})²/(p−1)²)
             × [∂_{e3} ln|L(1+ie32, χ₀)/L(1+ie31, χ₀)|² − i Σ_{q∤k} ln q ( (q^{ie12}−1)/((q^{1+ie23}−1)(q^{1+ie13}−1))
               + (q^{ie12}−1)/((q^{1+ie31}−1)(q^{1+ie32}−1)) + (1−q^{ie12})²/(q−2+q^{ie12}) (1/(q^{1−ie31}−1) + 1/(q^{1−ie23}−1)) ) ]
             + cyclic + c.c.

R3^diag is the six Q-terms (the triple of orbits p^m, p^n, p^{m+n}: "diagonal" for three points); R3^osc carries the
phase e^{2πi(r_i − r_j)} (period 1 in unfolded units) and is the six P-terms.

### 3.3 What a prime p | k changes, term by term (proved from 3.1; identical in 3.2)

At fixed L, adding p to the bad set S:
1. **Diagonal triple-orbit term**: the p-summand q_p(x, y) = −ln³p/((p^{1+x}−1)(p^{1+y}−1)) of Q is removed (all m, n).
2. **Off-diagonal amplitude P1**: multiplied by f_p(ε) = (p − p^{iε})/(p − 2 + p^{−iε}) = conj F_p(ε), BK's two-point
   factor: the off-diagonal Euler factor of p is frozen at its ε = 0 value (p/(p−1))², exactly as in the two-point function.
3. **The bracket P3**, three changes: −b1_p(a+b, a+c) (the p-summand of B1 removed), + ln p/(p^{1+a+c} − 1) and
   − ln p/(p^{1+c−b} − 1) (from L′/L(·, χ₀) = ζ′/ζ + Σ_{p|k} ln p/(p^{1+z} − 1)).
4. No other factor changes (the GUE skeleton, the phase e^{−Lx} at fixed L, the root number).
Because P1 · P3 is a product, the change is **not additive** over the bad primes (as for two points); the diagonal Q-part
is additive.

### 3.4 Checks (computed)

- B1 against the recipe (`check_B1_out.txt`): −∂_{β2} log of the local factor of the swapped ratio, computed as a
  θ-integral with 4096 points, equals CS 2008's B1 summand to 3.9·10⁻¹⁰ (18 random imaginary arguments, p = 2 … 101);
  the 2007 P(x, y) does not (misprint).
- GUE limit (`check_forms_out.txt` (1)): with ζ(1+x) → 1/x, ζ′/ζ → −1/x, A → 1, B1, Q → 0, the six terms give
  2 S12 S23 S31 (S = sin πr/πr) to 3.8·10⁻¹⁵.
- R3 = 0 at a coincidence (check (2)): T3(0, t, s) → −2T2(s) as t → 0, with T2 from the repo's frozen bklib, for
  S = ∅, {2}, {3,5}, {2,7,13}: differences 1.6·10⁻³ at t = 10⁻³ and 1.5·10⁻⁴ at t = 10⁻⁴ (linear in t, as it must be). This
  ties the normalisation, the sign of the phase and every bad-prime factor of T3 to the two-point law confirmed in rounds
  134 and 136.
- The two routes (check (3)): our ratios-form T3 against BK 2013's form, coded independently (mpmath ζ and L′/L at every
  point, Euler products and prime sums over p < 10⁵, no splines), 16 random triples, L = 9.4–10.7, S = ∅, {2}, {3,5},
  {2,7,13}: largest difference 8.5·10⁻⁸ (prime truncation). So the two heuristics give one function, with the bad primes.

## 4. κ3: computation

### 4.1 Quadrature (computed; `kf3_test_out*.txt`)

s ∈ [0, 8] in pieces of 0.25 with 8 Gauss nodes; t ∈ [0, s/2] in pieces of ≤ 0.5 with 8 nodes (9,216 nodes). The
P-part is evaluated as (arithmetic − its GUE limit), so the GUE skeleton 2SSS is integrated separately and exactly; the
remaining poles (order ≤ 2 at a coincidence) cancel between the six terms (CS prove the sum has none), and in double
precision the cancellation costs < 10⁻⁷ at the nodes used. Convergence: the GUE part is identical to 7 decimals on a grid
of 54,400 nodes (`kf3_test_out.txt`); the full arithmetic κ[3] (L = 9.4, S = {2}) agrees between the two grids to
5·10⁻¹⁰ (`kf3_grid_check_out.txt`); B1 over p ≤ 10⁴ against p ≤ 10⁵: ≤ 1·10⁻⁶. The oscillatory integrals are thus done
"period by period" in effect: 32 nodes per unit in s and ≥ 16 per unit in t, against the period 1 of e^{2πir}.

### 4.2 The diagonal part in closed form (proved; checked to 10⁻⁶)

The cube integral of the six Q-terms factorises (each is a product of e^{−iκ(m r_u + n r_v − (m+n) r_w)}, κ = 2π ln p/L):

    κ[3]_Q(ℓ) = −(6/π³) Σ_{p∤k} Σ_{m,n≥1} sin(mθ_p) sin(nθ_p) sin((m+n)θ_p) / (m n (m+n) p^{m+n}),   θ_p = πℓ ln p/L,

and a bad prime adds Δκ[3]_Q(ℓ; p, L) = +(6/π³) Σ_{m,n} (same summand); leading term (3/π³) sin²θ sin 2θ / p². This is
exactly the third cumulant of the counting function in the random-phase (diagonal) model, n − ℓ = Σ_p X_p with
X_p = (1/π) Im Σ_m χ(p)^m (e^{−imφ} − 1) e^{−imθ_p}/(m p^{m/2}), θ_p uniform and independent, φ = 2πℓ ln p/L:
E[X_p³] = (3/(4π³)) Im Σ y_{m}y_{n}ȳ_{m+n} = −(6/π³) Σ sin sin sin/(mn(m+n)p^{m+n}) (using sin(A+B) − sin A − sin B =
−4 sin(A/2) sin(B/2) sin((A+B)/2)). So the Q-term of the ratios conjecture is Berry's diagonal approximation one order up,
the three-point analogue of Berry's Σ²_diag = Σ_p E[X_p²]. (Derived here; the identification is elementary.)
The numerical cube integral of the Q-terms equals the closed form to 1·10⁻⁶ (`kf3_test_out2.txt`).

### 4.3 GUE (computed; checked by simulation)

κ3_GUE(ℓ) from 2 S S S and Σ²_GUE (`kappa3_tables_out.txt`), ℓ = 0.5 … 8:

    +0.04527 +0.01119 +0.01190 +0.00684 +0.00664 +0.00482 +0.00460 +0.00370 +0.00353 +0.00300 +0.00286 +0.00252 +0.00241 +0.00217 +0.00208 +0.00190

Small, positive, decaying roughly like 0.015/ℓ with a half-integer/integer sawtooth (odd cumulants of the sine-process
count tend to 0: the Basor–Widom/Fisher–Hartwig asymptotics of E e^{iλn} are even in λ beyond the linear term;
*literature*, not re-derived). CUE check (`cue_out_11.txt`, own code, 150 Haar matrices of U(800), periodic windows, step
0.01): the simulated mean agrees with the exact CUE(800) value at all 16 lengths, |z| ≤ 1.4 (se 0.0003–0.0019); CUE(800)
and the sine kernel agree to 10⁻⁵. The data estimator (windows from the first to the last point of a 750-point segment)
against the periodic one, paired per matrix: mean differences ≤ 6·10⁻⁴, all within 2 se (se 1–6·10⁻⁴), so the edge bias
of the estimator is below ≈5·10⁻⁴ (not resolved below that).

## 5. The prediction: κ3 and Δκ3 (computed; `kappa3_tables_out.txt`)

### 5.1 The level, k = 1 (no bad prime), ℓ = 0.5 … 8

    L = 9.4   κ3       +0.0379 +0.0038 -0.0073 -0.0086 -0.0184 -0.0135 -0.0190 -0.0119 -0.0148 -0.0065 -0.0074 +0.0002 +0.0003 +0.0023 +0.0073 +0.0043
              κ3 − GUE -0.0074 -0.0074 -0.0192 -0.0154 -0.0250 -0.0183 -0.0236 -0.0156 -0.0184 -0.0095 -0.0103 -0.0023 -0.0021 +0.0001 +0.0052 +0.0024
    L = 10.2  κ3       +0.0389 +0.0045 -0.0055 -0.0079 -0.0173 -0.0136 -0.0194 -0.0134 -0.0164 -0.0094 -0.0106 -0.0028 -0.0032 +0.0019 +0.0032 +0.0028
    L = 10.7  κ3       +0.0394 +0.0048 -0.0045 -0.0074 -0.0166 -0.0136 -0.0194 -0.0141 -0.0171 -0.0108 -0.0122 -0.0049 -0.0052 +0.0009 +0.0012 +0.0026
    GUE       κ3       +0.0453 +0.0112 +0.0119 +0.0068 +0.0066 +0.0048 +0.0046 +0.0037 +0.0035 +0.0030 +0.0029 +0.0025 +0.0024 +0.0022 +0.0021 +0.0019

Decomposition at L = 9.4 (κ3 − GUE = 3(Σ² − Σ²_GUE) + κ[3]_Q + (κ[3]_P − GUE)):
3(Σ² − Σ²_GUE) = −0.015 … −0.483 (the two-point arithmetic, tested in rounds 134–136); κ[3]_Q = −0.0016 … −0.0195 … +0.0092
(the diagonal triple-orbit term); the off-diagonal three-point part +0.009 … +0.476 cancels the first to within 0.005–0.02.
So the zeros' count is *less* skew than GUE's for 1.5 ≤ ℓ ≤ 5.5 (κ3 < 0 where GUE's is > 0), by 0.01–0.025: a
parameter-free arithmetic effect of the size of the diagonal term, not of the two-point term.

### 5.2 Per-prime change (S = ∅ otherwise): Δκ3 = Δκ[3]_Q (diagonal, closed form §4.2) + rest

    L = 9.4  p=2  Δκ3     -0.0021 -0.0012 +0.0013 +0.0044 +0.0082 +0.0107 +0.0139 +0.0116 +0.0166 +0.0070 +0.0145 +0.0007 +0.0056 -0.0020 -0.0079 -0.0044
                  diag(Q) +0.0003 +0.0019 +0.0048 +0.0085 +0.0120 +0.0147 +0.0163 +0.0166 +0.0157 +0.0136 +0.0105 +0.0067 +0.0025 -0.0019 -0.0062 -0.0101
             p=3  Δκ3     -0.0027 -0.0015 -0.0005 +0.0033 +0.0019 +0.0072 -0.0016 +0.0071 -0.0074 +0.0028 -0.0082 -0.0018 -0.0031 -0.0007 +0.0021 +0.0016
                  diag(Q) +0.0003 +0.0017 +0.0040 +0.0061 +0.0071 +0.0066 +0.0047 +0.0018 -0.0015 -0.0044 -0.0065 -0.0072 -0.0063 -0.0043 -0.0019 -0.0004
             p=5  Δκ3     -0.0029 -0.0017 -0.0022 +0.0019 -0.0028 +0.0035 -0.0037 +0.0027 +0.0001 +0.0027 +0.0014 -0.0007 -0.0021 -0.0019 +0.0004 -0.0020
    L = 10.2 p=2  Δκ3     -0.0018 -0.0011 +0.0009 +0.0036 +0.0072 +0.0099 +0.0127 +0.0125 +0.0160 +0.0098 +0.0159 +0.0039 +0.0109 -0.0009 -0.0000 -0.0028
             p=3  Δκ3     -0.0023 -0.0015 -0.0005 +0.0026 +0.0023 +0.0066 +0.0004 +0.0074 -0.0051 +0.0044 -0.0085 -0.0006 -0.0062 -0.0024 -0.0008 +0.0004
             p=5  Δκ3     -0.0025 -0.0016 -0.0019 +0.0014 -0.0020 +0.0031 -0.0035 +0.0024 -0.0016 +0.0023 +0.0016 +0.0013 -0.0006 -0.0019 -0.0015 -0.0013
    L = 10.7 p=2  Δκ3     -0.0017 -0.0011 +0.0007 +0.0032 +0.0066 +0.0093 +0.0121 +0.0127 +0.0156 +0.0111 +0.0162 +0.0060 +0.0129 +0.0004 +0.0043 -0.0022
             p=3  Δκ3     -0.0021 -0.0014 -0.0006 +0.0023 +0.0024 +0.0062 +0.0014 +0.0074 -0.0037 +0.0052 -0.0080 +0.0005 -0.0074 -0.0026 -0.0029 -0.0011
             p=5  Δκ3     -0.0023 -0.0016 -0.0018 +0.0011 -0.0016 +0.0030 -0.0032 +0.0024 -0.0023 +0.0019 +0.0011 +0.0019 +0.0005 -0.0011 -0.0018 -0.0015

(p = 7 and the parts of every row: `kappa3_tables_out.txt`; the same at L = 8.1 there.)

### 5.3 The structure of Δκ3 (in the form of the repo's ΔΣ²)

    Δκ3(ℓ; p, L, S) = Δκ3_diag(ℓ; p, L) + Δκ3_off(ℓ; p, L, S),
    Δκ3_diag = +(6/π³) Σ_{m,n≥1} sin(mθ) sin(nθ) sin((m+n)θ) / (m n (m+n) p^{m+n}),   θ = πℓ ln p/L      (proved, §4.2)
    Δκ3_off  = 3ΔΣ²(ℓ; p, L, S) + Δκ[3]_P(ℓ; p, L, S)                                                    (computed)

with ΔΣ² = ΔΣ²_diag + ΔΣ²_off of `bk/derivation.md` and Δκ[3]_P the cube integral of the change of the six P-terms (§3.3,
items 2–3). Findings (computed, all four L):
1. **The two-point drop does not reach the skewness.** 3ΔΣ² is −0.03 … −0.33 for p = 2; Δκ[3]_P is +0.03 … +0.33; their
   sum Δκ3_off is −0.007 … +0.006. A missing orbit changes the variance at order (ln²p)/p but the third cumulant only
   at order 1/p² (the diagonal term) plus a small off-diagonal remainder. This is what the random-phase picture says
   (n − ℓ = Σ_p X_p with nearly independent X_p: deleting X_p changes κ2 by −E[X_p²] and κ3 by −E[X_p³] only), and the
   ratios conjecture reproduces it through a cancellation between its two- and three-point parts.
2. **Diagonal dominates for p = 2 at 2 ≲ ℓ ≲ 6**: Δκ3_diag rises to +0.0166 at ℓ ≈ 4 (θ ≈ π/3 for the leading term
   (3/π³) sin²θ sin 2θ/p², = +0.0157), then changes sign near ℓ = L/(2 ln 2) ≈ 7 (L = 9.4). For p = 3 it peaks at +0.0071
   (ℓ ≈ 2.5) and turns negative beyond ℓ ≈ 4.3 (θ = π/2 for 3); for p = 5 it is ±0.0025.
3. **The remainder Δκ3_off** is ≈ −0.002 … −0.005 for ℓ ≲ 2.5 for every p (same sign and size as the "coupling" term
   3 Σ²_GUE′(ℓ) ΔΣ²_diag(ℓ; p): heuristic picture of BK 2013 — the window's content n = ℓ + δN + ξ with δN from the short
   orbits and ξ a GUE count of variance Σ²_GUE(ℓ + δN), so κ3 ⊃ 3Σ²_GUE′(ℓ) Var δN; *heuristic*, agrees in sign and to
   ±40 % in size for ℓ ≤ 2.5, not beyond), and for larger ℓ a period-1 sawtooth of amplitude ≈ 0.005 (half-integer ℓ up,
   whole ℓ down for p = 2 at ℓ ≥ 4.5; the inverse pattern for p = 3; weaker for p = 5) — the e^{2πir} off-diagonal term, like ΔΣ²_off.
4. **Not additive** over the bad primes (P1 is a product); predictions per L-function use the full bad set.
5. **Small ℓ.** As ℓ → 0, κ[3] = 2ℓ³ + O(ℓ⁵) for every k, so Δκ3 → 3ΔΣ², which is O(ℓ⁴) (bk/derivation.md §5.4).

### 5.4 Which alternatives the data can separate (per L-function predictions in `predictions_kappa3.json`)

- **GUE**: κ3 = κ3_GUE, Δκ3 = 0.
- **Q0** ("the conductor does not remove the triple-orbit term"): the full law with Q summed over all primes;
  differs from the full law by −Δκ3_diag per bad prime (−0.0166 for p = 2 at ℓ = 4).
- **two-point only** ("bad primes act through the two-point function only": κ[3] as for k = 1 at the same L, Σ² with the
  bad primes): Δκ3 = 3ΔΣ², i.e. −0.03 … −0.33 for p = 2. Rejected by anything (a strawman, but a sharp one: it is what a
  random-matrix ensemble with the BK two-point function and a GUE three-point function would give).
- **H_diag** ("diagonal only": κ3 of k = 1 at the same L plus Δκ3_diag): bad primes act on the skewness only through the
  triple-orbit term; it differs from the full law by the off-diagonal remainder (−0.007 … +0.006).
Expected values of every reading under each hypothesis: `registry_out.txt` (`registry.py`).

### 5.5 A numerical trap (found and fixed here)

B1(x, y) has poles where 1 − 2/3 + 3^{−1−x} = 0 (x = i(2j+1)π/ln 3), i.e. exactly at the zeros of A_3(x) in P1(x); the
product P1·B1 is finite (B1 = −∂ log A_p-type, so the pole is the log-derivative of A's only zero, the one noted in
`bk/derivation.md` §5.1). Interpolating B1 alone in L (to tabulate it once) broke the cancellation and gave errors up to
10² in κ[3] for some L-functions; taking the p = 3 summand out of the table and adding it exactly fixes it
(`interp_check.py` → `interp_check_out.txt`: cached against direct, ≤ 7·10⁻⁸ in the P-part, 6·10⁻⁷ in the Q-part, 7 cases incl. 3 | k and
large bad primes 89, 373, 383). When 3 | k the zero and the pole are both absent.

## 6. A by-product: κ3 is a completeness detector, and three used sets miss a pair of zeros

A missing zero lowers the count by 1 in a fraction ℓ/U of the windows (U = zeros in the band) and so changes κ3 by
≈ −(ℓ/U)(3Σ² + 1) ≈ −0.02 at ℓ = 8, five to seven times the per-L-function scatter; a missing pair by ≈ −(ℓ/U)(6Σ² + 8).
The noise run (`noise_out.txt`, residuals only) showed three such outliers, all with negative residuals growing with ℓ:

| set | L-function | κ3 residual at ℓ = 2, 4, 6, 8 | largest unfolded gap | N(t) − θ(t)/π at the band's top |
|---|---|---|---|---|
| C (600–1000, round 136) | D = −195 | −0.011, −0.029, −0.060, −0.084 | 3.22 at t = 984.4 | −2.09 at 1000 |
| Dn (100–600, round 136) | D = −939 | −0.005, −0.012, −0.010, −0.013 | 3.22 at t = 597.7 | −1.62 at 600 |
| A (complex, round 121/134) | q = 173 (m2), side −1 | −0.010, −0.041, −0.061, −0.065 | 4.30 at t = 488.5 | (not computed) |

(`gapcheck.py`; controls in C have N − θ/π within ±1.3 and gaps ≤ 2.5.) So each of these sets misses (most likely) a close
pair of zeros near the top of its band — the PARI subdivision problem found in round 136 for D = −283, here *inside* the
bands, where round 136's band-mean check (|mean| > 0.5) could not see them because the deficit covers only the last few
per cent of the band. Their effect on round 136's Σ² readings was not assessed here (a missing pair raises Σ²(8) of that
one L-function by roughly 4ℓ/U ≈ 0.05); the owner may want to re-run lfunzeros on these three windows with a finer step.
The registration below therefore requires a pointwise completeness check (N − θ/π on a grid *and* at the band's ends,
and the largest unfolded gap) before any statistic.

## 7. Noise (NOISE ONLY; `noise_clean_out.txt`)

Per L-function, κ3(ℓ) was computed on the sets used in rounds 120–136, regressed on [1, I(2|D), …, I(13|D)] (complex
set: 2, 3, 5, 7), and only the residuals were kept; no mean or coefficient was printed or compared with anything. With the
three incomplete sets of §6 dropped:

| sample | n | zeros/set | residual sd of κ3 at ℓ = 0.5, 1, 1.5, 2, 3, 4, 4.5, 6, 8 | n_eff of 16 lengths |
|---|---|---|---|---|
| R120 | 181 | 745 | 0.0036 0.0015 0.0033 0.0026 0.0031 0.0039 0.0044 0.0028 0.0030 | 13.5 |
| B | 184 | 813 | 0.0038 0.0017 0.0037 0.0028 0.0033 0.0037 0.0040 0.0035 0.0028 | 13.6 |
| Dn | 180 | 848 | 0.0036 0.0017 0.0038 0.0030 0.0036 0.0040 0.0046 0.0035 0.0033 | 13.5 |
| C | 180 | 654 | 0.0038 0.0019 0.0048 0.0036 0.0043 0.0046 0.0057 0.0040 0.0031 | 13.3 |
| A (complex) | 87 | 1172 (two sides) | 0.0031 0.0013 0.0033 0.0020 0.0030 0.0028 0.0035 0.0026 0.0023 | 12.5 |

The residuals are nearly uncorrelated across ℓ (|r| ≲ 0.3), with a sawtooth: half-integer ℓ noisier than whole ℓ. These
sds include the predicted between-L-function spread (L, other primes), so they are upper bounds on the measurement
noise. For comparison, the sd of the factorial cumulant κ[3] is 0.011–0.032 (it carries 3× the noise of Σ², whose sd is
0.004–0.009), and the GUE (CUE) sd of the same κ3 estimator on 750 points is 0.004–0.025, 3–8 times the zeros' — the
zeros' long-range count fluctuations are carried by a few short orbits and average out fast. Spacing statistics, for the
record (residual sd): var(s) 0.004, third moment of s 0.0025, lag-1 correlation of spacings 0.010–0.014.

## 8. Power (expected values from the frozen predictions, noise from §7; `power3_out.txt`)

Noise model: per L-function, the residual covariance across ℓ of a reference sample (Dn for Fn, C for Bp, own for the
used sets), scaled by zeros. Expected z = predicted difference / expected se; "combined" = GLS over the 16 lengths with
the residual covariance.

| reading | 546 pooled (R120+B+Dn) | Fn (fresh, 186) | Bp (fresh, 184) | C (181) |
|---|---|---|---|---|
| T1 level: se of the mean κ3(ℓ) | 0.00007–0.00021 | 0.00012–0.00035 | 0.00013–0.00041 | 0.00014–0.00043 |
| T1 full − GUE, z per ℓ (range) / combined | −121 … 0 / 332 | −66 … −2 / 198 | −60 … −2 / 187 | −60 … +1 / 175 |
| T1 full − Q0 (triple-orbit term kept), combined z | 107 | 65 | 60 | 54 |
| T2 γ_A (arithmetic amplitude of the level): se | 0.0030 | 0.0050 | 0.0054 | 0.0057 |
| T3 p = 2 coefficients: vs GUE, max cell z / combined | 46 / 104 | 24 / 68 | 22.5 / 61 | 21 / 54 |
| T3 p = 2 vs Q0 | 47 / 132 | 29 / 84 | 25 / 80 | 23 / 71 |
| T3 p = 3 vs GUE / vs Q0 (combined) | 38 / 60 | 25 / 38 | 21 / 29 | 20 / 27 |
| T3 p = 5 vs GUE / vs Q0 (combined; max cell vs GUE) | 16 / 19 (6.9) | 9.6 / 11.5 (4.2) | 8.4 / 10.0 (3.7) | 8.3 / 10.7 (3.6) |
| T4 γ_Q (triple-orbit amplitude, cells p = 2, 3, 5): se | 0.0077 | 0.0132 | 0.0149 | 0.0161 |
| T4 γ_off (off-diagonal bad-prime remainder, cells p = 2, 3, 5): se (H_diag at) | 0.0141 (71 se) | 0.0265 (38) | 0.0291 (34) | 0.0290 (34) |
| T4 γ_2pt: se (two-point-only at) | 0.0003 (3,173 se) | 0.0006 (1,747) | 0.0006 (1,618) | 0.0006 (1,559) |
| T5 κ[3] level vs GUE, z per ℓ | 61 … 750 | 32 … 565 | 30 … 348 | 33 … 351 |

So **every registered reading reaches far beyond 4 standard errors** on a single fresh set; the weakest is p = 5 per cell
(≈ 4 at its best length; ≈ 10 combined). The decisive bad-prime-specific reading, the amplitude of the triple-orbit
term in the contrasts (γ_Q), is measured to ±0.013 on Fn: the alternative "the conductor does not remove the
triple-orbit term" sits 76 se away, and the off-diagonal remainder (the bad-prime factors of P1 and P3 inside
κ3) has its own amplitude γ_off, ±0.027 on Fn, with "diagonal only" (H_diag) 38 se away. GUE is excluded at ≈ 200 se and "bad primes act through the two-point function only"
at ≈ 2,000 se. The binding limitation is not statistics but the theory's own finite-height accuracy: the level is read to
2–3·10⁻⁴, below the 4–7·10⁻⁴ to which the two-point level is known to hold (doubt 1). A level offset of 10⁻³ would fail
T1 at 3–8 se and move γ_A by ≈ 0.06; it would not move γ_Q (common to all L-functions) unless it depended on the bad
primes.

No alternative statistic is needed. For the record, the noise of the nearest three-point alternatives was measured
(§7): skewness-type statistics of the spacings (sd 0.0025 for the third moment, 0.010–0.014 for the lag-1 spacing
correlation per L-function); they were not developed because κ3 already has the power.

**Recommendation: run** — test Fn (all 186 fundamental discriminants 1000 < |D| ≤ 1300, heights 100–600, ≈ 25 min of
PARI) as the primary fresh test; Bp optional; the 546 + 181 sets on disk as test O (zeros known, κ3 never computed),
after repairing the three incomplete sets. Registration draft: `round_draft.txt`.

## 9. What is in print (literature; read as stated in §1)

- **ζ, three-point, with lower-order terms.** Formula: Bogomolny–Keating (semiclassical outline 1996; full form 2013,
  random-matrix-universality route) and Conrey–Snaith (ratios conjecture, 2007/2008); the two agree (CS 2007 say so; here
  checked to 10⁻⁷ including bad primes, §3.4). **Numerics: Conrey–Snaith 2007** compare their formula with the density
  of triples of the first 100,000 zeros (T ≈ 75,000): difference plot max 0.164, mean −0.00127, sd 0.03 (their
  normalisation), "probably due to the relatively small value of T". BK 2013 has no numerics. Leading-order (GUE)
  three-point: Hejhal 1994 (proved, restricted support), Rudnick–Sarnak 1996 (n-point, principal L-functions; not opened
  here — whether it shows numerical three-point data was not checked). Odlyzko's spacing and next-nearest-spacing
  distributions (1987–2001) depend on all n-point functions but are not a three-point test with arithmetic terms.
  Fazzari–Gerspach (2024, abstract) *use* a triple-correlation conjecture to get the third moment of log ζ; no numerics.
  Forrester–Mays (2015, abstract): finite-size corrections to two-point/spacing statistics of Odlyzko's zeros; no
  three-point statistic.
- **Dirichlet L-functions, three-point.** Nothing found: the two-point function with bad primes is BK 2013 (no numerics);
  two targeted searches for numerical triple correlation of Dirichlet L-zeros, or any conductor dependence of a
  three-point statistic, returned only ζ (CS 2007) and family statistics of low-lying zeros (one-level density, lower-
  order terms: e.g. Fiorilli–Miller) — a different object (average over a family near s = 1/2, not one L-function high up).
- **The third cumulant of the zero count.** Not found as a tested statistic, for ζ or for L-functions. Its decomposition
  κ3 = κ[3] + 3Σ² − 2ℓ and the GUE values are standard random-matrix material (the GUE values here are computed and
  CUE-checked, §4.3); the closed form of the diagonal part (§4.2), the per-prime decomposition (§5.3), the cancellation of
  3ΔΣ² inside κ3, and the conductor-dependence predictions are derived here.
- So, as far as found: **a test of the bad-prime ("missing orbit") factors of the three-point function, or of any
  three-point statistic of a Dirichlet L-function against its arithmetic lower-order terms, is not in print.** The search
  was two web queries plus the reference lists of CS 2007/2008 and BK 2013; it is not exhaustive (Rubinstein's thesis and
  Rudnick–Sarnak's figures not opened).

## 10. Doubts

1. **The law is a conjecture** (ratios recipe; BK's heuristic), with error O(T^{1/2+ε}) of unknown constant. At these
   heights the two-point level held to 4–7·10⁻⁴ (rounds 134, 136). For κ3 the analogous finite-height error is unknown;
   the pooled level is measured to ≈ 2–3·10⁻⁴ per length, so a 10⁻³ level offset would be visible and would not, by
   itself, contradict the arithmetic structure. The registration therefore separates the level (T1) from the structural
   readings (T2–T4) and states the theory's allowance separately.
2. **Normalisation and signs** are fixed by three independent checks (§3.4): exact GUE limit, R3 = 0 at coincidences
   against the repo's confirmed two-point law (this ties the phase convention e^{+iεL}·Φ_off(ε) of round 134 to the
   three-point terms), and the second route (BK 2013) to 10⁻⁷. The rule "p | k deleted" is derived in both routes, not
   checked against data.
3. **Convergence of the triple integrals.** The integrand is regular (no poles, CS 2008 Thm 5), but the six terms have
   poles that cancel; the GUE skeleton is integrated exactly and the remainder converges to 5·10⁻¹⁰ between grids. The
   pole of B1 at the zero of A_3 must be kept with it (§5.5) — a trap that a re-implementation must avoid.
4. **Truncations**: B1 over p ≤ 10⁴ (10⁻⁶), Q-part closed form over p ≤ 10⁶, terms p^{−(m+n)} < 10⁻¹⁷ dropped; Φ_off etc.
   from the repo's frozen grid (p < 10⁶); height average with 4 Gauss nodes in t (density weight), L-tables at step 0.1
   with cubic interpolation (checked ≤ 7·10⁻⁸).
5. **Estimator.** κ3 is the third central moment of overlapping window counts with the sample mean; the CUE paired test
   bounds the edge/estimator bias at ≲ 5·10⁻⁴ (not resolved below). The draft adds an edge-trimmed variant as a check.
6. **Completeness.** κ3 reacts to a missing zero ≈ 5–7 per-L-function sds at ℓ = 8 (§6); three used sets are incomplete.
   The draft requires a pointwise completeness check before any statistic, which is stricter than round 136's.
7. **Noise model.** Residual sds include predicted between-L-function spread, so the power figures are conservative; for
   Bp the reference is C (heights 600–1000, 654 zeros/set) scaled to 709 zeros.
8. **ℓ = 0.5** is nearly a function of Σ²(0.5) (correlation 1.00 of the residuals); it adds little three-point content.
9. **Real and complex characters** give the same law (only |χ|² enters); the complex sample was used for noise only.
10. **The coupling picture (§5.3, item 3)** is a heuristic reading of part of the remainder, not part of the prediction.
