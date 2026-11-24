# The within part of block transition tables of consecutive primes

> **In short.** For blocks of length H near x (λ = ln x), the within part of the block-homogeneity χ² of the table (p_n mod q, p_(n+k) mod q) is exactly 1 + (2/(C−1)²) Σ_j (1 − j/N) τ_j, where τ_j is the lag-j covariance of the within-projected transition increments (C = φ(q), N = H/λ). (1) For k = 1 the lag-1 term is fixed by triples of consecutive primes. At leading order it is Lemke Oliver and Soundararajan's repulsion between p_n and p_(n+2), and that gives **E[within]/dof = 1 − 1/((φ(q) − 1)λ)**. (2) For k = 2 the same term enters at lag 2 with weight 1/3. It is cancelled by the square of the k = 1 repeat bias B_q, which is why the k = 2 value is about 1: **1 − 1/(3(φ − 1)λ) + 2B_q²/(φ − 1)²**. (3) The next correction is a long-range tail of second order in 1/λ that grows with ln H. Its source is the Hardy–Littlewood pair correlations acting on the insertion responses of the table. Its slope is −V/((C − 1)²λ) per unit of ln H. The finite-λ Hardy–Littlewood evaluation (no parameter fitted) matches all 15 measured cells with |z| ≤ 1.4 (retrospective). The leading order alone misses q = 7 by 4σ.

Status labels: **proved** (algebra, no conjecture), **conjecture (HL)** (follows from the Hardy–Littlewood k-tuple conjecture, through Lemke Oliver–Soundararajan's (LOS) expansion where stated), **heuristic** (a modelling or truncation step), **computed** (a number from a script here), **retrospective** (compared with values known before the law was written).

Files (all in this folder): `law.py` → `law_output.txt` (assembly; retrospective table; predictions). `components.py` writes the Hardy–Littlewood components per height (`data/comp_*.json`). `hlcore.py` holds the singular series and the conditioned neighbour laws. `triple.py` holds the triple law and the exact lag-1 algebra. `chainlaw.py` holds the fine-class chain, the k = 2 covariances and the insertion responses. `tailmod.py` holds the long-range tail. `modelB.py` → `modelB_output.txt`. Validation on [1.5·10¹¹, 1.5·10¹¹ + 10⁸) only: `sieve_val.py`, `measure_val.py`, `tail_val.py` → `data/val_*.json`, `validation_report.py` → `validation_output.txt`, and `check_primes.py` → `check_primes_output.txt` (sympy/gmpy2 check of the sieve). `reg_numbers.py` → `reg_numbers_output.txt` holds the errors and powers of the registration. Nothing is imported from the repo. The measured cells are copied into `law.py` from `writeup/data/study/r135c.json` and `r135b_*.json`.

## 1. The statistic (as in r135b.py)

- Each window [x₀, x₀ + 10⁸) is cut into blocks β of length H.
- c_n is the class index of p_n mod q (C = φ(q) classes).
- X_ab,β = #{n : p_n ∈ β, c_n = a, c_(n+k) = b}, and N_a,β = Σ_b X_ab,β.
- P̂ is the window's pooled transition matrix.
- χ²_β = Σ_ab (X − N_a P̂_ab)²/(N_a P̂_ab).
- col_β = cᵀM⁺c, with c_b = Σ_a (X_ab − N_a P̂_ab) and M = Σ_a N_a (diag P̂_a − P̂_a P̂_aᵀ).
- within_β = χ²_β − col_β. It is reported divided by (C − 1)² and multiplied by NB/(NB − 1).

## 2. Exact reduction to lag covariances (proved, up to O(1/N) edge and ratio terms)

**2.1 Increments.** Let U_n(a, b) = 1[c_n = a, c_(n+k) = b] − 1[c_n = a] P_ab, with P the true transition matrix and π its stationary law. Write y_n(a, b) = U_n(a, b)/√(π_a P_ab) ∈ R^(C²) and Y = Σ_(n∈β) y_n.
- Then χ² = |Y|²/N up to replacing N_a by its mean, a ratio correction of O(1/N).
- Each row of y_n is orthogonal to √P_a, so Y lies in the range of R = ⊕_a (I − √P_a √P_aᵀ). This R is the row-multinomial covariance and has rank C(C − 1).
- The column sums are linear functionals VᵀY with V_(ab),b′ = δ_bb′ √(π_a P_ab), and VᵀRV = M/N.

Hence col = YᵀRV(VᵀRV)⁺VᵀRY/N and **within = YᵀΠ_W Y/N with Π_W = R − RV(VᵀRV)⁺VᵀR**, an orthogonal projection of rank (C − 1)².

**2.2 Expectation.** For a stationary sequence of N terms:

  **E[within]/(C − 1)² = 1 + (2/(C − 1)²) Σ_(j≥1) (1 − j/N) τ_j,  τ_j = E[y_nᵀ Π_W y_(n+j)] = tr(Π_W Γ_j).**   (2.1)

Here tr(Π_W Γ₀) = (C − 1)² when P is the true matrix. The finite-N corrections are the (1 − j/N) weights, which are exact, and O(1/N) terms from the random N_a and the block edges, which are not modelled (doubt 5). The factor NB/(NB − 1) of the statistic undoes the estimation of P̂.

**2.3 First-order chains give 1 (proved).** For k = 1 and a Markov chain on the classes, E[U_(n+j) | c_0..c_(n+j)] = 0. So the y_n are martingale differences, every τ_j with j ≥ 1 vanishes, and E[within]/dof = 1 (Anderson–Goodman). A departure from 1 measures the classes' memory beyond one step.

**2.4 The lag-1 term (proved).** Let D3(a, a′, b′) = P(c_n = a, c_(n+1) = a′, c_(n+2) = b′) and D2(a, a′) = π_a P_aa′. For k = 1, expanding the four products gives

  E[U_n(a, b) U_(n+1)(a′, b′)] = (δ_(a′b) − P_ab) Δ(a, a′, b′),  Δ(a, a′, b′) = D3(a, a′, b′) − D2(a, a′) P(b′ | a′).   (2.2)

So τ₁ is a fixed linear functional of the second-order Markov deviation Δ. Check (computed, `triple.py`): on the validation stretch, (2.2) applied to the measured triple frequencies gives τ₁/(C − 1)² = −0.02153, −0.00516, −0.00309 for q = 3, 5, 7. The directly measured lag-1 covariances are the same to 10⁻⁵.

**2.5 q = 3 in gap language (proved, up to O(δ²/N) with δ = N₁ − N/2, taking P₁₂ = P₂₁ = α).** With C = 2 the within direction is one-dimensional:

  within = (r₁N₂ − r₂N₁)²/(α(1 − α)N₁N₂N),  r₁ = X₁₁ − N₁P₁₁,  r₂ = X₂₁ − N₂P₂₁.

For N₁ = N/2 + δ this is (E − (1 − α)N)²/(Nα(1 − α)) + O(δ²/N). Here E = X₁₁ + X₂₂ counts the pairs p_(n+1) ≡ p_n (mod 3), that is, the gaps ≡ 0 (mod 6), and 1 − α = P its share. The covariance with the class contrast N₁ − N₂ vanishes by the reflection symmetry a → −a of the Hardy–Littlewood law. Write χ_n = ±1 for the class (the Legendre symbol mod 3), s_n = 1[g_n ≡ 0 mod 6] = (1 + χ_nχ_(n+1))/2 and b_j = E[χ_nχ_(n+j)]. Then

  Var(E − PN | class counts) = N P(1 − P) [1 + 2 Σ_j (1 − j/N) ρ_j],
  ρ₁ = (b₂ − b₁²)/(1 − b₁²)  (exactly, since χ² = 1),   ρ_j = (E[χ_nχ_(n+1)χ_(n+j)χ_(n+j+1)] − b₁²)/(1 − b₁²)  (j ≥ 2).   (2.3)

Read as pairs of gap events: the joint probability of two adjacent gaps ≡ 0 mod 6 is P(s_n = s_(n+1) = 1) = (1 + 2b₁ + b₂)/4, so Cov = (b₂ − b₁²)/4. It needs the law of three consecutive primes. The measured b₁ is the LOS repulsion between p_n and p_(n+1); it enters only through its square, and a chain has b₂ = b₁². **The operative quantity is b₂ − b₁²: the classes of p_n and p_(n+2) coincide less often than a chain allows.** On the validation stretch: b₁ = −0.08736, b₂ = −0.01373 (each ±0.0005), so ρ₁ = −0.02153 (computed).

## 3. Leading order from Hardy–Littlewood (conjecture HL, through LOS)

**3.1 The triple law at order 1/λ.** LOS's Main Conjecture (arXiv 1603.03720, p. 2–3; derivation in their §2 and §4), for r consecutive primes in classes a = (a₁, …, a_r), reads

  π(x; q, a) = li(x)/φ^r · [1 + c₁ lnln x/ln x + c₂/ln x + O((ln x)^(−7/4))],
  c₁(q; a) = (φ/2)((r − 1)/φ − #{i : a_i = a_(i+1)}),
  c₂(q; a) = Σ_i c₂(q; (a_i, a_(i+1))) + (φ/2) Σ_(j=1)^(r−2) (1/j) ((r − 1 − j)/φ − #{i : a_i = a_(i+j+1)}).

- c₁ is additive over adjacent pairs, so it is Markov at order lnλ/λ and gives no lag-1 covariance.
- The only non-Markov piece at order 1/λ is the term with j ≥ 1.
- For r = 3 it is (1/2)(1 − φδ_(a₁a₃)). Hence P(c | a, b) = P(c | b)[1 + (1/(2λ))(1 − φδ_ac)] + O(λ^(−7/4)).
- More generally, primes d indices apart interact through w_d(a, a′) = (1 − φδ_(aa′))/(2(d − 1)λ) (LOS Conjecture 1.4).

**Where it comes from** (LOS §4, D₀ term). The singular-series pair term S_(q,0)({0, h_(i+1) + … + h_j}) between non-adjacent primes enters through the sums S^k₀(q, 0; H) = Σ_(h≡0) h^k e^(−h/H) S_(q,0)(h) = −(φ/(2q)) Γ(k) H^k + O(H^(k−1/2)). This is the 1/h tail of the singular series on pairs in the same class. Summed over all index distances, the same tail gives the class-count rigidity −ln H/λ of round 135 (noise/derivation.md §3–4): b_j ≈ −1/(2(j − 1)λ) for q = 3. **So the lag-1 within deficit is the first link of the 1/h tail that makes class counts rigid. It is not the adjacent LOS bias, which only enters squared.**

**3.2 k = 1.** Take orthonormal mean-zero functions f_i, g_l on the classes. At leading order the within directions are f_i ⊗ g_l and each carries the mark m_n = f(c_n)[g(c_(n+1)) − (Pg)(c_n)], normalised by φ². With the law above,

  E[g(c_(n+2)) − (Pg)(c_(n+1)) | c_n, c_(n+1)] = (1/φ) Σ_c (1/(2λ))(1 − φδ_(c_n c)) g(c) = −g(c_n)/(2λ),

so E[m_n m_(n+1)] = −⟨f, g⟩²/(2λφ²) + O(λ^(−2) ln²λ). Using Σ_(i,l) ⟨f_i, g_l⟩² = C − 1 gives tr(Π_W Γ₁) = −(C − 1)/(2λ). For j ≥ 2 the four indices n, n + 1, n + j, n + j + 1 are distinct, and a first-order pair interaction cannot couple four mean-zero factors, so τ_j = O(λ^(−2)). Hence

  **E[within]/dof (k = 1) = 1 − 1/((φ(q) − 1)λ) + O(λ^(−2) ln²λ)**   (conjecture HL; LOS).   (3.1)

For q = 3: Var(#gaps ≡ 0 mod 6 | class counts) = N P(1 − P)(1 − 1/λ + …).

**3.3 k = 2.**
- Lag 1 involves the four distinct indices n, …, n + 3, so it is zero at first order.
- At lag 2 the shared class c_(n+2) leaves the pair (c_n, c_(n+4)): w₄ = (1 − φδ)/(6λ), so tr(Π_W Γ₂) = −(C − 1)/(6λ).
- The leading second-order term is the lag-1 pairing (n, n + 1)(n + 2, n + 3) of the Markov biases. It contributes (tr β̃/φ)² = B_q², where B_q = φ P(p_(n+1) ≡ p_n) − 1 is the repeat bias.

So

  **E[within]/dof (k = 2) = 1 − 1/(3(φ − 1)λ) + 2B_q²/(φ − 1)² + …,  B_q = [−((φ − 1)/2) ln λ + c₂(q; (a, a))]/λ**   (3.2)

with c₂(q; (a, a)) from LOS (1.1). For general k the non-Markov term is −1/((2k − 1)(φ − 1)λ), from the pair (n, n + 2k). **Why k = 2 gives about 1:** for q = 3 at λ = 27, −1/(3λ) = −0.0123 is cancelled by 2B₃² = +0.011. For q = 10 the values are −0.0041 and +0.0085.

**3.4 Numbers (computed, `law_output.txt` §1), at the mean height of the nine windows (λ = 27.02) and of the registration windows (λ = 22.69):**

| q | k = 1: 1 − 1/((φ − 1)λ) | k = 2: (3.2) |
|---|---|---|
| 3 | 0.9630 / 0.9559 | 0.9988 / 0.9998 |
| 5 (= 10) | 0.9877 / 0.9853 | 1.0044 / 1.0060 |
| 7 | 0.9926 / 0.9912 | 1.0045 / 1.0059 |

## 4. The next correction: finite-λ Hardy–Littlewood evaluation

The leading order misses (7, 1, 10⁴) by 4σ (0.9926 against 0.9859). The second-order terms are as large as the first for q = 5 and 7, so they are evaluated without expanding in 1/λ.

**4.1 Conditioned Hardy–Littlewood law, order 1 (conjecture HL + heuristic truncation).**
- Write Q = ∏_(p≤y) p (y = 13, so Q = 30030) and condition on x mod Q.
- By the Chinese remainder theorem the small primes are exact: a site x + t is either impossible or "rough".
- A rough site is prime with chance ρ₀ S_(>y)(t₁) S_(>y)(t₂)…, the Hardy–Littlewood pair factors with every earlier point of the pattern, where ρ₀ = (Q/φ(Q))/λ and S_(>y)(h) = C_(>y) ∏_(p>y, p|h) (p − 1)/(p − 2).
- The rough sites are taken as independent given the pattern. The interior pair corrections (order 2) and the triple-coincidence parts of the 3-point singular series are dropped.
- This is the LOS conditioning used in the repo's rounds 47–55.

Check (computed): the repeat biases B_q of this law are −0.08748, −0.22722, −0.32753 for q = 3, 5, 7. The primes give −0.08736, −0.22746, −0.32814 (λ = 25.73). With y = 7, 11 or 13 the values agree to 0.001.

**4.2 Lag 1 from the triple law (computed, `triple.py`).**
- For every rough class u mod Q, the joint law of the first and second neighbour is built with the pair factors of the second neighbour and of its interior sites to both x and x + t₁.
- Applying (2.2) gives τ₁/(C − 1)² = −0.02155, −0.00511, −0.00286 against −0.02153, −0.00516, −0.00309 measured (±0.0005, ±0.0002, ±0.0001).
- At y = 11 the values are −0.02112, −0.00505, −0.00285.
- A fine-class chain with pair factors only to the predecessor gives −0.0153 for q = 3. **About 30 % of the lag-1 term comes from the primes > 13 between p_n and p_(n+2)**, the S_(>y) part of the 1/h tail.
- Model B is LOS's first-order tensor on top of the actual P (`modelB.py`). It gives −0.0193, −0.0053, −0.0027. So the step from the closed form (−1/(2(C − 1)λ) = −0.0194, −0.0065, −0.0039) to the finite-λ values comes mainly from the non-uniform transition matrix for q = 5 and 7, and from finite-λ terms of the triple law for q = 3.

**4.3 Lags ≥ 2 (computed; heuristic: no S_(>y) memory beyond one step).**
- τ_j for j ≥ 2 comes from the Markov chain on the rough classes mod 30030, with kernel K(u → u + t) = the next-prime law of 4.1.
- The covariances are exact sums over powers of K (no simulation): τ_j = (1/φ(Q)) Σ_u Σ_(u′) K(u, u′)⟨F(u, u′), (K^(j−1)G)(u′)⟩, where G(u) = E[F | u].
- Validation, lag 2: +0.0015, +0.0006, +0.0004 against +0.0023, +0.0008, +0.0004 measured.

**4.4 The long-range tail (conjecture HL + heuristic linear response).**
- Inserting a prime at a rough site x ≡ u changes the table by y(c_L, c_x) + y(c_x, c_R) − y(c_L, c_R), where c_L and c_R are the classes of the neighbours. Its expectation r_u is computed from the per-class neighbour laws P_L(u, ·) and P_R(u, ·) of 4.1.
- To first order, the table is Ȳ + Σ_x r_x (1_x − E 1_x).
- The primes > y correlate as Cov(1_x, 1_(x+h)) = (Q/φ(Q))²(S_(>y)(h) − 1)/λ² for rough x and x + h. So

  W_tail(H) = (2/((C − 1)²λ)) (Q/φ(Q)) Σ_(h<H) (1 − h/H) w₂(h) A_W(h mod Q) (S_(>y)(h) − 1),  A_W(h) = mean_u ⟨Π_W r_u, Π_W r_(u+h)⟩,   (4.1)

  where w₂(h) = 1 − e^(−h/λ) removes adjacent pairs, which the triple law already holds (heuristic). A_W is computed by FFT over Z/Q, and the sum is exact over h < H.
- By §3 of noise/derivation.md, S_(>y) − 1 restricted to h ≡ 0 mod Q carries −(φ(Q)/Q) ln H, and the other classes carry only constants. So

  **dW/d ln H ≈ −V/((C − 1)²λ),  V = A_W(0) = mean_u |Π_W r_u|²**   (4.2)

  (for modes below H). V is proportional to ρ₀, so the tail is of order 1/λ² (times (Q/φ(Q))²), not 1/λ.
- The first-order argument of §3.2 explains why no 1/λ term appears: a uniform or class-wise density change projects onto the column part. V/(C − 1)² at λ = 27.12 is 0.036, 0.071, 0.093 for q = 3, 5, 7, so the slopes are −0.003, −0.006, −0.008 per decade.

**4.5 Short against long range (computed, λ = 27.12; contributions to E[within]/dof):**

| q | lag 1 (triple law) | lags ≥ 2 (chain) | tail at H = 10⁴ / 10⁵ / 10⁶ | total at 10⁴ |
|---|---|---|---|---|
| 3 | −0.0403 | +0.0021 | −0.0025 / −0.0052 / −0.0083 | 0.959 |
| 5 | −0.0097 | −0.0015 | −0.0045 / −0.0101 / −0.0160 | 0.984 |
| 7 | −0.0054 | −0.0027 | −0.0051 / −0.0124 / −0.0203 | 0.987 |

So **0.955 for q = 3 is the short-range lag-1 term**: the repulsion of p_n and p_(n+2), with the chain's lag 2 and the tail nearly cancelling at 10⁴. For q = 5 and 7 the long-range tail is already as large as the lag-1 term at 10⁴, and it makes them decline with ln H. For q = 3 the formula predicts a decline of −0.003 per decade, which is invisible under errors of ±0.014 and ±0.044.

**4.6 k = 2 (computed; heuristic).**
- Lag 1 comes from the chain: +0.0091, +0.0062, +0.0044 against +0.0098, +0.0060, +0.0044 measured.
- Lag 2 comes from the LOS leading term −1/(6(C − 1)λ): −0.0065, −0.0022, −0.0013 against −0.0055, −0.0023, −0.0016 measured. The chain alone gives −0.0007, −0.0013, −0.0010, because it has no memory beyond one step.
- The tail uses the k = 2 insertion response, which needs the second neighbours on each side.

**4.7 Validation of the tail (computed, `validation_output.txt`).**
- On the validation stretch, the lag-bin sums of τ_j/(C − 1)² for q = 7 over lags [11, 31), [31, 101), [101, 301) are −0.0019, −0.0029, −0.0034 (±0.0005, 0.0005, 0.0013). Chain plus tail gives −0.0013, −0.0016, −0.0018.
- The tail is real and of the predicted sign and order, but the model's tail is 50–70 % of the measured one at intermediate lags.
- For q = 5 the measured long-lag bins are larger still and noisy: [1001, 3001) gives −0.019 ± 0.006 against −0.002.

## 5. Retrospective table (the measured values were known; nine windows 2·10¹¹–10¹², `law_output.txt` §3)

err = the larger of the pooled block error and the window-scatter error. z uses err ⊕ the law's own error. The law's error combines truncation (y = 13 against y = 11), the w₂/w₃ choice, 40 % of the tail, and, for k = 2, 0.002 for the LOS lag-2 term. "Sieve chain" is the chain on rough classes mod 30030 with independent rough sites (S_(>y) ≡ 1) and no long-range part.

| cell (q, k, H) | measured ± err | law ± law err | z | leading order (z) | sieve chain | 1 (z) |
|---|---|---|---|---|---|---|
| (3, 1, 10⁴) | 0.9549 ± 0.0062 | 0.9592 ± 0.0014 | −0.7 | 0.9630 (−1.3) | 0.9736 | −7.3 |
| (3, 1, 10⁵) | 0.9562 ± 0.0144 | 0.9563 ± 0.0028 | −0.0 | 0.9630 (−0.5) | 0.9735 | −3.0 |
| (3, 1, 10⁶) | 0.9611 ± 0.0442 | 0.9532 ± 0.0044 | +0.2 | 0.9630 (−0.0) | 0.9735 | −0.9 |
| (5, 1, 10⁴) | 0.9844 ± 0.0019 | 0.9843 ± 0.0018 | +0.1 | 0.9877 (−1.7) | 0.9903 | −8.2 |
| (5, 1, 10⁵) | 0.9728 ± 0.0048 | 0.9785 ± 0.0044 | −0.9 | 0.9877 (−3.1) | 0.9902 | −5.7 |
| (5, 1, 10⁶) | 0.9486 ± 0.0152 | 0.9724 ± 0.0073 | −1.4 | 0.9877 (−2.6) | 0.9902 | −3.4 |
| (7, 1, 10⁴) | 0.9859 ± 0.0016 | 0.9868 ± 0.0021 | −0.4 | 0.9926 (−4.2) | 0.9923 | −8.8 |
| (7, 1, 10⁵) | 0.9798 ± 0.0029 | 0.9792 ± 0.0051 | +0.1 | 0.9926 (−4.4) | 0.9920 | −7.0 |
| (7, 1, 10⁶) | 0.9735 ± 0.0089 | 0.9713 ± 0.0087 | +0.2 | 0.9926 (−2.1) | 0.9920 | −3.0 |
| (3, 2, 10⁴) | 1.0047 ± 0.0047 | 1.0033 ± 0.0020 | +0.3 | 0.9988 (+1.3) | 1.0142 | +1.0 |
| (3, 2, 10⁵) | 1.0174 ± 0.0157 | 1.0023 ± 0.0021 | +1.0 | 0.9988 (+1.2) | 1.0142 | +1.1 |
| (3, 2, 10⁶) | 0.9555 ± 0.0552 | 1.0011 ± 0.0023 | −0.8 | 0.9988 (−0.8) | 1.0142 | −0.8 |
| (10, 2, 10⁴) | 1.0060 ± 0.0016 | 1.0063 ± 0.0020 | −0.1 | 1.0044 (+1.0) | 1.0089 | +3.8 |
| (10, 2, 10⁵) | 0.9980 ± 0.0049 | 1.0048 ± 0.0022 | −1.3 | 1.0044 (−1.3) | 1.0089 | −0.4 |
| (10, 2, 10⁶) | 0.9938 ± 0.0160 | 1.0032 ± 0.0026 | −0.6 | 1.0044 (−0.7) | 1.0089 | −0.4 |

Slopes per decade (10⁴→10⁵ / 10⁵→10⁶):

| table | law | measured |
|---|---|---|
| (3, 1) | −0.003 / −0.003 | +0.001 / +0.005 |
| (5, 1) | −0.006 / −0.006 | −0.012 / −0.024 |
| (7, 1) | −0.008 / −0.008 | −0.006 / −0.006 |

The law fits all 15 cells (|z| ≤ 1.4). The leading order alone fails q = 7 (−4.2σ and −4.4σ) and q = 5 at 10⁵ (−3.1σ). The sieve chain fails at 10⁴ for q = 3, 5, 7: (0.9549 − 0.9736)/0.0062 = −3.0, (0.9844 − 0.9903)/0.0019 = −3.1 and (0.9859 − 0.9923)/0.0016 = −4.0. q = 5 declines about twice as fast as the law (−1.4σ at 10⁶).

## 6. What is classical, what is in the repo, what is new

- **Classical / literature.**
  - The Hardy–Littlewood conjecture.
  - LOS's Main Conjecture with the non-adjacent c₂ term, and Conjecture 1.4. These are the source of the leading terms (3.1) and (3.2).
  - The restricted singular-series sums (LOS Prop. 2.1; Montgomery–Soundararajan).
  - Anderson–Goodman (first-order chains give χ²_(C(C−1))).
- **In the repo.**
  - The class-count rigidity and the column part (round 135, noise/derivation.md §3–5).
  - The class-response law (noise/derivation.md §6, eq. 6.1). The long-range tail (4.1) is its FO_HL part restricted to the within directions, evaluated here with exact finite-H pair sums and S_(>y) in place of a class-count covariance mod Q.
  - The conditioned expansion (rounds 47–55) and the k-th neighbour law (rounds 79–90), which already contain the triple and lag-k laws numerically.
  - Rounds 39–41 compute the long-range rigidity of gap-event counts per residue class mod 210 over whole windows (four-point sums {0, g, d, d + g} with d = 210k). That is the long-range mechanism of §4.4 for a different statistic; the interiors 210k apart are taken as independent there, so there is no short-range part. **Not in rounds 39–41:** the lag-1 triple term, which is the main effect for q = 3, and the within projection.
- **New here (as far as found):**
  - the exact reduction (2.1)–(2.2) of the within part to lag covariances and the second-order Markov deviation;
  - the closed forms (3.1) and (3.2), and the reading of the lag-1 term as the first link of the class-count 1/h tail;
  - the k = 2 cancellation;
  - the long-range tail of second order with slope −V/((C − 1)²λ);
  - the finite-λ evaluation and its component checks.

## 7. Prior art (Task 3)

**Opened in full text:** Lemke Oliver–Soundararajan, arXiv 1603.03720 v4 (PDF; Main Conjecture and Conjectures 1.1–1.6 on p. 2–4; §2 with Prop. 2.1 and (2.10)–(2.23); §4 sketch for r ≥ 3, "it is this last term which creates the additional bias against patterns with a non-immediate repetition").
- LOS state that the primes mod q are not Markovian, and that the k = 2 transition matrix is not the square of the k = 1 matrix.
- They do not consider block counts, variances, χ² statistics or within parts.

**Abstract pages only:**
- Goldston–Ledoan, arXiv 1111.3380: the count of consecutive primes at distance d. No variances or residue classes.
- Kuperberg, arXiv 2109.03767: odd moments of prime counts in short intervals. No gaps or residue classes.
- Holt, arXiv 2405.03540: sieve explanation of the mod-10 biases. No correlations or variances.

**Search snippets only:**
- Montgomery–Soundararajan 2004 (math/0409258), which gives the variance of prime counts.
- Funkhouser–Goldston–Ledoan, arXiv 1802.07609 (large gaps).
- Bloom–Kuperberg, arXiv 2312.09021.
- Ash–Beltis–Gross–Sinnott, Exp. Math. 20 (2011), on successive pairs of residues.
- Goodman, Ann. Math. Statist. 1958 (ψ² for Markov chains).
- No paper by "Leung" on this question was found.

No source found gives the variance of counts of gaps in a residue class or of a given size in short intervals, a χ² or within part of consecutive-prime transition tables, or a long-range tail of such counts.

**Verdict.**
- (3.1) is **KNOWN IN SUBSTANCE**: its content is LOS's non-adjacent c₂ term, read as a lag-1 covariance; the within-part formula itself was not found.
- (3.2), the k = 2 cancellation, the tail (4.1)–(4.2) and the finite-λ law: **NOT FOUND**.

## 8. Doubts

1. **Order-1 truncation.**
   - The conditioned law drops interior pair correlations and the triple-coincidence parts of the 3-point singular series.
   - q = 7's lag-1 term is 8 % short on the validation stretch: −0.00286 against −0.00309 ± 0.0001, which is 0.0005 in W.
2. **Lags ≥ 2 have no S_(>13) memory beyond one step.**
   - For k = 2 the lag-2 term was therefore replaced by LOS's leading term, a choice made after seeing the validation stretch (not the nine windows).
   - Measured against LOS at lag 2: q = 3 is −0.0055 against −0.0065, q = 7 is −0.0016 against −0.0013.
3. **The tail is a linear response.**
   - The weight w₂ is heuristic; w₃ changes it by up to 0.0002.
   - It is truncated at y = 13. y = 17 changes (7, 1) by +0.0008 at 10⁴ and −0.0017 at 10⁶ (λ = 25.7), and the response power keeps growing with y by about 0.02 per added prime, moving into higher modes.
   - On the validation stretch the model's tail is 50–70 % of the measured tail bins for q = 7 (2–2.7σ per bin) and less for q = 5.
   - The law's error carries ±40 % of the tail; the registration's power depends on that choice.
4. **q = 5 declines faster than the law:** −0.012 and −0.024 per decade against −0.006 (1.2σ and 1.4σ). The validation stretch shows large q = 5 lag bins beyond 300, at 2–3σ. This may be a single-stretch fluctuation or a missing piece.
5. **Unmodelled O(1/N) terms.** The ratio estimator (actual N_a in denominators) and block edges (p_(n+k) beyond the block) were not modelled. At H = 10⁴ (N ≈ 370) they could be about 10⁻³. The retrospective fit at 10⁴ to ±0.001 for q = 7 suggests they are small, but no bound was derived.
6. **The data were known.**
   - Every comparison in §5 is retrospective.
   - No parameter was fitted, but y = 13 (the repo's Q), w₂ and the k = 2 lag-2 choice were made by me while the nine-window numbers were in view.
   - The component checks used only the validation stretch.
7. **Error estimates.** The window-scatter errors come from nine values (about ±25 % each) and the registration errors are scaled from them. q = 7 at 10⁴ has window scatter 1.7 times the block error.
8. **The 1/λ dependence is modest.** Over λ = 27.0 → 22.7 it moves the k = 1 cells by 0.005–0.011, so the registration needs 24 windows and a pooled contrast reading.
9. **Gap language for q = 3** is exact only up to O(δ²/N), with δ = N₁ − N/2.
10. **LOS's expansions are asymptotic.** b₁ from LOS is −0.0775 against −0.0874 measured. The finite-λ law, not the expansion, is what the predictions use.
