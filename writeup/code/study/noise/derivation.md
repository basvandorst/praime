# Why the χ² of prime tables are low: class counts are rigid, and every table sees them through its responses

> **In short.** Many χ² values in the repo sit below their degrees of freedom. (1) From the Hardy–Littlewood pair conjecture, the counts N_a of primes in (x, x+H] in the reduced classes a mod q have Var/mean v_q = 1 − [ln(H/q) + Σ_{p|q} ln p/(p−1) − B]/λ, B = 1 − γ − ln 2π, and cross-covariance −Λ(m)/(φ(m)λ) times the mean, m = q/(b−a, q). The ln H rigidity sits on the diagonal only, so the class contrasts, which a χ² sees, are as rigid as the total. (2) A window statistic of any local prime pattern is, to first order, a linear function of the class counts mod Q plus a short-range remainder. So E χ² = E_M χ² − FO_M + FO_HL: a sieve model's value, with its class-count part replaced by the Hardy–Littlewood one. (3) This reproduces round 39's 0.847 (0.856), puts the last-digit tables at 25.9 / 27.7 / 31.1 on 36 (measured 25.4 / 24.5 / 20.0), and leaves round 13's train remainder and round 83's k = 4 value unexplained. Part (1) is classical in substance. Part (2) applied to prime tables was not found in the literature.

Status labels: **proved** (a mathematical step that holds without conjecture), **conjectural (HL)** (follows from the Hardy–Littlewood pair conjecture), **heuristic** (a modelling step), **computed** (a number from a script here), **retrospective** (compared with a value already known to the repo).

Files: [law.py](law.py) → [law_output.txt](law_output.txt); [check_T.py](check_T.py) → [check_T_output.txt](check_T_output.txt); [validate_cov.py](validate_cov.py) → [validate_cov.json](validate_cov.json); [simlaw.py](simlaw.py) → `sl_*.json`; [classresp.py](classresp.py) → `cr_*.json`; [rigid_check.py](rigid_check.py) → [rigid_3e11.json](rigid_3e11.json) (with sl_rigidcmp_3e11.json); [check_TQ.py](check_TQ.py); [granville_dk.py](granville_dk.py) (first look at the overlap inflation); [sieve.py](sieve.py). Registration draft: [round_draft.txt](round_draft.txt).

## 1. Setting

- Window (x, x+H], λ = ln x, modulus q, φ = φ(q). Reduced classes a mod q.
- N_a = #{primes p in the window, p ≡ a mod q}.
- μ₁ = (q/φ)/λ is the chance that an integer of a reduced class is prime (prime number theorem in progressions plus a slowly varying density). The mean is E N_a = (H/q)μ₁ = H/(φλ).

**Pair density (conjectural, HL).** For n ≡ a and n + h ≡ b, both reduced:

  P(n and n+h prime) = μ₁² S_q(h),   S_q(h) = ∏_{p∤q} (1 − ν_p(h)/p)(1 − 1/p)^{−2},  ν_p(h) = 1 if p | h, else 2.   (1.1)

This is the singular series without the primes that divide q (Lemke Oliver–Soundararajan's S_q). Conditioning on the classes mod q fixes the local factors at p | q, which is where the factor (q/φ)² in μ₁² comes from.

## 2. Averages of the singular series over a residue class (proved)

Ramanujan's expansion gives S_q(h) = Σ_{r ≥ 1, (r,q)=1} (μ(r)/φ(r))² c_r(h). For r > 1 coprime to q, the sum c_r(h) averages to 0 over h ≡ c (mod q), because h then runs over every residue mod r. So:

- **S_q has mean 1 on every residue class mod q.** (2.1)
- The full S has mean A_q(c) = ∏_{p|q, p|c} p/(p−1) · ∏_{p|q, p∤c} p(p−2)/(p−1)² on h ≡ c (mod q); it depends on c only through c mod rad(q). (2.2) For example: 3/2 and 3/4 mod 3; 2, 0, 2, 0 mod 4; 5/2 and 15/8 mod 10.

Checked numerically for q = 3, 4, 10, 30 to 10⁻⁵ ([check_T_output.txt](check_T_output.txt)).

## 3. The restricted Montgomery–Soundararajan sums

Define the two-sided sum

  T_q(c; H) = Σ_{0<|h|<H, h ≡ c (q)} (H − |h|)(S_q(h) − 1).

**Result (derived; agrees with Lemke Oliver–Soundararajan Prop. 2.1 and Friedlander–Goldston Prop. 3, which use other weights):**

  T_q(0; H) = −(φ/q) H [ln(H/q) + S_q + γ + ln 2π − 1] + H + O(H^{1/2+ε}),  S_q := Σ_{p|q} ln p/(p−1),   (3.1)

  T_q(c; H) = −H (φ(q)/φ(m)) Λ(m)/q + O(H^{1/2+ε}),  c ≢ 0, m = q/(c, q).   (3.2)

So **the H ln H term appears only on the class c ≡ 0**. For c ≢ 0 the sum is O(H). It is nonzero only when q/(c,q) is a prime power.

**Derivation (q even; the odd q are checked numerically below).**

*Step 1. A divisor sum for S_q (proved).* Since (p−1)/(p−2) = 1 + 1/(p−2),

  S_q(h) = C_q Σ_{d|h, (d,q)=1} μ²(d) g(d),  g(d) = ∏_{p|d} 1/(p−2),  C_q = ∏_{p∤q}(1 − (p−1)^{−2}).

*Step 2. A Dirichlet series (proved).* Write h = dm with m ≡ c d̄ (mod q). Then

  Σ_{h≥1, h≡c} S_q(h) h^{−s} = q^{−s} Σ_{e∈(Z/q)*} G_e(s) ζ(s, ⟨c ē⟩/q),

with G_e(s) = C_q Σ_{d≡e} μ²(d) g(d) d^{−s} and ⟨r⟩ ∈ (0, q]. Characters ψ mod q split G_e as G_e = φ⁻¹ Σ_ψ ψ̄(e) G_ψ, where G_ψ(s) = C_q ∏_{p∤q}(1 + ψ(p)p^{−s}/(p−2)) = C_q L(1+s, ψ) E_ψ(s), and E_ψ is holomorphic for Re s > −1/2.

*Step 3. The principal character near s = 0 (proved).*
- C_q E₀(0) = ∏_{p∤q} (1 − (p−1)^{−2})(1 + 1/(p−2))(1 − 1/p) = 1.
- (ln E₀)′(0) = Σ_{p∤q} ln p [−1/(p−1) + 1/(p−1)] = 0.
- Hence G_{ψ₀}(s) = (φ/q)[1/s + γ + S_q] + O(s).

*Step 4. The mean (proved).* At s = 1, the residue of the S_q series is (1/q)·G_{ψ₀}(1) = 1/q, which is that of Σ_{h≡c} h^{−s}. So the H² terms cancel; this is (2.1) again.

*Step 5. The Perron–Cesàro integral (derived).* Σ_{1≤h≤H, h≡c} (H − h) a_h = (2πi)⁻¹ ∫ D(s) H^{s+1}/(s(s+1)) ds. Shift the contour to Re s = −½ + ε; this needs polynomial growth of ζ(s, a) L(1+s, ψ) on vertical lines, which is standard but not written out here.

*Step 6. The residue at s = 0, for c ≡ 0 (derived).* Every ζ(s, ⟨cē⟩/q) is ζ(s). Only ψ₀ survives, and the subtracted "1" contributes q^{−s}ζ(s). Using ζ(s) = −½ − (s/2) ln 2π + …, q^{−s} = 1 − s ln q and H^{s+1}/(s+1) = H(1 + s(ln H − 1) + …), the residue is

  −(φ/(2q))H[ln H − 1 + ln 2π − ln q + γ + S_q] + H/2.

Doubling gives (3.1).

*Step 7. The residue at s = 0, for c ≢ 0 (derived).* Write g = (c, q) and m = q/g. Then cē runs over g·(Z/m)*, each value φ(q)/φ(m) times.
- Σ_e ζ(0, ⟨cē⟩/q) = Σ_e (½ − ⟨cē⟩/q) = 0, so there is no double pole.
- The principal character contributes (H/q) Σ_e ζ′(0, ⟨cē⟩/q), with ζ′(0, a) = ln Γ(a) − ½ ln 2π. By Gauss's multiplication formula and Möbius inversion, Σ_{(r,m)=1} ln Γ(r/m) = (φ(m)/2) ln 2π − Λ(m)/2. So this term is −(φ(q)/φ(m)) Λ(m) H/(2q).
- The non-principal ψ contribute φ⁻¹ G_ψ(0) Σ_e ψ̄(e)(½ − ⟨cē⟩/q) H. The subtracted part contributes −(½ − ⟨c⟩/q) H.

Both of the last two change sign under c → −c, because ⟨−r⟩ = q − ⟨r⟩. So in the two-sided sum they cancel. This is why the odd-character L(0,χ)L(1,χ) term and the B_q(v) term of LO–S Prop. 2.1 drop out here. Result (3.2).

*Step 8. Numerical check (computed, [check_T.py](check_T.py)).* S(h) was evaluated exactly for h ≤ 2·10⁶, and T_q(c; H) summed directly for q = 2, 3, 4, 5, 6, 7, 9, 10, 12, 30, 210, every class c, and H = 10⁵ and 10⁶.
- T/H agrees with (3.1)–(3.2) to within 0.001 at H = 10⁵ and within 3·10⁻⁴ at H = 10⁶. Examples: q = 210, c = 0: −1.1834 against −1.1829 (H = 10⁵), −1.7092 against −1.7092 (H = 10⁶); c = 105: −0.1584 against −0.1584 (H = 10⁶); q = 9, c = 3: −0.3662 against −0.3662 (H = 10⁶).
- Odd q follow the same formulas.
- An earlier run with C₂ truncated at 2·10⁶ (3·10⁻⁸ too large) showed a uniform shift of 0.03/q at H = 10⁶, the expected size 3·10⁻⁸·H/q. The exact C₂ removes it.

## 4. The covariance of class counts (conjectural, HL)

*Step 1.* The pairs (n, m) in the window with n ≡ a and m − n = h number (H − |h|)/q + O(1) for h ≡ b − a, and none otherwise.

*Step 2.* With (1.1),

  Cov(N_a, N_b) = δ_ab (H/q) μ₁(1 − μ₁) + (μ₁²/q) T_q(b − a; H).   (4.1)

*Step 3.* Divide by E N_a = Hμ₁/q and insert (3.1)–(3.2). The −μ₁ of the Bernoulli term cancels the +H of (3.1):

  Var(N_a)/E N_a = v_q(H) := 1 − [ln(H/q) + S_q − B]/λ,  B = 1 − γ − ln 2π = −1.4151,   (4.3)

  Cov(N_a, N_b)/E N_a = −Λ(m)/(φ(m) λ),  m = q/(b − a, q).   (4.4)

*Step 4. Consistency with the total (proved identity).* Σ_{b≠a} Λ(m_b)/φ(m_b) = ln q − S_q. Proof: count b by m. For p^k ∥ q, each j < k contributes ln p and j = k contributes (p−2) ln p/(p−1). The identity is checked for 11 moduli in law_output §1. Hence

  Var(N)/E N = v_q − (ln q − S_q)/λ = 1 − (ln H − B)/λ,

which is Montgomery–Soundararajan's variance for counts at constant density (Goldston 1990 / MS (47)).

*Step 5. The equidistribution χ² (fitted total) (conjectural, HL).* X = Σ_a (N_a − N/φ)²/(N/φ), and E X = [Σ_a Var N_a − Var(N)/φ]/E N_a. So

  E X/(φ − 1) = v_q + (ln q − S_q)/((φ − 1)λ).   (4.6)

*Step 6. The form seen from the zeros.* Since λ = ln x,

  v_q = [ln(xq/H) + B − S_q]/λ.

The ln q is the conductor's contribution: the "+ ln q" that the race variance gains in round 102. Summed over the φ classes and weighted by Λ, this is Goldston–Yildirim 1998 (1.13) in local form, h(ln(y/h) + ln q − S_q + B). Their constant integrated over [x, 2x] gives −(γ + ln(π/2)); this was checked by the literature agent, not by me.

**Validation on primes (computed; [validate_cov.py](validate_cov.py)).** Window [5·10¹¹, 5·10¹¹ + 10⁸): not a registered window, and not proposed below.
- H = 10⁵ (1000 windows): the contrast statistic (4.6) is measured / predicted 0.584/0.561 (q = 3), 0.580/0.572 (4), 0.574/0.580 (5, 10), 0.591/0.592 (7), 0.590/0.582 (12), 0.592/0.595 (30), 0.648/0.648 (210). Errors are 0.004–0.027; all within 1σ.
- The term S_q is needed: without it, q = 210 would be predicted at 0.721, against 0.648 ± 0.004 measured.
- The per-class Var/mean also fits, for example 0.645 against 0.645 for q = 210.
- H = 10⁶ (100 windows): the contrasts fit (0.488–0.582 against 0.475–0.563, within 1.8σ; the 1.8σ is q = 210). The total count is low, 0.286 against MS 0.435 ± 0.041. This is a single stretch of 10⁸; round 119 showed such window deviations are the zeros' terms for that stretch.
- The between-class part, which the χ² sees, is the part that fits.

**Regime.** For fixed q these formulas are an ensemble mean over x. For the full interval [1, x] the variance over residue classes has a limiting distribution governed by the zeros (Fiorilli 2015, Prop. 1.8), not a limit. For many short blocks (λq ≪ H ≪ x^{1−ε}) the χ² averages the zeros' terms over blocks, and the HL pair calculation gives the mean (MS / Goldston–Yildirim regime). The anomalous statistics of rounds 13, 83 and 85 use single long windows (H/x from 10⁻⁶ at 10¹⁵ to about 1 at train, where the constant-density form is outside the MS range; round 39's varying-density sums agree with it there to 0.01). There the realized class contrasts are one draw. For example, the four last-digit row totals at 10⁹ give a contrast χ²/dof of 0.04 against 0.12 predicted.

## 5. What a χ² sees: the between-class part is not near 1

A χ² with fitted or fixed margins removes the total count; only the contrasts remain.
- By (4.3)–(4.4), the H ln H rigidity lies on the diagonal of the class covariance. The cross-class terms are O(1/λ) and vanish unless q/(b−a, q) is a prime power.
- So **each class is as rigid as the total.** The contrast χ²/dof (4.6) differs from the total's 1 − (ln H − B)/λ only by (ln q − S_q)(1 + 1/(φ−1))/λ.
- For q = 210 in round 13's windows this gives 0.087, 0.319, 0.432, 0.572 (train, A, B, C), against 0.080, 0.386, 0.360, 0.454 measured (r13f).
- So the class-count spreads are explained, as round 36 found numerically. The closed form reproduces round 36's sums to 3 decimals: v = 0.0824, 0.316, 0.430, 0.570 against 0.0825, 0.3154, 0.4295, 0.5701.

This part of the "low-χ² family" is classical in substance. It is not an anomaly.

The rest of the family are **tables of events**: gap classes per residue, transitions between last digits, bands of D_k. These need the mechanism of section 6.

## 6. Tables of events: the class-response law

**6.1. Long-range covariance of local patterns (heuristic).**
- Let F_n be any local pattern indicator at n: "n prime, the next prime ≡ b", or "D_k(n) in band j".
- To first order in the long-range pair correlations, Cov(F_n, F_{n+d}) is a sum over pairs (point of F_n, point of F_{n+d}) of response × response × (S − 1). The response is +P(F)/P(prime) for points that must be prime and −P(F)/P(composite) for points that must be empty.
- Summed over n and d in a window, §3 applied at a modulus Q that resolves the patterns (Q = 30030 here) says:
  - the H ln H part comes only from point pairs in the same class mod Q;
  - the O(H) part comes from pairs whose difference c has Q/(c, Q) prime.
- These are exactly the class-count covariances of §4 at modulus Q. So the long-range part of every window statistic acts through the class counts N_c (c mod Q), with responses ψ_c = the change of the statistic per additional prime at a site of class c.

**6.2. Short-range model (proved within the model).**
- M_Q: every residue coprime to Q = 30030 is "prime" independently, with probability p = (Q/φ(Q))/λ. (It is a first-order chain on the classes mod Q with the sieve's own gap law.)
- In M_Q, the class counts are independent with Var = (1 − p)E N_c.
- By periodicity, the first-order Hoeffding projection of any window statistic is exactly its regression on the class counts N_c, with coefficients ψ_c. ψ_c is measured by insertion: the change of the statistic when a non-prime site of class c is made prime. That change does not depend on the site's own state, so averaging over non-prime sites is unbiased.

**6.3. The law (heuristic).** For a Pearson χ² whose expected counts are a calibrated law (times row totals where the table is row-conditioned):

  **E χ² = E_M χ² − FO_M + FO_HL,  FO_X = Σ_cells (1/E_cell) Σ_{c,c'} r_{cell,c} Cov_X(N_c, N_c') r_{cell,c'}**   (6.1)

- r is the response of the cell residual: count minus expected, the expected count responding through the row total or total.
- Cov_M = diag (1 − p)E N_c.
- Cov_HL = E N_c [v_Q(H) δ_cc' − Λ(m)/(φ(m)λ)(1 − δ_cc')], from (4.3)–(4.4) at modulus Q. For squarefree Q only c′ that differ from c in exactly one prime's residue enter, with −ln p/((p−1)λ).

In words: the primes are the sieve model with its class-count covariance replaced by the Hardy–Littlewood one.

Assumptions:
- (a) The primes' short-range structure is that of M_Q. The Hardy–Littlewood factors S_{>13} change transition probabilities by a few per cent; their effect on E χ² was not computed.
- (b) The first-order pair expansion holds at long range.
- (c) The part of the statistic orthogonal to the class counts is unchanged when the class counts are made rigid.

**Internal check of (c) (computed; [rigid_check.py](rigid_check.py), no primes).**
- Model R_Q is M_Q with every class count mod Q fixed in every block of 33 periods, so Cov = 0. The law predicts E χ² = E_M χ² − FO_M.
- Block-homogeneity statistics at λ = ln(3·10¹¹), 600 blocks:

| statistic | measured in R_Q | law | difference |
|---|---|---|---|
| mod-3 transitions, k = 1 | 0.973 ± 0.054 | 0.914 ± 0.074 | +0.6σ |
| mod-3 transitions, k = 2 | 0.986 ± 0.057 | 0.984 ± 0.071 | 0.0σ |
| mod-10 transitions, k = 1 | 8.40 ± 0.16 | 8.42 ± 0.20 | −0.1σ |
| mod-10 transitions, k = 2 | 9.02 ± 0.17 | 8.58 ± 0.18 | +1.7σ |
| D₄ bands | 23.12 ± 0.33 | 23.10 ± 0.36 | 0.0σ |

The decomposition holds within the models to about 0.2–0.4 in χ².

**6.4. Special cases.**

(i) *Equidistribution of classes.* The response is the identity, so (6.1) reduces to (4.6) (proved, by aggregation from Q to q; the identity of §4 Step 4 makes the levels consistent).

(ii) *Transition tables; the rigid-margin part (proved algebra + HL).*
- Rows are the class of p_n and columns the class of p_{n+k}, with expected N_a P(b|a).
- Up to O(k) boundary terms, the column sums of the residual table are c = (I − Pᵀ)δ, where δ is the class-count deviation. So C − 1 directions of the table are class contrasts.
- In z-space the row-multinomial null covariance is an orthogonal projection. Hence χ² = col + within exactly, with col = cᵀ M⁺ c and M = Σ_a N_a(diag P_a − P_a P_aᵀ).
- E[col] = tr(M⁺ (I−Pᵀ) Cov(δ) (I−P)). This is C − 1 for a chain on the same classes, (1 − p)(C − 1) for M_Q and v_q(C − 1) for the primes, each times a factor set by P (about 1.2 for the anti-persistent mod-3 table at k = 1).
- With rigid columns and nothing else, E χ² = R(C−1) − (C−1)(1 − ρ), i.e. **E χ²/dof ≈ 1 − (1 − ρ)/R** for a near-uniform R × R table.
- This is the classical (R−1)(C−1) of a contingency table with both margins fixed, reached because the primes fix the column margins. The rest of FO is the within-row composition of the fine classes, which (6.1) also makes rigid.

(iii) *Bands of D_k.*
- Σ_n D_k(n) telescopes to k × (span), so the mean direction is rigid.
- The first-order part of the band χ² in M_Q is FO_M = 1.28, 2.51, 3.60, 4.57 for k = 1–4 ([sl_band_1e9.json](sl_band_1e9.json)), essentially that one scale direction (about k, because overlapping windows share gaps; the class contrasts mod 30 explain nothing, [granville_dk.py](granville_dk.py)).
- The law removes it (FO_HL = 0.13–0.38), but M_Q's χ² is inflated by the overlap of the windows (22.4, 25.2, 27.7 for k = 2, 3, 4 on 19). So the law predicts **above** 19 for k ≥ 2.

(iv) *Residue × gap tables (round 13).* (6.1) evaluated with rows p mod 210 and columns the gap classes gives 0.856, 0.904, 0.914, 0.956, against round 39's four-point HL 0.847, 0.902, 0.925, 0.954. Two different computations agree.

## 7. Retrospective table (reported values were known; nothing here is a test)

All values come from [law_output.txt](law_output.txt), sections 3–7.

**Columns.**
- "nominal" is the null the repo used.
- M_Q is the sieve model, a first-order chain on the classes mod 30030 with no long-range rigidity.
- "law" is (4.6) or (6.1).
- The ± is the model's own scatter of a single window's value, scaled to the law. For round 13 it is round 13's null
  noise.

| statistic (round) | window | measured | nominal | M_Q | law | z (law) |
|---|---|---|---|---|---|---|
| prime counts over 48 classes mod 210, Var/mean (r13f) | train | 0.080 | 1 | 0.707 | 0.087 ± 0.018 | −0.4 |
| | A | 0.386 | 1 | 0.774 | 0.319 ± 0.066 | +1.0 |
| | B | 0.360 | 1 | 0.811 | 0.432 ± 0.089 | −0.8 |
| | C | 0.454 | 1 | 0.849 | 0.572 ± 0.118 | −1.0 |
| residue × gap table mod 210, χ²/dof, dof 401 (13) | train | 0.691 | 1 ± 0.058 | 0.992 | 0.856 (round 39: 0.847) | −2.8 |
| | A | 0.857 | 1 | 0.987 | 0.904 (0.902) | −0.8 |
| | B | 0.995 | 1 | 0.975 | 0.914 (0.925) | +1.4 |
| | C | 0.894 | 1 | 0.994 | 0.956 (0.954) | −1.1 |
| last-digit tables, Pearson χ² on 3 × 12 (85, 86) | 10⁹ | 25.4 | 36 | 33.6 ± 9.3 | 25.9 ± 7.1 | −0.1 |
| | 10¹² | 24.5 | 36 | 34.4 ± 9.5 | 27.7 ± 7.6 | −0.4 |
| | 10¹⁵ | 20.0 | 36 | 36.0 ± 10.9 | 31.1 ± 9.4 | −1.2 |
| same, Σz² of the 48 binomial z | 10⁹ / 10¹² / 10¹⁵ | 34.1 / 32.6 / 26.7 | 48 | – | ≈ 34.6 / 37.0 / 41.4 | (≈ 4/3 × law χ²; approximation) |
| block variance of the digit proportions (review 24) | all | 0.75–0.87 | 1 | – | 0.68–0.77 (10⁹), 0.77–0.80 (10¹²), 0.86–0.92 (10¹⁵) | block count unknown |
| bands of D_k, χ² on 19 (83) | k = 2 | 22.2 | 19 | 22.4 ± 7.7 | 20.2 ± 6.9 | +0.3 |
| | k = 3 | 16.7 | 19 | 25.2 ± 9.2 | 21.9 ± 8.0 | −0.6 |
| | k = 4 | 8.2 | 19 | 27.7 ± 9.9 | 23.5 ± 8.4 | −1.8 |

**The last-digit tables split exactly into the column part and the within part (§6.4 ii).**

| scale | measured column parts (k = 1, 2, 3) | law's column parts | measured within parts | law's within parts |
|---|---|---|---|---|
| 10⁹ | 0.16, 0.13, 0.14 | 0.44, 0.38, 0.37 | 8.62, 9.54, 6.79 | 7.51, 8.63, 8.58 |
| 10¹² | 0.34, 0.31, 0.28 | 0.81, 0.71, 0.71 | 7.53, 12.70, 3.31 | 8.24, 8.62, 8.64 |
| 10¹⁵ | 0.39, 0.34, 0.37 | 1.37, 1.24, 1.23 | 3.81, 10.26, 4.89 | 8.76, 8.98, 9.49 |

- The column part is 3 per table under the multinomial null; it is 0.1–0.4 measured. That realized class contrasts
  are small is classical (§4–5): one draw per window, about 0.3–0.4 of the ensemble value.
- The within parts total 24.95, 23.54 and 18.96, against 24.7, 25.5 and 27.2 by the law and 27 for a chain. Only
  10¹⁵ is low, by about 1σ.

**Reading of the table.**
- *Explained to the noise:*
  - the class-count spreads (classical in substance);
  - the last-digit tables at 10⁹ and 10¹², and their block errors (the law; 10¹⁵ is 1.2σ low);
  - round 13 at A, B, C.
- *Not explained:*
  - Round 13 at train: 0.691 against 0.856, −2.8σ. This is round 39's remainder; (6.1) reproduces round 39, so it
    adds nothing there.
  - Round 83, k = 4: 8.2 against 23.5 ± 8.4, −1.8σ. Here the nominal χ²₁₉ null was wrong in the other direction. The
    overlap of the D_k windows inflates the mean and the spread, and rigidity removes only the scale direction (about
    k(1 − v)). The law predicts more than 19 for k = 3, 4, so the low k = 4 value stays a lower-tail event: 0.6–2.4 %
    of the simulated block values (scaled or shifted by the law), against 1.5 % for χ²₁₉. k = 2 and k = 3 sit at the
    64 % and 29–32 % points. (k = 1, 26.2 on 16, is high, not low: it tests the calibration of the gap law, z = +2.0 against
    the law's 15.0 ± 5.6.)
- *Not reconstructed:*
  - MERSENNE §10, §45, §53 (prime bytes, lighthouse, pair clocks) and §28 (Goldbach comet). Whether they are
    class-count tables was not checked.
- The family also contains selection: it was assembled from the low values. Round 15 (0.99–1.09), round 82 (48.8/48),
  round 84 (50.6/48) and round 71 (25.2/30) are the same kind of statistic and are not low.

## 8. What is classical and what is ours

**Classical / known in substance.**
- Montgomery–Soundararajan's sum (MS 2004 (47), after Goldston 1990).
- The restricted sums with the log only at v ≡ 0 and the Λ(q/d)/φ(q/d) constant: Lemke Oliver–Soundararajan, arXiv 1603.03720 Prop. 2.1, read in full text by the literature agent. With their exponential weight their constants carry also B_q(v) and an odd-character L(0,χ)L(1,χ) term, which cancel in our two-sided Cesàro sum.
- Friedlander–Goldston Prop. 3, as quoted by Fiorilli; FG itself was not accessible.
- The total variance over all classes in short intervals: Goldston–Yildirim, Can. J. Math. 50 (1998) (1.13), read from page images by the literature agent; our (4.3) summed over a and weighted by Λ is its local form.
- The zero-side constant log q − γ − log 2π − S_q: Fiorilli 2015, Remark 1.9.
- The contingency-table fact that fixed margins give (R−1)(C−1), and Rao–Scott design effects for χ² under non-multinomial sampling (statistics, classical; recollection of the agent, not checked).

**Not found (ours as far as we searched).**
- The per-class and cross-class covariance (4.3)–(4.4) written as such: the closed form with −Λ(m)/(φ(m)λ). It follows from LO–S's constants with Cesàro weights, so it is "known in substance".
- The class-response law (6.1) for χ² statistics of prime tables, the M_Q decomposition with insertion responses, its rigid-model check, and its evaluation for the repo's tables. The agent found no source linking class-count rigidity to χ² tests of consecutive-prime tables; the Lemke Oliver–Soundararajan paper has no χ².

**Verdict on the final formula.**
- (4.3)/(4.6): KNOWN IN SUBSTANCE (GY 1998 (1.13) summed over classes; LO–S Prop. 2.1 for the restricted sums; FG Prop. 3).
- (6.1): NOT FOUND.
- In reading terms: LO–S Prop. 2.1 and MS (47) were read in full text (arXiv) by the agent. GY 1998 was read from page images (scan). FG and Anderson–Goodman were not accessible (FG quoted via Fiorilli; Anderson–Goodman abstract only). I did not read the papers myself; I checked the formulas numerically (§3).

## 9. Doubts

1. **The short-range model is not the primes.**
   - M_Q has the sieve to 13 but not the Hardy–Littlewood factors S_{>13}. Its gap law and transition probabilities
     differ from the primes' by a few per cent.
   - E_M χ² − FO_M, the part of every prediction that is not class counts, comes from M_Q. Its error for the primes
     was not estimated.
   - The exact HL chain (the repo's class_chain) would be the better short-range model; it was not run here.
2. **First-order pair expansion.**
   - Four-point terms of second order are dropped.
   - Round 39's exact four-point sums agree with (6.1) to 0.01 in χ²/dof (§7), which bounds the effect for that
     table only.
   - In the last-digit tables interior emptiness matters, and nothing checks it except the rigid-model test, which
     checks the decomposition, not the pair expansion.
3. **Assumption (c), checked only in models.** The rigid-model check (§6.3) agrees to 0.2–0.4 in χ², with one 1.7σ
   deviation (mod-10 transitions, k = 2: 9.02 against 8.58).
4. **One window, one draw.** The retrospective statistics come from single long windows (H/x up to about 1 at train).
   - The realized class contrasts there are one draw of a quantity that, for fixed q, has a distribution rather than a
     limit (Fiorilli).
   - The z-scores use the model's ensemble spread. The registration avoids this with 300 blocks.
5. **Finite H/Q.** The closed forms are asymptotic.
   - At Q = 30030 and H = 10⁶ (33 members per class) v_Q is within 0.001.
   - The cross-class T values are off by up to 0.02·H (Q/7: −0.069 against −0.062) ([check_TQ.py](check_TQ.py)),
     changing a covariance ratio by about 0.004.
   - Not propagated into FO_HL.
6. **Contour shift.** The O(H^{1/2+ε}) error in §3 rests on a standard contour shift that I did not write out. The
   numerical check (to 3·10⁻⁴ H at H = 10⁶) is the evidence.
7. **Validation stretch.** In the validation stretch at H = 10⁶ the total count was low (0.286 against 0.435, 3.6σ by
   the 100-window error). That error understates the slow zeros' term over a single 10⁸ stretch; the contrasts fitted.
   The registration's 1e6 blocks are in this regime, so S1's errors should be taken from the blocks, as registered.
8. **Block errors (review 24).** The number of blocks is not recorded, so that row is only a range. The per-cell
   variance ratio was approximated by the table's χ² ratio, scaling FO_HL with v_Q.
9. **Σz².** The law value is approximated as 4/3 of the law's Pearson χ² (cells with P ≈ 1/4); not computed cell by
   cell.
10. **Literature.**
    - The literature was read by a helper agent. LO–S and MS were read in full arXiv text, Goldston–Yildirim 1998 from
      decoded page images (a scan), Friedlander–Goldston only as quoted by Fiorilli and Keating–Rudnick, and
      Anderson–Goodman as an abstract only.
    - The agent flagged a possible misprint (−1 for +1) in GY 1998 Theorem 3's constant. This was not used.
    - I verified the formulas numerically, not against the papers' text.
11. **Prior art for (6.1).** "Not found" covers what we searched. Rao–Scott design effects (statistics) are the
    general framework that (6.1) instantiates, and hyperuniformity of the primes (Torquato et al.) is related in
    spirit. Neither was read.
12. **The registration's prediction error.**
    - It includes the simulation error of E_M − FO_M (0.004–0.010 in the ratio) but not the M_Q-versus-primes
      difference of doubt 1. That could move the q = 10 predictions by more than their stated 0.004.
    - The within part (≈ 0.96–1.0) is the reading most exposed to it. A failure of reading (3) alone with (2) passing
      should be read in that light.
