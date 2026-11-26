# Martin's Conjecture 1.10 for one free prime (levels P·q, P fixed): a proof without hypothesis

Check script: `martin_check.py` (own code; output `martin_check_out.txt`; **15 checks, 0 failed**, 25 s).
Cited notes of the project: **[ML]** = `writeup/code/study/lowzero/universal/murmur/theorem/murmurations_one_level.md`
(refereed three times), **[AL]** = `.../theorem/al/al_theorem.md` (refereed once; its Theorem A and Steps 1–4, 8 are used).

**Labels.** *proved*: argued here (or as a listed change of a step of [ML]/[AL], each change written out); *quoted*:
from the literature, with *read* (source text opened for this note), *read (secondary)* (statement read as restated
in a later source) or *not opened*; *checked*: an assertion of `martin_check.py`; **NOT closed**: a gap. Nothing in
this note is refereed.

## 0. The conjecture as printed (quoted; read in the LaTeX source of arXiv 2409.02338, file `Wqbias-corr.tex`, the corrected version after *Mathematika* 2025; numbering from its counters: thm, prop, cor, rem, conj share one counter within sections, so `conj:eps-murm` is 1.10)

> **Conjecture 1.10** (Murmurations on Atkin–Lehner eigenspaces). Fix k, r and 0 ≤ m < r. Fix primes p₁ < ⋯ < p_m.
> Let 𝒩 be the set of levels N = p₁⋯p_r ∈ ℕ^sqf_r such that p₁ < p₂ < ⋯ < p_r. Let ℱ be the family of weight k
> newforms with level in 𝒩. For ε = (ε₁, …, ε_r) ∈ {±1}^r, let ℱ^ε be the subset of f ∈ ℱ with Atkin–Lehner sign ε_i
> at p_i, and consider the averages
> A^ε_ℱ(ℓ, X; β) = (1/#ℱ^ε(X, βX)) Σ′_{f∈ℱ^ε(X,βX)} ℓ^(1−k/2) a_ℓ(f).
> These averages have murmurations which are scale invariant in ℓ/N. That is, as ℓ, X → ∞ such that ℓ/X → x,
> A^ε_ℱ(ℓ, X; β) → M^ε_ℱ(x; β) for a murmuration function M^ε_ℱ which is continuous on [0, ∞) × (1, ∞).

ℱ^ε(X, Y): forms of the family with level X ≤ N ≤ Y; Σ′: levels prime to ℓ (the denominator counts all levels);
ℓ prime. **Remark 1.11(3)** (quoted): for r = m + 1 "we expect the murmuration functions to be different on each
Atkin–Lehner eigenspace". **The sign at the free prime is fixed**: ε ∈ {±1}^r prescribes the sign at every p_i,
i ≤ r, and the figures for S_2(2q), S_4(2q) plot four colours (++, +−, −+, −−). So the ε_q-odd part does not cancel.
**Normalisation**: ℓ^(1−k/2)a_ℓ(f) = √ℓ·λ_f(ℓ), λ_f(ℓ) = a_ℓ(f)ℓ^(−(k−1)/2) (the normalisation of [ML], [AL]).

## 1. Result

**Setting.** k ≥ 2 even; P = p₁⋯p_m squarefree, m ≥ 0 (P = 1 allowed); q runs over primes q > p_m; N = Pq (this
is 𝒩 for r = m + 1). For Q | N let t_Q^(N)(ℓ) := ℓ^(1−k/2) tr(T_ℓW_Q | S_k^new(N)) = Σ_{f∈H_k(N)} w_Q(f)√ℓλ_f(ℓ).
Local data ([ML §2], [AL §1]): e_v(s) = v²/(v²−1) (v | s, s = 0 included), 1 − 1/((v−1)²(v+1)) (v ∤ s);
e_v^new(s) = −v²/(v²−1) (v ∤ s), −v/(v+1) (v | s); ν′(s) = Π_v e_v(s); κ(P) = P/φ(P). For Q₀ | P, M := P/Q₀:

    E(0) := ζ(2) Π_{v|Q₀}(1 − v⁻²) Π_{v|M}(−(v−1)/v),   E(s) := ν′(s) Π_{v|P} e_v(s)⁻¹ Π_{v|M} e_v^new(s)  (s ≥ 1),
    D^[Q₀,M](y) := (12κ(P)/((k−1)π)) [ (−1)^(k/2) E(0)√(y/M) − (1/M) Σ_{1≤s<2√(yM)} U_{k−2}(s/(2√(yM))) E(s)√(4yM − s²) ]
                   + 12κ(P)μ(M) y·[k = 2].                                                                    (1.1)

This is the q → ∞ limit of the density D_k^(Q₀q, M) of [AL (1.1)] (W_Q-density at one level N = Pq, y = ℓ/N).

**Theorem.** Fix k, P, β > 1, Y > 0 and ε = (ε₁, …, ε_m, ε_q) ∈ {±1}^(m+1). Put

    M^ε(x; β) := ε_q Σ_{Q₀|P} ε(Q₀) ∫₁^β D^[Q₀,P/Q₀](x/t) t dt / ∫₁^β t dt      (ε(Q₀) = Π_{p_i|Q₀} ε_i).

Then, **with no hypothesis**, uniformly for primes ℓ ≤ YX with ℓ ∤ P, as X → ∞,

    A^ε(ℓ; X, βX) = M^ε(ℓ/X; β) + O_{k,P,β,Y}(1/log X)    (ineffective constant, through Siegel–Walfisz).

With the prime-density weight w_X(t) := log X / log(tX/P) inserted in both integrals the error is O((log X)^(−A))
for every A. M^ε is continuous on [0, ∞) × (1, ∞) and M^ε(0; β) = 0.

**Corollary 1 (Conjecture 1.10 for r = m + 1; proved).** For every m ≥ 0, every fixed p₁ < ⋯ < p_m, every even k and
every ε ∈ {±1}^(m+1), A^ε(ℓ, X; β) → M^ε(x; β) as X → ∞ and ℓ/X → x ∈ [0, ∞) (ℓ prime; ℓ → ∞ not even needed).
The functions: M^(ε′, −ε_q) = −M^(ε′, ε_q) exactly; for m = 1, P = p₁: **M^(ε₁,ε_q) = ε_q(M_q + ε₁M_N)**, four
functions in two opposite pairs (M_q from Q = q, M_N from Q = N), as Remark 1.11(3) expects. For P = 2, k = 2:

    D_q(y) = π√(2y) + (12/π) Σ_{s<2√(2y)} c(s)ν′(s)√(8y − s²) − 24y,   c(s) = 2 (s odd), ½ (s even);
    D_N(y) = −3π√y − (24/π) Σ_{r<2√y} c₂(r)ν′(r)√(4y − r²) + 24y,      c₂(r) = 3/2 (r odd), 3/4 (r even)

(the q → ∞ densities of [AL §1], recomputed from (1.1): checked F1, max deviation 4.6·10⁻⁴ from [AL]'s table).

**Corollary 2 (proved by Steps 4–7 alone).** Martin's Conjecture 1.8 for the pairs (N, Q) = (Pq, Q₀q), Q₀ | P, in his
normalisation A^𝒬 = (1/#ℱ(X,βX)^(ℓ)) Σ′_N Σ_f w_Q(f)√(N/Q)ℓ^(1−k/2)a_ℓ(f): A^(Q₀q) → √(P/Q₀)·∫₁^β D^[Q₀,P/Q₀](x/t)t dt/∫₁^β t dt,
**without smoothing**. This covers Type III pairs with the free prime in Q (Q₀ ≠ P), and Type I with Q prime and M = P
fixed (Q₀ = 1; Martin's "variant where Q is restricted to ℕ^sqf_r", r = 1). The pairs with q ∤ Q (weight √(N/Q) ≍ √X)
are **NOT closed** (Section 5).

## 2. Proof

Throughout I := {q prime : X ≤ Pq ≤ βX, q ∤ P, q ≠ ℓ}; X ≥ X₀(k, P, β, Y), so q > p_m, q² > 4YXP ≥ 4ℓP and
#I ≍ X/log X. Constants depend on k, P, β, Y, A.

**Step 1. Eigenspaces and the trace formula (quoted; checked A1, A2).** On S_k^new(N), N squarefree, the W_v (v | N)
are commuting involutions commuting with T_ℓ (ℓ ∤ N), newforms are common eigenvectors (Atkin–Lehner 1970, not
opened), so Σ_{f∈H^ε(N)}√ℓλ_f(ℓ) = 2^(−m−1) Σ_{Q|N} ε(Q) t_Q^(N)(ℓ). *Martin, Proposition 5.5* (read): for N = QM
squarefree, ℓ ∤ N prime,

    t_Q^(N)(ℓ) = −½ Σ_{s∈ℤ, s²Q<4ℓ} U_{k−2}(s√Q/(2√ℓ)) ξ_{Δ_s}(M) H(4ℓQ − s²Q²) − δ_{N=1}ℓ^(1−k/2) + μ(M)(ℓ + 1)[k = 2],

Δ_s = s²Q² − 4ℓQ, ξ_Δ multiplicative on squarefree M, ξ_Δ(v) = (Δ/v) − 1 (v² ∤ Δ), vH(|Δ|/v²)/H(|Δ|) − 1 (v² | Δ);
|ξ_Δ(v)| ≤ 2 [AL Step 8]. *Checked:* A1 exact equality with PARI trace forms, tr(T_ℓW_Q) = (−1)^(ω(Q))Q^(1−k/2)a_{Qℓ}(trace
form of the new space), 2,120 triples (18 levels 2q, 3q, 6q; k = 2, 4; ℓ ≤ 47; every Q | N); A2 eigenspace dimensions
and eigenspace traces 2^(−r)Σ_Q ε(Q)t_Q against PARI newforms with their Atkin–Lehner eigenvalues, 1,520 cases. So

    A^ε = Σ_{q∈I} Σ_{Q|Pq} ε(Q) t_Q^(Pq)(ℓ) / Σ_{X≤Pq≤βX} 2^(m+1)|H^ε(Pq)|.                                     (2.1)

**Step 2. Split (proved).** Q = Q₀q ("q-terms", Step 4–7) or Q = Q₀ ("fixed terms", Step 3), Q₀ | P; ε(Q₀q) = ε(Q₀)ε_q.

**Step 3. The fixed terms vanish in the limit: Lemma T′ (proved, no hypothesis).**
(a) *Exact splitting (proved; checked A3, 4,126 cases).* Since |Δ_s| ≤ 4ℓQ₀ < q², ξ_{Δ_s}(Mq) = ξ_{Δ_s}(M)((Δ_s/q) − 1),
and μ(Mq) = −μ(M). Hence t_{Q₀}^(Pq)(ℓ) = Leg_{Q₀}(q) − t_{Q₀}^(P)(ℓ) − δ_{P=1}ℓ^(1−k/2) with

    Leg_{Q₀}(q) := −½ Σ_{s²Q₀<4ℓ} U_{k−2}(s√Q₀/(2√ℓ)) ξ_{Δ_s}(M) (Δ_s/q) H(|Δ_s|),    Δ_s = s²Q₀² − 4ℓQ₀ (independent of q).

(b) *Fixed level (quoted: Deligne, not opened).* |t_{Q₀}^(P)(ℓ)| ≤ 2√ℓ·dim S_k^new(P); summed over q ∈ I: ≪ X^(3/2).
(For P = 2, 3 and k = 2 it is 0: checked A4, 1,716 cases.)
(c) **Lemma T′.** *Σ_{q∈I} Leg_{Q₀}(q) ≪_δ X^(7/4+δ) for every δ > 0, uniformly for ℓ ≤ YX.* Proof:
(i) Write Δ_s = d_s f_s² with d_s a negative fundamental discriminant. q ∤ f_s (q² > |Δ_s|), so (Δ_s/q) = χ_{d_s}(q),
χ_d := (d/·) the primitive real character of conductor |d| ≤ Z := 4ℓQ₀ ≤ 4YPX. Put S(d) := Σ_{q∈I} χ_d(q).
(ii) *Multiplicity (proved; checked D1).* If d_s = d for s ≥ 0, then (Q₀s)² + |d|f_s² = 4ℓQ₀, so α = Q₀s + f_s√d ∈ O_K,
K = ℚ(√d), has norm 4ℓQ₀, and distinct s give distinct α. The number of elements of norm n in O_K is ≤ 6·#{ideals of
norm n} = 6Σ_{e|n}χ_d(e) ≤ 6τ(n) (quoted classical, not opened). So each d occurs for ≤ F := 6τ(4ℓQ₀) ≤ 48τ(P)
values s ≥ 0 (checked D1: max fibre 3 over 12 primes ℓ, Q₀ = 1, 2, 3, 6).
(iii) *Quadratic large sieve (quoted: Heath-Brown, Acta Arith. 72 (1995), Theorem 1 and Corollary 2; original not
opened; read (secondary) in arXiv 2505.09637 v2, eq. (10) and its restatement of Corollary 2):*
Σ_{χ∈S(Z)} |Σ_{n≤N′} a_n χ(n)|² ≪_δ (ZN′)^δ (Z + N′) Σ_{n₁n₂=□}|a_{n₁}a_{n₂}|, S(Z) the primitive quadratic characters of
conductor ≤ Z. With a_n = 1[n ∈ I] (primes, so n₁n₂ = □ iff n₁ = n₂), N′ = βX/P: Σ_d |S(d)|² ≪ X^(1+δ)·#I ≪ X^(2+δ).
(Equivalently from Theorem 1 over odd squarefree moduli: write |d| = 2^e m, m odd squarefree, e ∈ {0, 2, 3}; by
reciprocity χ_d(q) = (q/m)·λ with λ = ±1 depending only on q mod 8, m mod 4 and e; split q into 4 classes mod 8.)
(iv) H(|Δ_s|) ≪ √ℓ log²X log log X ([ML 6.1]); |U_{k−2}| ≤ k − 1; |ξ| ≤ 2^(ω(P)); #s ≪ √ℓ. Cauchy–Schwarz:
|Σ_q Leg| ≤ Σ_s c_s|S(d_s)| ≤ (Σ_s c_s²)^(1/2)(2F Σ_d |S(d)|²)^(1/2) ≪ (ℓ^(3/2)log⁵X)^(1/2) X^(1+δ) ≪ X^(7/4+2δ). ∎
The term with d small (d = −3, −4, …; S(d) as large as #I) is inside the mean square; d is never a square (d < 0).
Since the denominator of (2.1) is ≍ X²/log X (Step 6), **the fixed terms contribute O(X^(−1/4+δ)) to A^ε, without
GRH** (in [AL] the same statement in the other order of averaging needed GRH). Numerically the fixed part decays faster:
rms 0.255 at X = 2500 and 0.121 at X = 10⁴ (ratio 0.48, i.e. ≈ X^(−1/2); not proved).

**Step 4. Local Fubini (Lemma F; proved; checked B1–B3).** For the q-terms the discriminant is
D(ℓ, q) := Q₀q(4ℓ − s²Q₀q), quadratic in q at fixed ℓ, and D(λℓ, λq) = λ²D(ℓ, q).
*Lemma F.* Let v be a prime and Φ: ℤ_v → ℂ bounded measurable with Φ(c²z) = Φ(z) for all c ∈ ℤ_v^×. Then for every
ℓ₀, q₀ ∈ ℤ_v^×: ∫_{ℤ_v^×} Φ(D(ℓ₀, q))dq = ∫∫_{(ℤ_v^×)²} Φ(D) dℓ dq = ∫_{ℤ_v^×} Φ(D(ℓ, q₀))dℓ (Haar probability).
*Proof.* f(ℓ) := ∫Φ(D(ℓ, q))dq satisfies f(λℓ) = ∫Φ(D(λℓ, λq′))dq′ = ∫Φ(λ²D(ℓ, q′))dq′ = f(ℓ) (substitute q = λq′),
so f is constant on ℤ_v^×, equal to its mean, the joint mean; the same for the mean over ℓ. ∎
*Application.* The local factor F_v(D) (v ∤ P), the new local factor G_v(D) = ξ_Δ(v)F_v(D) (v | M) [AL Step 2(ii)], and
every coefficient T_v(i, j; D)·c_v of the (m, n)-expansion ([ML Step 2(iv)], [AL Step 2(iii)]) depend on D only
through its valuation and the class of its unit part modulo squares (for v = 2: c² ≡ 1 mod 8), so they are invariant
under D ↦ c²D. Hence **the mean over the free prime q at fixed ℓ of each local datum equals the mean over ℓ at a fixed
level**, computed in [ML Lemma 4] and [AL Step 3]: e_v(s) for v ∤ P, e_v^new(s) for v | M, and exactly 1 for v | Q₀
(v ∥ D for odd v | Q₀; v₂(Δ₀) ∈ {2, 3} with 2 ∤ f₀ for 2 | Q₀, as in [ML Lemma 3], whose proof uses only that ℓ, q
are odd and prime to Q₀). The joint local density over pairs (ℓ, q) is the product measure, and both orders of
averaging give its mean. *Checked:* B1 exact equality of the distributions of D mod v^J up to unit squares over q for
every unit ℓ and over ℓ at fixed q (v = 2, 3, 5, 7; 40,628 histogram comparisons); B2 q-means of F_v and G_v at fixed ℓ
inside rigorous enclosures of e_v(s), e_v^new(s) (440 cases, width ≤ 3.9·10⁻³); B3 F_v = 1 at v | Q₀ (158,290 cases).
The local factor at v = ℓ (a large prime) is 1 + O(1/ℓ); at v = q it is 1 (q ∥ D).

**Step 5. Key lemma in the free prime (proved by the listed changes of [ML Step 6] and [AL Step 4]).**
*For 0 ≤ s ≤ 2√(YP), G of bounded variation on [1, β] with G(t) = 0 when s²Q₀tX/P ≥ 4ℓ, ‖G‖ := sup|G| + Var G,
uniformly for primes ℓ ≤ YX, ℓ ∤ P:*

    Σ_{q∈I} G(Pq/X) ξ_{Δ_s(q)}(M) L*(D_s(q)) = E(s)(1 + O(ℓ⁻²)) Σ_{q∈I} G(Pq/X) + O(‖G‖ X (log X)^(−A)),

D_s(q) := Q₀q(4ℓ − s²Q₀q), Δ_s(q) = −D_s(q), L*(D) := πH(D)/√D = Σ_{m∈𝓜(D)}m⁻¹L(1, ξ_{−D/m²}) ([ML Lemma 2]).
Parameters as in [ML 6]: M_B = X^(7/8) (truncation), m₀, L, Q₀′ powers of log X. Note: if s ≥ 1 occurs then
ℓ > s²Q₀q/4 ≥ X/(4P) > M_B, so every prime ≤ M_B is prime to ℓ.
- [6.0] PNT for q ∈ I; the O(1) primes q | P or q = ℓ are removed at cost ≪ X log³X.
- [6.1] m ∈ 𝓜(D) is prime to Q₀q (Step 4), so m² | D forces m² | 4ℓ − s²Q₀q: for s = 0, m ∈ {1, 2}; for s ≥ 1 an odd
  v | (m, s) would divide 4ℓ, i.e. v = ℓ > s, impossible; so q lies in one class mod m_odd² (and mod a divisor of the 2-part): #{q ∈ I : m² | D}
  ≤ 4|I|/m² + 1. The bound |L(1, ξ)| ≪ log X log log X is unchanged.
- [6.2] Burgess (r = 2) for every q: conductor ≤ D_s(q) ≤ 4ℓQ₀q ≤ 4YβX²; error X^(−7/16)·X^(3/8+ε) = X^(−1/16+ε) per q.
- [6.3] T(m, n; q) := 1[m ∈ 𝓜(D)]·(−D/m² / n) is periodic in q mod 4nm² (times 8M² for the factors c_v, v | M).
  With d := Π{v odd : v_v(n) odd, v ∤ msP} and n″ := n/d♮: T = W(q)·σ·χ_d(q(4ℓ − s²Q₀q)), σ = (−Q₀/d) = ±1, W periodic
  mod K = 4n″m²·8M², |W| ≤ C(P). **The new ingredient** is the expansion of x ↦ χ_v(x(a − cx)), a = 4ℓ, c = s²Q₀,
  in characters ψ mod v (v | d odd, v ∤ ac): the coefficient is ψ̄(a/c)χ_v(c)J(χ_vψ̄, χ_v)/(v−1) for ψ ∉ {ψ₀, χ_v}
  (substitute x = (a/c)y; |J| = √v: Jacobi sums, quoted classical), and −χ_v(−c)/(v−1), −χ_v(a)/(v−1) for ψ₀, χ_v
  (complete quadratic sums with discriminant a² ≢ 0). So |â(ψ)| ≤ √d/φ(d), |â(ψ₀)| = 1/φ(d), exactly as in [ML 6.3(c)].
  Degenerate cases: v | s is excluded from d (kept in n″), and the quadratic has a double root mod v only if v | a, i.e.
  v ∈ {2, ℓ}: 2 ∤ d, ℓ > M_B ≥ d. For s = 0, χ_d(4ℓQ₀q) = χ_d(4ℓQ₀)χ_d(q) exactly (if ℓ | d, T ≡ 0). q ∤ d since q > M_B.
  Checked B4: 7,668 coefficients, v ≤ 17.
- [6.4] the per-n cut is unchanged (n″ ≤ L·rad(msP)·n₂², K ≤ (log X)^(C₃)).
- [6.5] Siegel–Walfisz for q in classes mod 4nm² ≤ (log X)^(B₂): main term e(m, n)Σ_{q∈I}G, e(m, n) the mean of
  T(m, n; ·) over unit classes; error ≪ ‖G‖X exp(−c√log X) per pair. The only ineffective step.
- [6.6] large sieve with a_q := G(Pq/X)1[q ∈ I] on an interval of length ≍ X, d ≤ M_B: piece
  ≪ ‖G‖X(log X)^(C)[Q₀′^(−1/2) + (M_B/X)^(1/2)], (M_B/X)^(1/2) = X^(−1/16).
- [6.7] by CRT e(m, n) = Π_v e_v(i, j) with e_v(i, j) the q-mean of T_v(i, j; D(ℓ, q)) (v ≠ ℓ). By Lemma F each
  e_v(i, j) equals the coefficient of the one-level setting, so the bounds |e_v(i, j)| ≤ 8v^(1−2i), the closed forms,
  the local sums Σ e_v(i, j)v^(−i−j) ∈ {e_v(s), e_v^new(s), 1} and the weighted absolute convergence of [ML 6.7] hold
  verbatim; pairs with ℓ | mn lie beyond the truncation (n ≤ M_B < ℓ, m ≤ m₀) when s ≥ 1; for s = 0 and small ℓ the
  local factor at ℓ is 1 instead of e_ℓ(0) = 1 + O(ℓ⁻²). Sum: E(s)(1 + O(ℓ⁻²)) + O((log X)^(−A)). ∎
*Checked C1* (empirical, not a proof): at ℓ = 20011, the means over all primes q in the ranges of Section 4 of
ξ_Δ(M)L*(D_s(q)) for P = 2, 3, Q ∈ {q, Pq}, s = 0, 1, 2 agree with E(s): 12 cases, max |z| 0.96.

**Step 6. The count (proved from [AL Step 8] and the dimension formula).** |H_k(N)| = (k−1)φ(N)/12 + O_k(2^ω(N))
(Martin, J. Number Theory 2005, quoted as in [ML Step 8], not opened) and |H^ε(N)| = 2^(−ω(N))|H_k(N)| + O(2^ω(N)√N log N)
([AL Step 8], from Prop. 5.5 at ℓ = 1; checked A2). Hence the denominator of (2.1) is
((k−1)φ(P)/12)Σ_{X≤Pq≤βX}(q − 1) + O(X^(3/2)log X) ≍ X²/log X. Normalisation: (2.1) and t_Q use √ℓλ_f(ℓ), as Martin.

**Step 7. Assembly (proved).** With H(D) = (√D/π)L*(D), Step 1 gives for Q = Q₀q

    t_{Q₀q}(ℓ) = −(1/2π) Σ_s w_s U_{k−2}(s√(Q₀q)/(2√ℓ)) √(D_s(q)) ξ_{Δ_s}(M) L*(D_s(q)) + μ(M)(ℓ + 1)[k = 2]   (w₀ = 1, w_s = 2).

Apply Step 5 with G_s(t) := X⁻¹ U_{k−2}(s√(Q₀tX/P)/(2√ℓ)) √(Q₀(tX/P)(4ℓ − s²Q₀tX/P))₊ (continuous, ‖G_s‖ ≪ 1) for each
of the O(1) values s: Σ_{q∈I} t_{Q₀q}(ℓ) = Σ_{q∈I} τ(q) + O(X²(log X)^(−A)), where τ(q) is the same expression with
ξL* replaced by E(s). The algebra of [AL Step 5] (√(D_s) = Q₀q√(4yM − s²) at ℓ = yPq; (−1)^(k/2) = −U_{k−2}(0);
12κ(N)/((k−1)N)·(k−1)φ(N)/12 = 1/N) gives τ(q) = ((k−1)φ(P)(q−1)/12)·D^[Q₀,M](ℓ/(Pq)) + O(1), using
κ(Pq) = κ(P)(1 + 1/(q−1)) and that E^(Q₀q,M)(s) of [AL] differs from E(s) by the factor e_q(s)⁻¹(1 − q⁻²)^[s=0] = 1 + O(q⁻²).
With Step 3 and Step 6:

    A^ε = ε_q Σ_{Q₀|P} ε(Q₀) Σ_{q∈I}(q−1)D^[Q₀,P/Q₀](ℓ/(Pq)) / Σ_{q∈I}(q−1) + O((log X)^(1−A) + X^(−1/4+δ))   (rename A).

*Weights → integral.* D^[Q₀,M] is continuous and of bounded variation on [0, Y]. By the prime number theorem with
error O(x exp(−c√log x)) and partial summation, Σ_{q∈I}(q−1)f(ℓ/(Pq)) / Σ_{q∈I}(q−1) =
∫₁^β f(x_ℓ/t) t w_X(t) dt / ∫₁^β t w_X(t) dt + O(exp(−c√log X)), x_ℓ = ℓ/X, w_X(t) = log X/log(tX/P) = 1 + O(1/log X).
Dropping w_X costs O(1/log X). ∎ (Theorem, error as stated.)

**Step 8. Uniformity in ℓ, continuity, x = 0 (proved).** Every bound above is uniform for primes ℓ ≤ YX, ℓ ∤ P (for
ℓ < X/(4P) only s = 0 occurs in Step 5, and the s = 0 case of Step 5 needs no lower bound on ℓ). D^[Q₀,M] is
continuous on [0, ∞) with D(0) = 0 (each s-term vanishes at y = s²/(4M); the s = 0 term is ∝ √y; the Eisenstein
term ∝ y) and locally bounded, so by dominated convergence M^ε(x; β) is continuous on [0, ∞) × (1, ∞) and M^ε(0; β) = 0.
Hence ℓ/X → x implies A^ε → M^ε(x; β), including x = 0 (checked F2: |M^(++)(10⁻⁴; 2)| = 0.040 < |M(10⁻³; 2)| = 0.128). ∎

## 3. Inputs

| input | where | status |
|---|---|---|
| Martin, Conj. 1.10, Rem. 1.11, Conj. 1.8, Thm 1.9, Prop. 5.5, §6.3 | §0, Step 1 | quoted, read (LaTeX source); Prop. 5.5 checked A1–A2 |
| Atkin–Lehner theory: commuting involutions on the new space, projector, a_f(v) = −v^(k/2−1)w_v(f) | Step 1 | quoted, not opened; checked A1, A2 |
| class-number formula for orders, L*(D) | Steps 3, 5 | quoted classical via [ML Lemma 2] (Cox, not opened); checked in [ML] A1–A3 |
| ξ_Δ, G_v, e_v^new; key lemma at one level for W_Q | Steps 4, 5 | proved in [AL Steps 2–4] (refereed once) |
| [ML] Step 6 machinery (6.0–6.7) | Step 5 | proved in [ML] (refereed three times); changes listed here, not refereed |
| Lemma F (homogeneity ⇒ local Fubini) | Step 4 | proved; checked B1–B3 |
| Jacobi sums of modulus √v; complete quadratic character sums | Step 5 [6.3] | quoted classical, not opened; checked B4 |
| Burgess (r = 2, any modulus) | Step 5 [6.2] | quoted as in [ML 6.2] (read (secondary)) |
| Siegel–Walfisz; PNT with de la Vallée Poussin error | Steps 5, 7 | quoted (not opened); ineffective |
| multiplicative large sieve | Step 5 [6.6] | quoted as in [ML 6.6] (read (secondary)) |
| Heath-Brown's quadratic large sieve (Acta Arith. 72 (1995), Thm 1, Cor. 2) | Step 3 | quoted; **read (secondary)** in arXiv 2505.09637 v2; original not opened |
| elements of norm n in O_K ≤ 6τ(n) | Step 3(ii) | quoted classical, not opened; checked D1 |
| Deligne's bound | Step 3(b) | quoted, not opened (for P = 2, 3, k = 2 the term is 0: A4) |
| dimension of S_k^new(N), O(2^ω) | Step 6 | quoted (Martin 2005, not opened) |

## 4. Numerics (own code; P = 2, k = 2, X = 10⁴, β = 2: all 560 levels 2q with 5003 ≤ q ≤ 9973, N ≤ 19946)

Data: A^ε from (2.1), traces by Prop. 5.5 with PARI class numbers, exact eigenspace sizes. Finite-X model:
t_Q → τ (Step 7, exact q, exact eigenspace sizes in the denominator), fixed terms → 0. Limit: M^ε(x; 2) (midpoint
rule). Errors: three interleaved groups of levels, se = sd/√3, pooled per ℓ over the four eigenspaces. Decision rule
written into the script before the run: χ²(own) ≤ 100 on 48 terms; every alternative ≥ own + 100.

| x = ℓ/X | 0.151 | 0.300 | 0.451 | 0.601 | 0.801 | 1.001 | 1.250 | 1.501 | 1.801 | 2.100 | 2.501 | 2.901 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| (++): A | −1.404 | 0.987 | −0.008 | −1.306 | −1.361 | −0.775 | −0.834 | −1.144 | −1.197 | −0.705 | 0.809 | 0.505 |
| (++): model_X | −1.420 | 1.050 | 0.186 | −1.372 | −1.543 | −0.891 | −0.944 | −1.020 | −1.158 | −1.042 | 0.641 | 0.729 |
| (++): M^ε(x; 2) | −1.366 | 1.007 | 0.191 | −1.317 | −1.495 | −0.863 | −0.891 | −0.982 | −1.138 | −1.009 | 0.593 | 0.718 |
| (−+): A | −0.277 | 0.348 | 1.184 | 1.976 | 1.437 | 0.464 | −0.792 | −1.363 | −0.848 | 0.703 | 0.939 | 1.050 |
| (−+): model_X | −0.227 | 0.342 | 1.096 | 1.702 | 1.214 | 0.365 | −0.736 | −1.029 | −0.554 | 0.625 | 1.177 | 0.683 |
| (−+): M^ε(x; 2) | −0.247 | 0.379 | 1.216 | 1.897 | 1.359 | 0.430 | −0.796 | −1.161 | −0.636 | 0.711 | 1.338 | 0.742 |
| se (pooled) | 0.090 | 0.110 | 0.131 | 0.069 | 0.190 | 0.258 | 0.239 | 0.369 | 0.261 | 0.383 | 0.491 | 0.520 |

(+−) and (−−): all 48 values in `martin_check_out.txt` (their data are close to −(++), −(−+); the differences are the ε_q-even finite-X terms).
χ² on 48 terms: **own finite-X model 39.2**, limit M^ε 18.2; W_N-only model 880.5; zero 3,632; the other eigenspaces'
models 3,404–14,400 (E1, E2 pass). rms(A − model_X) 0.164, rms(A − M^ε) 0.152, rms se 0.298 (the three-group errors
overstate the scatter). The fixed terms (Q = 1, 2): rms 0.123 against model rms 1.03 (E3). Finite-X effects visible:
eigenspace sizes differ from dim/4 by −4, +4, +11, −11 % (++, +−, −+, −−), and model_X − M^ε reaches 0.23; the data do
not distinguish model_X from M^ε at this X. Class numbers above 5·10⁵ come from PARI's class-group routine, which by its
documentation (recalled, not opened) may assume GRH: a caveat for the numerics only.
**Martin's figure (S_2(2q), 3000 < q < 6000, X = 6000; illustration, read by eye from `pq-profile3000.png`).**
M^ε(x; 2)/4 at ℓ = 250, 800, 2000, 4000 for (++, −−, +−, −+): (−0.21, −0.24, 0.21, 0.24), (−0.36, 0.03, 0.36, −0.03),
(0.27, −0.13, −0.27, 0.13), (−0.37, −0.45, 0.37, 0.45); the figure shows ≈ (−0.2, −0.25, 0.2, 0.25), (−0.36, 0.05,
0.36, −0.08), (0.25, −0.1, −0.28, 0.12), (−0.38, −0.4, 0.4, 0.5). So the plotted quantity is A^ε/4 = 2^(−r)A^ε, i.e.
(1/#ℱ)Σ_{f∈ℱ^ε}, not A^ε as defined in Conjecture 1.10. Support (read): Martin, arXiv 2505.01093, §"Local root
numbers": "Adding all 4 colors corresponds to looking at murmurations with no root number."

## 5. Relation to Martin's paper; what is and is not claimed

**What Martin proves** (read): Theorem 1.9: Type I (N = MQ, M fixed, Q over all squarefree Q prime to M), only for
x < 1/(4M) − ε, where only the s = 0 term of Prop. 5.5 is present: A^𝒬 → c√x + δ_{k=2}d (for k = 2 the Eisenstein term
δ_{k=2}μ(M)σ(ℓ) gives a term linear in x; as printed "d" is a constant: [AL §1, item 4], not rechecked here);
Corollary 6.3 (A⁺ + A⁻ → 0 for all squarefree levels); Lemma 6.2. **No case of Conjecture 1.10 is proved there**; §6.3
relates it heuristically to Conjecture 1.8. As printed, §6.3 sums over index sets {1, …, m} ⊂ I; his own criterion
(c_{I,X} ↛ 0 needs N/Q_I bounded) selects I ⊇ {m+1, …, r}, which is what the proof gives (Q ∋ q). Our reading of a
probable slip, not confirmed with the author.
**Added here (not found in print: Martin's papers 2409.02338 and 2505.01093 read; one web search on 26 Nov 2026 found
no proof; arXiv 2506.01640 seen as a search snippet only):** the Theorem (Conjecture 1.10 for r = m + 1, every m ≥ 0,
every even k, every ε, every x ≥ 0, no hypothesis, explicit M^ε); Corollary 2 (Conjecture 1.8 for the pairs
(Pq, Q₀q), unsmoothed); Lemma T′ without GRH. **Not claimed:** the trace formula (Martin), the one-level densities and
local factors ([ML], [AL]; the local recipe is in substance in Cardi, arXiv 2607.08982 §5, per [ML §7]), the
unconditional method (Zubrilina's trace route; the large-sieve treatment as in [ML] and arXiv 2603.25564), the
quadratic large sieve (Heath-Brown). For P = 1 (prime levels, root numbers) the statement may be in print in substance;
not searched.

**NOT closed.**
1. r ≥ m + 2 (two or more free primes, e.g. Martin's S_4(pq) with p, q both varying): not treated; the q-terms then
   need sums over two free primes of class numbers, and the fixed-prime structure of Step 3 is lost.
2. Conjecture 1.8 for pairs with q ∤ Q (Type II/III with weight √(N/Q) ≍ √X, e.g. √N·W₁): Lemma T′ gives only
   O(X^(−1/4+δ)) relative, times √X: no conclusion. Numerically the fixed part is ≈ X^(−1/2) (ratio 0.48 from X = 2500
   to 10⁴), consistent with a √X-normalised limit, which is not identified.
3. Rate: O(1/log X) against M^ε comes only from the prime density; O((log X)^(−A)) against the weighted M_X^ε; no power
   saving (Siegel–Walfisz in [6.5]); the constant is ineffective.
4. Short level intervals (β → 1, a "murmuration density" in Martin's Remark 1.11(1)): not written out. Step 5 holds for
   q-intervals of relative length h ≥ (log X)^(−B) by the same proof; Step 3 then gives O(X^(−1/4+δ)h^(−1/2)), enough for
   h ≥ X^(−1/2+3δ); neither is written out in detail.
5. Non-squarefree levels; weight k ≥ 4 numerics (the proof covers every even k; no numerics for k ≥ 4 here).
6. Nothing here is refereed; Step 5 is a list of changes to refereed steps, not a re-derivation in full.

## 6. Corrections after the referee report ([referee/referee_report.md](referee/referee_report.md); no failing step found)

1. Step 6.1: replace "#{q ∈ I : m² | D} ≤ 4|I|/m² + 1" by the bound over integers, ≤ βX/(P m²) + 1; the tail becomes ≪ X log²X loglog X / m₀², which is still admissible.
2. For s ≥ 1 the local factor at ℓ is exactly (1 − 1/ℓ)^(−1), since Δ_s ≡ (sQ₀q)² mod ℓ. The relative error is O(1/X), not O(ℓ^(−2)); the Burgess tail absorbs it.
3. The definition of d ("v ∤ msP") is empty for s = 0; that case is the separate s = 0 term.
4. If ℓ | P, then A = 0 = M(0) (not stated before).

Still open after the report: Heath-Brown's 1995 paper was read only in a secondary source (arXiv 2505.09637 v2, eq. (10) and Corollary 2); two or more free primes (r ≥ m + 2); a power saving.
