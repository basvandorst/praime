# Murmurations for Atkin–Lehner operators W_Q and on Atkin–Lehner eigenspaces at one squarefree level

Extension of `murmurations_one_level.md` (refereed three times; cited below as **[ML]**, its steps as [ML 6.2] etc.) from
the full involution W_N to W_Q, Q | N, and to the Atkin–Lehner eigenspaces of levels N = p₁q (p₁ fixed, q → ∞), which is
the setting of a conjecture of K. Martin. Check script: `al_check.py` (own code; output `al_check_out.txt`; **18 checks,
0 failed**). Blind test frozen, not run: `registration.txt`, `al_blind_test.py`, `predictions_blind.json`.

**Labels.** *proved*: argued here in full or as a listed change of a step of [ML]; *quoted*: from the literature, with
*read* (source text opened for this note), *read (secondary)*, or *not opened*; *checked*: an assertion of `al_check.py`;
**NOT closed**: a gap. "GRH" = GRH for Dirichlet L-functions. Nothing here is refereed.

## 0. Martin's statements (quoted; read in the LaTeX source of arXiv 2409.02338, file `Wqbias-corr.tex`; numbers computed from its counters: `conj`, `thm`, `rem` share the `lem` counter numbered within sections)

Notation there: ℱ(X, Y) = newforms of the family with level X ≤ N ≤ Y; Σ′ = only levels coprime to ℓ; w_Q(f) the
W_Q-eigenvalue; normalisation ℓ^(1−k/2)a_ℓ(f) = √ℓ·λ_f(ℓ) (the same as ours).

> **Conjecture 1.10** (Murmurations on Atkin–Lehner eigenspaces). Fix k, r and 0 ≤ m < r. Fix primes p₁ < ⋯ < p_m. Let 𝒩
> be the set of levels N = p₁⋯p_r ∈ ℕ^sqf_r such that p₁ < p₂ < ⋯ < p_r. Let ℱ be the family of weight k newforms with
> level in 𝒩. For ε = (ε₁, …, ε_r) ∈ {±1}^r, let ℱ^ε be the subset of f ∈ ℱ with Atkin–Lehner sign ε_i at p_i, and consider
> the averages A^ε_ℱ(ℓ, X; β) = (1/#ℱ^ε(X, βX)) Σ′_{f∈ℱ^ε(X,βX)} ℓ^(1−k/2) a_ℓ(f). These averages have murmurations which are
> scale invariant in ℓ/N. That is, as ℓ, X → ∞ such that ℓ/X → x, A^ε_ℱ(ℓ, X; β) → M^ε_ℱ(x; β) for a murmuration function
> M^ε_ℱ which is continuous on [0, ∞) × (1, ∞).

> **Remark 1.11(3).** […] the trace of T_ℓ on ℱ^ε(N) is a linear combination of the traces of T_ℓW_Q on ℱ(N) where one sums
> over all Q | N. Correspondingly, we expect the murmuration functions to be different on each Atkin–Lehner eigenspace when
> r = m+1, i.e., when all but one prime is fixed in the level as in Figures [S_2(2q), S_4(2q)].

Figure "Murmurations for AL-eigenspaces on S_2(2q)" (file `pq-profile3000.png`, commented caption: 3000 < q < 6000; blue
++, green −−, red +−, orange −+, first sign at 2). Also quoted: **Conjecture 1.8** (murmurations of
A^𝒬 = (1/#ℱ)Σ′_N Σ_f w_Q(f)√(N/Q) ℓ^(1−k/2)a_ℓ(f) for Types I (N = QM, M fixed), II (Q fixed), III), **Theorem 1.9**
(Type I, x < 1/(4M) − ε: "A^𝒬 → c√x + δ_{k=2}d … for some constants c, d"), **Proposition 5.5** (trace formula, Step 1
below), and §6.3: tr(T_ℓ | S_k^new(N)^ε) = 2^(−r)Σ_{Q|N} ε(Q) tr(T_ℓW_Q) (citing his refdim paper, Prop. 3.2, not opened).
**Martin's averages are over levels N ∈ [X, βX] at a fixed prime ℓ; the theorems below are at one level N, averaged over
primes p in windows. This is the main difference (Section 4).**

## 1. Results

**Setting.** N squarefree, k ≥ 2 even, Q | N, M := N/Q, H_k(N) the newforms of S_k^new(Γ₀(N)), p ∤ N prime,

    t_Q(p) := p^(1−k/2) tr(T_p W_Q | S_k^new(N)) = Σ_{f∈H_k(N)} w_Q(f) √p λ_f(p)      (t_1 = plain trace, t_N = W_N-trace).

Windows J = (y₀N, (y₀+h)N] and ⟨X⟩_J as in [ML §1] (weights log p, p ∤ N). κ(N) = N/φ(N), π(Q) = Π_{ℓ|Q}(1 + 1/ℓ),
e_ℓ(s), ν′(s) = Π_ℓ e_ℓ(s), u_ℓ as in [ML (2.1)]. New local factor at ℓ | M:

    e_ℓ^new(s) := −ℓ²/(ℓ² − 1)  (ℓ ∤ s),     −ℓ/(ℓ + 1)  (ℓ | s, s = 0 included).

**The W_Q-density.** For s ≥ 1 put E^(Q,M)(s) := ν′(s) Π_{ℓ|N} e_ℓ(s)^(−1) Π_{ℓ|M} e_ℓ^new(s), and
E^(Q,M)(0) := ζ(2) Π_{ℓ|Q}(1 − ℓ^(−2)) Π_{ℓ|M}(−(ℓ−1)/ℓ). With y = p/N:

    D_k^(Q,M)(y) := (12κ(N)/((k−1)π)) [ (−1)^(k/2) E(0) √(y/M) − (1/M) Σ_{1≤s<2√(yM)} U_{k−2}(s/(2√(yM))) E(s) √(4yM − s²) ]
                    + 12 κ(N) μ(M) y·[k = 2].                                                                 (1.1)

The s = 0 term equals (−1)^(k/2+ω(M)) (2π/(k−1)) π(Q) √(y/M) (checked D2). The natural variable is x = p/Q = yM: the
s-th term starts at x = s²/4. **For M = 1, D^(N,1) = (−1)^(k/2) M_k^(N) of [ML (2.1)]** (checked D1, 1,600 points, 6·10⁻¹⁴).

**Theorem A (W_Q at one level).** Fix k, Y, A, M₀. (a) *No hypothesis*: for every B ≥ 0, uniformly for squarefree N ≥ N₀
with H_k(N) ≠ ∅, Q | N with M = N/Q ≤ M₀, 0 ≤ y₀ ≤ Y and (log N)^(−B) ≤ h ≤ 1,
⟨t_Q / |H_k(N)|⟩_J = (1/h)∫_{y₀}^{y₀+h} D_k^(Q,M)(y) dy + O((log N)^(−A)) (ineffective constant, as in [ML]).
(b) *Under GRH* the same for N^(−1/2+δ) ≤ h ≤ 1.

**Theorem B (eigenspaces, N = mq).** Fix k, Y, A and m squarefree; N = mq, q prime → ∞, H_k(N) ≠ ∅.
(i) *No hypothesis* (windows as in A(a)). For a sign vector ε on the primes of m,
Σ_{f: w_ℓ(f)=ε_ℓ ∀ℓ|m} w_q(f)√pλ_f(p) = 2^(−ω(m)) Σ_{Q₀|m} ε(Q₀) t_{Q₀q}(p), and
⟨Σ_{Q₀|m} ε(Q₀) t_{Q₀q}/|H_k(N)|⟩_J = (1/h)∫ Σ_{Q₀|m} ε(Q₀) D_k^(Q₀q, m/Q₀) + O((log N)^(−A)).
For m = p₁: **⟨(t_q + ε₁t_N)/|H_k(N)|⟩_J → (1/h)∫(D^(q,p₁) + ε₁ D^(N,1))** (the part of the eigenspace sums that is odd in the
sign at q).
(ii) *Under GRH* (Lemma T, Step 7), for every Q₀ | m, δ > 0 and N^(−1/50+δ) ≤ h ≤ 1: ⟨t_{Q₀}/|H_k(N)|⟩_J = O(N^(−1/50)h^(−1) + N^(−1/2)).
Hence for every sign vector ε on the primes of N, the eigenspace mean E^ε(p) := |H^ε|^(−1) Σ_{f∈H^ε} √p λ_f(p) satisfies

    ⟨E^ε⟩_J = (1/h)∫ Σ_{Q₀|m} ε(Q₀q) D_k^(Q₀q, m/Q₀)(y) dy + o(1);   for N = p₁q:  E^(ε₁,ε_q) → ε_q D^(q,p₁) + ε₁ε_q D^(N,1).

(iii) *No hypothesis*: |H^ε| = 2^(−ω(N))|H_k(N)| + O_k(2^(ω(N))√N log N) (Step 8).

**Densities for N = 2q and N = 3q (q → ∞; proved from (1.1)).** With c₂(r) = 3/2, 3/4 (r odd, even) as in [ML], and

    D_q^[p₁](y) := (−1)^(k/2+1)(2π/(k−1))√(y/p₁) + (12/((k−1)π(p₁−1))) Σ_{1≤s<2√(p₁y)} U_{k−2}(s/(2√(p₁y))) c_{p₁}(s) ν′(s) √(4p₁y − s²)
                   − (12p₁/(p₁−1)) y [k = 2],    c_{p₁}(s) = p₁(p₁−1)/(p₁²−p₁−1) (p₁ ∤ s),  (p₁−1)/p₁ (p₁ | s)

(c₂ = 2, 1/2; c₃ = 6/5, 2/3), and D_N = (−1)^(k/2) M_k^(p₁q). For **N = 2q, k = 2**:

    D_q(y) = π√(2y) + (12/π) Σ_{s<2√(2y)} c(s)ν′(s)√(8y − s²) − 24y,   c(s) = 2 (s odd), 1/2 (s even);
    D_N(y) = −3π√y − (24/π) Σ_{r<2√y} c₂(r)ν′(r)√(4y − r²) + 24y;
    (+,+) = D_q + D_N = (√2 − 3)π√y + (12/π)[Σ_s c(s)ν′(s)√(8y−s²) − 2Σ_r c₂(r)ν′(r)√(4y−r²)],   (+,−) = −(+,+),
    (−,+) = D_q − D_N = (√2 + 3)π√y + (12/π)[Σ_s c(s)ν′(s)√(8y−s²) + 2Σ_r c₂(r)ν′(r)√(4y−r²)] − 48y,   (−,−) = −(−,+)

(signs (ε₂, ε_q)). Values (k = 2, q → ∞; checked D3 printout):

| y | 0.05 | 0.1 | 0.2 | 0.5 | 1 | 1.5 | 2 | 3 |
|---|---|---|---|---|---|---|---|---|
| 2q: D_q | −0.207 | −0.995 | 0.827 | −0.719 | −2.425 | 1.440 | −0.458 | −2.774 |
| 2q: D_N | −0.907 | −0.580 | 0.585 | −1.713 | 2.366 | −1.273 | 1.924 | −2.291 |
| 2q: (+,+) | −1.114 | −1.575 | 1.412 | −2.432 | −0.059 | 0.167 | 1.466 | −5.065 |
| 2q: (−,+) | 0.701 | −0.415 | 0.242 | 0.994 | −4.791 | 2.714 | −2.382 | −0.483 |
| 3q: D_q | −0.089 | −0.022 | −0.310 | 0.705 | −0.094 | 0.612 | 0.116 | −1.012 |
| 3q: (+,+) | −1.062 | −0.872 | −0.456 | 0.022 | 3.017 | −1.688 | −0.716 | −1.769 |
| 3q: (−,+) | 0.884 | 0.827 | −0.163 | 1.388 | −3.205 | 2.912 | 0.948 | −0.255 |

**What this says about Martin's conjecture (r = 2, m = 1).**
1. The four eigenspace densities are ±D_q ± D_N: **four different functions, as Remark 1.11(3) expects**, in two pairs of
   opposite functions: flipping the sign at the large prime q negates the density (exact for the odd part, Theorem B(i),
   no hypothesis; for each eigenspace given Lemma T, i.e. under GRH). Flipping the sign at p₁ exchanges D_q + D_N and D_q − D_N.
2. The only possible limit function. Weighting levels N = tX ∈ [X, βX] by their dimension (∝ t), the windowed theorem gives
   M^ε(x; β) = ∫₁^β D^ε(x/t) t dt / ∫₁^β t dt. If Conjecture 1.10 holds (pointwise in ℓ, continuous in x), its limit
   must be this function (average the conjectured limit over ℓ in windows). **Conjecture 1.10 itself (a fixed ℓ, average
   over levels) is NOT proved here** (Section 4).
3. Figure `pq-profile3000` (X = 6000, β = 2), read by eye against M^ε/4 (checked D4 printout): ℓ = 250: (++, −−, +−, −+)
   predicted −0.21, −0.24, 0.21, 0.24, figure ≈ −0.2, −0.25, 0.2, 0.25; ℓ = 800: −0.36, 0.03, 0.36, −0.03, figure ≈ −0.36,
   0.05, 0.36, −0.08; ℓ = 2000: 0.27, −0.13, −0.27, 0.13, figure ≈ 0.25, −0.1, −0.28, 0.12; ℓ = 4000: −0.37, −0.45, 0.37,
   0.45, figure ≈ −0.38, −0.4, 0.4, 0.5. Shapes, signs and the ± pairing agree; **the plotted values are 1/4 = 2^(−r) of
   the eigenspace averages A^ε as defined in the conjecture**. Martin's §6.3 approximation A^ε ≈ 2^(−r)Σ_I ε(Q_I)c_I A^(Q_I)
   carries the same 2^(−r), i.e. a division by #ℱ instead of #ℱ^ε; so the figure appears to use that normalisation
   (inference, not stated in the paper as read).
4. Theorem 1.9 (Type I, x < 1/(4M)): (1.1) gives there c√x plus, for k = 2, the term 12κ(N)μ(M)y, which is **linear** in
   x (it comes from δ_{k=2}μ(M)σ(ℓ) in Prop. 5.5). As printed, Theorem 1.9 has "δ_{k=2}d" with d a constant; our
   computation gives d·x. (Likely a slip in the statement; not checked against the arXiv PDF.)

**The scale question.** In t_Q the discriminants are 4pQ − s²Q², s² < 4p/Q, so the variable is p/Q. For N = p₁q:
Q = N: x = y; Q = q: x = p₁y (a murmuration at p ≍ N, first r-term at y = 1/(4p₁)); Q = p₁ (and Q = 1): x = p/p₁ ≍ N, an
unbounded number (≍ √N) of terms. At p ≍ N these terms are not a murmuration: exactly (Step 6)
t_{p₁}(p) = Leg_{p₁}(p) − t^(level p₁)_{p₁}(p), where the second term is O(√p) by Deligne (0 for p₁ ∈ {2, 3, 5, 7, 13},
k = 2) and Leg is a sum of ≍√N class numbers twisted by the Legendre symbol mod q; its window mean is o(|H_k(N)|) under
GRH (Lemma T) and numerically ≈ N^(−1/2)|H_k(N)| (Section 3). The scale p ≍ p₁ (bounded p) is Martin's small-ℓ regime
(his §5.1), not studied here.

## 2. Proof

**Step 1. Trace formula on the new space (quoted: Martin, Prop. 5.5, read; checked C1–C3).** For N = QM squarefree and a
prime p ∤ N:

    t_Q(p) = −½ Σ_{s∈ℤ, s²Q<4p} U_{k−2}(s√Q/(2√p)) ξ_{Δ_s}(M) H(4pQ − s²Q²) + μ(M)(p + 1)·[k = 2],   Δ_s := s²Q² − 4pQ,

ξ_Δ multiplicative on squarefree M, ξ_Δ(ℓ) = (Δ/ℓ) − 1 if ℓ² ∤ Δ, ℓH(|Δ|/ℓ²)/H(|Δ|) − 1 if ℓ² | Δ (H(n) = 0 for
n ≡ 1, 2 mod 4). Martin derives it from Skoruppa–Zagier (2.7) (not opened) and the trace of W_q on an old class π_g^N,
σ₀(M/M₀)a_ℓ(g)w_q(g) (Atkin–Lehner (5.1)–(5.2); read as restated by Martin, §3), followed by Möbius inversion with
μ∗μ. **Old forms do not cancel for Q ≠ N** (an old class contributes σ₀(M/M₀)a_ℓ(g)w_Q(g)); the inversion removes them
exactly, and this is where the factor ξ_Δ(M) = Π_{ℓ|M}((1 + ξ⁰_Δ(ℓ)) − 2) comes from. So **no gap from old forms**: the
theorems use only the new-space formula. *Checked:* C1 exact equality with PARI trace forms
(tr(T_pW_Q) = (−1)^(ω(Q))Q^(1−k/2)a_{Qp}(trace form)) for N ∈ {6, 10, 14, 15, 21, 30, 42}, k = 2, 4, all Q | N, primes
p ≤ 150: 2,344 triples; C2 against PARI's Atkin–Lehner matrices (mfatkininit), independent of a_f(q) = −q^(k/2−1)w_q(f):
588 triples; C3 tr W_Q (p = 1) and the eigenspace traces against PARI newform orbits: 340 cases.

**Step 2. Local structure (proved; checked A2, A3).** (i) Lemma 3 of [ML] holds with Q in place of N for
D = Q(4p − s²Q): every prime ℓ | Q divides Δ₀ exactly, F_ℓ = 1, 𝓜(D) is prime to Q, T(m, n; p) = 0 if (n, Q) > 1 (same
proof; A3: 73,656 cases). (ii) For ℓ | M write −D = Δ₀f₀², e := v_ℓ(f₀), χ := (Δ₀/ℓ). Since ℓH(|Δ|/ℓ²)/H(|Δ|) =
F_ℓ(Δ/ℓ²)/F_ℓ(Δ) (the other local factors of L*(D) are unchanged), ξ_Δ(ℓ)F_ℓ(Δ) = F_ℓ(Δ/ℓ²) − F_ℓ(Δ) or ((Δ/ℓ) − 1)F_ℓ(Δ),
and both equal

    G_ℓ(Δ) := ℓ^(−e)(χ − 1)/(1 − χ/ℓ)      (0 if ℓ splits, −2ℓ^(−e)/(1 + 1/ℓ) if inert, −ℓ^(−e) if ramified).

This reproduces Martin's explicit ξ (A2: 8,000 cases with class numbers, D ≤ 4000, ℓ = 2, 3, 5, 7). (iii) In the
(m, n)-expansion of [ML Step 2(iv)], ξ_Δ(ℓ)H(D) = Σ_{m∈𝓜(D)} h′(−D/m²) c_ℓ(m, D) with c_ℓ = ℓ·1[ℓ | m] − 1 (ℓ² | D) or
(Δ/ℓ) − 1 (ℓ² ∤ D), because H(|Δ|/ℓ²) collects the terms m ∈ 𝓜(D) with ℓ | m. So T_M(m, n; p) := T(m, n; p)Π_{ℓ|M}c_ℓ
is periodic in p modulo nm²·8M², bounded by Π_{ℓ|M} max(ℓ − 1, 2), and its local sums at ℓ | M equal G_ℓ.

**Step 3. Local means at ℓ | M (Lemma B; proved; checked B1).** ∫_{ℤ_ℓ^×} G_ℓ(Q(s²Q − 4p)) dp = e_ℓ^new(s).
*Proof.* ℓ odd, ℓ ∤ s: Δ = Qu′ with u′ := s²Q − 4p uniform on ℤ_ℓ minus the disc s²Q + ℓℤ_ℓ (density ℓ/(ℓ−1));
multiplication by the unit Q preserves Haar measure. Over all of ℤ_ℓ (x = 1/ℓ): v(Δ) = 2a with χ = ±1 equally (mean
−x^a/(1 + x)), v = 2a + 1 with χ = 0 (G = −x^a): I = −(1−x)Σ_a[x^(3a)/(1+x) + x^(3a+1)] = −1/(1 + x). On the removed disc
Δ ≡ (sQ)² is a unit square, χ = 1, G = 0. So the mean is (ℓ/(ℓ−1))·(−ℓ/(ℓ+1)) = −ℓ²/(ℓ²−1). ℓ | s: Δ ≡ −4pQ is a unit,
χ = ±1 equally: mean −1/(1+x) = −ℓ/(ℓ+1). ℓ = 2 (Q odd): s odd: Δ ≡ 5 (mod 8), χ = −1, G = −4/3. s = 2s′, s′ even:
Δ = 4w, w uniform on ℤ₂^×: G = 0, −2/3, −1 for w ≡ 1 (8), 5 (8), 3 (4): mean −2/3. s′ odd: Δ = 8w″, w″ uniform on ℤ₂:
G = −1 (w″ odd); for w″ = 2w‴, G = F(4w‴) − F(16w‴) = F(4w‴)/2 − 1 with mean (4/3)/2 − 1 = −1/3 ([ML Lemma 4]: mean of
F₂(4w′) is 4/3): mean −2/3. ∎ *Checked B1:* exact enumeration of p mod ℓ^J with enclosure (ℓ = 2, 3, 5, 7; nine Q per ℓ,
including Q ≡ 1, 3, 5, 7 mod 8; s = 0…8): 324 cases, width ≤ 2·10⁻³. Consistency: e^new + 2e_ℓ is the mean of the
root-count factor (1 + ξ⁰)F on the full space (ℓ/(ℓ−1) resp. ℓ(ℓ−2)/(ℓ−1)²), as "new = full − 2·old" requires.

**Step 4. Key lemma for (Q, M), M ≤ M₀ (proved by listed changes of [ML Steps 5–7]).** For 0 ≤ s < 2√((Y+1)M₀) and G
as in [ML Step 5], Σ_{p∈J, p∤N} log p G(p/N) ξ_{Δ_s(p)}(M) L*(D_s(p)) = E^(Q,M)(s)∫_J G(u/N)du + O(‖G‖hN(log N)^(−A))
(a) without hypothesis for h ≥ (log N)^(−B), (b) under GRH for h ≥ N^(−1/2+δ); constants depend on k, Y, A, B, M₀.
Changes: [6.0–6.2] unchanged (|T_M| ≤ C(M₀)|T|; Burgess per p; conductor |Δ₀| ≤ 4(Y+1)NQ ≤ 4(Y+1)N²). [6.3(a)] period
nm²·8M² instead of nm²; T_M = 0 unless (mn, Q) = 1; the moduli are prime to Q but may share primes with M (harmless:
p ∤ M). [6.3(b)] D = Qu with u = 4p − s²Q, the sign (−Q/q′); all primes of M are kept in n″ (the factors c_ℓ depend on
p only M-adically). [6.4] n″ ≤ L·rad(msM)·n₂², so K ≤ (log N)^(C₃)·C(M₀). [6.5, 6.6] unchanged. [6.7] the local sum at
ℓ | M is Σ_{i,j} e^M_ℓ(i, j)ℓ^(−i−j) = ∫G_ℓ = e_ℓ^new(s) (dominated convergence: |T^M_ℓ| ≤ ℓ + 1 and
P(ℓ^(2i) | D) ≤ 8ℓ^(1−2i)); at ℓ | Q it is 1 (Step 2(i)); at ℓ ∤ N it is e_ℓ(s) ([ML Lemma 4] uses only that the level is an
ℓ-adic unit; here Q is). Weighted absolute convergence: finitely many ℓ | M, each with a bounded factor. [7.1–7.4] F_ℓ is
replaced by G_ℓ at ℓ | M, whose expansion Σ_{j≥e}(χ^(j−e+1) − χ^(j−e))ℓ^(−j) has coefficients ≤ 2 depending on D mod
ℓ^(2j+2) (mod 2^(2j+6) at 2). **The unconditional steps carry over** (they never used the sign or the size of the local
factor beyond boundedness).

**Step 5. Assembly of Theorem A (proved; checked D1, D2).** H(D) = (√D/π)L*(D) [ML Lemma 2] turns Step 1 into
t_Q(p) = −(1/(2π))Σ_s U(…)√(4pQ − s²Q²)·ξL*(D_s) + μ(M)(p+1)[k=2]; with p = yN, √(pQ) = N√(y/M),
√(4pQ − s²Q²) = (N/M)√(4yM − s²), |H_k(N)| = (k−1)φ(N)/12 + O(2^ω) and PNT as in [ML Step 8], Step 4 gives (1.1). The
s = 0 term: κ(N)φ(M)/M = κ(Q) and κ(Q)Π_{ℓ|Q}(1 − ℓ^(−2)) = π(Q). ∎

**Step 6. Eigenspaces and the exact splitting (proved; checked C3, C4).** On S_k^new(N) the W_ℓ (ℓ | N) are commuting
involutions commuting with T_p (Atkin–Lehner 1970, not opened; C3), so the projector on the ε-eigenspace is
Π_{ℓ|N}(1 + ε_ℓW_ℓ)/2 = 2^(−ω(N))Σ_{Q|N} ε(Q)W_Q; restricting the product to ℓ | m and multiplying by W_q gives the identity
of Theorem B(i), and Theorem A for Q = Q₀q (M = m/Q₀ ≤ m) proves B(i). For Q₀ | m and q² > 4(Y+1)mQ₀N/q (so q² ∤ Δ_s),
ξ_Δ(q) = (Δ/q) − 1; splitting off the "−1" and using Step 1 at level m:

    t_{Q₀}(p) = Leg_{Q₀}(p) − t^(m)_{Q₀}(p),   Leg_{Q₀}(p) := −½ Σ_s U_{k−2}(s√Q₀/(2√p)) ξ_{Δ_s}(m/Q₀) (Δ_s/q) H(4pQ₀ − s²Q₀²),

t^(m)_{Q₀} the same trace on S_k^new(m) (the Eisenstein terms cancel: μ(Mq) = −μ(M)). Deligne (quoted, not opened):
|t^(m)_{Q₀}(p)| ≤ 2 dim S_k^new(m)·√p. C4: 1,424 exact cases.

**Step 7. Lemma T (proved given GRH; written out here, not refereed).** *Under GRH, for N = mq, Q₀ | m, 0 ≤ y₀ ≤ Y, h ≤ 1:
Σ_{p∈J} log p·Leg_{Q₀}(p) ≪_{k,Y,m} N^(2−1/50).* Since Σ_{p∈J}log p·|H_k(N)| ≍ hN², this gives B(ii).
*Proof.* Write Leg = Σ_s G_s(p)χ_q(Δ_s)·ξL*(D_s(p)) with G_s(p) := −U(…)√(4pQ₀ − s²Q₀²)/π (s ≥ 1 doubled), D_s ≤ 4(Y+1)mN.
(1) [ML 7.1–7.2] under GRH for the real characters χ_{Δ₀} (conductor ≤ 4(Y+1)mN): ξL*(D) = Σ_{n≤L, P⁺(n)≤z} β_D(n)/n
+ O((log N)^(−A−2)) uniformly, L = N^(1/300), log z = log L/((A+3) log log N); β_D(n) depends on D mod K(n) ≤ 16n⁴, hence
on (p mod K, s mod K); z < q, so the factor at q is not touched. (2) Smooth the edge: replace √u₊ by g(u) = √u (u ≥ ηN),
u/√(ηN) (0 ≤ u ≤ ηN), η = N^(−1/24); this changes ≪ η√N + 1 values of s per p by ≤ √(ηN) each: total ≪ η^(3/2)N²log²N
= N^(2−1/16)log²N. (3) For a class (b, c) mod K: 1[p ≡ b] = φ(K)^(−1)Σ_η η̄(b)η(p), and for q ∤ s (s < q)
χ_q(4p − s²Q₀) = Σ_{ψ mod q} â_s(ψ)ψ(p) with â_s(ψ) = χ_q(Q₀)ψ̄(Q₀/4)ψ̄²(s)â*(ψ), â*(ψ) = φ(q)^(−1)Σ_{y≠0}χ_q(y−1)ψ̄(y)
(substitute p = s²Q₀y/4); |â*(ψ)| = √q/φ(q) (Jacobi sum) for ψ ∉ {ψ₀, χ_q}, = 1/φ(q) for ψ ∈ {ψ₀, χ_q} (these two
contribute ≪ (1/q)·hN·N ≪ N). (4) For ψ ∉ {ψ₀, χ_q} and every η, ψη is non-principal mod qK; GRH (θ(x, χ) ≪ √x log²(qKx),
Davenport ch. 20, not opened) and partial summation: F(s) := Σ_{p∈J}log p (ψη)(p)G_s(p) ≪ ‖G_s‖√N log²N and
F(s) − F(s+K) ≪ ‖G_s − G_{s+K}‖√N log²N, with ‖G_s‖ ≪ √N and ‖G_s − G_{s+K}‖ ≪ K/√η after step (2). (5) Abel summation
over s ≡ c (mod K) and Burgess (r = 2, quoted as in [ML 6.2]) for the primitive character ψ̄² mod q on the progression
(length ≤ 2√((Y+1)N)/K ≥ q^(3/8+ε)): |Σ_{s≡c} ψ̄²(s)F(s)| ≪ (√N/K)^(1/2)q^(3/16+ε)·√N log²N·√N η^(−1/2). (6) Summing over
ψ (Σ|â*| ≤ √q), η (weight φ(K)^(−1)), the K·φ(K) classes and n ≤ L (Σ_n n^(−1)K(n)^(3/2) ≪ L⁶):
≪ N^(1/2+1/4+3/16+1+ε) L⁶ η^(−1/2) = N^(94/48+6/300+ε) ≤ N^(2−1/50). The s = 0 term is one character sum (≪ N^(1+ε)). ∎
The GRH families: characters mod qK, K ≤ 16N^(1/75) (all characters of these moduli), and χ_{Δ₀}, conductor ≤ 4(Y+1)mN.
The saving N^(−1/50) is far from the observed N^(−1/2) (Section 3).

**Step 8. Eigenspace dimensions (proved on Prop. 5.5; checked C3).** |H^ε| = 2^(−ω)Σ_Q ε(Q) tr(W_Q | S_k^new(N)); for Q > 1
Prop. 5.5 with ℓ = 1 has s²Q ≤ 4 only, |ξ_Δ(ℓ)| ≤ 2 (from G_ℓ/F_ℓ), so |tr W_Q| ≪_k 2^(ω(M))H(4Q) + 1 ≪ 2^(ω)√Q log Q.
Then E^ε = (|H_k(N)|/(2^ω|H^ε|))·Σ_Q ε(Q)t_Q/|H_k(N)| with ratio 1 + O(4^ω N^(−1/2)log N). With B(i) and Lemma T this
proves B(ii). ∎

## 3. Numerics (discovery, declared; weight 2; every p ≤ 4N, p ∤ N; 16 windows of y = p/N; log p weights)

Levels 2·5003, 2·5009, 2·5011 and 3·3331, 3·3343, 3·3347 (seen before in the project as W_N discovery levels). Data
from Step 1 (vectorised; equal to the exact formula on 484 values, E0), model = (1.1) evaluated at each prime with the
exact level constants and exact dimensions. Errors: one scale per type c_T × the single-value standard error, c_T from
the scatter of the residuals between the three levels (c_T = 0.67 (2q), 0.74 (3q); 0.54, 0.74 for the odd parts).

| type | eigenspaces, χ² on 192 terms: own | W_N only (no W_q term) | other eigenspace densities | odd parts (t_q ± t_N)/dim, 96 terms: own / W_N only |
|---|---|---|---|---|
| 2q | **254.4** | 1,714 | 9,065–9,572 | **87.1** / 1,522 |
| 3q | **221.0** | 688 | 8,175–9,566 | **73.7** / 364 |

The W_q-term is seen at both types (Δχ² 1,460 and 467 against the W_N-only model). The symmetric parts
(t₁ ± t_{p₁})/dim (expected → 0, Lemma T) have window-mean rms 0.435 (2q), 0.315 (3q) at N ≈ 10⁴ and 1.455, 1.202 at
N ≈ 10³ (levels 2·{503, 509, 521}, 3·{337, 347, 349}): ratios 0.30, 0.26, against 0.32 for N^(−1/2). They are **not pure
noise**: correlation between neighbouring levels −0.22 / 0.59 (2q, signs +/− at p₁) and 0.54 / 0.77 (3q). A secondary
term of relative size ≈ N^(−1/2) is plausible (Martin's √N-normalised W₁ murmurations, Conjecture 1.8(2), are of this
kind); it is not identified here. It is the main reason why χ²_own exceeds the number of terms.

## 4. Status

| step | status |
|---|---|
| Martin's Conjecture 1.10, Remark 1.11, Theorem 1.9, Prop. 5.5, §6.3 | quoted, read (LaTeX source) |
| Step 1 trace formula for T_pW_Q on the new space, old forms removed by μ∗μ | quoted (Martin Prop. 5.5, read; SZ (2.7) and Atkin–Lehner not opened); checked C1–C3 exactly |
| Step 2 primes of Q in Δ₀; G_ℓ = ξF; the c_ℓ(m, D) expansion | proved; checked A2, A3 |
| Step 3 local means e^new_ℓ(s) | proved; checked B1 |
| Step 4 key lemma for M ≤ M₀ | proved (listed changes of [ML 6–7]); not refereed |
| Theorem A (a) no hypothesis, (b) GRH | proved on the inputs of [ML §4] plus Prop. 5.5; checked D1, D2 |
| Theorem B(i), odd part in the sign at q, no hypothesis | proved |
| Step 6 splitting t_{Q₀} = Leg − t^(m) | proved; checked C4 |
| Lemma T (symmetric parts → 0) under GRH | proved (written out, not refereed) |
| **Lemma T without GRH** | **NOT closed** |
| Theorem B(ii), each eigenspace | proved given GRH |
| Conjecture 1.10 as stated (fixed ℓ, average over levels N = p₁q, q prime) | **NOT closed** |
| rate N^(−1/2) and the systematic secondary term of the symmetric parts | **NOT closed** (numerics only) |

**Not closed, precisely.** (1) *Lemma T unconditionally.* The odd part needs nothing; the even part needs
Σ_{p∈J} log p Σ_s G_s(p)(Δ_s(p)/q)ξL*(D_s(p)) = o(N·hN). Without GRH, steps 6.1, 6.2, 6.6 of [ML] carry over (the factor
(Δ_s/q) can be absorbed into the large-sieve sequence a_p), but 6.5 then needs, for moduli P ≤ (log N)^B and uniformly in
the shift c = s²Q₀ (s up to 2√((Y+1)N/Q₀)), a bound Σ_{p∈J, p≡b (P)} log p·((4p − c)/q) = o(hN/φ(P)) for a modulus q
comparable to the length hN — a Karatsuba-type estimate for characters of shifted primes (known for lengths ≥ q^(3/4+ε) by
Vinogradov's method; quoted from memory, not opened, uniformity in P and c not checked); and 6.3–6.4 must be redone for
shifts s up to √N (primes of s give characters of p, not of the shifted prime). (2) *Martin's order of averaging.* At a
fixed ℓ, the W_q- and W_N-terms need sums over primes q of H(q(4ℓ − s²q)) (a quadratic polynomial in q); the W₁- and
W_{p₁}-terms need Σ_q (D/q) over primes q for ≍√ℓ fixed discriminants D, which a quadratic large sieve over D (Heath-Brown
1995; not opened) should control on average. Neither is done. (3) Levels with two or more large primes (r > m + 1) and
non-squarefree levels: not treated. (4) Weight k ≥ 4: theorems cover it; no numerics here.

**Literature (prior-art status).** Gate of the project (26 Nov 2026) and this reading: Martin 2409.02338 has the trace
formula, the conjectures, Theorem 1.9 (s = 0 term only, levels averaged) and the plots; no density for the eigenspaces,
no one-level statement. Not found in print (in the sources opened for [ML] and here; no new literature search was made):
the W_Q-density (1.1) with the factors e^new_ℓ, Theorem A at one level, Theorem B and the identification of Martin's four
eigenspace functions as ±D_q ± D_N.

## 5. Blind test (frozen, not run)

Levels 2q for the ten primes q > 10000 other than 10007 (20018 … 20198) and 3q for the ten primes q > 6666.7 other than
6673 (20037 … 20283); weight 2; p ≤ 4N; 16 windows; four eigenspaces; own density against W_N-only and the other three
eigenspace densities; odd parts against W_N-only. Decisions, hashes, commands: `registration.txt`.
