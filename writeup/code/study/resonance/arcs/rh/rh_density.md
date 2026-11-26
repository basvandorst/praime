# Gaussian primes in a residue class in sectors of width x^(−θ): a second proof, by zero density, that uses Coleman's zero-free region and nothing else from unread sources

Labels: **proved** (complete argument written here), **proved in rh_proof.md** (lemma of
[rh_proof.md](../rh_proof.md), refereed in [referee_rh.md](../referee/referee_rh.md), cited by number and used inside
its hypotheses), **quoted (checked)** (statement read in the source), **quoted (secondary)** (statement known only
through an opened paper that quotes it; verbatim or not is said), **classical** (textbook fact, with reference, not
re-opened), **computed** (`rh_density_check.py`, output `rh_density_check_out.txt`), **heuristic**. Nothing here is
called new.

## Verdict

**Yes: among unread sources the theorem now rests on Coleman's Theorem 2 alone**, the zero-free region as quoted word
for word by Sun–Zhang (arXiv 2008.09677). Theorem 1 of rh_proof.md, which is Hypothesis R_h of arcs_theorem.md §2 in the
exact form used there (norm ≤ x, sharp sector, fixed anchor, uniform in the anchor and in the width), is proved below
from

- (i) Coleman's Theorem 2 (zero-free region of Vinogradov–Korobov type, with the modulus) — **quoted (secondary)**;
- (ii) a pointwise bound P_κ: Wu's Theorem 1.1 (arXiv 1604.08551) with ϑ = 7/64 (Blomer–Brumley, arXiv 1003.0559,
  Theorem 1) — both **quoted (checked)**; or Michel–Venkatesh (read by the earlier authors);
- (iii) Lemmas 1, 2, 3, 5, 6, 7, §2 and Steps 1–4 of §5 of rh_proof.md — **proved in rh_proof.md**;
- (iv) textbook facts (Hecke's functional equation with the factor Γ(s + |k|/2), the Hadamard product, Stirling,
  L(1 + it) ≠ 0, the prime number theorem mod 4, polynomial growth in vertical strips) — **classical**, with the
  uniformity in the infinity type k derived here (Lemmas S, ZC, PF), not asserted.

What is **no longer used**: Coleman's Theorem 1 (the growth bound near σ = 1, known only from a zbMATH review), the
bounds for 1/L inside the region (derived in rh_proof.md by an argument not written out), the Heath-Brown identity and
Lemmas 4, 8, 9 of rh_proof.md. The zero-density estimate (Theorem D) uses no zero-free region at all; the zero-free
region enters only in §5.

Three things remain that are not "read in the source": the textbook facts (iv) were not re-opened (the Γ-factor and
|W| = 1 are checked numerically, §10); the *generality* of Coleman's statement (every Grössencharacter mod 𝔣, with
V = max(‖m‖, |t|)) is that of Sun–Zhang's definitions, not checked in Coleman's text (Doubt 1); and the exponent
3/10 needs Söhne's bound, which is still secondary.

**Exponent.** θ < 1/(2 + 8κ) for every κ ∈ (0, 1/4) with P_κ (rh_proof.md needed κ ≥ 1/8; this route does not).
Read sources give κ = 103/512 + ε and **θ < 64/231 = 0.27706**. With Söhne's bound (secondary) κ = 1/6 + ε and
**θ < 3/10**. Michel–Venkatesh give θ < 1/(4 − 8δ). All are > 1/4, which is all R_h needs.

| step | content | status |
|---|---|---|
| §2 | classical facts (C1)–(C7) and where their uniformity in k comes from; Lemma S (Γ, ψ with explicit constants), Lemma ZC (zeros in a unit box), Lemma PF (partial fractions) | (C1)–(C7) **classical**; Lemmas S, ZC, PF **proved** from them; computed (§B, §C) |
| §3 | zero-free region for the family, uniform in k; the exceptional zero | **proved** from (i) [**quoted (secondary)**] + (C4) + compactness |
| §4 | zero density N_q(σ; M, T₀) ≤ C T₀ M^{(2+8κ)(1−σ)}(log M)^{27}, all 1/2 ≤ σ ≤ 1 (Theorem D) | **proved** from P_κ, Lemmas 5, 6, 7 of rh_proof.md, (C1)–(C4), (C6), (C7) ((C4) only through Lemma ZC); no zero-free region |
| §5 | explicit formula with the log-Gaussian weight; the bound for the zero sum is o(x) when θ(2 + 8κ) < 1 and ≥ x when θ(2 + 8κ) > 1 | **proved** from Theorem D, §3, Lemma 3 of rh_proof.md |
| §6 | removal of the smoothings | **proved in rh_proof.md** (§5, Steps 1–4; hypotheses checked: the same) |
| §7 | Theorem 1′ with explicit θ | **proved**; exponent arithmetic computed (§A) |

## 1. Statement

Notation as in rh_proof.md §1–§2: O = ℤ[i], K = ℚ(i), q ∈ O∖{0} fixed, Φ(q) = |(O/q)^*|, Λ(z) = Λ((z)),
e_k(α) = (α/|α|)^k, λ = e_4 on ideals. For a character χ of (O/q)^*, k_χ ∈ {0, 1, 2, 3} with χ(i) = i^{−k_χ};
for k ≡ k_χ (mod 4), ξ_{χ,k}((α)) := χ(α)e_k(α) is a Grössencharacter of modulus q and infinity type k, and
ξ_{χ, k_χ + 4m} = ξ_χ λ^m, m ∈ ℤ (rh_proof.md §2, proved). L(s, ξ) := Σ_{(𝔞,q)=1} ξ(𝔞)N𝔞^{−s}. If ξ ≠ 1, ξ* denotes
the primitive character inducing ξ: conductor 𝔣 | q, the same infinity type k, and
L(s, ξ) = L(s, ξ*)·Π_{𝔭|q, 𝔭∤𝔣}(1 − ξ*(𝔭)N𝔭^{−s}); the zeros of L(s, ξ) with Re s > 0 are those of L(s, ξ*).

**The family.** For M ≥ 1: 𝓕_q(M) := {ξ_χ λ^m : χ mod q, |m| ≤ M} ∖ {1}. It has ≤ Φ(q)(2M + 1) members, each
nontrivial, with |k| ≤ 4M + 3, and it contains every ξ_{χ,k} ≠ 1 with |k| ≤ 4M − 3.

**(Z) Coleman's Theorem 2, as quoted verbatim by Sun–Zhang** (arXiv 2008.09677, §2, "Theorem [coleman, Theorem 2]";
LaTeX read): "Let μλ^**m** mod 𝔣 be a Hecke character and L(σ + it, μλ^**m**) the associated Hecke L-function.
Denote V = max(‖**m**‖, t). For a constant A > 0 depending only on K and sufficiently large V, we have
L(σ + it, μλ^**m**) ≠ 0 for σ ≥ 1 − A/max(log N𝔣, (log V)^{2/3}(log log V)^{1/3}) apart from a possible exceptional
zero when **m** = **0** and μ is a real character mod 𝔣." In Sun–Zhang's §1.1 (read) λ₁, …, λ_d (d = n − 1 = 1 here)
is a basis of the infinite-order Hecke characters **mod 𝔣** (characters of I_𝔣 given on (α) ∈ P_𝔣 by the infinity
type, "we need λ(ε) = 1 for all units ε ≡ 1 mod 𝔣") and μ is a character of I_𝔣/P_𝔣; "each Hecke character mod 𝔣 can be
written as μλ^**m**". For K = ℚ(i) this covers every Grössencharacter of conductor dividing q: if u₀(𝔣) ∈ {1, 2, 4} is
the positive generator of {u ∈ ℤ : ε^u = 1 for all units ε ≡ 1 mod 𝔣}, then ξ* = μλ_𝔣^m with λ_𝔣 of type u₀(𝔣) and
|m| = |k|/u₀(𝔣) ≤ |k| (μ := ξ*λ_𝔣^{−m} is trivial on P_𝔣). Status: **quoted (secondary)**; also paraphrased, with the
modulus entering as log q in the max, in X. Li, arXiv 2111.05403 (Duke Math. J. 171 (2022)), §"Improved zero free
region type bound" (read). Sun–Zhang apply it only with 𝔣 = (1).

**Hypothesis P_κ** (as in rh_proof.md §1). For every q there is C_q such that every primitive Grössencharacter ξ of
conductor dividing q and infinity type k satisfies |L(1/2 + it, ξ)| ≤ C_q(2 + |k| + |t|)^{2κ} for all t ∈ ℝ.

**Theorem 1′.** Assume P_κ with 0 < κ < 1/4. Let θ < 1/(2 + 8κ) and ζ ∈ (O/q)^*. Then, with (Z) as the only quoted
input about zeros (besides the classical non-vanishing (C4) on Re s = 1), as x → ∞,

    S(x; φ, ω) := Σ_{z ≡ ζ (mod q), N z ≤ x, arg z ∈ (φ, φ+ω]} Λ(z) = (2/π)·ω·x/Φ(q)·(1 + o(1)),

uniformly in φ ∈ ℝ and ω ∈ [x^{−θ}, 2π].

**Corollary.** P_κ holds with κ = 103/512 + ε (§7), so Theorem 1′ holds for every θ < 64/231; given Söhne's bound
(secondary), for every θ < 3/10. Since 64/231 > 1/4, this is Hypothesis R_h of arcs_theorem.md §2, by the same
reduction as rh_proof.md §1 (Corollary): for fixed φ, c > 0 the width c·x^{−1/4} is ≥ x^{−θ} for x ≥ x₀(c), and
(φ − ω, φ] = (φ', φ' + ω] with φ' = φ − ω is covered by the uniformity in φ.

## 2. Classical facts, and where their uniformity in the infinity type comes from

Throughout, ξ* is a nontrivial primitive Grössencharacter of K of conductor 𝔣 and infinity type k, and
A_𝔣 := √(|d_K| N𝔣)/(2π) = √(N𝔣)/π.

- **(C1) Functional equation** (Hecke 1920; Neukirch, *Algebraic Number Theory*, Ch. VII §8; Iwaniec–Kowalski,
  *Analytic Number Theory*, §3.8 — classical, not re-opened). Λ(s, ξ*) := A_𝔣^s Γ(s + |k|/2) L(s, ξ*) is entire and
  Λ(s, ξ*) = W(ξ*) Λ(1 − s, ξ̄*), |W(ξ*)| = 1. Also Λ(s, ξ̄*) = conj Λ(s̄, ξ*). *Uniformity:* an identity for each ξ*;
  k enters only through Γ(s + |k|/2) and A_𝔣 does not depend on k, so no constant has to be uniform. Computed (§C):
  |W| = 1 to 10⁻¹², W independent of s and of the theta split to 10⁻¹⁴, for k = 4, 8, −12 (q = 1) and k = −1, 1, 2, 3, 5, 6, 9
  (conductor 2 + i); negative controls (Γ(s + |k|), a pair (χ, k) with k ≢ k_χ, an imprimitive character) fail.
- **(C2) Hadamard product** (Λ entire of order 1; Davenport, *Multiplicative Number Theory*, Ch. 12; IK Thm 5.6 —
  classical). Λ(s, ξ*) = e^{a + bs} Π_ρ (1 − s/ρ)e^{s/ρ} over the zeros of Λ (the nontrivial zeros of L, 0 ≤ β ≤ 1).
  *Uniformity:* only identities derived from it are used (Re b = −Σ Re 1/ρ, and differences in which b cancels); no
  constant from (C2) enters a bound.
- **(C3) Γ and ψ: the shift, Binet/Stirling, Cahen–Mellin** (Whittaker–Watson Ch. XII, Binet §12.31–12.33 —
  classical). In the form used:
  (i) Γ is meromorphic on ℂ without zeros, its poles are exactly 0, −1, −2, …, all simple, Γ(1) = 1 (so the residue
  at 0 is 1), and Γ(z + 1) = zΓ(z), ψ(z + 1) = ψ(z) + 1/z for z ∉ {0, −1, −2, …};
  (ii) for Re z > 0 (principal branch of log z, the branch of log Γ that is real on (0, ∞)):
  log Γ(z) = (z − 1/2)log z − z + (1/2)log 2π + μ(z), μ(z) = ∫₀^∞ g(t)e^{−zt}dt/t, and
  ψ(z) = Γ'/Γ(z) = log z − 1/(2z) − ∫₀^∞ g(t)e^{−zt}dt, with g(t) = 1/(e^t − 1) − 1/t + 1/2, 0 ≤ g(t) ≤ t/12;
  (iii) e^{−y} = (1/2πi)∫_{(c)} Γ(w)y^{−w}dw for y > 0 and c > 0 (Mellin inversion of Euler's integral for Γ).
  Formula (ii) is used only on Re z > 0; on Re z ≤ 0, Γ and ψ enter only through the shifts in (i) (Lemma S(c), the
  bound (Γ) in Lemma DI, §4.5, Lemma EF).
  *Uniformity:* statements about the single functions Γ, ψ on half-planes; k enters as the shift z = s + |k|/2.
- **(C4) Non-vanishing on Re s = 1** (Hecke; e.g. Lang, *Algebraic Number Theory*, Ch. XV; Narkiewicz, *Elementary and
  Analytic Theory of Algebraic Numbers*, Ch. 7 — classical). L(1 + it, ξ*) ≠ 0 (3–4–1 inequality
  |ζ_K³ L(σ+it, ξ*)⁴ L(σ+2it, ξ*²)| ≥ 1 for σ > 1, and for real ξ* of type 0 the factorisation ζ_K L(s, ξ*) = ζ_{K'}).
  Used only qualitatively: Λ(0, ξ*) = W Λ(1, ξ̄*) ≠ 0, and for finitely many characters in a compact region (§3).
- **(C5) Prime number theorem mod 4** (Davenport Ch. 18 and Ch. 20 — classical): ψ(y) = y + O(y e^{−c√log y}),
  ψ(y, χ_{−4}) = O(y e^{−c√log y}). With Σ_{N𝔞=n}Λ(𝔞) = Λ(n)(1 + χ_{−4}(n)) (computed for n ≤ 3000, §E) this is the
  prime ideal theorem for K, input (d) of rh_proof.md. Two fixed L-functions: no uniformity question.
- **(C6) Polynomial growth in vertical strips** (Phragmén–Lindelöf / convexity; Rademacher, Math. Z. 72 (1959);
  IK §5.1 — classical): for each fixed ξ*, L(σ + it, ξ*) ≪_{ξ*} (1 + |t|)^{C} on 1/2 ≤ σ ≤ 3. Used only to justify
  contour shifts (Lemma DI), qualitatively; no constant from it enters a bound.
- **(C7)** For σ > 1: |L'/L(σ + it, ξ)| ≤ −ζ_K'/ζ_K(σ) (coefficients Λ(𝔞)ξ(𝔞), |ξ| ≤ 1). Uniform trivially.
  Computed: −ζ_K'/ζ_K(2) = 0.4809, −ζ_K'/ζ_K(5/4) = 3.2769.
- Jensen's formula is not needed.

**Lemma S (proved from (C3)).** (a) |ψ(z) − log z| ≤ 4/3 for Re z ≥ 1/2. (b) |Γ(a + iv)| ≤ 3e^{−|v|} for 1/2 ≤ a ≤ 1,
v ∈ ℝ. (c) For Re z ≥ −1/4 and |z| ≥ 1/4: |ψ(z)| ≤ log(|z| + 3) + 8.
Proof. (a) |ψ(z) − log z| ≤ 1/(2|z|) + ∫₀^∞ (t/12)e^{−t Re z}dt = 1/(2|z|) + 1/(12(Re z)²) ≤ 1 + 1/3.
(b) |μ(z)| ≤ ∫₀^∞ e^{−t Re z}dt/12 = 1/(12 Re z) ≤ 1/6; v·arg(a + iv) = |v| arctan(|v|/a) ≥ π|v|/2 − a; so
log|Γ(a+iv)| ≤ (a − 1/2)log|a + iv| − π|v|/2 + (1/2)log 2π + 1/6 ≤ (1/2)log(1 + |v|) − π|v|/2 + 1.086, and
(1 + |v|)^{1/2}e^{−(π/2 − 1)|v|} ≤ 1 (its log has derivative 1/(2(1+v)) − 0.5708 < 0); e^{1.086} < 3.
(c) If Re z ≥ 1/2: (a) and |log z| ≤ |log|z|| + π/2 ≤ log(|z| + 3) + π/2 (|z| ≥ 1/2), total ≤ log(|z| + 3) + 2.91.
Otherwise ψ(z) = ψ(z + 1) − 1/z ((C3)(i)) with Re(z + 1) ≥ 3/4 and |1/z| ≤ 4, total ≤ log(|z| + 3) + log(4/3) + 2.91 + 4. ∎
(Computed, §B: sup |ψ − log| = 1.27 on a grid; sup (|ψ(z)| − log(|z| + 3)) = 3.05 ≤ 8 on a grid; sup |Γ(a+iv)|e^{|v|}
= 1.97; 0 ≤ g ≤ t/12 at 60 digits; |μ| ≤ 1/(12 Re z).)

**Lemma ZC (zeros in a unit box; proved from (C1)–(C4), (C7), Lemma S).** For every T ∈ ℝ,

    #{ρ : |γ − T| ≤ 1} ≤ 1.25 log N𝔣 + 2.5 log(|k| + |T| + 2) + 3.35     (with multiplicity).

Proof. Λ(0) ≠ 0 by (C1), (C4). Differentiating the functional equation, Λ'/Λ(s, ξ*) = −Λ'/Λ(1 − s, ξ̄*); with (C2)
for ξ* and for ξ̄* (whose zeros are the ρ̄, and b(ξ̄*) = conj b(ξ*)) at s = 0:
b = −b̄ − Σ_ρ (1/(1 − ρ̄) + 1/ρ̄), the sum absolutely convergent. The map ρ ↦ 1 − ρ̄ permutes the zeros (functional
equation and conjugation), and Σ Re 1/ρ = Σ β/|ρ|² converges absolutely; so 2 Re b = −2 Σ_ρ Re 1/ρ and
Re Λ'/Λ(s) = Σ_ρ Re 1/(s − ρ). At s₁ = 2 + iT: Λ'/Λ(s₁) = log A_𝔣 + ψ(s₁ + |k|/2) + L'/L(s₁), so by (C7) and Lemma S(a)
Σ_ρ Re 1/(s₁ − ρ) ≤ (1/2)log N𝔣 − log π + log(|k| + |T| + 2) + 4/3 + 0.481. Every term is ≥ 0, and for |γ − T| ≤ 1,
Re 1/(s₁ − ρ) = (2 − β)/((2 − β)² + (T − γ)²) ≥ 2/5 (minimum of a/(a² + 1) on 1 ≤ a ≤ 2). Multiply by 5/2;
(5/2)(0.481 − log π + 4/3) = 1.67 ≤ 3.35. ∎ *Uniformity in k:* only through log(|k| + |T| + 2), from Stirling at
Re z ≥ 2. (Computed: the constants, §B; on the critical line of L(s, λ) (q = 1, k = 4) and of L(s, ξ_{χ₁,1})
(conductor 2 + i) the sign changes in [0, 14] are 3 and 6, at most 0.18 of the bound in every window, §C.)

**Lemma PF (partial fractions; proved).** There is an absolute constant C_PF such that for s = σ + it with
−1/4 ≤ σ ≤ 2 and (σ ≥ 1/2 or |t| ≥ 1):

    | L'/L(s, ξ*) − Σ_{ρ: |γ − t| < 1} 1/(s − ρ) | ≤ C_PF log(N𝔣(|k| + |t| + 3)).

Proof. From (C2): L'/L(s) = b − log A_𝔣 − ψ(s + |k|/2) + Σ_ρ(1/(s − ρ) + 1/ρ). Subtract the same at s₁ = 2 + it:
L'/L(s) − L'/L(s₁) = ψ(s₁ + |k|/2) − ψ(s + |k|/2) + Σ_ρ (1/(s − ρ) − 1/(s₁ − ρ)). Here |L'/L(s₁)| ≤ 0.481 (C7);
both ψ-terms are ≤ log(|k| + |t| + 6) + 8 by Lemma S(c) (the hypothesis σ ≥ 1/2 or |t| ≥ 1 gives |z| ≥ 1/2 for
z = s + |k|/2, so z stays away from the poles of ψ); for |γ − t| ≥ 1, |1/(s − ρ) − 1/(s₁ − ρ)| = |2 − σ|/(|s − ρ||s₁ − ρ|) ≤ (9/4)/(γ − t)², and the
zeros with n ≤ |γ − t| < n + 1 are ≤ 2·ZC(|t| + n + 1), so they contribute
≤ (9/2)Σ_n ZC(|t| + n + 1)/n² ≪ log(N𝔣(|k| + |t| + 3)); for |γ − t| < 1, |1/(s₁ − ρ)| ≤ 1 and there are ≤ ZC(t)
of them. ∎ *Uniformity in k:* Lemma ZC and Lemma S only.

Consequences used in §5 (proved). Write ZC(T) := 1.25 log N𝔣 + 2.5 log(|k| + |T| + 2) + 3.35 for the bound of
Lemma ZC; it depends on T only through |T| and increases with |T|.
(a) For each n ≥ 1 and each sign ε ∈ {+1, −1} there is T ∈ [n, n + 1] such that t := εT satisfies |t − γ| ≥ 1/(2Z_n + 2)
for every zero, Z_n := ZC(n + 1/2); write T_n := T for ε = +1 and T_n' := T for ε = −1. For these t and −1/4 ≤ σ ≤ 2,

    |L'/L(σ + it, ξ*)| ≤ Z_n(2Z_n + 2) + 2(ZC(n − 1/2) + ZC(n + 3/2)) + C_PF log(N𝔣(|k| + n + 4))
                      ≤ (60 + C_PF)·log²(N𝔣(|k| + n + 4)).

Proof. Take ε = +1; for ε = −1 replace every ordinate window [a, b] below by [−b, −a] (ZC depends on |T| only, so the
counts are the same). The zeros with γ ∈ [n − 1/2, n + 3/2] satisfy |γ − (n + 1/2)| ≤ 1, so there are at most Z_n of
them (with multiplicity). Put r := 1/(2Z_n + 2). The open intervals (γ − r, γ + r) around their ordinates have total
length ≤ 2rZ_n = Z_n/(Z_n + 1) < 1, so some T ∈ [n, n + 1] lies in none of them: |T − γ| ≥ r for these zeros; every
other zero has γ < n − 1/2 or γ > n + 3/2, hence |T − γ| > 1/2 ≥ r. Apply Lemma PF at s = σ + iT (|T| ≥ 1). The zeros
with |γ − T| < 1 have γ ∈ (n − 1, n + 2) and are of two kinds. (α) γ ∈ [n − 1/2, n + 3/2]: at most Z_n of them, each with
|1/(s − ρ)| ≤ 1/|T − γ| ≤ 2Z_n + 2. (β) γ ∈ (n − 1, n − 1/2) ∪ (n + 3/2, n + 2): such zeros can occur, since the window
(T − 1, T + 1) leaves [n − 1/2, n + 3/2] unless T = n + 1/2. The first ones satisfy |γ − (n − 1/2)| < 1/2, the second
ones |γ − (n + 3/2)| < 1/2, so there are at most ZC(n − 1/2) + ZC(n + 3/2) of them; each has |T − γ| > 1/2 (T ∈ [n, n + 1]),
so |1/(s − ρ)| < 2. Lemma PF adds C_PF log(N𝔣(|k| + T + 3)) ≤ C_PF log(N𝔣(|k| + n + 4)). This is the first line. For the
second, put ℓ := log(N𝔣(|k| + n + 4)) ≥ log 5 > 1.6. Every ZC(·) above is ≤ ZC(n + 3/2) ≤ 2.5ℓ + 3.35 ≤ 4.6ℓ (use
1.25 log N𝔣 ≤ 2.5 log N𝔣 and 3.35 ≤ 2.1ℓ), so the first line is ≤ 2Z² + 6Z + C_PF ℓ with Z := 4.6ℓ, that is
≤ (42.32 + 27.6/1.6 + C_PF/1.6)ℓ² ≤ (60 + C_PF)ℓ². ∎
(b) |L'/L(1, ξ*)| ≤ Σ_{|γ|<1} 1/|1 − ρ| + C_PF log(N𝔣(|k| + 3)).

## 3. The zero-free region for the family (proved from (Z), (C4) and compactness)

**Lemma ZFR.** For fixed q there are c₃ = c₃(q) > 0 and V₃(q) such that for V* ≥ V₃: every zero ρ = β + iγ of
L(s, ξ*), ξ* primitive of conductor 𝔣 | q with |k| ≤ V*, and |γ| ≤ V*, satisfies

    β ≤ 1 − c₃ (log V*)^{−2/3} (log log V*)^{−1/3}.

Proof. Let A₀ and V₀ be the constant and the threshold of (Z) (if V₀ depends on 𝔣, take the largest over the
finitely many 𝔣 | q). With ξ* = μλ_𝔣^m (§1), Coleman's V is
V = max(|m|, |γ|) ≤ V*. If the quoted "t" means t ≥ 0 only, zeros with γ < 0 are handled by conjugation:
L(s̄, ξ̄*) = conj L(s, ξ*) and ξ̄* = μ̄λ_𝔣^{−m} has the same |m|.
(a) V ≥ V₁ := max(V₀, V₂), V₂ = V₂(q) such that (log V)^{2/3}(log log V)^{1/3} ≥ log N q for V ≥ V₂. Then the max in (Z)
is (log V)^{2/3}(log log V)^{1/3} and β < 1 − A₀(log V)^{−2/3}(log log V)^{−1/3} ≤ 1 − A₀(log V*)^{−2/3}(log log V*)^{−1/3}
(the function decreases for V ≥ e^e). No exceptional zero occurs here: it is real (γ = 0) with m = 0, so V = 0 < V₁.
(b) V < V₁. Then |k| ≤ 4V₁ and |γ| < V₁: finitely many characters (𝔣 | q, |k| ≤ 4V₁), and the zeros of each in the
compact set {1/2 ≤ β ≤ 1, |γ| ≤ V₁} are finitely many (zeros of a nonzero entire function); each has β < 1 by (C4)
(and the Euler product for β > 1). Let 1 − c₂(q) be the largest of these β; c₂(q) > 0. Any exceptional zero is among
them. Take c₃ := A₀ and V₃ so large that A₀(log V₃)^{−2/3}(log log V₃)^{−1/3} ≤ c₂. ∎

**The exceptional zero.** For fixed q it costs nothing in an asymptotic without rate: it belongs to one of the
finitely many real characters of type 0 mod q, it is a fixed real number β₁ < 1, and it is absorbed in c₂(q). Exactly
what it would cost: its term in S(x; φ, ω) is −(4/Φ(q))χ̄(ζ)ψ̂_±(0)Ṽ(β₁) = −(2/π)(ω x/Φ(q))·χ̄(ζ)x^{β₁−1}/β₁·(1 + o(1)),
a relative error ≍ x^{−(1−β₁)}: o(1) for fixed q; an explicit rate would need a numerical bound for β₁ (finitely many
fixed L-functions); a version uniform in q would meet the Siegel-zero problem. No rate is claimed here.

**Where it is used.** With V* := 4M + T₀ + 6 ≤ x and η := c₃(log x)^{−2/3}(log log x)^{−1/3}, every zero of every
L(s, ξ), ξ ∈ 𝓕_q(M), and of L(s, ξ̄*), with |γ| ≤ T₀, has β ≤ 1 − η (x ≥ x₀(q)). This is the only appearance of (Z).

## 4. The zero-density estimate

**Theorem D (proved).** Assume P_κ, 0 < κ < 1/4, A := 2 + 8κ. Let N_q(σ; M, T₀) be the number of zeros
ρ = β + iγ (with multiplicity) of the L(s, ξ), ξ ∈ 𝓕_q(M), with β ≥ σ and |γ| ≤ T₀. There are C = C(q, κ) and
M₀ = M₀(q, κ) such that for M ≥ M₀, 2 ≤ T₀ ≤ M and 1/2 ≤ σ ≤ 1:

    N_q(σ; M, T₀) ≤ C · (T₀ + 4) · M^{A(1−σ)} · (log M)^{27}.

The power 27 of log is crude and not optimised; T₀ enters linearly (it will be (log x)³). The principal character is
not in 𝓕_q(M) (its zeros, those of ζ_K, are never needed; its term is the main term, §6). No zero-free region is used.

**Where it is needed, and where it must be sharp.** In §5 only σ ∈ [σ₁, 1 − η] matter, σ₁ := 1 − 1/A
(= 167/231 = 0.7229 for κ = 103/512, 7/10 for κ = 1/6): below σ₁ the trivial bound N ≪ T₀ M log M is at least as good
as M^{A(1−σ)} (because A(1 − σ) ≥ 1), and above 1 − η there are no zeros (§3). On [σ₁, 1 − η] the exponent A(1 − σ)
is needed with the same A everywhere: the condition of §5 is sup_σ (exponent/(1 − σ))·θ < 1, and the bound gives the
ratio exactly A on all of [σ₁, 1). (Computed, §A: with X = M, Y = M^{1+4κ} the sup over σ of min(trivial,
detection)/(1 − σ) equals 2 + 8κ to 10⁻⁹ for κ = 0.02 … 0.25, and no (X, Y) = (M^a, M^b) on a grid
a ∈ [0.5, 1.5], b ∈ [1, 2.5] does better.)

**4.1 Trivial range (proved).** Each ξ ∈ 𝓕_q(M) has ≤ (T₀ + 2)Z₁ zeros with |γ| ≤ T₀, where
Z₁ := 1.25 log Nq + 2.5 log(4M + T₀ + 6) + 3.35 bounds every window count of Lemma ZC that occurs (|k| ≤ 4M + 3,
|T| ≤ T₀ + 1). So N_q(σ) ≤ Φ(q)(2M + 1)(T₀ + 2)Z₁ ≤ C T₀ M log M ≤ C T₀ M^{A(1−σ)} log M for σ ≤ σ₁. It remains to
treat σ ∈ [σ₁, 1], where σ − 1/2 ≥ σ₁ − 1/2 = 2κ/(1 + 4κ) > 0 (computed exactly, §A).

**4.2 Mollifier and the detection identity.** For X ≥ 1 let M_X(s, ξ) := Σ_{N𝔞≤X} μ(𝔞)ξ(𝔞)N𝔞^{−s} and
c_X(𝔞) := Σ_{𝔡|𝔞, N𝔡≤X} μ(𝔡). Then L(s, ξ)M_X(s, ξ) = Σ_𝔞 c_X(𝔞)ξ(𝔞)N𝔞^{−s} for Re s > 1 (complete multiplicativity of
ξ, with ξ = 0 off ideals prime to q), c_X(O) = 1, c_X(𝔞) = 0 for 1 < N𝔞 ≤ X (Σ_{𝔡|𝔞} μ(𝔡) = 0), and |c_X(𝔞)| ≤ d(𝔞).

**Lemma DI (proved).** Let ξ ∈ 𝓕_q(M), ρ = β + iγ a zero of L(s, ξ) with β > 1/2, X ≥ 1, Y > 0. Then

    e^{−1/Y} + Σ_{N𝔞>X} c_X(𝔞)ξ(𝔞)N𝔞^{−ρ}e^{−N𝔞/Y}
        = (1/2π) ∫_ℝ L(1/2 + i(γ+v), ξ) M_X(1/2 + i(γ+v), ξ) Γ(1/2 − β + iv) Y^{1/2−β+iv} dv =: I(ρ).

Proof. e^{−y} = (1/2πi)∫_{(2)} Γ(w)y^{−w}dw ((C3)(iii)); on Re w = 2 the series Σ|c_X(𝔞)|N𝔞^{−β−2} converges, so the
left side is (1/2πi)∫_{(2)} L(ρ + w, ξ)M_X(ρ + w, ξ)Γ(w)Y^w dw. Move to Re w = 1/2 − β ∈ [−1/2, 0): L(·, ξ) is entire
(ξ ≠ 1, (C1)), M_X is entire, Γ has in −1/2 ≤ Re w ≤ 2 only the pole w = 0 ((C3)(i)), whose residue L(ρ)M_X(ρ) = 0.
The horizontal segments vanish, and the integral on Re w = 1/2 − β converges absolutely, by (C6) and the bound

    (Γ)  |Γ(a + iv)| ≤ 10|v|^{3/2}e^{−π|v|/2}   for −1/2 ≤ a ≤ 2, |v| ≥ 1.

Proof of (Γ). Let w = a + iv. (1/2 ≤ a ≤ 2.) By (C3)(ii), log|Γ(w)| = Re log Γ(w) = (a − 1/2)log|w| − v·arg w − a
+ (1/2)log 2π + Re μ(w), with |μ(w)| ≤ ∫₀^∞ e^{−at}dt/12 = 1/(12a) ≤ 1/6 (g(t)/t ≤ 1/12) and
v·arg w = |v|arctan(|v|/a) = |v|(π/2 − arctan(a/|v|)) ≥ π|v|/2 − a. Hence

    |Γ(a + iv)| ≤ √(2π) e^{1/(12a)} |w|^{a−1/2} e^{−π|v|/2}   (1/2 ≤ a ≤ 2, v ∈ ℝ),

and, for |v| ≥ 1, |w|^{a−1/2} ≤ |w|^{3/2} ≤ (5v²)^{3/4} (|w| ≥ |v| ≥ 1, |w|² ≤ 4 + v² ≤ 5v²), so
|Γ(w)| ≤ √(2π)e^{1/6}5^{3/4}|v|^{3/2}e^{−π|v|/2} ≤ 9.91|v|^{3/2}e^{−π|v|/2}. (−1/2 ≤ a < 1/2.) Binet's formula is not
available here; by (C3)(i), Γ(w) = Γ(w + 1)/w with Re(w + 1) = a + 1 ∈ [1/2, 3/2). The display above for w + 1, with
0 ≤ a + 1/2 < 1 and 1 ≤ |w + 1| ≤ (√13/2)|v| (|w + 1|² ≤ 9/4 + v² ≤ (13/4)v²), gives |Γ(w + 1)| ≤ √(2π)e^{1/6}(√13/2)|v|
e^{−π|v|/2}; dividing by |w| ≥ |v| gives |Γ(w)| ≤ √(2π)e^{1/6}(√13/2)e^{−π|v|/2} ≤ 5.34e^{−π|v|/2}. Both cases give (Γ).
Use of (Γ). On a horizontal segment Im w = ±V (V ≥ 1, 1/2 − β ≤ Re w ≤ 2) one has 1/2 ≤ Re(ρ + w) ≤ 3, so
|L(ρ + w, ξ)| ≪_ξ (2 + |γ| + V)^C by (C6) for ξ* and the finitely many Euler factors |1 − ξ*(𝔭)N𝔭^{−ρ−w}| ≤ 2;
|M_X(ρ + w, ξ)| ≤ Σ_{N𝔞≤X}N𝔞^{−1/2} and |Y^w| ≤ Y². With (Γ) the integrand is ≪ (2 + |γ| + V)^C V^{3/2}e^{−πV/2} → 0;
on the line Re w = 1/2 − β the same bounds hold for |v| ≥ 1, and Γ is continuous on |v| ≤ 1 there (Re w ≠ 0). ∎
(Computed, §D: the identity in the form "sum = L M + line Re w = −1/4" holds to 10⁻⁸ for L(s, χ_{−4}) at its first zero
1/2 + 6.0209489046976i and at s₀ = 0.8 + 3i, X = 10, Y = 30.)

**4.3 Classification (proved).** Let Y ≥ 100, Y₁ := Y(log Y)², X ≤ Y. Tail: the facts (T1)–(T3) below give

    𝒯(Y) := Σ_{N𝔞>Y₁} d(𝔞)N𝔞^{−1/2}e^{−N𝔞/Y} ≤ 4.7·10⁻⁵ ≤ 10⁻⁴   for every real Y ≥ 100,

and e^{−1/Y} ≥ 1 − 1/Y ≥ 0.99. Here d(𝔞) is the number of ideals dividing 𝔞, r₀(n) := #{𝔞 : N𝔞 = n}, and
d_j := 1 ∗ ⋯ ∗ 1 (j factors, Dirichlet convolution), so d₁ = 1, d₂ = d, d₄ = d ∗ d.
- (T1) Σ_{n≤u} d₄(n) ≤ u(1 + log u)³ for real u ≥ 1. Proof: by induction on j, Σ_{n≤u} d_j(n) ≤ u(1 + log u)^{j−1} for
  u ≥ 1. For j = 1 the sum is ⌊u⌋ ≤ u. If it holds for j, then, since d_{j+1} = d_j ∗ 1 and u/m ≥ 1 for m ≤ u,
  Σ_{n≤u} d_{j+1}(n) = Σ_{m≤u} Σ_{l≤u/m} d_j(l) ≤ Σ_{m≤u} (u/m)(1 + log(u/m))^{j−1} ≤ u(1 + log u)^{j−1} Σ_{m≤u} 1/m
  ≤ u(1 + log u)^j, because Σ_{m≤u} 1/m ≤ 1 + ∫₁^u dt/t = 1 + log u (1/m ≤ ∫_{m−1}^m dt/t for m ≥ 2). Take j = 4.
- (T2) Σ_{N𝔞=n} d(𝔞) = (r₀ ∗ r₀)(n) ≤ d₄(n), hence Σ_{N𝔞≤u} d(𝔞) ≤ Σ_{n≤u} d₄(n) ≤ u(1 + log u)³ for u ≥ 1.
  Proof: by unique factorisation of ideals, d(𝔞) is the number of pairs (𝔟, 𝔠) with 𝔟𝔠 = 𝔞, and N(𝔟𝔠) = N𝔟·N𝔠, so
  Σ_{N𝔞=n} d(𝔞) = #{(𝔟, 𝔠) : N𝔟·N𝔠 = n} = Σ_{ab=n} r₀(a)r₀(b). Next r₀ ≤ d: by the splitting law of ℤ[i] (classical:
  (2) = (1 + i)², a prime p ≡ 1 (mod 4) is 𝔭𝔭̄ with 𝔭 ≠ 𝔭̄ of norm p, a prime p ≡ 3 (mod 4) is inert, of norm p²) and
  unique factorisation, r₀ is multiplicative with r₀(2^e) = 1, r₀(p^e) = e + 1 for p ≡ 1 (mod 4) (the ideals 𝔭^a𝔭̄^{e−a},
  0 ≤ a ≤ e), and r₀(p^e) = 1 or 0 for p ≡ 3 (mod 4) as e is even or odd. So r₀(p^e) ≤ e + 1 = d(p^e), r₀(n) ≤ d(n), and,
  as r₀ ≥ 0, (r₀ ∗ r₀)(n) ≤ (d ∗ d)(n) = d₄(n). The last claim is (T1) summed.
- (T3) Dyadic blocks. Put L := log Y ≥ log 100 > 4.6, so Y₁/Y = L² and Y₁ = L²e^L. Split N𝔞 > Y₁ into the blocks
  (2^jY₁, 2^{j+1}Y₁], j ≥ 0. On block j, N𝔞^{−1/2} < (2^jY₁)^{−1/2}, e^{−N𝔞/Y} < e^{−2^jL²}, and Σ d(𝔞) over the block is
  ≤ 2^{j+1}Y₁(1 + log(2^{j+1}Y₁))³ by (T2). Hence

      𝒯(Y) ≤ Σ_{j≥0} t_j(L),   t_j(L) := 2^{1+j/2} L e^{L/2} (1 + (j + 1)log 2 + L + 2 log L)³ e^{−2^j L²}.

  Each t_j decreases on L ≥ 4: d/dL log t_j(L) = 1/2 + 1/L + 3(1 + 2/L)/(1 + (j + 1)log 2 + L + 2 log L) − 2^{j+1}L
  ≤ 1/2 + 1/4 + 9/10 − 8 < 0 (the denominator is ≥ 1 + L ≥ 5). The ratios are tiny: for L ≥ log 100,
  t_{j+1}/t_j = √2·(1 + log 2/(1 + (j + 1)log 2 + L + 2 log L))³·e^{−2^j L²} ≤ √2(1 + log 2/(1 + L))³e^{−L²}
  ≤ √2(1.124)³e^{−21.2} ≤ 1.3·10⁻⁹. So 𝒯(Y) ≤ t_0(log Y)(1 + 2·10⁻⁹) ≤ t_0(log 100)(1 + 2·10⁻⁹) for every Y ≥ 100, and
  t_0(log 100) = 2·10·L(1 + log 2 + L + 2 log L)³e^{−L²} with L = log 100 = 4.6052: 1 + log 2 + L + 2 log L = 9.3527,
  L² = 21.2076, t_0 = 92.10·818.1·6.161·10⁻¹⁰ = 4.642·10⁻⁵. Therefore 𝒯(Y) ≤ 4.7·10⁻⁵.

Levels U_j := 2^j X, 0 ≤ j ≤ J, J minimal with 2^{J+1}X ≥ Y₁ (J + 1 = ⌈log₂(Y₁/X)⌉ < log₂Y₁ + 1 ≤ 2 log Y₁, as X ≥ 1 and
log Y₁ > 7), pieces
P_j(ρ) := Σ_{U_j<N𝔞≤min(2U_j, Y₁)} c_X(𝔞)ξ(𝔞)N𝔞^{−ρ}e^{−N𝔞/Y}. By Lemma DI, Σ_j|P_j(ρ)| + |I(ρ)| ≥ 0.98. Call ρ
**type II** if |I(ρ)| ≥ 1/2; otherwise Σ_j|P_j(ρ)| ≥ 0.48 and ρ is **type I at level U_j** for the least j with
|P_j(ρ)| ≥ δ_I := 1/(5 log Y₁) (0.48/(J + 1) ≥ δ_I).

**4.4 Type I zeros (proved; Lemmas 5, 6, 7 of rh_proof.md).** Fix χ and a level U = U_j. Removing the dependence on β
(partial summation): put S_u(t; ξ) := Σ_{U<N𝔞≤min(u,Y₁)} c_X(𝔞)ξ(𝔞)e^{−N𝔞/Y}N𝔞^{−σ−it} for U ≤ u ≤ 2U. Since
P_j(ρ) = Σ a_n g(n) with g(u) = u^{σ−β}, 0 < g ≤ 1, |g'(u)| ≤ (β − σ)/u ≤ 1/(2u), Abel summation gives
|P_j(ρ)| ≤ |S_{2U}(γ)| + (1/2)∫_U^{2U}|S_u(γ)|du/u, hence, by Cauchy–Schwarz,

    |P_j(ρ)|² ≤ 2|S_{2U}(γ)|² + (log 2/2)∫_U^{2U}|S_u(γ)|² du/u.

Summing over the type-I zeros at level U of all ξ = ξ_χλ^m, |m| ≤ M (the u-integral commutes with the sum):
#𝒵_I(U)·δ_I² ≤ (2 + (log 2)²/2)·sup_{U≤u≤2U} Σ_ξ Σ_ρ |S_u(γ_ρ; ξ)|². For fixed u, S_u is a Dirichlet polynomial in t
with coefficients independent of ρ, so Lemma 7 of rh_proof.md (|F(t₀)|² ≤ ∫_{|t−t₀|≤1/2}(|F|² + 2|F||F'|), F C¹: a
finite Dirichlet polynomial is) and the window count give
Σ_ρ |S_u(γ_ρ; ξ)|² ≤ Z₁ ∫_{|t|≤T₀+1/2} (|S_u|² + 2|S_u||S_u'|) dt, and with Cauchy–Schwarz over (ξ, t),
Σ_ξ ∫(…) ≤ S₀ + 2√(S₀S₁), S₀ := Σ_ξ ∫_{|t|≤T₀+1/2}|S_u|², S₁ the same for S_u' (coefficients times −i log N𝔞).
**Lemma 5 of rh_proof.md** applies with: coefficients c(𝔞) := c_X(𝔞)ξ_χ(𝔞)e^{−N𝔞/Y}N𝔞^{−σ}·1_{U<N𝔞≤min(u,Y₁)}
(arbitrary complex coefficients are allowed; the fixed character ξ_χ sits in the coefficients, |ξ_χ| ≤ 1), frequencies
λ^m(𝔞)N𝔞^{−it} with |m| ≤ M (both signs, as in the lemma), support bound N := 2U, height T := T₀ + 1/2. Since
D(𝔞) = cont(𝔞)√(2U/N𝔞) ≤ √2 cont(𝔞) on the support, |c(𝔞)|² ≤ d(𝔞)²U^{1−2σ}/N𝔞 (σ ≥ 1/2, N𝔞 > U) and
T + √(2π)cont(𝔞) ≤ (T₀ + 4)cont(𝔞):

    S₀ ≤ e²√π(√πM + πU)(T₀ + 4) U^{1−2σ} L₆(2U),   L₆(Y) := Σ_{N𝔟≤Y} d(𝔟)² cont(𝔟)/N𝔟 ≤ C (log Y)^{20}

by **Lemma 6 of rh_proof.md** with g = 2 (g⁴ + g² = 20); S₁ ≤ (log 2U)²·(the same), so S₀ + 2√(S₀S₁) ≤ (1 + 2 log 2U)·(bound
for S₀). The derivative and the cut 1_{N𝔞≤min(u,Y₁)} do not leave the lemma's hypotheses (coefficients on N𝔞 ≤ N).
With X := M we have U ≥ M, so √πM + πU ≤ 5U, and summing over the χ (factor Φ(q)) and the levels (U ≤ Y₁,
J + 1 ≤ 2 log Y₁):

    #(type I) ≤ C(q) · 25(log Y₁)² · Z₁ · (1 + 2 log 2Y₁) · (T₀ + 4) · L₆(2Y₁) · (J + 1) · Y₁^{2−2σ},

and Y₁^{2−2σ} ≤ Y^{2−2σ}(log Y)² for σ ≥ 1/2. (Computed, §D: Lemma 5 holds for these vectors — c_X(𝔞)ξ_{χ₁}(𝔞)e^{−N𝔞/Y}
N𝔞^{−σ} on (U, u], conductor 2 + i, and their log-derivatives — with LHS/bound ≤ 0.056; partial summation and Gallagher
with multiplicity on 300 resp. 20 random cases.)

**4.5 Type II zeros (proved; P_κ, Lemmas 5, 6 of rh_proof.md).** For σ ≤ β < 1, the shift (C3)(i) and Lemma S(b) give
|Γ(1/2 − β + iv)| = |Γ(3/2 − β + iv)|/|1/2 − β + iv| ≤ 3e^{−|v|}/(σ − 1/2), and Y^{1/2−β} ≤ Y^{1/2−σ}. By P_κ for ξ*
and |1 − ξ*(𝔭)N𝔭^{−1/2−iu}| ≤ 2: |L(1/2 + iu, ξ)| ≤ 2^{ω(q)}C_q(2 + |k| + |u|)^{2κ}
≤ 2^{ω(q)}C_q(4M + T₀ + 5)^{2κ}(1 + |v|)^{2κ} for u = γ + v, |γ| ≤ T₀, |k| ≤ 4M + 3. With Cauchy–Schwarz and
∫(1 + |v|)^{4κ}e^{−|v|}dv ≤ ∫(1 + |v|)e^{−|v|}dv = 4 (κ ≤ 1/4), |I(ρ)| ≥ 1/2 implies

    ∫ |M_X(1/2 + i(γ+v), ξ)|² e^{−|v|} dv ≥ [π(σ − 1/2) Y^{σ−1/2} / (6·2^{ω(q)}C_q (4M + T₀ + 5)^{2κ})]².

Sum over the type-II zeros of ξ: Σ_ρ ∫|M_X(1/2 + i(γ_ρ + v))|²e^{−|v|}dv = ∫|M_X(1/2 + it, ξ)|² w_ξ(t)dt with
w_ξ(t) = Σ_ρ e^{−|t−γ_ρ|}. Unit windows give w_ξ(t) ≤ (2 + 2/(e − 1))Z₁ ≤ 3.2 Z₁ for all t, and, all |γ_ρ| ≤ T₀,
w_ξ(t) ≤ (T₀ + 2)Z₁e^{−1}e^{−|t|/2} for |t| ≥ 2T₀ + 2. Mean values: split M_X into the pieces on (X/2^{i+1}, X/2^i],
0 ≤ i ≤ log₂X (Cauchy–Schwarz, factor ≤ log₂X + 2); on a piece with top N', **Lemma 5** (coefficients
μ(𝔞)ξ_χ(𝔞)N𝔞^{−1/2}, support N', D(𝔞) ≤ √2 cont(𝔞)) and **Lemma 6** (g = 1) give
Σ_ξ ∫_{−T}^{T}|M_X(1/2 + it, ξ)|²dt ≤ Φ(q)(log₂X + 2)e²√π(√πM + πX/2)(T + 3)L₁(X), L₁(X) := Σ_{N𝔞≤X}cont(𝔞)/N𝔞
≤ C(log X)². For the tail weight use ∫G(t)e^{−|t|/2}dt = ∫₀^∞ (1/2)e^{−T/2}(∫_{−T}^{T}G)dT and
∫₀^∞ (1/2)e^{−T/2}(T + 3)dT = 5. Altogether Σ_ξ ∫|M_X|²w_ξ ≤ 10(T₀ + 2)·Z₁·Φ(q)(log₂X + 2)e²√π(√πM + πX/2)L₁(X), so

    #(type II) ≤ C(q, κ) (σ − 1/2)^{−2} · Y^{1−2σ} (4M + T₀ + 5)^{4κ} · (M + X) · (T₀ + 2) · Z₁ (log X)³.

(Computed, §D: Lemma 5 holds for the mollifier pieces with LHS/bound ≤ 0.104.)

**4.6 Choice of X and Y; the exponent; B (proved).** Take X := M and Y := M^{1+4κ}. Then (4M + T₀ + 5)^{4κ}
≤ 8M^{4κ} (T₀ ≤ M), M + X = 2M, and M^{1+4κ}Y^{1−2σ} = M^{(1+4κ)(2−2σ)} = M^{A(1−σ)}; and Y^{2−2σ} = M^{A(1−σ)}. Both
zero types give exactly A = 2 + 8κ (as the referee's sketch predicted; computed for κ = 1/6, 103/512, §A). The
factors of log: type I 2 (δ_I^{−2}) + 1 (Z₁) + 1 (1 + 2 log 2Y₁) + 20 (L₆) + 1 (J + 1) + 2 (Y₁/Y) = 27; type II
1 (Z₁) + 3 = 4. Since (σ − 1/2)^{−2} ≤ ((1 + 4κ)/(2κ))² on [σ₁, 1], Theorem D follows on [σ₁, 1]; with §4.1, on
[1/2, 1]. C(q, κ) depends on q through Φ(q), 4^{ω(q)}, C_q², log Nq and on κ through ((1 + 4κ)/(2κ))². ∎

**4.7 Edge cases (proved).** σ near 1/2: the detection degenerates (the factor (σ − 1/2)^{−2}), but there the trivial
bound is used (σ ≤ σ₁, and σ₁ − 1/2 = 2κ/(1 + 4κ) > 0 for every κ > 0). Principal character: excluded from the
family; L(s, ξ) is entire for every member, which Lemma DI needs (for ζ_K the pole would add the residue
M_X(1)Γ(1 − ρ)Y^{1−ρ}, not small for β near 1/2). Dependence on T₀: linear, because the t-aspect of Lemma 5 is not
exploited; harmless for T₀ = (log x)³, fatal for T₀ = x^{ε} — this is why the log-Gaussian norm weight of rh_proof.md
(Lemma 3) is kept. Zeros with β = 1 do not exist (C4) and would not break Lemma DI anyway (Re w = −1/2).

## 5. The explicit formula and the zero sum

Weight: V(u) := P(Z ≥ H log(u/x)), H := (log x)², with Ṽ(s) = x^s e^{s²/(2H²)}/s (**Lemma 3 of rh_proof.md**; the
formula continues meromorphically, and V(u) = (1/2πi)∫_{(c)}Ṽ(s)u^{−s}ds on every line c > 0 by Lemma 3 and Cauchy,
Ṽ having Gaussian decay in |Im s|). For ξ ∈ 𝓕_q(M) put Ψ(ξ) := Σ_𝔞 Λ(𝔞)ξ(𝔞)V(N𝔞) (the sums of rh_proof.md §5, Step 2).

**Imprimitive to primitive (proved).** |Ψ(ξ) − Ψ(ξ*)| ≤ Σ_{𝔭|q}Σ_{j≥1} log N𝔭·V(N𝔭^j) ≤ 3ω(q)·log Nq·log x
(V ≤ 1; V(u) ≤ e^{−H²/2} once log(u/x) ≥ 1).

**Lemma EF (proved).** For ξ* nontrivial primitive (conductor 𝔣 | q, type k), x ≥ 3, H ≥ 1:

    Ψ(ξ*) = − Σ_ρ Ṽ(ρ) + R(ξ*),
    |R(ξ*)| ≤ log x + 2|log A_𝔣| + 2 log(|k| + 4) + 20 + |L'/L(1, ξ̄*)| + C x^{−1/4} H log(N𝔣(|k| + H + 3)),

the sum over the nontrivial zeros with multiplicity, absolutely convergent.
Proof. Lemma 3 and Fubini give Ψ(ξ*) = (1/2πi)∫_{(2)}Ṽ(s)(−L'/L)(s, ξ*)ds. Truncate at Im s = T_n and
Im s = −T_n' (consequence (a) of Lemma PF) and move to Re s = −1/4. On these horizontal segments |Ṽ(s)| =
x^σe^{(σ² − T²)/2H²}/|s| ≤ x²e^{(4 − n²)/2H²} (T ∈ {T_n, T_n'} ⊂ [n, n + 1], |s| ≥ 1), and the second line of consequence (a)
bounds |L'/L|; the segments have length 9/4, so the two horizontal integrals are together
≤ (9/2)(60 + C_PF)·x²e^{(4 − n²)/2H²}·log²(N𝔣(|k| + n + 4)) → 0 as n → ∞.
Poles in −1/4 < Re s < 2: the nontrivial zeros ρ (0 < β < 1: Euler product, (C4), and Λ(it) = WΛ(1 − it) ≠ 0), with
residue −Ṽ(ρ) each; s = 0. Trivial zeros −|k|/2 − n lie in Re s > −1/4 only for k = 0, n = 0. From (C1):
L'/L(s, ξ*) = −2 log A_𝔣 − ψ(s + |k|/2) − ψ(1 − s + |k|/2) − L'/L(1 − s, ξ̄*). If k ≠ 0, L(0, ξ*) ≠ 0 and the residue
at 0 is −L'/L(0, ξ*) = 2 log A_𝔣 + ψ(|k|/2) + ψ(1 + |k|/2) + L'/L(1, ξ̄*) (each ψ-term ≤ log(|k| + 4) + 8 by
Lemma S(c), |k|/2 ≥ 1/2). If k = 0, L has a simple zero at 0 (Λ(0) ≠ 0,
Γ has a simple pole), −L'/L(s) = −1/s + c₀ + O(s) with c₀ = 2 log A_𝔣 − 2γ_E + L'/L(1, ξ̄*), Ṽ(s) = 1/s + log x + O(s),
residue c₀ − log x. On Re s = −1/4: |Ṽ| ≤ 4.2 x^{−1/4}e^{−t²/2H²}, and |L'/L(−1/4 + it)| ≤ 2|log A_𝔣| + 2 log(|k| + |t| + 5)
+ 16 + 3.28 (Lemma S(c): |−1/4 + |k|/2 + it| ≥ 1/4; (C7) at 5/4). ∎ *Uniformity in k:* Lemma S, Lemma PF and (C7) only.

**Proposition 5.2 (proved).** Assume (Z) and P_κ, 0 < κ < 1/4, A = 2 + 8κ, let 0 < λ < 1/A, M := ⌈x^λ⌉,
H := (log x)², T₀ := H log x = (log x)³, δ' := 1 − λA > 0. Then, for some c = c(q, κ, λ) > 0,

    Σ_{ξ∈𝓕_q(M)} |Ψ(ξ)| ≪_{q,κ,λ} x · exp(−c (log x)^{1/3}(log log x)^{−1/3}).

Proof. By the two displays above, Σ_ξ|Ψ(ξ)| ≤ Σ_ξ Σ_ρ |Ṽ(ρ)| + Σ_ξ(|R(ξ*)| + 3ω(q)log Nq log x), ρ over the zeros of
L(s, ξ*). Note |Ṽ(ρ)| = x^β e^{(β² − γ²)/2H²}/|ρ| ≤ e^{1/2H²}x^β e^{−γ²/2H²}/|ρ|.
- R-terms: by consequence (b) of Lemma PF and §3 (zeros of L(s, ξ̄*) with |γ| < 1 have 1 − β ≥ η),
  |L'/L(1, ξ̄*)| ≤ Z₁/η + C log x ≪ (log x)²; so Σ_ξ |R| ≪ M(log x)² = o(x).
- Zeros with β ≤ 1/2: for |γ| < 1, 1 − ρ̄ is a zero of L(s, ξ̄*) with real part ≥ 1/2, so β ≥ η and |ρ| ≥ η; for
  |γ| ≥ 1, |ρ| ≥ |γ|. With Lemma ZC, Σ_ρ e^{−γ²/2H²}/|ρ| ≤ Z₁/η + Σ_{n≥1} 2 ZC(n + 1)e^{−n²/2H²}/n ≪ (log x)³. Over the
  family: ≪ M x^{1/2}(log x)³ ≤ x^{1/2+λ}(log x)⁴ = o(x), since λ < 1/A < 1/2.
- Zeros with β > 1/2, |γ| > T₀: ≤ 3x Σ_{n≥T₀} 2 ZC(n + 1)e^{−n²/2H²} ≪ x(log x)² H e^{−(log x)²/2}, times ≤ 3M members:
  negligible.
- Zeros with β > 1/2, |γ| ≤ T₀: |Ṽ(ρ)| ≤ 2e^{1/2H²}x^β (|ρ| > 1/2). For each such ρ, x^β = x^{1/2} + log x ∫_{1/2}^{β}x^σdσ;
  summing (a finite sum) and using N(σ) := N_q(σ; M, T₀) = 0 for σ > 1 − η (§3):

      𝒮 := Σ x^β = x^{1/2}N(1/2⁺) + log x ∫_{1/2}^{1−η} x^σ N(σ) dσ.

  On [1/2, σ₁] use N(σ) ≤ N(1/2) ≪ T₀ M log M (§4.1): the contribution is ≪ x^{σ₁ + λ}(log x)⁵ = x^{1 − (1/A − λ)}(log x)⁵.
  On [σ₁, 1 − η] use Theorem D and M^A ≤ 2^A x^{1−δ'}: x^σ N(σ) ≤ C(T₀ + 4)(log M)^{27}·x·2^A x^{−δ'(1−σ)}, so

      log x ∫_{σ₁}^{1−η} x^σ N(σ) dσ ≤ C' x (log x)^{30} x^{−δ'η}/δ' = C' x (log x)^{30} exp(−δ'c₃(log x)^{1/3}(log log x)^{−1/3})/δ'.

  Every other term is O(x^{1−δ''}). ∎

Conversely, if λA > 1 the bound x^σ N(σ) ≤ C x(log x)^{30}(M^A/x)^{1−σ} exceeds x on all of [σ₁, 1 − η]: the method then
gives nothing, so θ(2 + 8κ) < 1 is exactly the condition under which this argument works (not a statement about the
truth of the asymptotic).

**Blocks in 1 − σ.** The integral above is the same as cutting 1 − β into blocks [j/log x, (j + 1)/log x): on a block
x^β ≤ x e^{−j}. The blocks that meet [σ₁, 1 − η] have j/log x ≤ 1 − σ₁ = 1/A; for them σ_j := 1 − (j + 1)/log x
≥ σ₁ − 1/log x ≥ 1/2 (x large: log x ≥ (1 + 4κ)/(2κ)), so Theorem D applies at σ_j, and with
M^{A(1−σ_j)} = (M^A)^{(j+1)/log x} ≤ (2^A x^{1−δ'})^{(j+1)/log x} the count is
≤ N(σ_j) ≤ C(T₀ + 4)(log M)^{27}·2^{A(j+1)/log x}·e^{(1−δ')(j+1)}. Here A(j + 1)/log x ≤ 1 + A/log x ≤ 2, so
2^{A(j+1)/log x} ≤ 4 (log x ≥ 4 > A). With T₀ + 4 ≤ 2(log x)³ and log M ≤ log x (M ≤ 2x^λ, λ < 1/2, log x ≥ 4), the block contributes
≤ x e^{−j}·8C(log x)^{30}e^{(1−δ')(j+1)} ≤ 8eC(log x)^{30}x e^{−δ'j}. A zero with 1 − β ≥ η lies in a block with
(j + 1)/log x > η, so j ≥ j₀ := ⌊η log x⌋ > η log x − 1, and, with 1 − e^{−δ'} ≥ δ'/(1 + δ') ≥ δ'/2 (δ' ≤ 1),
Σ_{j≥j₀} e^{−δ'j} = e^{−δ'j₀}/(1 − e^{−δ'}) ≤ (2e/δ')e^{−δ'η log x}: the blocks sum to
≤ 16e²C x(log x)^{30}e^{−δ'η log x}/δ', the bound of the integral version. A **dyadic** cut
1 − β ∈ (2^{−i−1}, 2^{−i}] is too coarse: there x^β ≤ x^{1 − 2^{−i−1}} while the count is only bounded by
N(1 − 2^{−i}) ≪ M^{A 2^{−i}}, and the product x·(M^{2A}/x)^{2^{−i−1}} needs M^{2A} < x, i.e. θ < 1/(2A): half the
exponent is lost. So the decomposition must have width O(1/log x) in 1 − σ (equivalently, the Stieltjes integral).

**Why the classical region c/log V does not suffice (the computation).** With a region 1 − β ≥ η₀ := c/log V*,
V* ≍ x^λ, one has η₀ log x ≍ c/λ, a constant, and the last display becomes ≍ x(log x)^{30}e^{−cδ'/λ}/δ': not o(x),
by a factor (log x)^{30}. Even a log-free density estimate (B = 0, no T₀) would leave x e^{−cδ'/λ}/δ' = O(x): an
upper bound of the right order, not an asymptotic (compare the ≍ of Duke's Theorem 3.1 reported in rh_proof.md §10;
whether this is the reason there was not checked).
What is needed is η(x) log x − (30/δ')log log x → ∞; since log V* ≈ λ log x, a region of width
η(V) ≥ C log log V/log V gives η log x ≈ (C/λ)log log x, so any C > 30λ/δ' would do; Coleman's width (log V)^{−2/3}(log log V)^{−1/3} gives η log x ≍ (log x)^{1/3}(log log x)^{−1/3}, far more
than needed. (Illustration with placeholder constants c = δ' = 0.1, B = 30, §G: log₁₀ of bound/x at log x = 10²⁴ is
−1.1·10⁵ for the VK width and +720 for the classical width — heuristic, not evidence.) The zero-free region enters at
three places; only the first needs this strength: (a) N(σ) = 0 for σ > 1 − η; (b) |ρ| ≥ η for zeros near 0;
(c) |L'/L(1, ξ̄*)| ≪ (log x)² (any width ≥ (log x)^{−C} would do for (b), (c)).

## 6. Removal of the smoothings (proved in rh_proof.md §5; hypotheses checked)

Steps 1–4 of rh_proof.md §5 are used verbatim, with the same objects: Δ_a := ωx^{−ε}, r := ⌈4/ε⌉, K₁ := ω^{−1}x^{2ε},
H := (log x)², and the height (log x)³ (T₁ there, T₀ here). Hypotheses: Lemma 1 needs 0 < Δ_a < ω ≤ π and r ≥ 2 (ω > π
is split into halves, as there); Lemma 2 needs r₁ ≥ 2 (x large); Lemma 3 needs x, H > 0. What the steps need from the
exponent is only ε > 0 and θ + 2ε < 1/2: Step 1 (sandwich ψ₋ ≤ 1_arc ≤ ψ₊ and the norm edge,
Err ≪ (log x)(ωx/H + √x) = o(ωx) since ω ≥ x^{−θ}, θ < 1/2; the referee's correction — Err over supp ψ₊, width
≤ 2ω — gives the same bound); Step 2 (orthogonality of χ, Fourier expansion of ψ_± — absolutely convergent since
|ψ̂_±(k)| ≪ |k|^{−r−1} — and the unit sum of §2; the classes not prime to q contribute exactly 0); Step 3 (main term from
(C5): the (χ₀, 0) term is (2/π)ωx/Φ(q)·(1 + o(1)), since Δ_a = o(ω)); Step 4 (the |k| > K₁ tail is ≪ x(2r)^r x^{−εr}
= o(ωx)). Lemmas 4, 8, 9, Step 5, the Main Lemma and §6–§7 of rh_proof.md are **not** used. In place of rh_proof.md's
(★) and the Main Lemma, the remaining terms are bounded by Proposition 5.2:

    (★★)  |(4/Φ) Σ_χ χ̄(ζ) Σ_{(χ,k)≠(χ₀,0), |k|≤K₁} ψ̂_±(k) Ψ(ξ_{χ,k})| ≤ (4/Φ)(ω/π) Σ_{ξ∈𝓕_q(M)} |Ψ(ξ)|,  M := ⌈x^{θ+2ε}⌉,

because |ψ̂_±(k)| ≤ ω/π, each (χ, k) with |k| ≤ K₁ ≤ x^{θ+2ε} is ξ_χλ^m with |m| ≤ (K₁ + 3)/4 ≤ M, and distinct
(χ, k) give distinct members. The right side does not depend on φ, and on ω only through the factor ω.

## 7. Proof of Theorem 1′ and the explicit exponents

Given θ < 1/(2 + 8κ), choose ε > 0 with **(E)** (θ + 2ε)(2 + 8κ) < 1 (then θ + 2ε < 1/2). Proposition 5.2 with
λ = θ + 2ε makes (★★) o(ωx) uniformly; Steps 1–4 give S_± = (2/π)ωx/Φ(q)(1 + o(1)) and S = (2/π)ωx/Φ(q)(1 + o(1)),
uniformly in φ and ω ∈ [x^{−θ}, 2π] (every error term of Steps 1–4 is uniform, rh_proof.md and referee §3(b); the
zero sum is independent of φ, ω). ∎ (Computed, §A: for θ = 0.2770 the choice ε = 1.86·10⁻⁵ satisfies (E).)

**Explicit P_κ.** Wu, Theorem 1.1 (read, arXiv 1604.08551, LaTeX): "Let χ be a Hecke character of F with analytic
conductor C(χ). We have |L(1/2, χ)| ≪_{F,ε} C(χ)^{1/4 − (1−2θ)/16 + ε}, where θ is any constant towards the
Ramanujan–Petersson conjecture"; C(χ) = C(χ_∞)C(χ_fin), for varying Hecke characters, so uniform in the infinity type
and (via χ|·|^{it}) in t. For ξ*|·|^{it}, C ≍ N𝔣(2 + |k|/2 + |t|)² (the complex place has L-factor Γ_ℂ(s + |k|/2 + it);
rh_proof.md §3(e)). Blomer–Brumley, Theorem 1 (read, arXiv 1003.0559v4, LaTeX): "Hypothesis H_n(δ_n) holds with
δ₂ = 7/64", H_n(δ): for any number field, any cuspidal π on GL_n(𝔸) with unitary central character and any place v,
m(π, v) ≤ δ. Hence κ = 1/4 − (1 − 14/64)/16 + ε = 103/512 + ε, A = 231/64 + 8ε, **θ < 64/231 = 0.27706**. (Computed,
exact fractions, §A: also ϑ = 1/9 (Kim–Shahidi, quoted in Blomer–Brumley) → 18/65; ϑ = 3/10 (Luo–Rudnick–Sarnak) → 5/19;
ϑ = 0 → 2/7; any ϑ < 1/2 gives θ > 1/4.) **Söhne** (J. Number Theory 66 (1997), main theorem, quoted verbatim in
Wu's introduction, read): L(1/2 + it, χ) ≪_{ε,F} (C_∞(χ|·|^{it})N𝔣)^{1/6+ε} + N(𝔣₀)^{1/2+ε} + N(𝔣/𝔣₀)^{1/4+ε}, any 𝔣₀ | 𝔣;
with 𝔣₀ = O and fixed q, κ = 1/6 + ε and **θ < 3/10** — quoted (secondary). **Michel–Venkatesh**: κ = 1/4 − δ,
θ < 1/(4 − 8δ), δ not explicit — quoted (checked) in rh_proof.md §3(e), not re-read here.

## 8. What the theorem now rests on

| input | where used | status |
|---|---|---|
| §2 of rh_proof.md (the characters ξ_{χ,k}, unit sums) | everywhere | proved in rh_proof.md (refereed) |
| Lemma 1 (angular cut-offs), Lemma 2 (polar boxes) | §6, Steps 1, 2, 4 | proved in rh_proof.md (refereed) |
| Lemma 3 (log-Gaussian weight, Ṽ) | §6 Step 1; Lemma EF | proved in rh_proof.md (refereed) |
| Lemma 5 (mean value with content weight, all coefficients, both signs of m) | Theorem D, types I and II | proved in rh_proof.md (refereed); used with arbitrary coefficients, support bound N, any T > 0 — inside its hypotheses |
| Lemma 6 (divisor moments with content weight), g = 1, 2 | Theorem D | proved in rh_proof.md (refereed) |
| Lemma 7 (Gallagher) | Theorem D, type I | proved in rh_proof.md (refereed); F a finite Dirichlet polynomial |
| Steps 1–4 of rh_proof.md §5 | §6 | proved in rh_proof.md (refereed); hypotheses the same (θ + 2ε < 1/2) |
| Lemmas S, ZC, PF, DI, ZFR (deduction), EF; Theorem D; Proposition 5.2; (★★); Theorem 1′ | §2–§7 | **proved here** |
| (C1) Hecke's functional equation, factor Γ(s + \|k\|/2), \|W\| = 1 | Lemmas ZC, PF, EF; Lemma DI (L(s, ξ) entire), hence Theorem D | classical (Hecke; Neukirch VII §8; IK §3.8), not re-opened; computed check (§C) |
| (C2) Hadamard product | Lemmas ZC, PF | classical (Davenport Ch. 12; IK Thm 5.6), not re-opened; identities only |
| (C3) Γ and ψ: the shift Γ(z + 1) = zΓ(z) and the poles, Binet/Stirling on Re z > 0, Cahen–Mellin | Lemma S; Lemma DI (bound (Γ) on −1/2 ≤ Re w ≤ 2, through the shift) and §4.5, hence Theorem D; Lemma EF | classical (Whittaker–Watson Ch. XII), not re-opened; constants computed (§B, §H) |
| (C4) L(1 + it, ξ) ≠ 0 | Lemma ZC (Λ(0) ≠ 0), hence **Theorem D** (the window counts Z₁ in §4.1, §4.4, §4.5) and Lemma PF; Lemma ZFR (compact part), Lemma EF | classical (Hecke; Lang Ch. XV; Narkiewicz Ch. 7), not re-opened |
| splitting law of ℤ[i] ((2) = (1 + i)², p ≡ 1 (mod 4) split, p ≡ 3 (mod 4) inert) | (T2) in §4.3, hence Theorem D | classical, not re-opened; r₀ = 1 ∗ χ_{−4} and r₀ ≤ d computed (§E, §H) |
| (C5) prime number theorem for ζ and L(s, χ_{−4}) | main term (Step 3) | classical (Davenport Ch. 18, 20), not re-opened |
| (C6) polynomial growth in vertical strips, per character | Lemma DI (contour shift) | classical (Rademacher 1959; IK §5.1), not re-opened; qualitative only |
| P_κ: Wu, Theorem 1.1 | Theorem D, type II | quoted (checked): statement and the definition of C(χ) read (LaTeX) by the referee and here |
| ϑ = 7/64 over every number field: Blomer–Brumley, Theorem 1 | value of κ | quoted (checked): read here (LaTeX: introduction, Definition H_n(δ), Theorem 1) |
| P_κ: Michel–Venkatesh (alternative, δ not explicit) | — | quoted (checked) in rh_proof.md and by the referee; not re-read here |
| Söhne's bound (κ = 1/6 + ε) | only for θ < 3/10 | **quoted (secondary)**: verbatim in Wu's introduction (read), zbMATH review (referee) |
| **(Z) Coleman, Mathematika 37 (1990), Theorem 2** | Lemma ZFR | **quoted (secondary)**: verbatim in Sun–Zhang (read); paraphrased with the modulus in X. Li, arXiv 2111.05403 (read); described in the zbMATH review (Kaczorowski); the generality (all Grössencharacters mod 𝔣) is Sun–Zhang's definition |
| Coleman's Theorem 1 (growth), bounds for 1/L | — | **no longer used** |

The only secondary quotation left that Theorem 1′ with θ < 64/231 (hence R_h) depends on is Coleman's Theorem 2 as
quoted verbatim by Sun–Zhang. The textbook facts (C1)–(C6) and the splitting law of ℤ[i] are "classical, not re-opened"; that is a separate
category from "secondary", and their uniformity in k is not quoted but derived in §2.

## 9. Generality (remark, without proof; heuristic)

Steps that use K = ℚ(i):
- **Units of order 4** (rh_proof.md §2): the class-and-sector indicator expands into ξ_{χ,k} with k ≡ k_χ (mod 4),
  λ = e_4, the factor 4/Φ(q) and the main term (2/π)ωx/Φ(q). For ℚ(√−3): units of order 6, k ≡ k_χ (mod 6),
  λ = e_6, factor 6/Φ(q), main term (3/π)ωx/Φ(q) in the normalisation by elements. For other imaginary quadratic
  fields: units ±1, k ≡ k_χ (mod 2).
- **Class number 1**: elements ↔ ideals up to units. For h_K > 1 the element condition z ≡ ζ (mod q) with arg z in a
  sector is detected on principal ideals; one needs, besides χ, the characters of the class group (ray class group
  mod q), and the family becomes {ψ λ^m}, ψ running over ray class characters mod q extended by Grössencharacters of
  each type — Hecke's general setting, in which (Z) as quoted, Wu's and Michel–Venkatesh's theorems and (C1)–(C6) are
  all stated.
- **The lattice ℤ[i] in Lemma 5** (spacing): Im(α'·conj(εβ')) ∈ ℤ∖{0} gives |Δ arg| ≥ 1/N. In ℤ[ω] the same
  quantity lies in (√3/2)ℤ∖{0}: the same lemma with another constant. In a general O_K the angle of a non-principal
  ideal is defined through a Grössencharacter; for 𝔞, 𝔟 in the same class, 𝔞𝔟̄ is principal, so the spacing argument
  should go through class by class, with constants depending on h_K and |d_K|. cont(𝔞) (largest rational integer
  dividing 𝔞) and Lemma 6 generalise directly. Lemma 2 (polar boxes) uses ℤ²: replace by the lattice O of covolume
  √|d_K|/2.
- **Field-independent parts**: Lemma 3, Lemma 7, Lemmas S, ZC, PF, DI, Theorem D's structure, the explicit formula
  (for imaginary quadratic K the same Γ(s + |k|/2) with A_𝔣 = √(|d_K|N𝔣)/(2π)); the constants c_K = −ζ_K'/ζ_K(2) and the
  ZC constant change. Wu's exponent does not depend on the field, so the same θ < 1/(2 + 8κ), κ = 103/512, is expected
  for every imaginary quadratic field and every fixed modulus — not proved here.

## 10. Numerical verification (`rh_density_check.py`; output `rh_density_check_out.txt`: 122 checks, 0 failed, 145 s on the shared machine)

- **§A exponent.** For κ ∈ {0.02, 0.06, 1/8, 1/6, 0.2, 103/512, 0.24, 1/4}: with X = M, Y = M^{1+4κ} the sup over σ of
  min(trivial, detection)/(1 − σ) equals 2 + 8κ to 10⁻⁹; the minimum over (X, Y) = (M^a, M^b), a ∈ [0.5, 1.5],
  b ∈ [1, 2.5] (powering allowed for short type-I levels) is never below 2 + 8κ, attained at a = 1, b = 1 + 4κ. Exact
  fractions: κ = 103/512, A = 231/64, θ* = 64/231, σ₁ = 167/231; 18/65, 5/19, 2/7, 3/10; σ₁ − 1/2 = 2κ/(1 + 4κ).
- **§B constants.** Binet 0 ≤ g ≤ t/12; |ψ − log| ≤ 1.27 (≤ 4/3); |ψ(z)| − log(|z| + 3) ≤ 3.05 (≤ 8) for Re z ≥ −1/4,
  |z| ≥ 1/4; |μ| ≤ 1/(12 Re z); |Γ(a + iv)|e^{|v|} ≤ 1.97 (≤ 3);
  Re 1/(s₁ − ρ) ≥ 0.4; −ζ_K'/ζ_K(2) = 0.4809, (5/4): 3.2769; ZC constant 1.67 ≤ 3.35; ∫(1 + |v|)e^{−|v|} = 4;
  ∫(1/2)e^{−T/2}(T + 3) = 5; 2 + 2/(e − 1) = 3.164; 2 + (log 2)²/2 = 2.240; the dyadic tail bound ≤ 4.6·10⁻⁵ on 10² ≤ Y ≤ 10¹²
  (proved for all Y ≥ 100 in §4.3 (T3); tested further in §H).
- **§C functional equation and zeros.** Root numbers from the theta split at two split points and three values of s:
  |W| = 1 to 10⁻¹² and W constant to ≤ 7.5·10⁻¹⁵ for 10 characters (q = 1: k = 4, 8, −12, W = 1; conductor 2 + i: k = −1, 1, 2, 3,
  5, 6, 9; e.g. W = 0.6 + 0.8i for k = 2). Negative controls fail as they must (Γ(s + |k|); k ≢ k_χ; an imprimitive
  character with N𝔣 = 5). W^{−1/2}Λ(1/2 + it) is real (relative imaginary part < 10⁻⁸); sign changes in [0, 14]: 3
  (L(s, λ), first at 7.53) and 6 (conductor 2 + i, k = 1, first at 3.62), ≤ 0.18 of the Lemma ZC bound per window.
- **§D detection.** The detection identity for L(s, χ_{−4}) at its first zero and at 0.8 + 3i (to 10⁻⁸); the numeric
  mollified tail 5.5·10⁻⁶; partial summation (300 random cases, max ratio 0.998 — the inequality is nearly sharp);
  Gallagher with multiplicity (20 clustered sets); Lemma 5 for the type-I vectors (≤ 0.056) and the mollifier pieces
  (≤ 0.104), conductor 2 + i, χ₁ in the coefficients, M ∈ {1, 10, 100}.
- **§E ideal arithmetic.** r₀(n) = Σ_{d|n}χ_{−4}(d); Σ_{N𝔞=n}d(𝔞) = (r₀ ∗ r₀)(n) ≤ d₄(n) (n ≤ 3000);
  Σ_{n≤u}d₄(n) ≤ u(1 + log u)³ (u ≤ 10⁵); Σ_{N𝔞=n}Λ(𝔞) = Λ(n)(1 + χ_{−4}(n)) (n ≤ 3000). (The first three are
  proved in §4.3 (T1), (T2); here and in §H they are tests.)
- **§F sanity only (not evidence for θ).** All 2,759,572 Gaussian prime-power elements with N z ≤ 1.04·10⁷:
  ψ_K(10⁷)/10⁷ = 0.99987; 48 sectors (6 anchors, widths x^{−1/4}, 2x^{−1/4}, 4x^{−1/4}, x^{−0.27}, both orientations):
  Λ-sum/((2/π)ωx) between 0.9880 and 1.0141, mean 1.0006; 960 sector–class pairs (4 classes mod 2 + i, 16 mod 5):
  z-scores mean +0.016, rms 0.645, max |z| 2.25 (under-dispersed, as in rh_proof.md §9). With the weight V
  (H = (log x)²): max over the nontrivial family |k| ≤ 200 of |Ψ(ξ)| is 4.4–5.5·√x = 1.4–1.7·10⁻³x; the element sums
  vanish for k ≢ k_χ (mod 4).
- **§G** illustration of VK vs classical width (heuristic, placeholder constants).
- **§H the repaired inequalities** (tests of proved statements, not their proofs). Consequence (a) of Lemma PF on 6600
  synthetic ordinate sets (both signs, window counts in place of ZC): T_n exists, the new bound holds (max ratio 0.82);
  in 1249 cases a zero with |γ − T_n| < 1 lies outside [n − 1/2, n + 3/2], and in 279 cases the old bound Z_n(2Z_n + 2)
  is exceeded (negative control {n − 0.6}: 1.67 against 0); the second line (≤ 47.3ℓ² on a grid of N𝔣, k, n).
  (T1) for d₁, …, d₄ and Σ_{m≤u}1/m ≤ 1 + log u at every integer u ≤ 10⁶ (hence every real u ∈ [1, 10⁶ + 1)); d ∗ d = d₄;
  (T2): r₀ = 1 ∗ χ_{−4}, r₀ ≤ d, (r₀ ∗ r₀) ≤ d₄ and Σ_{N𝔞≤u}d(𝔞) ≤ u(1 + log u)³ for n, u ≤ 10⁶, the prime-power values
  of r₀, and Σ_{N𝔞=n}d(𝔞) = (r₀ ∗ r₀)(n) by factorising every ideal of norm ≤ 22000. (T3): t_j decreasing (j ≤ 24,
  L ∈ [4, 200]), ratios ≤ 1.08·10⁻⁹ ≤ 1.3·10⁻⁹, t₀(log 100) = 4.642·10⁻⁵, Σ_j t_j(log Y) ≤ 4.7·10⁻⁵ on Y ∈ [100, 10³⁰⁰];
  the true tail (7.6·10⁻⁹ at Y = 100) is below the bound for Y = 100, 150, 200, 300. (Γ) on −1/2 ≤ Re w ≤ 2,
  1 ≤ |Im w| ≤ 500 (51 × 368 points): max ratio 0.355; the two intermediate bounds 0.959 and 0.474; the shift to 10⁻²⁰;
  |Γ(2 + iv)|e^{π|v|/2}/|v|^{3/2} → √(2π) (the power 3/2 is needed at Re w = 2); Cahen–Mellin to 10⁻¹². Blocks of §5:
  2^{A(j+1)/log x} ≤ 4 and the geometric sum on grids.

## 11. Sources: what was tried, what opened

**None of the primary sources opened; none was found to contain a residue-class sector theorem.**

| source | attempt | status |
|---|---|---|
| Coleman, A zero-free region for the Hecke L-functions, Mathematika 37 (1990) 287–304, doi 10.1112/S0025579300013000 | Wiley pdf and abstract pages (onlinelibrary and londmathsoc subdomain): 403; doi.org → Wiley; Cambridge Core landing page: opened (only the first sentences as "abstract": "Let K be an algebraic number field of degree n and discriminant d. Let K(1), …, K(n) be the embeddings …"), its PDF link redirects to the paywall; Semantic Scholar: "CLOSED", abstract elided; OpenAlex: closed, no repository copy; web search: no PDF; zbMATH review (Zbl 0721.11050, Kaczorowski) read via the API | **could not open**; abstract or snippet only |
| Söhne, An upper bound for Hecke zeta-functions with Groessencharacters, J. Number Theory 66 (1997) 225–250, doi 10.1006/jnth.1997.2167 | Elsevier API metadata: openArchiveArticle = true (open archive); ScienceDirect /pdf and /pdfft: 403 to automated requests; full-text API: 401 (key required); doi.org → linkinghub → ScienceDirect | **could not open** (open access in principle; the statement is known from Wu's introduction, verbatim) |
| Koval'chik, Density theorems and the distribution of primes in sectors and progressions, Dokl. Akad. Nauk SSSR 219 (1974) 31–34 | Math-Net.Ru (dan38636 and getFT): connection timed out (curl, 40 s) and refused (ECONNREFUSED); geodesic.mathdoc.fr record opened (bibliography only, language ru, full-text link to Math-Net.Ru); zbMATH: no entry found under the author | **could not open** |
| Koval'chik, Litovsk. Mat. Sb. 15 (1975) 133–151 | via Math-Net.Ru (unreachable) | could not open |
| Ricci, Local distribution of primes, PhD thesis, Michigan 1976 | web search (no ProQuest record found) | **could not open** (not located) |
| Maknys, On the distance between consecutive prime ideal numbers in sectors, Acta Math. Hungar. 42 (1983) 131–138, doi 10.1007/BF01960557 | Springer landing page opened: paywalled; reference list read (cites Koval'chik 1974, Maknys 1975–1977, Kaufman 1978, Heath-Brown 1977, Huxley 1972, Johnson 1979, Jutila 1969); zbMATH review by Heath-Brown exists, text not served by the API | **could not open**; reference list only |
| Maknys, Zeros of Hecke Z-functions and the distribution of primes of an imaginary quadratic field, Lith. Math. J. 15 (1975) 140–149 (Zbl 0329.12005); also Litov. Mat. Sb. 15 (1975) 173–184, 16 (1976) 173–180, and "Refinement of the remainder term … of an imaginary quadratic field in sectors", Lith. Math. J. 17 (1977) (Zbl 0372.12006); Large value estimates for Hecke polynomials (1992, Zbl 0768.11047) | zbMATH API: titles; "contents unavailable due to conflicting licenses" or no review text; Jutila's review of the 1992 paper: "large value estimates for Dirichlet polynomials involving Hecke Grössencharacters in an imaginary quadratic field … analogous to … Dirichlet L-functions" | **could not open**; titles and one review only |
| Kaufman, On the short equations of A. F. Lavrik, Zap. Nauchn. Sem. LOMI 76 (1978) 124–158 | Math-Net.Ru unreachable | **could not open** |
| Blomer–Brumley, On the Ramanujan conjecture over number fields, arXiv 1003.0559v4 (Ann. of Math. 174 (2011), from memory) | arXiv source | **read the relevant section** (introduction, Definition H_n(δ), Theorem 1) |
| Wu, Burgess-like subconvexity for GL₁, arXiv 1604.08551 | referee's download (byte copy) | **read the relevant section** (Theorem 1.1, Corollary, introduction with Söhne's bound) |
| Sun–Zhang, arXiv 2008.09677 | referee's download | **read the relevant section** (§1.1 Hecke characters mod 𝔣; §2 zero-free region, Coleman's Theorem 2 verbatim; zero density; explicit formula) |
| X. Li, Prime values of a sparse polynomial sequence, arXiv 2111.05403 | arXiv source | **read the relevant section** (§"Improved zero free region type bound": Coleman's Theorems 1 and 2 paraphrased, with modulus q via log q) |

Downloads (each in its own folder of this directory): `dl_coleman/` (Cambridge landing page), `dl_sohne/` (Elsevier
API metadata), `dl_koval/` (empty), `dl_maknys/` (Springer landing page), `dl_bb/` (Blomer–Brumley source),
`dl_2111/` (Li source).

## 12. Doubts and limits

1. **(Z) is secondary, and so is its generality.** The proof needs the region for *every* Grössencharacter of
   conductor dividing q, including those with χ(i) ≠ 1 (infinity types k ≢ 0 mod 4), which the class condition
   z ≡ ζ (mod q) on elements cannot avoid. Sun–Zhang's definitions (λ a basis of the infinite-order characters
   *mod 𝔣*) give this; Li's paraphrase (characters θλ with θ of modulus q), Duke's Urbialis statement (χλ^m mod q,
   c₂ > 5/7, for the zero-free part) and the zbMATH review ("Größencharacters mod 𝔣") agree. If Coleman's Theorem 2
   were stated only for λ of modulus (1) times ray class characters, the characters with k ≢ 0 (mod 4) would not be
   covered. Not checkable without the original.
2. **Textbook facts not re-opened.** (C1)–(C6). Their k-uniformity is derived (Lemmas S, ZC, PF); (C1)'s Γ-factor is
   computed (§C). Theorem/section numbers of the textbooks are from memory.
3. **No rate.** c₂(q) (compactness, including a possible exceptional zero) is not computed; the theorem is an
   asymptotic with o(1). Proposition 5.2 has the rate exp(−c(log x)^{1/3}(log log x)^{−1/3}), but Steps 1 and 3
   contribute relative errors O(1/log x) and O(H^{−2}).
4. **Theorem D is not sharp near σ = 1.** A sharper density estimate near σ = 1 (Halász–Montgomery or Huxley-type large
   values in the hybrid family, or a Vinogradov-type density near 1) would raise θ; not attempted. The power 27 of
   log is crude.
5. **Exponent via Wu needs ϑ for GL₂ over ℚ(i)** (Blomer–Brumley, read); Wu's "θ any constant towards the
   Ramanujan–Petersson conjecture" is taken to mean H₂(θ) over F (Wu does not define it further in the text read), as the
   referee did.
6. Numerics test identities, inequalities and the main term at x = 10⁷; they cannot test θ > 1/4.

## Files

- `rh_density.md` (this file), `rh_density_check.py`, `rh_density_check_out.txt`
- downloads: `dl_coleman/`, `dl_sohne/`, `dl_koval/`, `dl_maknys/`, `dl_bb/`, `dl_2111/`; prototypes: `proto/`


## Referee report (26 November 2026)

An independent referee (own code in [referee2/](referee2/)) re-derived Theorem D by the classical detection argument before
reading the proof and obtained the same exponent (2 + 8κ)(1 − σ) for both classes of zeros, with X = M and Y = M^(1+4κ).
**Verdict: proved, given Coleman's Theorem 2 as quoted (with the generality of Doubt 1) and the pointwise bound P_κ; no
failing step.** Checked: only two coefficient vectors occur (the mollifier, and 1∗μ_(≤X) with |c| ≤ d), no squaring or
powering; the sum over m is inside the mean value for both classes; zeros are counted with multiplicity through Gallagher's
inequality; the zeros with |γ| > T₀ cost e^(−(log x)²/2) by the zero count alone; the block width 1/log x is needed and
suffices; the exceptional zero affects only the characters of type 0 and costs a relative x^(β₁−1). Lemma ZC tested on PARI
Hecke characters with k = 400, 1000, 401, −999 (largest window count 0.22 of the bound; zero counts 70, 86, 83, 99 against
main terms 70.1, 85.5, 83.7, 99.1); Lemma EF tested at k = 12 and 100 (agreement 1.4·10⁻⁷ and 1.1·10⁻⁶).
The referee listed five slips, none affecting the theorem; they are repaired in the text above. Not opened by the
referee: Coleman's text and the textbook facts.

**Corrections after the referee report.** (1) Consequence (a) of Lemma PF (lines 167–187) now chooses a height for
each sign (T_n and T_n') and counts the zeros with |γ − T| < 1 that lie outside [n − 1/2, n + 3/2]: at most
ZC(n − 1/2) + ZC(n + 3/2) of them, each at distance > 1/2, so the bound gains the term 2(ZC(n − 1/2) + ZC(n + 3/2)) and a
second line ≤ (60 + C_PF)log²(N𝔣(|k| + n + 4)); the horizontal integrals of Lemma EF are redone with it (lines
394–398), with the same conclusion. (2) §4.3 (lines 276–307) proves the three facts that were labelled computed: (T1)
Σ_{n≤u}d₄(n) ≤ u(1 + log u)³ for real u ≥ 1, by induction on the number of factors from Σ_{m≤u}1/m ≤ 1 + log u; (T2)
Σ_{N𝔞=n}d(𝔞) = (r₀ ∗ r₀)(n) ≤ d₄(n), from r₀ ≤ d (splitting law of ℤ[i], values at prime powers); (T3)
𝒯(Y) ≤ 4.7·10⁻⁵ for every Y ≥ 100, by dyadic blocks, monotonicity of each block bound in log Y, a ratio ≤ 1.3·10⁻⁹
between blocks and the value 4.642·10⁻⁵ at Y = 100; the constant 10⁻⁴ used downstream is unchanged, and J + 1 ≤ 2 log Y₁
is derived (lines 306–307). (3) (C3) (lines 103–113) is stated in the form used: poles and shift relations of Γ and ψ,
Binet's formulas on Re z > 0 only, and Cahen–Mellin; the proof of Lemma DI (lines 251–272) proves the bound (Γ),
|Γ(a + iv)| ≤ 10|v|^{3/2}e^{−π|v|/2} for −1/2 ≤ a ≤ 2, |v| ≥ 1, through Γ(w) = Γ(w + 1)/w when a < 1/2, and writes out
the horizontal segments; the uses of the shift are marked at lines 135 and 344. (4) The block paragraph of §5 (lines
438–447) has 2^{A(j+1)/log x} in place of 2^A, bounds it by 4 on the blocks that meet [σ₁, 1 − η] (where Theorem D
applies at σ_j = 1 − (j + 1)/log x ≥ 1/2 for x large), and redoes the sum: ≤ 16e²C x(log x)^{30}e^{−δ'η log x}/δ', as in
the integral version. (5) (C4) is listed among the inputs of Theorem D in the table of §8 (line 517) and in the table of
the Verdict (line 43); the §8 rows of (C1) and (C3) (lines 514, 516) now also name Lemma DI, and a row for the splitting
law of ℤ[i], used in (T2), is added (line 518; also line 529). §10 (lines 555, 564–565, 576–577, 585–596) and section H
of rh_density_check.py test the repaired inequalities (122 checks, 0 failed).
No statement of a theorem or lemma changed except consequence (a) of Lemma PF (the added term, the second line, both
signs); (C3) is restated as a list of classical facts.
