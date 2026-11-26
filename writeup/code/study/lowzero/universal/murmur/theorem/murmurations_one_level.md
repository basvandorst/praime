# Murmurations at one squarefree level, natural weights, without hypothesis for windows of relative length ≥ (log N)^(−B)

This note merges three refereed notes of the project into one statement and one proof:
`mprime_theorem.md` (Theorem 2: one prime level, under GRH for Dirichlet L-functions; steps T1–T6),
`mprime_unconditional.md` (Theorem U: one prime level, no hypothesis, windows ≥ (log N)^(−B); steps U1–U6) and
`level_theorem.md` (Theorem 5: one squarefree level, under GRH; steps S1–S7), all in
`writeup/code/study/lowzero/universal/murmur/theorem/`. The unconditional route is carried out here for every squarefree
level (Section 3, steps 3 and 6, written out in full); the other steps are restated so that the note can be read alone,
and every slip listed in the three referee reports is repaired in the text (list in Section 5). Check script:
`one_level_check.py` (own code; output `one_level_check_out.txt`; **32 checks, 0 failed**).

**Labels.** **proved**: the argument is written out here. **quoted**: a statement taken from the literature with its
reference; *read* means that the source text was opened by the author of one of the three earlier notes or by their
referees, *read (secondary)* that the statement was read as restated in a later source, *not opened* that it is cited
from memory or from an earlier note; no source was reopened for this merge (no web). **checked**: a numbered assertion of
`one_level_check.py`. **NOT closed**: a gap. "GRH" always means GRH for Dirichlet L-functions; in the project it is
fact 2047, here it is treated as a hypothesis and used only in Theorem part (b).

## 1. The main theorem

**Setting.** N ≥ 2 squarefree, k ≥ 2 even, H_k(N) the newforms in S_k^new(Γ₀(N)), normalised by a_f(1) = 1,
λ_f(n) := a_f(n)n^(−(k−1)/2), ε_f the root number of L(s, f). For a prime p ∤ N put

    E_N(p) := |H_k(N)|^(−1) Σ_{f∈H_k(N)} ε_f √p λ_f(p)        (natural weights; for k = 2: (1/dim)Σ_f ε_f a_f(p)).

For y₀ ≥ 0, h > 0 let J := (y₀N, (y₀+h)N] and ⟨X⟩_J := Σ_{p∈J, p∤N} log p·X(p) / Σ_{p∈J, p∤N} log p.
The density M_k^(N) is defined in Section 2 (formula (2.1)).

**Theorem.** Fix k ≥ 2 even, Y ≥ 1, A ≥ 1.

(a) *(no hypothesis)* For every B ≥ 0, uniformly for squarefree N with H_k(N) ≠ ∅, 0 ≤ y₀ ≤ Y and (log N)^(−B) ≤ h ≤ 1,

    ⟨E_N⟩_J = (1/h) ∫_{y₀}^{y₀+h} M_k^(N)(y) dy + O_{k,Y,A,B}((log N)^(−A)).

The implied constant is ineffective (only through the Siegel–Walfisz theorem, step 6.5).

(b) *(under GRH for Dirichlet L-functions)* For every δ ∈ (0, ½) the same holds uniformly for N^(−1/2+δ) ≤ h ≤ 1, with
constant O_{k,Y,A,δ}. GRH is used for exactly two families: the real primitive characters χ_{Δ₀} of conductor
≤ 4(Y+1)N² (Δ₀ the fundamental discriminant of −N(4p − r²N)), and all characters to moduli ≤ 16N^(4δ/5).
No hypothesis on L(s, f) or L(s, Sym² f) is used in (a) or (b).

**Remarks (proved).** (1) H_k(N) ≠ ∅ is needed: S_2^new(N) = 0 for N = 6, 10, 22 (checked C4). By the dimension bound
(step 8) H_k(N) ≠ ∅ for every squarefree N ≥ N₀(k). (2) For N bounded the statement is trivial, so "uniformly in N" is a
statement about N → ∞ through squarefree integers with any number of prime factors. (3) Consequences, for N in the same
range of h: (i) ⟨E_N⟩_J = M_k^(N)(y₀) + O_{k,Y}(κ(N)h^(1/2)) + O((log N)^(−A)), since M_k^(N) is Hölder-½ on [0, Y+1] with
constant O_{k,Y}(κ(N)) (each r-term is 2√y·sin((k−1)θ_r(y)), cos θ_r = r/(2√y)); (ii) for Φ of bounded variation with
support in [0, Y+1] and ∫Φ > 0, Σ_{p∤N} log p Φ(p/N)E_N(p) / Σ_{p∤N} log p Φ(p/N) = ∫ΦM_k^(N)/∫Φ + O_Φ((log N)^(−A)) (the
key lemma below is stated for any weight of bounded variation, so no approximation by step functions is needed); (iii) with
weight 1 in place of log p and y₀ ≥ η > 0 fixed, the window mean changes by O_{k,Y,B,η}(κ(N)/log N) (apply (ii) with the
weight log N/log(vN) = 1 + O_η(1/log N) on the window). **Nothing is claimed for a single prime p.**

## 2. The density M_k^(N): local factors and three special cases

For every prime ℓ (ℓ = 2 included) and r ≥ 0 let u_ℓ := 1/((ℓ−1)²(ℓ+1)) and

    e_ℓ(r) := ℓ²/(ℓ² − 1) if ℓ | r (r = 0 included),   e_ℓ(r) := 1 − u_ℓ if ℓ ∤ r        (e₂(r) = 4/3 or 2/3),
    ν′(r) := Π_ℓ e_ℓ(r),   ν_N(r) := Π_{ℓ∤N} e_ℓ(r) = ν′(r)·Π_{q|N} e_q(r)^(−1)   (r ≥ 1),
    κ(N) := N/φ(N),   π(N) := Π_{q|N}(1 + 1/q) = κ(N)·Π_{q|N}(1 − q^(−2)),
    M_k^(N)(y) := (1/(k−1))·[2π π(N)√y + (−1)^(k/2−1)(12/π) κ(N) Σ_{1≤r<2√y} ν_N(r) U_{k−2}(r/(2√y)) √(4y − r²)]
                  − 12 κ(N) y·[k = 2].                                                                    (2.1)

U_{k−2} is the Chebyshev polynomial of the second kind; each r-term vanishes at y = r²/4, so M_k^(N) is continuous on
[0, ∞). ν′(1) = C₀ := Π_ℓ(1 − u_ℓ) = 0.615132657 (checked D1).

**Where the factors come from (proved in Section 3).** e_ℓ(r) is the mean over primes p of the local factor at ℓ of the
class number H(4pN − r²N²) when ℓ ∤ N (Lemma 4); at every prime q | N that mean is exactly 1 (Lemma 3), which removes
e_q(r) from ν′(r); κ(N) is the ratio N/(12·dim/(k−1)) of the normalisations; for r = 0 the Euler product is
ζ(2)Π_{q|N}(1 − q^(−2)) and (12/π)κ(N)ζ(2)Π(1 − q^(−2)) = 2π π(N) (checked D2). This local recipe (the factors e_ℓ(r)
and the factor 1 at q | N) is in print as a step inside a theorem averaged over levels (Cardi, Section 7); what is
proved here is the statement at one fixed level.

**How M_k^(N) depends on N (proved; repaired).** M_k^(N) depends on N only through the set P(N) of its prime factors.
Let S ⊂ P(N) contain every prime factor q ≤ 2√(Y+1) and put ε := Σ_{q|N, q∉S} 1/q. For q ∉ S one has q ≥ 3, q > r for every
r that occurs (r < 2√(Y+1)), −log(1 − 1/q) ≤ 1.5/q, log(1 + 1/q) ≤ 1/q and 0 ≤ −log(1 − u_q) ≤ 1/q. Hence
κ(N)/κ(S) ∈ [1, e^(1.5ε)], π(N)/π(S) ∈ [1, e^ε], ν_N(r)/ν_S(r) ∈ [1, e^ε] (checked D6 for all squarefree N ≤ 2·10⁵), and

    M_k^(N)(y) = M_k^(S)(y) + O_{k,Y}(κ(S)·ε·e^(3ε))   on [0, Y],

M_k^(S) being (2.1) for the level Π_{q∈S} q. ε is not bounded in general (for N = 5·7·11·…·127, ε = 1.03; checked D6),
so the error is e^(O(ε))-type; for ε ≤ 1 it is O(κ(S)ε). In particular N = mq with m fixed and q → ∞ gives M^(N) → M^(m)
uniformly on [0, Y]. Always 1 ≤ κ(N) ≪ log log N and C₀ ≤ ν_N(r) ≤ ζ(2), so M_k^(N) = O_{k,Y}(log log N).

**Special case 1: N = q prime (proved; checked D3).** ν_q(r) = ν′(r)/(1 − u_q) = ν′(r)(1 + O(q^(−3))) for q > r,
κ(q) = 1 + 1/(q−1), π(q) = 1 + 1/q, so M_k^(q) = M_k + O_Y(1/q) with

    M_k(y) := (1/(k−1))·[2π√y + (−1)^(k/2−1)(12/π) Σ_{1≤r<2√y} ν′(r) U_{k−2}(r/(2√y)) √(4y − r²)] − 12y·[k = 2],

the density of the prime-level theorems (for k = 2: M_2(y) = 2π√y − 12y + (12/π)Σ_r ν′(r)√(4y − r²)). Max |M^(q) − M_k| on
y ≤ 6 at q = 1000003: 3.4·10⁻⁶ (≤ 10/q) (k = 2, 4). Part (a) of the Theorem for N prime is Theorem U, part (b) is Theorem 2.

**Special case 2: N = 2q, q → ∞ (proved; checked D4).** e₂(r)^(−1) = 3/2 (r odd), 3/4 (r even), κ(2q) → 2, π(2q) → 3/2:

    M_k^(2q)(y) = (1/(k−1))·[3π√y + (−1)^(k/2−1)(24/π) Σ_r ν′(r)c₂(r) U_{k−2}(r/(2√y)) √(4y − r²)] − 24y·[k = 2] + O_Y(1/q),
    c₂(r) := 3/2 (r odd), 3/4 (r even).

The same recipe gives N = 3q (κ → 3/2, π → 4/3, c₃ = 16/15 for 3 ∤ r, 8/9 for 3 | r) and N = 6q (κ → 3, π → 2, c₂c₃).
Values for k = 2 (q → ∞; checked D4 printout):

| y | 0.1 | 0.2 | 0.5 | 1 | 2 |
|---|---|---|---|---|---|
| N = q (prime) | 0.7869 | 0.4099 | 0.7925 | −1.6471 | 0.5008 |
| N = 2q | 0.5804 | −0.5851 | 1.7132 | −2.3662 | −1.9239 |
| N = 3q | 0.8492 | 0.1466 | 0.6833 | −3.1109 | 0.8318 |
| N = 6q | 0.3738 | −1.5801 | 2.1640 | −3.8991 | −1.8326 |
| average over squarefree N (Special case 3) | 0.7829 | 0.1093 | 1.0077 | −2.3039 | 0.0040 |

A factor 2 in N changes the density by O(1) (y = 0.2: 0.41 for q against −0.59 for 2q and −1.58 for 6q). The term −12κy
and the line r = 1 (κν_N(1) = 0.615, 1.845, 0.984, 2.953 for q, 2q, 3q, 6q) carry most of the difference.

**Special case 3 = Proposition (the average over squarefree levels; proved; checked D5).** Fix k and y > 0. As X → ∞,

    Σ^□_{N≤X} φ(N) M_k^(N)(y) / Σ^□_{N≤X} φ(N) = M_zub,k(y) + O_{k,y}(X^(−1/2)),
    M_zub,k(y) := (β/(k−1))√y + (α(−1)^(k/2−1)/(k−1)) Σ_{1≤r<2√y} ν(r) U_{k−2}(r/(2√y)) √(4y − r²) − γy·[k = 2],

Σ^□ over squarefree N, with β = 2πΠ_ℓ (ℓ³+ℓ²−1)/(ℓ(ℓ²+ℓ−1)), γ = 12Π_ℓ ℓ(ℓ+1)/(ℓ²+ℓ−1), α = 2πΠ_ℓ (ℓ⁴−2ℓ²−ℓ+1)/(ℓ⁴−2ℓ²+ℓ),
ν(r) = Π_{ℓ|r}(1 + ℓ²/(ℓ⁴−2ℓ²−ℓ+1)): α = 3.965411, β = 7.862553, γ = 17.034753 (Euler products over ℓ < 10⁶; checked D5). These are the constants of Zubrilina's density (arXiv 2310.07681,
Thm 1, as restated with U_{k−2} inside the r-sum in arXiv 2609.29631, Thm 1.1; quoted, read in `level_theorem.md`; the
referee of that note matched them to her source to 7 digits). The same holds for N ∈ (X, X + X′] with X′ ≥ X^(1/2+ε),
error O(X^(1/2)/X′). Since |H_k(N)| = (k−1)φ(N)/12 + O_k(2^(ω(N))), the weight φ(N) is the dimension weight.
*Proof* in step 9. *Reading.* A prime ℓ divides a random squarefree N with probability 1/(ℓ+1); at ℓ | N the local factor
is 1, at ℓ ∤ N it is e_ℓ(r); her constants are, prime by prime, the mixtures (ℓe_ℓ(r) + 1)/(ℓ+1), normalised by the mean of
φ(N)/N. (This identity is in print in substance for α·ν(r) and β, as a local Fubini identity inside Cardi's averaged
theorem; γ = 12⟨N/φ(N)⟩ is a one-line addition: Section 7.) Logically this is an identity between densities; with the uniformity of the Theorem it gives the double average
over (N, p), not Zubrilina's theorem itself (she fixes p/N in a range and averages over N). It **suggests** why families of
levels with different small prime factors have different densities; it does not explain the footnote of 2310.07681, which
concerns non-squarefree levels.

## 3. Proof

**Plan.** Step 1 turns E_N(p) into class numbers H(4pN − r²N²) (trace formula); step 2 writes class numbers as L(1)-values
of quadratic characters; step 3 shows that at a squarefree level every prime q | N divides the fundamental discriminant
exactly and contributes the local factor 1; step 4 computes the other local means; step 5 states the key lemma (mean over
primes p of these L-values); step 6 proves it with no hypothesis (Burgess, Siegel–Walfisz, large sieve), step 7 under GRH;
step 8 assembles; step 9 proves the Proposition. Throughout D_r(p) := 4pN − r²N² = N·u_r(p), u_r(p) := 4p − r²N, for
0 ≤ r < 2√(Y+1) and primes p ∤ N with u_r(p) > 0; D_r(p) ≤ 4(Y+1)N² for p ≤ (Y+1)N.

### Step 1. The trace formula on the new space (quoted; old forms proved; checked C1–C3)

**Lemma 1.** For N squarefree, k ≥ 2 even and a prime p ∤ N,

    Σ_{f∈H_k(N)} ε_f √p λ_f(p) = ½H(4pN) + (−1)^(k/2−1) Σ_{1≤r<2√(p/N)} U_{k−2}(r√N/(2√p)) H(4pN − r²N²) − δ_{k,2}(p + 1),

H the Hurwitz class number. *Source (quoted):* Skoruppa–Zagier, Invent. Math. 94 (1988), §2 (5), (7) (not opened), the
trace of T_p∘W_N on S_k(Γ₀(N)); stated in this form in arXiv 2609.29631, Lemma (Skoruppa–Zagier trace formula) (read),
which notes that 2310.07681 §2 prints U_{k−2} outside the r-sum. r²N = 4p is impossible (p ∤ N, p prime).
*Root numbers (quoted):* ε_f = (−1)^(k/2) w_N(f), w_N = Π_{q|N} w_q the Atkin–Lehner eigenvalue, and a_f(q) = −q^(k/2−1)w_q(f)
for q | N (Atkin–Lehner 1970, not opened); so ε_f a_f(p) = (−1)^(k/2+ω(N)) a_f(Np)/N^(k/2−1), which is how the check
compares with PARI's trace forms. *Old forms contribute 0 (proved):* for a newform g of level M | N, M < N, the space
V_g = span{g(dz) : d | N/M} is T_p-stable (scalar a_g(p)), and W_N maps g(dz) to a nonzero multiple of g(((N/M)/d)z)
(Atkin–Lehner 1970, Lemma 26, not opened). N/M is squarefree and > 1, so d ≠ (N/M)/d: W_N permutes 2^(ω(N/M)) lines with no
fixed line and Tr(T_pW_N | V_g) = 0. Hence the trace on S_k(Γ₀(N)) equals the trace on the new space.
**Checked:** C1 exact rational equality with PARI trace forms of the new space at N = 23, 46, 69, 78 (prime, 2q, 3q, 6q),
k = 2 (all p ≤ 500), k = 4 (p ≤ 250), k = 6 (p ≤ 120), p = 2 included for odd N: 372, 204 and 112 pairs (N, p); C2 old forms give 0 (full
space = new space, 95 traces); C3 root numbers from the functional equation (lfunrootres) = (−1)^(k/2)w_N and
a_f(q) = −q^(k/2−1)w_q (11 rational newforms). Earlier: 1,386 pairs (level_theorem C1–C3) and 1,644 triples with
weights 2 to 12 and up to three prime factors (its referee).

### Step 2. Class numbers as L-values (quoted classical; checked A1–A3)

For D > 0 with D ≡ 0, 3 (mod 4) write −D = Δ₀f₀² with Δ₀ a fundamental discriminant, χ := (Δ₀/·), and let
𝓜(D) := {m ≥ 1 : m² | D, −D/m² ≡ 0, 1 (mod 4)}. For a discriminant Δ < 0 let ξ_Δ := (Δ/·) (Kronecker symbol), a character
mod |Δ| induced by χ_{Δ₀}: ξ_Δ(n) = χ_{Δ₀}(n)·1[(n, f) = 1] for Δ = Δ₀f².

**Lemma 2.** (i) H(D) = Σ_{m∈𝓜(D)} h′(−D/m²), h′(Δ) := 2h(Δ)/w(Δ) (H counts all classes of forms of discriminant −D,
each primitive class of discriminant −D/m² once, weighted 2/w). (ii) h′(Δ) = √|Δ|·L(1, ξ_Δ)/π for every discriminant Δ < 0
(class number formula for orders; Cox, *Primes of the form x² + ny²*, Thm 7.24, Cor. 7.28: quoted, not opened). Hence

    H(D) = (√D/π)·L*(D),   L*(D) := Σ_{m∈𝓜(D)} m^(−1) L(1, ξ_{−D/m²})   (a finite sum).

(iii) *(proved)* L*(D) = L(1, χ)·Π_{ℓ^a∥f₀} Q_ℓ with Q_ℓ := (1 − χ(ℓ)/ℓ)Σ_{j<a}ℓ^(−j) + ℓ^(−a), and
**1 ≤ Q_ℓ < (ℓ+1)/(ℓ−1)** (Q_ℓ = 1 if χ(ℓ) = 1; (ℓ − ℓ^(−a))/(ℓ−1) if χ(ℓ) = 0; (ℓ + 1 − 2ℓ^(−a))/(ℓ−1) if χ(ℓ) = −1; for
example ℓ = 2, χ(2) = −1, a = 2 gives 5/2, so the range [1 − 1/ℓ, 1 + 2/ℓ] stated in an earlier draft was wrong). The local
factor F_ℓ(D) := Q_ℓ/(1 − χ(ℓ)/ℓ) = Σ_{j<a}ℓ^(−j) + ℓ^(−a)/(1 − χ(ℓ)/ℓ) satisfies 0 < F_ℓ(D) ≤ ℓ/(ℓ−1)
(F_ℓ = (1 − χ(ℓ)/ℓ)^(−1) for ℓ ∤ f₀, i.e. a = 0).
(iv) *(proved)* Expanding each L(1, ξ) as a Dirichlet series, L*(D) = Σ_m m^(−1) Σ_n T(m, n)/n with
T(m, n) := 1[m ∈ 𝓜(D)]·ξ_{−D/m²}(n), and pointwise Σ_{i,j≥0} T_ℓ(i, j; D)ℓ^(−i−j) = F_ℓ(D), where T = Π_ℓ T_ℓ(v_ℓ(m), v_ℓ(n); D)
with T_ℓ(i, j; D) := 1[ℓ^(2i) | D, and for ℓ = 2: −D/4^i ≡ 0, 1 (4)]·(−Dℓ^(−2i)/ℓ)^j depending on D only ℓ-adically (since
(ab²/ℓ^j) = (a/ℓ^j) for ℓ ∤ b, also for ℓ = 2). (Proof of the pointwise identity: −D/ℓ^(2i) is a discriminant iff i ≤ a, and
(−D/ℓ^(2i)/ℓ) is 0 for i < a and χ(ℓ) for i = a.)
**Checked:** A1 own count of reduced forms = PARI's H(D), D ≤ 800; A2 (i) exactly for all 1,500 D ≤ 3000 with own counts of
primitive forms; A3 (ii) by L(1, ξ) = −|Δ|^(−1)Σ_a ξ(a)ψ(a/|Δ|) (ψ = digamma) for all 947 Δ = −D/m², D ≤ 1200: max |difference| 1.4·10⁻¹⁴.

### Step 3. The discriminants of a squarefree level (proved; checked A4, A5)

**Lemma 3.** Let N be squarefree, p ∤ N prime, r ≥ 0 with u := 4p − r²N > 0, D := Nu, −D = Δ₀f₀².
(i) For every prime q | N: q | Δ₀ and q ∤ f₀; hence χ(q) = 0, a = 0 and **F_q(D) = 1**. In particular N | Δ₀ and
|Δ₀| ≥ N. For even N, v₂(Δ₀) ∈ {2, 3}.
(ii) Every m ∈ 𝓜(D) is prime to N.
(iii) For m ∈ 𝓜(D) and n with (n, N) > 1: ξ_{−D/m²}(n) = 0. So T(m, n) = 0 unless (mn, N) = 1.
*Proof.* (i) q odd: u ≡ 4p ≢ 0 (mod q), so v_q(D) = 1, which is odd: q divides the squarefree kernel, hence Δ₀, and q ∤ f₀.
q = 2, N = 2N′ (N′ odd): −D = 4N′(r²N′ − 2p). If r is odd, w := N′(r²N′ − 2p) is odd and w ≡ N′²r² − 2pN′ ≡ 1 − 2 ≡ 3 (mod 4);
with w = w₀g² (w₀ squarefree, g odd) w₀ ≡ 3 (mod 4), so Δ₀ = 4w₀, v₂(Δ₀) = 2, 2 ∤ f₀. If r = 2r′ (r = 0 included),
−D = 8N′(2r′²N′ − p) = 8·odd, so Δ₀ = 8·(odd squarefree), v₂(Δ₀) = 3, 2 ∤ f₀. In all cases q | Δ₀ gives χ(q) = 0, and a = 0, so
F_q = Σ_{j<0} + q⁰/(1 − 0) = 1. (ii) m² | D and q | m would need q² | D: impossible for odd q (v_q(D) = 1); for q = 2,
v₂(D) ∈ {2, 3}, so 4 ∤ m, and for 2 ∥ m the number −D/m² ≡ −D/4·(odd square)^(−1) is ≡ 3 or 2 (mod 4), not a
discriminant; so no even m lies in 𝓜(D). (iii) q | (n, N) gives q | Δ₀ | −D/m²
(by (i), (ii)), so the Kronecker symbol vanishes at q. ∎
*Consequence.* The conductor |Δ₀| of χ lies in [N, 4(Y+1)N²]. Only the upper bound is used (Burgess, step 6.2); the lower
bound is recorded because it means that none of the characters χ_{Δ₀} has small conductor. **Checked:** A4 on 28,671 cases
(all squarefree N ≤ 130 and ten levels with 3 to 6 prime factors up to 39270; all p < 400, p ∤ N; all r): N | Δ₀, q ∤ f₀,
F_q = 1 at q | N, 𝓜(D) prime to N, v₂(Δ₀) ∈ {2, 3} (both occur); A5 T(m, n; p) = 0 for (n, N) > 1 (12,406 cases).

### Step 4. Local means (proved; checked B1)

For a prime ℓ let E_ℓ^(N)(r) := ∫_{ℤ_ℓ^×} F_ℓ(D_r(p)) dp (Haar probability on the units; p ∤ N is the only condition, and
p prime ∤ ℓ is an ℓ-adic unit).

**Lemma 4.** E_ℓ^(N)(r) = 1 if ℓ | N, and E_ℓ^(N)(r) = e_ℓ(r) if ℓ ∤ N. Hence
E^(N)(r) := Π_ℓ E_ℓ^(N)(r) = ν_N(r) (r ≥ 1) and E^(N)(0) = ζ(2)Π_{q|N}(1 − q^(−2)).
*Proof.* ℓ | N: Lemma 3(i). ℓ ∤ N uses only that N is an ℓ-adic unit. Put x := 1/ℓ, S_a := (1 − x^a)/(1 − x).
- ℓ odd, ℓ | r (r = 0 included): −D ≡ −4pN (mod ℓ) is a unit; χ(ℓ) = (−pN/ℓ) is ±1 with probability ½ each and a = 0:
  E = ½[(1 − x)^(−1) + (1 + x)^(−1)] = ℓ²/(ℓ² − 1).
- ℓ odd, ℓ ∤ r: −D = N·u′ with u′ := r²N − 4p Haar-uniform on ℤ_ℓ ∖ (r²N + ℓℤ_ℓ) (density ℓ/(ℓ−1)); multiplying by the unit N
  maps Haar measure on ℤ_ℓ to itself. Over all of ℤ_ℓ (P(v = 2a) = (1−x)x^(2a), χ = ±1 equally; P(v = 2a+1) = (1−x)x^(2a+1),
  χ = 0) the mean of F_ℓ is Σ_a (1−x)x^(2a)[S_a + x^a/(1−x²)] + (1−x)x^(2a+1)[S_a + x^a] = 1/(1 − x²). On the removed disc −D is
  a unit square (≡ N²r²), F_ℓ = (1 − x)^(−1), measure x. So E = (ℓ/(ℓ−1))[ℓ²/(ℓ²−1) − 1/(ℓ−1)] = ℓ(ℓ²−ℓ−1)/((ℓ−1)²(ℓ+1)) = 1 − u_ℓ.
- ℓ = 2 ∤ N, r odd: −D ≡ N²r² − 4pN ≡ 1 − 4 ≡ 5 (mod 8), fundamental at 2 with χ(2) = −1: F₂ = 2/3 = 1 − u₂.
  r = 2r′: −D = 4Δ′, Δ′ := r′²N² − pN. r′ even: Δ′ is uniform on ℤ₂^×; F₂ = 1 (Δ′ ≡ 3 mod 4), 2 (Δ′ ≡ 1 mod 8), 4/3
  (Δ′ ≡ 5 mod 8): mean 4/3. r′ odd: Δ′ = 2w with w uniform on ℤ₂, −D = 8w; F₂(8·odd) = 1, and F₂(4Δ₁) = 1 + F₂(Δ₁)/2 for a
  discriminant Δ₁, so the mean is ½ + ½(1 + I/2) with I := mean of F₂(4w′) over w′ uniform on ℤ₂, I = 7/6 + I/8, I = 4/3:
  mean 4/3 = 2²/(2² − 1). ∎
**Checked B1:** exact enumeration of p mod ℓ^J with a rigorous enclosure of the undetermined residues, ℓ = 2, 3, 5, 7
(J = 18, 11, 7, 6), levels 19, 23, 9973, 10009 (prime, all classes mod 8), 26, 9998 (2q), 21, 9987 (3q), 42, 9978 (6q),
35, 2310 (5 | N), r = 0…10: 528 cases, every prediction inside its enclosure, width ≤ 1.9·10⁻⁶. Earlier: 676 cases (level_theorem B1), 396 cases with other code (its referee).

### Step 5. The key lemma

**Key lemma.** Let 0 ≤ r < 2√(Y+1) and G of bounded variation on [0, Y+1], ‖G‖ := sup|G| + Var G, with G(v) = 0 for
v ≤ r²/4 when r ≥ 1. Uniformly in squarefree N, 0 ≤ y₀ ≤ Y and G,

    Σ_{p∈J, p∤N} log p·G(p/N)·L*(D_r(p)) = E^(N)(r)·∫_J G(u/N) du + O(‖G‖·hN·(log N)^(−A)),

(a) with no hypothesis for (log N)^(−B) ≤ h ≤ 1 (constant O_{Y,A,B}, ineffective); (b) under GRH for N^(−1/2+δ) ≤ h ≤ 1.

### Step 6. Proof of the key lemma without hypothesis (U2–U6 for squarefree N; proved)

Parameters: M := N^(7/8), m₀ := (log N)^(A+3), L := (log N)^(2A+6), Q₀ := (log N)^(B₁) with B₁ = B₁(A) fixed in 6.6.
All errors below are O(‖G‖hN(log N)^(−A)) and none of the constants depends on N or on ω(N); N enters only through log N,
through ω(N) ≤ 2 log N and through τ(D) ≤ N^(o(1)).

**6.0 Two elementary facts.** (a) Σ_{p∈J} log p ≪ hN (PNT with de la Vallée Poussin's error, Step 8, since
h ≥ (log N)^(−B); under GRH in case (b)) and Σ_{p∈J} log²p ≤ (hN + 1)log²((Y+1)N) (trivial). (b) The primes p | N, at most ω(N) ≤ 2 log N of them, may be added to or removed from any sum below at cost
O(‖G‖·log⁴N·log log N) for the sum with L*(D_r(p)) (|L*(D)| ≤ Σ_{m≤√D}m^(−1)·max|L(1, ξ)| ≪ log²N log log N by 6.1), and
O(‖G‖log²N) for a sum with a single term T(m, n; p) (|T| ≤ 1); both are O(‖G‖hN(log N)^(−A)).

**6.1 Large m (U3; proved).** |L(1, ξ_Δ)| ≤ L(1, χ)Π_{ℓ|f}(1 + 1/ℓ) ≤ (log|Δ₀| + 2)·σ(f)/f ≪ log N·log log N
(L(1, χ) ≤ log q + 2 by partial summation: |Σ_{n≤q}χ(n)/n| ≤ 1 + log q and |Σ_{n>q}χ(n)/n| ≤ q·(1/q)). The number of p ∈ J
with m² | 4p − r²N is ≤ 4hN/m² + 1 (one class mod m²/4 at most). Hence the terms m > m₀ of L*(D_r(p)), summed over p with
log p·|G|, are ≪ ‖G‖ log²N log log N·(hN/m₀² + log N) = O(‖G‖hN(log N)^(−A)).

**6.2 Truncation at M for every p (U2; proved from Burgess; repaired).** *Quoted (Burgess, r = 2, any modulus):* for χ
primitive mod q > 1, ε > 0 and every interval, |Σ_{a<n≤a+H} χ(n)| ≤ C(ε)H^(1/2)q^(3/16+ε) (Burgess, Proc. LMS 13 (1963), not
opened; read (secondary) in Treviño, arXiv 2511.17778, §1). *Partial summation (proved):* for Z > 0,
|Σ_{n>Z} χ(n)/n| = |∫_{max(Z,1)}^∞ S_Z(t)t^(−2)dt| ≤ 2C(ε) max(Z, 1)^(−1/2) q^(3/16+ε), S_Z(t) := Σ_{Z<n≤t}χ(n).
For Δ = Δ₀f² < 0, Möbius gives Σ_{n>M} ξ_Δ(n)/n = Σ_{d|f} μ(d)χ(d)d^(−1) Σ_{n′>M/d} χ(n′)/n′, and
d^(−1)max(M/d, 1)^(−1/2) ≤ d^(−1/2)M^(−1/2) for every d ≥ 1. So

    |L(1, ξ_Δ) − Σ_{n≤M} ξ_Δ(n)/n| ≤ 2C(ε) τ(f) M^(−1/2) |Δ₀|^(3/16+ε) ≪_{Y,ε} N^(−7/16+3/8+3ε) = N^(−1/16+3ε),

using |Δ₀| ≤ 4(Y+1)N² and τ(f) ≤ N^ε. (Earlier draft: N^(−1/32+o(1)) by an unnecessary split at d = N^(1/32); Burgess applies
to every divisor.) χ is never principal (Δ₀ < 0). This holds for **every** p and every m ≤ m₀; summed with 1/m and with
log p·|G| it costs O(‖G‖hN N^(−1/16+4ε)). Pólya–Vinogradov would need M > N, beyond hN, where the large sieve below fails;
Burgess avoids conductors comparable to the length. **Checked E3** at the composite levels 9998 and 9978: rms per-p truncation error at M = N^(1/2), N^(3/4), N^(7/8): 0.114, 0.039, 0.022 (r = 0) and 0.059, 0.019, 0.010 (r = 1) at N = 9998; 0.077, 0.025, 0.015 and 0.060, 0.019, 0.012 at N = 9978; the mean over the ~1030 primes in (N, 2N] of the series cut at N^(7/8) is 1.2253 against E^(N)(0) = 1.2337 (z −0.48) and 0.9143 against ν_N(1) = 0.9227 (z −0.79) at 9998, 1.0944 against 1.0966 (z −0.22) and 0.9838 against 0.9842 (z −0.05) at 9978.

**6.3 The terms as characters of the shifted prime (U4; proved; checked A5, A6).** Write
L*(D_r(p)) = Σ_{m≤m₀} m^(−1) Σ_{n≤M} T(m, n; p)/n + (errors of 6.1, 6.2), T(m, n; p) := T(m, n) of Step 2 at D = D_r(p).
(a) *Periodicity.* (Δ/n) depends only on Δ mod 4n ((Δ/ℓ) on Δ mod ℓ for odd ℓ, (Δ/2) on Δ mod 8). 1[m ∈ 𝓜(D)] depends on
Nu mod 4m², and Δ_m := −Nu/m² mod 4n on Nu mod 4nm², i.e. on 4p mod 4nm². So **T(m, n; ·) is periodic mod nm²**
(not only mod 2nm² as stated earlier; checked A5, 39,984 tests). By Lemma 3(ii)(iii), T ≡ 0 unless (mn, N) = 1, so
**every period nm² that occurs is prime to N**.
(b) *Factorisation, r ≥ 1.* Let q := Π{ℓ odd : v_ℓ(n) odd, ℓ ∤ mrN} (odd squarefree), q♮ := Π_{ℓ|q} ℓ^(v_ℓ(n)) = q·s² (s | q^∞),
n″ := n/q♮. Kronecker symbols are completely multiplicative in the lower entry, (Δ/ℓ^(2b+1)) = (Δ/ℓ), and Δ_m m² = −Nu with
(m, q) = (N, q) = 1, so (Δ_m/q) = (−N/q)(u/q):

    T(m, n; p) = W(p)·χ_q(4p − r²N),   χ_q := (·/q) real primitive mod q,   W(p) := (−N/q)·T(m, n″; p),

W periodic mod K := n″m², (K, q) = 1, |W| ≤ 1; (−N/q) = ±1 because (q, N) = 1. *r = 0:* q := Π{ℓ odd : v_ℓ(n) odd, ℓ ∤ mN};
then χ_q(4p) = χ_q(p) and T = (−N/q)χ_q(p)T(m, n″; p). For squarefree N the only change against the prime case is that all
primes of N are kept out of q and, by (a), out of n″ (T = 0 otherwise). Checked A5 (39,984 random cases, squarefree N < 400
and 2310, 4290, 39270, 9978, 9998).
(c) *Expansion in characters of p (quoted classical + proved; checked A6).* For p ∤ q, χ_q(4p − c) = Σ_{ψ mod q} â(ψ)ψ(p),
c := r²N, â multiplicative over ℓ | q (CRT). For ℓ ∤ c: â_ℓ(ψ₀) = −(−c/ℓ)/(ℓ−1) (complete sum of a Legendre symbol over a
line is 0, minus the missing class x ≡ 0); for ψ ≠ ψ₀, (ℓ−1)â_ℓ(ψ) is a Jacobi sum of the Legendre symbol and ψ̄, of modulus
√ℓ, or 1 when ψ is the Legendre symbol (|J(χ₁, χ₂)| = √ℓ for χ₁χ₂ ≠ ψ₀: quoted classical). So **|â(ψ₀)| = 1/φ(q),
|â(ψ)| ≤ √q/φ(q)**. For r = 0, â = 1[ψ = χ_q] exactly, for every p.

**6.4 The per-p cut (U5; proved).** Write n = 2^v n₁ n₂² with n₁ odd squarefree, n₂ odd. Let 𝒩 := {n ≤ M : 2^v ≤ L,
n₂ ≤ L}. Since |T| ≤ 1, dropping n ∉ 𝒩 costs per p at most (1 + log m₀)(Σ_{2^v>L}2^(−v) + Σ_{n₂>L}n₂^(−2))·2(1 + log M)
≪ log log N·log N/L, in total O(‖G‖hN(log N)^(−2A−4)). For n ∈ 𝒩 with T ≢ 0, n″ consists of 2^v ≤ L, of the full
ℓ-parts of n for odd ℓ | mr (ℓ | N cannot occur, Lemma 3(iii)), each ℓ^(v_ℓ) ≤ ℓ·(ℓ^(⌊v_ℓ/2⌋))², and of even-exponent parts;
so n″ ≤ L·rad(mr)·n₂² ≤ L³m₀·2√(Y+1) (rad(m) in place of rad(mr) for r = 0), and K = n″m² ≤ L³m₀³·2√(Y+1) =: (log N)^(C₃),
C₃ = C₃(A).
The cut depends on n only (not on p), so it can be made before averaging; this is what keeps the moduli of 6.5
polylogarithmic.

**6.5 Small conductors q ≤ Q₀ (U5(ii); proved from Siegel–Walfisz).** *Quoted:* for every B′ there is c_{B′} > 0
(ineffective) with θ(x; P, b) = x/φ(P) + O_{B′}(x exp(−c_{B′}√log x)) uniformly for P ≤ (log x)^{B′}, (b, P) = 1 (Siegel 1935,
Walfisz 1936: not opened; read (secondary) in the earlier note). Here n = q♮n″ ≤ Q₀L²·n″ (s ≤ n₂ ≤ L), so T(m, n; ·) has period
P = nm² ≤ (log N)^(B₂), B₂ = B₂(A), and (P, N) = 1. Splitting p into classes b mod P and integrating by parts against
θ(x; P, b) − x/φ(P) (trivial bound 2x for x ≤ √N): Σ_{p∈J, p∤N} log p G(p/N)T(m, n; p) = e(m, n)∫_J G(u/N)du +
O(P‖G‖N exp(−c√log N)), e(m, n) := φ(P)^(−1)Σ_{b∈(ℤ/P)^×}T(m, n; b) (6.0(b) for p | N). There are (log N)^(O(1)) such (m, n),
each with weight 1/(mn) ≤ 1; the total error is O(‖G‖hN(log N)^(−A)) because h ≥ (log N)^(−B). **This is the only step that
needs h ≥ (log N)^(−B), and the only ineffective one.**

**6.6 Conductors Q₀ < q ≤ M (U5(iii); proved from the large sieve).** *Quoted:* Σ_{d≤R}(d/φ(d))Σ*_{ψ mod d}|Σ_{n∈I}a_nψ(n)|²
≤ (|I| + R²)Σ|a_n|² for an interval I (Montgomery–Vaughan 1973, Selberg; read as (1.1) in Friedlander–Iwaniec, arXiv
2409.19634, in the earlier note). Fix (m, n″) (outer weight 1/(mn″)), r ≥ 1, and a dyadic block q ∈ (Q, 2Q]. The sum over s
costs at most ζ(2) (T does not depend on s). With a_p := log p·G(p/N)·1[p ∤ N] and S_b(ψ) := Σ_{p∈J, p≡b (K)} a_pψ(p):

    Σ_{q∼Q} q^(−1)|Σ_p a_p W(p) χ_q(4p − c)| ≤ Σ_{b∈(ℤ/K)^×} [ Σ_{q∼Q} |S_b(ψ₀)|/(qφ(q)) + Σ_{1<d≤2Q} w_d Σ*_{ψ mod d} |S_b(ψ)| ],

where each ψ ≠ ψ₀ mod q is replaced by the primitive ψ* mod d | q inducing it, and the primes p | K (at most log K ≪ log log N
of them, all ≤ (log N)^(C₃)) are removed first at trivial cost. **No correction for p | q is needed** (repair):
for r ≥ 1, G(p/N) = 0 for p ≤ r²N/4, and q ≤ M < N/4, so p ∤ q on the support; for r = 0 the expansion is exact. Here
w_d := Σ_{q∼Q, d|q} q^(−1/2)/φ(q) ≤ c₀ Q^(−1/2)/φ(d) with **c₀ := 2ζ(2)ζ(3)/ζ(6) < 3.9** (proved: q = dq′ and q/φ(q) ≤
(d/φ(d))(q′/φ(q′)) give q^(−1/2)/φ(q) ≤ (d/φ(d))(q′/φ(q′))Q^(−3/2); Σ_{q′≤X} q′/φ(q′) = Σ_e (μ²(e)/φ(e))⌊X/e⌋ ≤
(ζ(2)ζ(3)/ζ(6))X; checked A7: max φ(d)Q^(1/2)w_d = 1.41 for Q ≤ 2^14). Cauchy–Schwarz over (d, ψ) (at most φ(d) primitive ψ
mod d, Σ_{d≤2Q}1/d ≪ log Q) and the large sieve with R = 2Q, |I| ≤ hN + 1:
Σ_d w_d Σ*|S_b| ≪ Q^(−1/2)(log Q)^(1/2)((hN + 4Q²)Σ_{p≡b}a_p²)^(1/2); Cauchy–Schwarz over b gives φ(K)^(1/2); 6.0(a) gives
Σa_p² ≪ ‖G‖²hN log²N. The principal part is ≤ Σ_{q∼Q} 1/(qφ(q))·Σ_p|a_p| ≪ ‖G‖hN/Q. Block total:

    ≪ ‖G‖·[hN/Q + φ(K)^(1/2) log^(3/2)N·(hN Q^(−1/2) + (hN)^(1/2)Q^(1/2))].

For r = 0, â = 1[ψ = χ_q] and Cauchy–Schwarz over q (Σ_{q∼Q}q^(−2) ≤ Q^(−1)) gives the same bound. Summing over dyadic
Q₀ ≤ Q ≤ M and over (m, n″) (Σ 1/(mn″) ≪ (log log N)², φ(K)^(1/2) ≤ (log N)^(C₃/2)):

    ≪ ‖G‖hN (log N)^(C₃/2+2)·[Q₀^(−1/2) + (M/(hN))^(1/2)],   (M/(hN))^(1/2) ≤ N^(−1/16)(log N)^(B/2).

Choose **B₁ := 2A + C₃ + 4**: the piece is O(‖G‖hN(log N)^(−A)). (The sequence a_p does not depend on q; W depends on q
only through the sign (−N/q), which is pulled out in absolute value.) q > M does not occur (n ≤ M).

**6.7 The main term (U6; proved; checked B2).** By Step 2(iv) and CRT, e(m, n) = Π_ℓ e_ℓ(v_ℓ(m), v_ℓ(n)) with
e_ℓ(i, j) := ∫_{ℤ_ℓ^×} T_ℓ(i, j; D_r(p)) dp (T_ℓ(i, j; ·) depends on p mod ℓ^(2i+j) for ℓ odd and mod 2^(2i+j) for j ≥ 1,
which divides nm²). Since |T_ℓ| ≤ 1 and Σℓ^(−i−j) < ∞, dominated convergence and Step 2(iv) give
**Σ_{i,j} e_ℓ(i, j)ℓ^(−i−j) = E_ℓ^(N)(r)** (Lemma 4). Explicitly:
- ℓ | N: **e_ℓ(i, j) = [i = j = 0]** (Lemma 3: ℓ² ∤ D, (Δ/ℓ) = 0); the local factor is exactly 1.
- ℓ ∤ 2rN: e_ℓ(0, 0) = 1, e_ℓ(0, j) = −1/(ℓ−1) (j odd), (ℓ−2)/(ℓ−1) (j even ≥ 2), e_ℓ(i, 0) = ℓ^(1−2i)/(ℓ−1) (i ≥ 1),
  |e_ℓ(i, j)| ≤ e_ℓ(i, 0). For r ≥ 1, e_ℓ(0, odd) ≠ 0: non-principal n carry part of the main term (the principal-n sum
  alone is ≈ 1.64 for every r, against E(1) = 0.615: checked in the earlier note, B3–B4).
- ℓ | r, ℓ ∤ 2N: e_ℓ(0, j) = [j even], e_ℓ(i, j) = 0 for i ≥ 1.
- For every ℓ and i ≥ 1: |e_ℓ(i, j)| ≤ P(ℓ^(2i) | D_r(p)) ≤ 8ℓ^(1−2i) (one class of p mod ℓ^(2i) for odd ℓ ∤ N, none for
  ℓ | N; for ℓ = 2 ∤ N one class mod 4^(i−1); for 2 | N, v₂(D) ≤ 3).
*Weighted absolute convergence (proved).* With the weight ω(m, n) := m·n₂^(1/2)·2^(v₂(n)/2)·q(n)^(1/2) (q(n) the q of 6.3),
Σ_{m,n} |e(m, n)|ω(m, n)/(mn) = Π_ℓ(Σ_{i,j}|e_ℓ(i, j)|ℓ^(−i−j)ω_ℓ(i, j)) ≤ C(Y) uniformly in N: the local factor is exactly 1
for ℓ | N, 1 + O(ℓ^(−3/2)) for ℓ ∤ 2rN (term (0,1): ℓ^(−1/2)/(ℓ−1); (0,2): < ℓ^(−3/2); (1, j): O(ℓ^(−2)); the weight q(n)^(1/2) enters only for i = 0, j odd), and bounded for the
finitely many ℓ | 2r. Rankin's trick then bounds the left-out pairs: m > m₀ by C(Y)/m₀; 2^v > L or n₂ > L by C(Y)L^(−1/2);
q(n) > Q₀ by C(Y)Q₀^(−1/2); and n > M with 2^v ≤ L, n₂ ≤ L forces q(n) ≥ n/(L³·rad(mr)) > M/(L³m₀2√(Y+1)) > Q₀, so these
pairs are among the previous ones. All are
O((log N)^(−A)), so

    Σ_{m≤m₀, n∈𝒩, q(n)≤Q₀} e(m, n)/(mn) = E^(N)(r) + O((log N)^(−A)).

**Checked B2** (exact e_ℓ(i, j), ℓ = 3, 5, 7, levels 19, 9973, 26, 21, 42, 35, 2310, r = 0…6): 147 cases; every closed form above exact, every local sum within its rigorous tail of E_ℓ^(N)(r).

**6.8 Conclusion.** 6.1 + 6.2 + 6.4 reduce the left side of the key lemma, up to admissible errors, to
Σ_{m≤m₀}Σ_{n∈𝒩} (mn)^(−1)Σ_{p∈J, p∤N}log p G(p/N)T(m, n; p); 6.5 evaluates the pairs with q(n) ≤ Q₀ as e(m, n)∫_J G; 6.6 bounds
those with q(n) > Q₀; 6.7 sums the main terms to E^(N)(r)∫_J G(u/N)du. ∎ (Key lemma (a).)

**What changed for squarefree N, and that it closes (summary of steps 3 and 6).** (1) The conductor: N | Δ₀ including at 2
(v₂(Δ₀) ∈ {2, 3} for even N), so |Δ₀| ∈ [N, 4(Y+1)N²]; Burgess needs only the upper bound. (2) 𝓜(D) is prime to N and
T(m, n; p) = 0 when (n, N) > 1: every modulus P = nm² and K = n″m² is prime to N, and the primes of N stay out of q and n″.
(3) (−N/q) is still a sign. (4) The local factor at every q | N is exactly 1 (e_q(i, j) = [i = j = 0]), so the main term is
ν_N(r), resp. ζ(2)Π_{q|N}(1 − q^(−2)). (5) p ∤ N removes at most 2 log N primes. (6) No constant depends on ω(N); the factor
1 at q | N only helps the absolute convergence. **Nothing fails: part (a) holds for every squarefree N.**

### Step 7. Proof of the key lemma under GRH (T3, T4, T6 for squarefree N; proved given GRH)

Parameters: L := N^(δ/5), log z := log L/((A+1) log log N).
**7.1 (T3).** *Quoted consequence of GRH:* θ(x, χ) := Σ_{ℓ≤x}χ(ℓ)log ℓ ≪ √x log²(qx) for χ ≠ χ₀ mod q, and
θ(x; q, b) = x/φ(q) + O(√x log²(qx)) for (b, q) = 1 (Davenport, *Multiplicative Number Theory*, ch. 20, not opened). Then
Σ_ℓ χ(ℓ)ℓ^(−s) converges on Re s > ½, Σ_{ℓ>z}χ(ℓ)/ℓ ≪ z^(−1/2)log²(qz), and
L(1, χ) = Π_{ℓ≤z}(1 − χ(ℓ)/ℓ)^(−1)exp(O(z^(−1/2)log²(qz) + 1/z)). With ω(f₀) ≤ log D and Q_ℓ = 1 + O(1/ℓ) (Lemma 2(iii)),
Π_{ℓ|f₀, ℓ>z}Q_ℓ = exp(O(log D/z)). With P_z(D) := Π_{ℓ≤z}F_ℓ(D) and **P_z(D) ≪ log z** (0 < F_ℓ ≤ (1 − 1/ℓ)^(−1), Mertens;
the earlier (log z)³ was true but loose), uniformly for D ≤ 4(Y+1)N²:
L*(D) = P_z(D) + O(log z·[z^(−1/2)log²(Nz) + log N/z]).
**7.2 (T4).** F_ℓ(D) = Σ_j c_{ℓ,j}(D)ℓ^(−j), c_{ℓ,j} = 1 (j < a), χ(ℓ)^(j−a) (j ≥ a); |c| ≤ 1; c_{ℓ,j}(D) depends only on D mod
ℓ^(2j+2) (ℓ odd), mod 2^(2j+6) (ℓ = 2). So P_z(D) = Σ_{P⁺(n)≤z} β_D(n)/n with β_D(n) := Π_{ℓ^j∥n}c_{ℓ,j}(D), depending on D mod
Q(n) ≤ 16n⁴. Rankin (σ = 1/log z): Σ_{n>L, P⁺(n)≤z}1/n ≪ log z·exp(−log L/log z). For ℓ | N, c_{ℓ,j} = [j = 0].
**7.3 (T6).** For n ≤ L the class of p mod Q(n) fixes D_r(p) mod Q(n), so
Σ_{p∈J, p∤N}log p G(p/N)β_{D_r(p)}(n) = e_r(n)∫_J G(u/N)du + O(‖G‖Q(n)√N log²N + (L + log²N)sup|G|),
e_r(n) := φ(Q)^(−1)Σ_{b∈(ℤ/Q)^×}β_{D_r(b)}(n) = Π_{ℓ^j∥n}∫c_{ℓ,j}(D_r(p))dp (CRT), and Σ_{P⁺(n)≤z}e_r(n)/n = Π_{ℓ≤z}E_ℓ^(N)(r).
**7.4 Assembly of (b).** 7.1 for each p (error ≪ exp(−c log N/(A log log N)) per p); 7.2 keeps n ≤ L (Rankin
≪ (log N)^(−A)); 7.3 for each n ≤ L (total ≪ ‖G‖L⁴√N log³N ≤ ‖G‖hN N^(−δ/5)log³N since hN ≥ N^(1/2+δ)); the main terms
give Π_{ℓ≤z}E_ℓ^(N)(r) = E^(N)(r)(1 + O(1/z)) since E_ℓ^(N)(r) = 1 + O(ℓ^(−2)) for ℓ > z > 2√(Y+1) (1 exactly for ℓ | N).
The moduli are Q(n) ≤ 16N^(4δ/5) and the conductors of χ_{Δ₀} are ≤ 4(Y+1)N²: the two GRH families of the Theorem. The
errors depend on N only through log N and ω(N) ≤ 2 log N. ∎ (Key lemma (b).) (The referee of `mprime_theorem.md` notes that
with log z ≍ √(log N) both errors become exp(−c√log N).)

### Step 8. Assembly (proved; dimension quoted; checked C4)

*Dimension (quoted):* |H_k(N)| = (k−1)φ(N)/12 + O_k(2^(ω(N))) for squarefree N (Martin, J. Number Theory 2005, Thm 1, not
opened; used in this form by arXiv 2609.29631, read). Only O(2^ω) is used: the sharper |dev| ≤ 2^(ω(N)−1), true for
squarefree N ≤ 3000, fails at N = 701569 = 11·23·47·59 (dev 25/3 > 8; checked C4, which also checks |dev| ≤ 2^(ω+1) for
squarefree N ≤ 3000, k = 2, 4: max |dev|/2^ω = 0.500 (k = 2), 0.250 (k = 4)). Hence 12|H_k(N)|/((k−1)φ(N)) = 1 + O(N^(−1+ε)).
*PNT (quoted, not opened):* θ(x) = x + O(x exp(−c√log x)) (de la Vallée Poussin); so Σ_{p∈J, p∤N}log p = hN(1 + O((log N)^(−A)))
and Σ_{p∈J}log p·p/N = N∫_{y₀}^{y₀+h}y dy + O(N exp(−c√log N)) for h ≥ (log N)^(−B) (under GRH: error √N log²N).
With y_p := p/N, Lemma 2 gives ½H(4pN) = (N√y_p/π)L*(D₀(p)) and H(D_r(p)) = (N√(4y_p − r²)/π)L*(D_r(p)), so Lemma 1 gives

    E_N(p) = (1 + O(N^(−1+ε)))·{(12κ(N)/((k−1)π))[√y_p L*(D₀) + (−1)^(k/2−1)Σ_r U_{k−2}(r/(2√y_p))√(4y_p − r²)L*(D_r)]
             − δ_{k,2}·12κ(N)(y_p + 1/N)}.

Apply the key lemma with G₀(v) := √v and G_r(v) := U_{k−2}(r/(2√v))√(4v − r²) = 2√v sin((k−1)θ_r(v)) for v > r²/4
(cos θ_r = r/(2√v)), G_r := 0 otherwise: continuous, ‖G_r‖ ≪_k √(Y+1). With (12/π)κ(N)E^(N)(0) = 2ππ(N) and E^(N)(r) = ν_N(r)
(Lemma 4) the window mean is (1/h)∫M_k^(N) with error O(κ(N)(log N)^(−A)) = O((log N)^(−A+1)); rename A. ∎ (Theorem.)

### Step 9. Proof of the Proposition (proved; checked D5)

*Mean-value lemma (proved).* Let g(N) = Π_{q|N}c(q) on squarefree N with c bounded and c(ℓ) − 1 = O(1/ℓ). Then
Σ_n μ²(n)g(n)n^(−s) = Π_ℓ(1 + c(ℓ)ℓ^(−s)) = (ζ(s)/ζ(2s))·Π_ℓ(1 + (c(ℓ)−1)ℓ^(−s)/(1 + ℓ^(−s))), so μ²g = μ² ∗ b with
b(ℓ^j) = (c(ℓ)−1)(−1)^(j−1) and Σ_d|b(d)|d^(−1/2) = Π_ℓ(1 + O(ℓ^(−3/2))) < ∞. With Σ_{m≤x}μ²(m) = 6x/π² + O(√x):
Σ^□_{N≤X}g(N) = (6X/π²)G(1) + O(√X), G(1) = Π_ℓ (ℓ + c(ℓ))/(ℓ + 1); by partial summation Σ^□_{N≤X}N g(N) = (3X²/π²)G(1) + O(X^(3/2));
over (X, X + X′] the main terms are (6X′/π²)G(1) and (6XX′/π²)G(1)(1 + O(X′/X)) while the errors stay O(√X) and
O(X^(3/2)), hence the condition X′ ≥ X^(1/2+ε).
*Application.* φ(N)M_k^(N)(y) = (N/(k−1))[2π√y Π_{q|N}(1 − q^(−2)) + (−1)^(k/2−1)(12/π)Σ_r ν′(r)Π_{q|N}e_q(r)^(−1)U√(4y − r²)]
− 12Ny[k = 2], and Σφ(N) = Σ N·Π(1 − 1/q). The functions c needed are 1 − ℓ^(−1), 1 − ℓ^(−2), e_ℓ(r)^(−1) (bounded;
c − 1 = O(ℓ^(−3)) for ℓ ∤ r): all admissible. *Euler factors (exact, prime by prime; checked D5 for every ℓ < 1000):*
G_{1−1/ℓ}: (ℓ²+ℓ−1)/(ℓ(ℓ+1)), whose inverse product is γ/12; G_{1−ℓ^−2}/G_{1−1/ℓ}: (ℓ³+ℓ²−1)/(ℓ(ℓ²+ℓ−1)) = the factor of β/2π;
ν′(r)G_{e^−1}/G_{1−1/ℓ}: ℓ(ℓe_ℓ(r)+1)/(ℓ²+ℓ−1), and with (12/π)/(2π) = Π(1 − ℓ^(−2)) the factor for ℓ ∤ r is
(ℓ⁴−2ℓ²−ℓ+1)/(ℓ⁴−2ℓ²+ℓ) = the factor of α/2π, the ratio for ℓ | r is 1 + ℓ²/(ℓ⁴−2ℓ²−ℓ+1) = her ν; the mixture
(ℓe_ℓ(r)+1)/(ℓ+1) is her local B (ℓ ∤ r) and A (ℓ | r, r = 0 included). ∎ **Brute force (D5):** over squarefree N ≤ 3·10⁵, the φ-weighted means of M_2^(N)(y) at y = 0.1, 0.5, 1, 2, 3.3 are 0.7829, 1.0077, −2.3039, 0.0041, 0.1453 against M_zub,2 = 0.7829, 1.0077, −2.3039, 0.0040, 0.1453.

## 4. Inputs

| input | where | status |
|---|---|---|
| Skoruppa–Zagier trace formula of T_p∘W_N | Lemma 1 | quoted (SZ 1988 not opened; read in 2609.29631 and 2310.07681 §2); checked C1 |
| W_N on old spaces (Atkin–Lehner 1970, Lemma 26); old forms give trace 0 | Lemma 1 | lemma quoted, not opened; vanishing proved; checked C2 |
| ε_f = (−1)^(k/2)w_N(f), a_f(q) = −q^(k/2−1)w_q(f) | Lemma 1 | quoted (Atkin–Lehner), not opened; checked C3 |
| H(D) = Σ_{m∈𝓜(D)}h′(−D/m²); h′(Δ) = √|Δ|L(1, ξ_Δ)/π | Lemma 2 | quoted classical (Cox Thm 7.24, Cor. 7.28, not opened); checked A1–A3 |
| N \| Δ₀, 𝓜(D) prime to N, F_q = 1 at q \| N | Lemma 3 | proved; checked A4, A5 |
| local means e_ℓ(r) | Lemma 4 | proved; checked B1 |
| Burgess, r = 2, any modulus | 6.2 | quoted (Burgess 1963 not opened; read (secondary): Treviño arXiv 2511.17778 §1); effective |
| L(1, χ) ≤ log q + 2; Möbius for ξ_Δ; partial summation; Kronecker periodicity and multiplicativity | 6.1–6.3 | proved; checked A5 |
| Jacobi sums of modulus √ℓ | 6.3(c) | quoted classical, not opened; checked A6 |
| Siegel–Walfisz, moduli ≤ (log x)^B | 6.5 | quoted (not opened; read (secondary)); **ineffective** |
| multiplicative large sieve, constant \|I\| + R² | 6.6 | quoted (read as (1.1) of Friedlander–Iwaniec arXiv 2409.19634); effective |
| Σ_{n≤X} n/φ(n) ≤ (ζ(2)ζ(3)/ζ(6))X; w_d bound | 6.6 | proved; checked A7 |
| main-term coefficients e_ℓ(i, j), weighted absolute convergence | 6.7 | proved; checked B2 |
| PNT with de la Vallée Poussin error | Step 8 | quoted, not opened; effective |
| dimension of the new space, O(2^ω) | Step 8 | quoted (Martin 2005 Thm 1, not opened); checked C4 |
| GRH for Dirichlet L-functions ⇒ θ(x, χ), θ(x; q, b) with error √x log²(qx) | Step 7, part (b) only | hypothesis (fact 2047 of the project); consequence quoted from Davenport ch. 20, not opened |
| mean of μ², Rankin's trick, CRT, dominated convergence | Steps 6.7, 9 | proved |
| Zubrilina's constants α, β, γ, ν(r) | Proposition | quoted (2310.07681 Thm 1, 2609.29631 Thm 1.1: read in the earlier note); matched to 7 digits by its referee |

## 5. Slips of the three earlier notes and where they are repaired

From the referee of `mprime_theorem.md` (already repaired there in place; carried over correctly):
1. Range of Q_ℓ: 1 ≤ Q_ℓ < (ℓ+1)/(ℓ−1), not [1 − 1/ℓ, 1 + 2/ℓ] — Lemma 2(iii).
2. Primes p ≤ ℓ in the Petersson step (Theorem 1 of that note) cost O(τ(ℓ)ℓ√N log N) — that step (harmonic weights with
   Hecke twists) is not used in this note; the Theorem is proved by the trace route only.
3. The reason given for the low χ² values is unproved — Section 6 states it as measured, not explained.
4. P_z(D) ≪ log z (not (log z)³) — 7.1; E(r) = ν′(r)(1 + O(N^(−3))) at prime level (not O(N^(−2)); also in
   `mprime_unconditional.md`) — Special case 1.
From the referee of `mprime_unconditional.md`:
5. Truncation error N^(−1/16+o(1)), not N^(−1/32+o(1)); Burgess applies to every divisor d | f — 6.2.
6. N divides the fundamental discriminant, so the conductor is never small — Lemma 3(i), now for squarefree N and at 2.
7. The period of T(m, n; ·) is nm², not 2nm² — 6.3(a), checked A5.
8. The term "p | q" in the large sieve is empty for r ≥ 1 (and absent for r = 0) — 6.6.
From the referee of `level_theorem.md`:
9. The statement needs H_k(N) ≠ ∅ (N = 6, 10, 22 have no newforms of weight 2) — Theorem, Remark (1), checked C4.
10. "ε ≤ 1" is false in general — Section 2: ratios in [1, e^(O(ε))], error O(κ(S)εe^(3ε)), checked D6.
11. "explains the footnote" of 2310.07681 → "suggests" (her footnote concerns non-squarefree levels) — Section 2, Special case 3.
12. The dimension deviation 2^(ω−1) fails beyond N = 3000 (N = 701569: 25/3 against 8); only O(2^ω) is used — Step 8, C4.
Made explicit in this merge (not slips of a statement, but steps the earlier notes left implicit): the explicit constant
c₀ of w_d (6.6); the Rankin bounds of the left-out main-term pairs with weight ω(m, n) and the choice L = (log N)^(2A+6)
instead of (log N)^(A+3) (6.7); the bound P(ℓ^(2i) | D) ≤ 8ℓ^(1−2i) for all ℓ (6.7); the key lemma for arbitrary weights of
bounded variation, which gives consequence (ii) directly.

## 6. Numerical evidence

**Round 160 (blind, prime levels; STUDY.md).** All 110 prime levels 2300–3200, weight 2, every p ≤ 12N, 48 bins: χ² 17.1
against M_2 and 1,342 against Zubrilina's average density (ratio 78); rms of data − density 0.058 against 0.409; y ≤ 6:
7.0 against 875; y > 6: 10.1 against 468. Every registered decision passed. The χ² of 17 on 48 bins is too low for errors
of the right size: the registered error (the scatter of single values in a bin) overstates the error of a bin mean, so
the comparison of the two densities decides.
**Round 163 (blind, squarefree levels; STUDY.md).** 40 levels never computed before, ten of each type q, 2q, 3q, 6q
(N = 30,034 to 30,522), weight 2, every p ≤ 4N, 16 windows; errors frozen from levels near 10⁴ as the scatter between
levels (predictions sha256 2f67dc3a…). χ² on 160 terms for the own density: 47.7 (q), 68.6 (2q), 33.3 (3q), 49.9 (6q); for
the other three type densities 2,910–17,658, 3,607–9,099, 1,673–9,307, 2,062–6,762; for Zubrilina's average 734, 3,829,
1,263, 4,581; rms(data − own) 0.103, 0.163, 0.140, 0.178. Every registered decision passed (own χ² ≤ 256; every alternative
worse by ≥ 100). The own χ² are low because the frozen errors from N ≈ 10⁴ are too large at N ≈ 3·10⁴; the separation
decides. Disclosed in the registration: a timing run had summed class numbers at one test level and printed only the time.
(Not recomputed here; no level in 30030–30530 was touched.)
**The low χ² (repair 3).** Not proved. Measured in `level_theorem.md`: the scatter of window means between levels is about
0.53 times the single-value standard error (20 squarefree levels near 10⁴). A heuristic (not proved): most of the variance of
L*(D(p)) comes from the local factors at small ℓ, i.e. from p mod small moduli, which average out over a window faster than
independent samples would.
**The referees' own tests.** *Of `mprime_theorem.md`:* local means by enumeration with other code (2/3, 4/3 at ℓ = 2; 15/16,
9/8 at 3; 95/96, 25/24 at 5; F_N = 1); T1 exact for every p ≤ 4N at k = 2 (N = 167, 173, 197), 4 (29, 31), 6 (17), 8 (19),
14 (11), 12 (13, 17: old forms give 0); one level N = 30011, all p ≤ 4N, 16 windows: χ² 3.2 against M_2 and 70.5 against the
squarefree average (weight 4: 1.5 and 388.8); N = 30029: 4.0 and 69.5 (1.8 and 410.2). *Of `mprime_unconditional.md`:* the
class number formula for orders on 300 discriminants (2·10⁻¹⁴); the factorisation of 6.3 in 23,705 random cases; the
main-term double series 1.64493, 0.61488, 1.22888, 0.73798 against ζ(2), ν′(1), ν′(2), ν′(3) = 1.64493, 0.61513, 1.23027,
0.73816 (the gap is the tail in m); at N = 20011, all 1,942 primes in (N, 2N]: means by conductor block q = 1: +0.8201
(predicted +0.8193), 1 < q ≤ 32: −0.2080 (−0.2075), every block above 32 at most 0.003. *Of `level_theorem.md`:* trace
formula exact at 1,644 triples (18 levels, up to three prime factors, weights 2–12, old forms included); local factors in
396 cases; Zubrilina's constants from her source equal the φ-weighted averages to 7 digits; levels 2·10007, 3·6673,
6·3343, 5·4003: χ² 17.3, 21.5, 17.7, 14.1 for the own density against 435–1,623 for the prime-level density (72 for 5q)
and 91–1,063 for the squarefree average.
**This note's script (section E; levels 9973, 2·4999, 3·3329, 6·1663, not used before; weight 2; all p ≤ 4N; 16 windows of
width ¼; errors = single-value standard errors of the window means, which overstate the error; rule fixed before the run:
own χ² ≤ 32, every other density ≥ own + 9):**

| level | type | dim | primes p ≤ 4N | χ² own density | Zubrilina's average | type q | type 2q | type 3q | type 6q |
|---|---|---|---|---|---|---|---|---|---|
| 9973 | q | 830 | 4,194 | **3.7** | 26.0 | — | 311.7 | 85.3 | 588.9 |
| 9998 = 2·4999 | 2q | 417 | 4,201 | **6.3** | 118.3 | 174.2 | — | 251.9 | 103.7 |
| 9987 = 3·3329 | 3q | 555 | 4,196 | **6.1** | 59.0 | 76.0 | 347.4 | — | 430.0 |
| 9978 = 6·1663 | 6q | 277 | 4,193 | **6.6** | 218.8 | 305.3 | 90.9 | 284.9 | — |

Every decision of the rule passes; the weakest separation is the prime level against Zubrilina's average (Δχ² 22.3 with
these conservative errors; 78 times in round 160 with 110 levels).

Class-number means over the primes (key lemma main term at one composite level): log-weighted means of L*(D_r(p)) = πH(D_r)/√D_r over the primes p ≤ 4N against E^(N)(0), ν_N(1), ν_N(2), ν_N(3): max |z| 0.84 over the 16 cases (for example 2·4999, r = 1: 0.9230 against ν_N(1) = 0.9227, while the prime-level value is ν′(1) = 0.6151; 6·1663, r = 0: 1.0983 against ζ(2)(3/4)(8/9) = 1.0966).

## 7. Relation to the literature

What follows is taken from the project's prior-art gate of 26 November 2026 and from the three earlier notes; the
sources were opened there, not for this merge.
- **Zubrilina, arXiv 2310.07681, and the correction arXiv 2609.29631** ("On Murmurations"): murmuration densities with
  natural weights averaged over squarefree levels N in an interval, at fixed p/N; 2609.29631 corrects the place of U_{k−2}
  in the trace formula. She fixes the prime and sums over levels; the statements here are not implicit in her proof.
  The Proposition (Special case 3) recovers her density as the dimension-weighted average of the M^(N).
- **Cardi, arXiv 2607.08982 v2, §5** (averaged theorem over squarefree N ≍ X with weight k ≥ X^ε): as a step inside the
  average it has, for each fixed squarefree N, the prime-mean of the local data of N(r²N − 4p) as an Euler product with
  local sums 1 − 1/((ℓ−1)(ℓ²−1)) for ℓ ∤ Nr, (1 − ℓ^(−2))^(−1) for ℓ | r, ℓ ∤ N, 2/3 or 4/3 at 2, and exactly 1 at every
  ℓ | N; a proposition g_r(N) = C_r·Π_{ℓ|N}a_{ℓ,r}^(−1) (our ν_N(r)); and a local Fubini identity (primes first or levels
  first, P(ℓ | N) = 1/(ℓ+1)). It has no statement for one fixed level and fixed weight, no factors N/φ(N), Π(1 + 1/q), no
  term −12(N/φ(N))y, no densities for 2q, 3q, 6q. So the **local recipe of the squarefree theorem is known in substance**,
  and the **Proposition is known in substance for α·ν(r) and β** (Cardi); γ = 12⟨N/φ(N)⟩ is a one-line addition.
  (Cardi arXiv 2609.00777 was seen as an abstract only in the earlier notes.)
- **Tomczak, arXiv 2606.08353 (June 2026), §1, Theorem 1.2(i):** levels ℓ^(2a), a ≥ 2 fixed, ℓ → ∞ prime, natural weights,
  log-weighted mean over primes P with P/ℓ^(2a) in a fixed interval, one level at a time, under GRH for Dirichlet
  L-functions, error O(ℓ^(−a/5+ε) + ℓ^(−1)); his density M_(∞,k) equals the prime-level M_k of Special case 1 (checked by
  the gate numerically to 6·10⁻⁸ at k = 2, 4, 8); his Remark 4 notes Zubrilina's shape with other coefficients; he says odd
  exponents "should hold". So **the prime-level density is known in substance** for this neighbouring family.
- **arXiv 2603.25564 v2** ("Murmurations in the depth aspect"): levels ℓ^(2a+1), ℓ fixed, a → ∞, one level at a time,
  natural weights, unconditional (a Bombieri–Vinogradov variant for square moduli, a mean-value bound for real characters;
  in v2 the tail is bounded by Heath-Brown's fourth moment, after Tomczak's Remark 1 doubted the tail step of v1). Its
  density is not the squarefree recipe on the set {ℓ} (extra terms from a Ramanujan sum) and tends to M_k as ℓ → ∞. So the
  **method "one level, mean over the primes" is in print, unconditionally, for prime-power levels.**
- **Cowan, arXiv 2408.12723:** a fixed prime level under GRH with **harmonic** weights (Thm 1.4, from Miller–Montague).
- **Martin, arXiv 2409.02338, §1.3:** conjectures murmurations for levels p₁q (p₁ fixed) per Atkin–Lehner eigenspace,
  averaged over q.
- **Local means are classical:** ν′(r) = (π/2)·c_r with c_r the average Lang–Trotter constant (David–Pappalardi 1999,
  class numbers H(r² − 4p) averaged over primes p); Fouvry–Murty 1996 for r = 0. Averages of L(1, χ) over shifted primes by
  the large sieve are classical.
**Claimed as not found in print** (and nothing else): (1) the theorem at a single prime level N (exponent 1, where the
term −12y comes from the Eisenstein part), in particular part (a) without hypothesis for windows ≥ (log N)^(−B) as N → ∞
through primes, and part (b) for windows down to N^(−1/2+δ); (2) the theorem at a single fixed squarefree level with fixed
weight, parts (a) and (b), with the factors κ(N), π(N) and the term −12κ(N)y; (3) the explicit densities and the blind
numerics for 2q, 3q, 6q (round 163). **Not claimed:** the method (Zubrilina's trace route with the roles of N and p
exchanged; the unconditional large-sieve treatment as in 2603.25564), the local recipe (Cardi), the averaging identity for
α·ν(r) and β (Cardi), the prime-level density as a function (Tomczak), the local means (David–Pappalardi, Fouvry–Murty).

## 8. Open

- **Shorter windows without GRH: NOT closed.** 6.2–6.6 hold for h ≥ N^(−1/8+δ) (M = N^(7/8)); with Burgess r = 3 (any
  modulus, M = N^(2/3+δ)) for h ≥ N^(−1/3+2δ); with Burgess for all r on real characters (odd part of Δ₀ squarefree, the
  2-part split mod 8; a sketch, not written out) for h ≥ N^(−1/2+2δ). The obstacle is 6.5 alone: primes in progressions to
  moduli ≤ (log N)^B in intervals of length hN ≤ N^(1−η) (Huxley-type results; not opened).
- **Effectivity: NOT closed.** The only ineffective input is Siegel–Walfisz in 6.5 (a possible Siegel zero of a real
  character of conductor ≤ (log N)^(B₂)). Page's theorem (one exceptional modulus ≤ (log N)^(2−ε)) does not suffice as
  set up, because the large-sieve losses force Q₀ far above (log N)². The rate (log N)^(−A) is the same as under GRH in
  this note; a power saving in (a) would need one in 6.5.
- **Non-squarefree levels: NOT treated.** Old forms with N/M a square contribute (2310.07681 §2, read); q² | D_r(p) for
  q² | N, so the local factor at q is a different computation (compare 2603.25564 and Tomczak for ℓ^(2a+1), ℓ^(2a)).
- **Atkin–Lehner signs separately: NOT closed.** Σ_f w_Q(f)a_f(p) for Q ∥ N, Q ≠ N is a trace of T_pW_Q: old forms no
  longer cancel (inversion over Q | M | N with weights 2^(ω(N/M)) needed), and the trace formula has class numbers of
  r²Q² − 4pQ with a local factor counting solutions mod N/Q (not opened); the natural variable is p/Q. Martin's conjecture
  for p₁q (averaged over q) is of this kind.
- **Single primes p: nothing claimed.** The Theorem is an average over p in a window of length ≥ N(log N)^(−B)
  (≥ N^(1/2+δ) under GRH); for a single p, E_N(p) is a finite sum of class numbers with no averaging.
- **The Petersson route.** Harmonic → natural weights as a direct estimate needs GRH for L(s, Sym² f) (ILS Lemma 2.12) and
  stays NOT closed; it is not needed for the Theorem.

## 9. The check script

`one_level_check.py` (stand-alone; sympy, numpy, scipy, gmpy2, cypari2; run as
`(ulimit -v 1200000; timeout 600 nice python -I -B one_level_check.py)`, about 80 s). Sections: A class numbers, the
discriminants of a squarefree level, the factorisation and the large-sieve weights (A1–A7); B local means and the
coefficients e_ℓ(i, j) by exact enumeration (B1–B2); C trace formula at N = 23, 46, 69, 78 against PARI for k = 2, 4, 6, old
forms, root numbers, dimensions (C1–C4); D Euler products, special cases, the Proposition, the repaired ε-comparison
(D1–D6); E one level of each type below 10⁴ (E0–E3). Final line of `one_level_check_out.txt`: **32 checks, 0 failed**.
