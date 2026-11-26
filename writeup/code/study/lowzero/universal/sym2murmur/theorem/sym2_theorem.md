# The symmetric-square family mean at one prime level, p ≍ N²: what is proved given fact 2047, and where the proof stops

Folder: `scratchpad/sym2_thm/` (nothing in the repository was changed). Check script `sym2_check.py`, output
`sym2_check_out.txt`. Template: `writeup/code/study/lowzero/universal/murmur/theorem/mprime_theorem.md` (cited by its step
labels T1–T6 and "key lemma"; this note inherits the referee's verdict on those steps). Labels: **proved** (written out
here), **quoted** (reference; "read" = I opened the text, "not opened" = cited from memory or via the template),
**checked** (a named assertion of the script), **NOT closed** (a gap).

## 0. Verdict

**STOPS at step (d), at a named hypothesis (Hypothesis T, §5.4).** Given fact 2047 (GRH for Dirichlet L-functions only),
everything else is proved: the exact trace identity and the reduction t = ±(2p − mN²) (steps a, b); the class-number part
with the local factors ν(m) (Theorem A, windows of relative length h ≥ N^{−1+δ}); and, inside step (d), the *head* of the
twisted sum R: after expanding L(1, χ_{t²−4p²}) into its Dirichlet coefficients n ≤ L = N^{δ/10}, keeping the local
factor at N exactly, Poisson summation in t (complete sums mod N·Q(n), Weil) and GRH in p give exactly the claimed smooth
term 24[k=2]y^{3/2} (Proposition H, windows h ≥ N^{−1/2+δ}). What is missing is the *tail*: the twisted sum of
τ_L(D) := L*(D) − (head) over the family D = 4p² − t² must average to o(N^{−2}); pointwise GRH gives L^{−1/2+ε}, the head
method reaches only n ≤ N^{1/16}, and the quadratic large sieve gives O(h^{−1/2}N^ε). Given the proved parts,
**Hypothesis T is equivalent to the statement** (Theorem C): it is an arithmetic reformulation of the claim (class
numbers and real characters only, no modular forms), not a weaker assumption. For windows inside (0, 1/4), where the
class-number part is empty, the statement is literally the hypothesis R of the review.

## 1. Statements

**Setting.** N ≥ 5 prime, k ∈ K := {2, 4, 6, 8, 10, 14} (S_k(SL₂(ℤ)) = 0, so S_k(Γ₀(N)) consists of newforms; template
step (e)), H_k(N) its Hecke eigenbasis, d := |H_k(N)|, ρ_N := (k−1)N/(12d) = 1 + O_k(1/N) (template F3). λ_f(n) =
a_f(n)n^{−(k−1)/2}; for p ≠ N, λ_f(p²) = λ_f(p)² − 1 and Σ_f λ_f(p²) = Tr T_{p²}/p^{k−1}. Statistic
S_N(p) := d^{−1}Σ_f √p λ_f(p²). Window J := (y₀N², (y₀+h)N²], ⟨X⟩_J := Σ_{p∈J, p≠N} log p·X(p)/Σ_{p∈J, p≠N} log p, y_p := p/N².
ν(m) := Π_ℓ ν_ℓ(v_ℓ(m)) over all primes, x := 1/ℓ, S_a := Σ_{j<a} x^j, u_ℓ := 1/((ℓ−1)²(ℓ+1)):

    ν_ℓ(0) = 1 − u_ℓ,   ν_ℓ(2a) = S_a + x^a/(1 − x²) (a ≥ 1),   ν_ℓ(2a+1) = S_a + x^a      (= derivation.md §0, check C3)
    M_k(y) = 24[k=2]y^{3/2} − 6√y/(k−1) − (12/(π(k−1)√y))·Σ_{1≤m<4y} ν(m)√(m(4y−m))·U_{k−2}(1 − m/(2y)).

H(D) is the Hurwitz class number, L*(D) := πH(D)/√D = L(1, χ_{Δ₀})Π_{ℓ|f₀}Q_ℓ, F_ℓ(D), β_D(n), Q(n) ≤ 16n⁴ as in template
T2, T4. For N ∤ D, F_N(D) = (1 − χ_N(−D)/N)^{−1}. D_t := 4p² − t², D_m(p) := m(4p − mN²).

**Exact split (step a–b).** S_N = S^I + S^H + S^{CN} + S^R with
S^I = ρ_N/√p, S^H = −(√p/(2d))(1 + 2p^{1−k}), S^{CN} = −(12ρ_N/(π(k−1)))Σ_{1≤m<4y_p} G_m(y_p)L*(D_m(p)),
G_m(y) := y^{−1/2}√(m(4y−m))U_{k−2}(1 − m/(2y)) for y > m/4 (0 otherwise), and the **twisted part**

    S^R(p) := −(1/(2√p d))·Σ_{|t|<2p} U_{k−2}(t/2p)·H(D_t)·χ_N(t² − 4p²)
            = −(1/(2π√p d))·Σ_t W_p(t)χ_N(−D_t)L*(D_t),    W_p(t) := U_{k−2}(t/2p)√(4p² − t²).

Head and tail (L ≥ 1): head_L(D) := F_N(D)Σ_{n≤L, N∤n} β_D(n)/n, τ_L(D) := L*(D) − head_L(D); S^R = S^{R,head}_L + S^{R,tail}_L
by inserting head_L, τ_L for L*.

**Theorem A (class-number part; proved given fact 2047).** Fix k ∈ K, Y ≥ 1, δ ∈ (0, 1), A ≥ 1. Uniformly for N prime,
0 ≤ y₀ ≤ Y, N^{−1+δ} ≤ h ≤ 1:  ⟨S_N − S^R⟩_J = (1/h)∫_{y₀}^{y₀+h}(M_k(y) − 24[k=2]y^{3/2})dy + O_{k,Y,δ,A}((log N)^{−A}).

**Proposition H (head of the twisted sum; proved given fact 2047).** Fix k ∈ K, Y, δ ∈ (0, ½); L := N^{δ/10}. Uniformly
for 0 ≤ y₀ ≤ Y, N^{−1/2+δ} ≤ h ≤ 1:  ⟨S^{R,head}_L⟩_J = 24[k=2]·(1/h)∫_{y₀}^{y₀+h} y^{3/2}dy + O_{k,Y,δ}(N^{−δ/25}).

**Hypothesis T(δ) (NOT closed).** With L = N^{δ/10}: ⟨S^{R,tail}_L⟩_J = o(1) uniformly for 0 ≤ y₀ ≤ Y, N^{−1/2+δ} ≤ h ≤ 1.

**Theorem C (proved given fact 2047).** Under the hypotheses of Proposition H,
⟨S_N⟩_J − (1/h)∫_{y₀}^{y₀+h} M_k(y)dy = ⟨S^{R,tail}_L⟩_J + O(N^{−δ/25}). Hence the murmuration statement
"⟨S_N⟩_J → (1/h)∫M_k" holds for these windows **if and only if** Hypothesis T(δ) holds. For y₀ + h ≤ ¼ it is equivalent
to ⟨S^R⟩_J = 24[k=2](1/h)∫y^{3/2} + o(1) (hypothesis R). Nothing is claimed for single primes.

## 2. Step (a): the exact trace identity (quoted, read; checked A1, A2, A6)

**Source (read).** Kaplan–Petrow, *Elliptic curves over a finite field and the trace formula*, arXiv:1510.03980, Theorem
"Eichler–Selberg trace formula for Γ₀(N), χ" (label ESTFKL; they quote it from Knightly–Li): for k ≥ 2, χ(−1) = (−1)^k,
Tr(T_q | S_k(Γ₀(N), χ)) = ((k−1)/12)ψ(N)χ(q^{1/2})q^{k/2−1} − ½Σ_{t²<4q}U_{k−2}(t, q)Σ_{m²|t²−4q}h_w((t²−4q)/m²)μ_χ(t, m, q)
− ½Σ_{b|q}min(b, q/b)^{k−1}Σ_τ φ((τ, N/τ))χ(y_τ) + δ_{k,2}1_{χ=1}Σ_{c|q, (N, q/c)=1}c, with U_{k−2}(t, q) =
q^{k/2−1}U_{k−2}(t/(2√q)) and μ_χ(t, m, q) = (ψ(N)/ψ(N/(N, m)))·#{c ∈ (ℤ/N)^× lifting to a root of c² − tc + q mod N(N, m)}.

**Specialisation (proved).** χ = 1, N prime, (q, N) = 1: for N ∤ m the count of roots of c² − tc + q mod N is
1 + χ_N(t² − 4q); for N | m (then N² | t² − 4q) exactly c ≡ t/2 lifts, μ = N + 1; τ ∈ {1, N} gives the factor 2 in the
hyperbolic term; the last term is σ(q). At level 1 all μ = 1. Grouping Σ_{m}h_w(D/m²) by N | m:
Σ_m h_w μ = H(D)(1 + χ_N(−D)) + N·H(D/N²)[N² | D] (when N² | D, χ_N(−D) = 0). For q = p² and k ∈ K (Tr on S_k(SL₂ℤ) = 0),
subtracting the level-1 formula: t with N ∥ D_t or N | D_t, N² ∤ D_t have equal μ at both levels and cancel, and

    Tr T_{p²} = (k−1)Np^{k−2}/12 − (2 + p^{k−1})/2 + A2^χ + A2^{N²},
    A2^χ = −½Σ_{|t|<2p}U_{k−2}(t, p²)H(D_t)χ_N(t² − 4p²),   A2^{N²} = −(N/2)Σ_{|t|<2p, N²|D_t}U_{k−2}(t, p²)H(D_t/N²).

**Checked:** A1 (own implementation = PARI `mftraceform` for N = 37, 43, k = 2, 4, all n ≤ 600 prime to N, exact);
A2 (the decomposition = PARI Tr T_{p²} for N = 11, 23, k = 2, 4, all primes p ≤ 150, including 25 + 3 primes with m-terms);
A6 (dim = Tr T₁ = PARI `mfdim`); A4–A5 (the decomposition = an independent Gram route Σ_f a_f(p)² − p^{k−1}d for every prime
p ≤ 1500 at N = 37, 43, 139, 211; the Gram route also exactly over ℚ).

## 3. Step (b): which t carry the main terms (proved; checked A3)

**Lemma B.** N odd prime, p ≠ N prime. N² | D_t = (2p−t)(2p+t) with |t| < 2p iff t = ±(2p − mN²), 1 ≤ m < 4p/N²; then
D_t/N² = m(4p − mN²) = D_m(p); the 2·#{m} values of t are distinct; U_{k−2}(±(2p − mN²), p²) = p^{k−2}U_{k−2}(1 − m/(2y_p)).
*Proof.* (2p−t) + (2p+t) = 4p is prime to N, so N divides at most one factor and N² divides it: 2p ∓ t = mN², m ≥ 1,
and 2p ± t = 4p − mN² > 0. Distinctness: 2p − mN² = −(2p − m′N²) would give 4p = (m + m′)N²; t = 0 would give N | 2p.
U_{k−2} is even and homogeneous: U_{k−2}(t, p²) = p^{k−2}U_{k−2}(t/2p). ∎
So A2^{N²} = −NΣ_m p^{k−2}U_{k−2}(1 − m/(2y_p))H(D_m(p)); with √D_m = N√(m(4y_p − m)) and division by p^{k−1}d/√p this is
S^{CN}; the identity and hyperbolic terms give S^I and S^H (the "2" in 2 + p^{k−1} gives the factor 1 + 2p^{1−k});
A2^χ gives S^R. **The main term of the class-number part comes from the 2·⌊4y⌋ values t = ±(2p − mN²); every other t
with N ∤ D_t sits in S^R; the t with N | D_t, N² ∤ D_t contribute nothing.** For m < 4(Y+1) < N, D_m(p) ≡ 4mp (mod N), so
N ∤ D_m(p), and −D_m ≡ 0, 1 (mod 4) (D_t ≡ 0, 3 (4), N² ≡ 1 (8)).

## 4. Step (c): the class-number part (Theorem A)

**Lemma L (local means; proved; checked C1, C2).** For a prime ℓ ≠ N, p Haar-uniform on ℤ_ℓ^× and e := v_ℓ(m):
E_ℓ(m) := ∫F_ℓ(D_m(p))dp = ν_ℓ(e). For ℓ = N: E_N = N²/(N² − 1).
*Proof, ℓ odd.* u := 4p − mN² is Haar-uniform on ℤ_ℓ ∖ (−mN² + ℓℤ_ℓ) with density ℓ/(ℓ−1), and −D = −mu.
(e = 0) The mean of F_ℓ(−mu) over all of ℤ_ℓ is 1/(1 − x²): it depends only on the law of (v_ℓ(u), quadratic character of
the unit part), which a unit factor −m does not change (template T5, case ℓ ∤ r). On the removed disk −D ≡ m²N² is a
non-zero square: F = (1 − x)^{−1}. So E_ℓ = (ℓ/(ℓ−1))[ℓ²/(ℓ²−1) − 1/(ℓ−1)] = ℓ(ℓ²−ℓ−1)/((ℓ−1)²(ℓ+1)) = 1 − u_ℓ.
(e ≥ 1) The removed disk is ℓℤ_ℓ, so u is a uniform unit, v_ℓ(D) = e exactly and the unit part of −D is uniform. e = 2a:
ℓ^a ∥ f₀, χ(ℓ) = ±1 with probability ½ each: F = S_a + x^a·½[(1−x)^{−1} + (1+x)^{−1}] = S_a + x^a/(1−x²). e = 2a+1:
ℓ^a ∥ f₀, ℓ | Δ₀, χ(ℓ) = 0: F = S_a + x^a.
*ℓ = 2* (m = 2^e m′, N² ≡ 1 (8), p odd). For −D = 2^V w (w a uniform odd unit) the mean of F₂ is G(2b) = 2 − (4/3)2^{−b}
(w ≡ 1 (4): a = b, χ = ±1; w ≡ 3 (4): a = b−1, χ = 0) and G(2b+1) = 2 − 2^{1−b} (a = b−1, χ = 0). e = 0: −D ≡ 5 (8), F₂ = 2/3
= 1 − u₂. e = 1: −D = 4w, w = m′²N² − 2pm′ ≡ 3 (4): F₂ = 1. e ≥ 3: −D = 2^{e+2}·(uniform unit): G(e+2) = 2 − 2^{−(e−1)/2}
(e odd) or 2 − (2/3)2^{−e/2} (e even), which is S_a + x^a resp. S_a + x^a/(1−x²). e = 2: −D = 2^{4+s}·(uniform unit) with
P(s) = 2^{−s}: Σ_s 2^{−s}G(4+s) = 2 − 2/7 − 1/21 = 5/3 = S₁ + x/(1−x²).
*ℓ = N:* N ∤ D_m, χ(N) = (−mp/N) = ±1 equally often: ½[(1−1/N)^{−1} + (1+1/N)^{−1}] = N²/(N²−1). ∎
Hence E(m) := Π_{ℓ≠N}E_ℓ(m) = ν(m)/(1 − u_N) = ν(m)(1 + O(N^{−3})). (C1: 216 cases, ℓ = 2, 3, 5, 7, e ≤ 5, exact
enumeration of p mod ℓ^J with a rigorous enclosure of width ≤ 10^{−5}; C3: the closed form equals the formula printed in
derivation.md, exactly, ℓ < 60, e ≤ 12; C6: ν(1..16) and C₀ = Π(1 − u_ℓ) = 0.615133 as printed there.)

**Lemma K (proved given fact 2047; = the template's key lemma with D_r(p) replaced by D_m(p)).** For 1 ≤ m < 4(Y+1),
uniformly for windows J ⊂ (0, (Y+1)N²] of length hN² ≥ N^{1+δ} and G of bounded variation on [0, Y+1]:
Σ_{p∈J} log p·G(p/N²)L*(D_m(p)) = E(m)∫_J G(u/N²)du + O_A(‖G‖hN²(log N)^{−A}).
*Proof.* The template's proof (T3, T4, T6, Rankin) with four changes, each checked against what the template uses:
(i) D_m(p) = 4mp − m²N² is linear in p, so β_{D_m(p)}(n) depends only on p mod Q(n) (T4, T6 unchanged);
(ii) T3 is applied to D ≤ max_m m(4p − mN²) ≤ 4p²/N² ≤ 4(Y+1)²N² (real primitive characters of conductor ≤ 4(Y+1)²N²);
(iii) T6 with x ≤ (Y+1)N²: error per class O(N log²N), total O(‖G‖L⁵N log²N) for n ≤ L; with L := N^{δ/10} this is
≤ ‖G‖hN²·N^{−δ/2}log²N; (iv) the local means are Lemma L, and E_ℓ(m) = 1 − u_ℓ = 1 + O(ℓ^{−3}) for ℓ > z > 4(Y+1). The local
factor at N (1 + O(1/N) per discriminant) lies beyond z and is covered by the GRH error of T3, as in the template. ∎

**Proof of Theorem A.** By (a)–(b), S_N − S^R = S^I + S^H + S^{CN}. ⟨S^I⟩_J ≤ ρ_N⟨p^{−1/2}⟩_J ≪ (N√h)^{−1} ≤ N^{−1/2}.
⟨S^H⟩_J = −(6ρ_N/(k−1))⟨√y_p⟩_J + O(1/N) and ⟨√y_p⟩_J = (1/h)∫√y + O(log²N/(hN)) (T6 with q = 1). Lemma K with
G = G_m (continuous, vanishing at y = m/4, ‖G_m‖ ≪_k √(Y+1)) for each of the O(Y) values of m, E(m) = ν(m)(1 + O(N^{−3}))
and ρ_N = 1 + O(1/N) give ⟨S^{CN}⟩_J = −(12/(π(k−1)))Σ_m ν(m)(1/h)∫G_m + O((log N)^{−A}). Summing: the integrand is
−6√y/(k−1) − (12/(π(k−1)√y))Σ_m ν(m)√(m(4y−m))U_{k−2}(1 − m/(2y)) = M_k(y) − 24[k=2]y^{3/2}. ∎
**Checked:** D1 (one level N = 10007: means of πH(D_m)/√D_m over 2500 random primes p ∈ (2N², 3N²] against E(m), m = 1..8:
|z| ≤ 1.35); E4 (N = 211, k = 2: the exact class-number term over its ν-model, pooled, = 1.0000 ± 0.0006; N = 139, k = 4:
0.989 ± 0.014).

## 5. Step (d): the twisted sum R

### 5.1 The local model (Lemma M; proved; checked C4, C5)
(i) For every p ∈ (ℤ/N)^×: N^{−1}Σ_{t mod N, N∤D_t} χ_N(−D_t)F_N(D_t) = −2/(N² − 1).
*Proof.* Σ_t χ_N(t² − 4p²) = −1 (a quadratic with non-zero discriminant), with two zeros, so (N−3)/2 values χ = +1 and
(N−1)/2 values χ = −1; χ/(1 − χ/N) gives N^{−1}[(N−3)/2·N/(N−1) − (N−1)/2·N/(N+1)] = −4/(2(N²−1)). ∎
(ii) For ℓ ∤ 2pN, E′_ℓ := ∫_{ℤ_ℓ}F_ℓ(4p² − t²)dt = 1; also E′₂ = 1 (p odd).
*Proof.* ℓ odd: with probability 1 − 2/ℓ, −D is a unit and χ(ℓ) = ((t²−4p²)/ℓ) takes the value +1 for (ℓ−3)/2 and −1 for
(ℓ−1)/2 residues; with probability 2/ℓ exactly one of 2p ∓ t is divisible by ℓ, v_ℓ(D) = v ≥ 1 has law (1−x)x^{v−1} and a
uniform unit part, giving G_ℓ(2a) = S_a + x^a/(1−x²), G_ℓ(2a+1) = S_a + x^a. So E′_ℓ = x[((ℓ−3)/2)/(1−x) + ((ℓ−1)/2)/(1+x)]
+ 2x(1−x)Σ_{v≥1}x^{v−1}G_ℓ(v), four geometric series; the resulting rational function of x equals 1 identically (sympy,
C5). ℓ = 2: t odd (prob ½): −D ≡ 5 (8), F = 2/3; t = 2s, s even (¼): −D = 4w, w ≡ 3 (4), F = 1; s odd (¼):
−D = 2^{3+r}·(uniform unit), P(r) = 2^{1−r} (r ≥ 2), mean Σ_{c≥1}[2^{1−2c}(2 − 2^{−c}) + 2^{−2c}(2 − 2^{−c}/3)] = 5/3;
total ⅓ + ¼ + 5/12 = 1. ∎ (C5 also checks (ii) by enumeration with enclosure, ℓ = 2, 3, 5, 7.)
(iii) Consequence: replacing L*(D_t) by the product of its local means gives
S^R_loc(p) = −(1/(2π√p d))·(−2/(N²−1))·Ŵ_p(0) = [k=2]·2p^{3/2}/((N²−1)d) = 24[k=2]ρ_N y_p^{3/2}·N²/(N²−1), with
Ŵ_p(0) = ∫W_p = 4p²∫_{−1}^{1}U_{k−2}(u)√(1−u²)du = 2πp²[k=2]. This is the derivation's step 2.4 main term, now with the
exact constant −2/(N²−1) (the −1/N of χ_N and the +1/N of its local factor at N cancel to this order).

### 5.2 Proposition H: the head is the main term (proved given fact 2047)
Let n ≤ L = N^{δ/10} (so N ∤ n, and n < p except for O(L) primes), q := N·Q(n) ≤ 16NL⁴,
φ_{p,n}(t) := χ_N(−D_t)F_N(D_t)β_{D_t}(n) (0 if N | D_t). Then S^{R,head}_L(p) = −(1/(2π√p d))Σ_{n≤L}n^{−1}Σ_t W_p(t)φ_{p,n}(t).
**H1 (periodicity; proved).** By T2/T4, φ_{p,n} depends only on t mod q, and on p only through p mod q; |φ| ≤ 2. By CRT
φ = φ^{(N)}⊗φ^{(Q)} with φ^{(N)}(x) = (χ + χ²/N)/(1 − N^{−2}), χ = χ_N(x² − 4p²) (χ² = 1 − 1_{x≡±2p}), and |φ^{(Q)}| ≤ 1.
**H2 (Poisson; quoted standard).** W_p is continuous, compactly supported, of bounded variation, and Ŵ_p ∈ L¹ (H4), so
Σ_t W_p(t)φ(t) = Σ_{j∈ℤ}φ̂(−j)Ŵ_p(j/q), φ̂(j) := q^{−1}Σ_{x mod q}φ(x)e(−jx/q) = φ̂^{(N)}(jQ̄)φ̂^{(Q)}(jN̄) (Q̄Q ≡ 1 (N),
N̄N ≡ 1 (Q)), |φ̂^{(Q)}| ≤ 1.
**H3 (N-adic coefficients; proved + quoted).** φ̂^{(N)}(0) = −2/(N²−1) (Lemma M(i)). For i ≢ 0 (N),
|φ̂^{(N)}(i)| ≤ (1 − N^{−2})^{−1}(2/√N + 2/N²) ≤ 3/√N by Weil's bound |Σ_x χ_N(x² − c)e(ix/N)| ≤ 2√N (N ∤ c) (Weil; e.g.
Iwaniec–Kowalski Thm 11.23, not opened; checked F2, F3). Averaged over b ∈ (ℤ/N)^× (p ≡ b), for i ≢ 0 (proved by the
substitution x = 2bs and Σ_b e(2bi(1∓s)/N) = N[s ≡ ±1] − 1):
avg_b φ̂^{(N)}_b(i)e(±2bi/N) = [1/(N(N−1)) − (N−2)/(N²(N−1))]/(1 − N^{−2}) = 2/((N−1)(N² − 1)) (checked F4: N³|avg| = 2.020).
**H4 (Bessel; quoted, checked F1).** Ŵ_p(ξ) = 4p²ĝ(2pξ), ĝ(η) = ∫_{−1}^{1}U_{k−2}(u)√(1−u²)e(−ηu)du =
(k−1)(−1)^{k/2−1}J_{k−1}(2πη)/(2η) (Gegenbauer's integral, DLMF §18.17, not opened). With |J_ν(x)| ≪_ν x^{−1/2}:
|Ŵ_p(ξ)| ≪_k p^{1/2}|ξ|^{−3/2}; and for x ≥ 1, J_ν(x) = Re(e^{ix}h_ν(x)) with |h_ν| ≪ x^{−1/2}, |h′_ν| ≪ x^{−3/2} (Hankel
expansion, DLMF §10.17, not opened). So for p ≥ q, Ŵ_p(j/q) = Σ_± e(±2pj/q)A^±_j(p), sup|A| + Var_J A ≪ (Y+1)^{1/2}N(q/|j|)^{3/2}.
**H5 (j = 0; proved).** φ̂(0) = (−2/(N²−1))e′(n), e′(n) := mean over t mod Q(n) of β_{D_t}(n) = Π_{ℓ^j∥n}E′_ℓ[c_{ℓ,j}],
independent of p (t ↦ 2pt multiplies D by the unit square 4p²). |e′(ℓ)| ≤ 1/ℓ + 2/ℓ², |e′(ℓ^j)| ≤ 1, so Σ_n|e′(n)|n^{−σ} < ∞
for σ > ½ and Σ_{n≤L}e′(n)/n = Π_{ℓ≠N,p}E′_ℓ + O(L^{−1/2+ε}) = 1 + O(L^{−1/2+ε}) (Lemma M(ii); e′(n) := 0 for p | n).
**H6 (j ≠ 0, N | j; proved).** |φ̂(−j)| ≤ 2/(N²−1) and Σ_{j′≠0}|Ŵ_p(j′/Q)| ≪ p^{1/2}Q^{3/2}: contribution ≪ Q^{3/2}N^{−3} per
prime, no average needed (this includes the multiples of q, i.e. the Euler–Maclaurin corrections).
**H7 (j ≢ 0 (N): average over p; proved given fact 2047).** For p ≥ q the summand is Ψ^±_j(p mod q)Γ^±_j(p) with
Ψ^±_j(b) := φ̂_{b,n}(−j)e(±2bj/q) (well defined mod q), Γ^±_j := A^±_j/(2π√p d), ‖Γ_j‖ ≪ N^{−1}(q/|j|)^{3/2}. GRH for all
χ mod q (template T6 quote, Davenport ch. 20, not opened): Σ_{p∈J}log p χ(p)Γ(p) = δ_{χ₀}∫_JΓ + O(‖Γ‖N√(Y+1)log²N); with
Ψ = Σ_χΨ̂(χ)χ on units and Σ_χ|Ψ̂(χ)| ≤ φ(q)^{1/2}sup|Ψ| ≤ 3√Q (Parseval, H3), and the CRT mean |Ψ̄| ≤ 2/((N−1)(N²−1)) (H3):
Σ_{p∈J, p≥q}log p·Ψ·Γ ≪ N^{−3}‖Γ‖hN² + ‖Γ‖√Q·N log²N. Summed over j: ≪ N^{−5/2}Q^{3/2}hN² + N^{3/2}Q²log²N. Primes p < q:
Σ_j|Ŵ_p(j/q)| ≪ pq, so they contribute ≪ Σ_{p<q}log p·√p Q/√N ≪ NQ^{5/2}.
**H8 (collect; proved).** Summing H5–H7 over n ≤ L with weight 1/n (Q(n) ≤ 16n⁴) and dividing by Σ_{p∈J}log p ≍ hN²:
error ≪ L⁶N^{−5/2} + L⁸log²N/(h√N) + L^{10}/(hN) + L^{−1/2+ε} ≪ N^{−δ/25} for h ≥ N^{−1/2+δ}. The j = 0 term is S^R_loc(p)
(1 + O(L^{−1/2+ε})), and ⟨y_p^{3/2}⟩_J = (1/h)∫y^{3/2} + O(log²N/(hN)). ∎

### 5.3 Proof of Theorem C (proved)
S_N = S^I + S^H + S^{CN} + S^{R,head}_L + S^{R,tail}_L exactly; Theorem A and Proposition H evaluate all but the last
(h ≥ N^{−1/2+δ} ≥ N^{−1+δ}); the sum of their main terms is (1/h)∫M_k. ∎

### 5.4 Hypothesis T: where exactly the proof stops (NOT closed)
Needed: Σ_{p∈J}log p Σ_t W_p(t)χ_N(−D_t)τ_L(D_t) = o(hN⁴), over a family of ≍ hN⁴/log N pairs (p, t) of weight W ≍ N².
Equivalently the twisted tail must average to o(N^{−2}). What the available tools give (bounds computed here, used in no
proof):
- **Pointwise GRH** (Perron + Lindelöf on GRH; quoted standard): τ_L(D) ≪ L^{−1/2+ε}N^ε, so ⟨S^{R,tail}_L⟩ ≪ y^{3/2}N²L^{−1/2+ε}:
  needs L ≥ N^{4+ε}. But the head method (H7) controls only moduli with L⁸ ≪ h√N, i.e. n ≤ N^{1/16}: the range
  N^{1/16} < n < N^{4+ε} is not covered by either.
- **Quadratic large sieve** (Heath-Brown 1995; not opened) on a dyadic block n ~ X after Cauchy–Schwarz in n:
  ⟨·⟩ ≪ h^{−1/2}(1 + N²X^{−1/2})N^ε, never o(1): it keeps the cancellation in D but loses the one in n.
- **GRH for each L(s, Sym²f)** (not available; IK Thm 5.15, not opened) gives only ⟨S_N⟩_J ≪ (log N)²/h, the right size
  but no asymptotic: the statement is a family statement, not a per-form one.
- The missing input is therefore a twisted first moment of L(1, χ_{t²−4p²}) over the sparse family {(p, t)} with
  square-root cancellation jointly in the discriminant and in the Dirichlet variable n, in the range N^{1/16} < n < N^{4+ε}.
  (Heuristic, not proved: after Poisson in t these dual terms become Salié/Kloosterman-type sums to moduli nN, i.e. the
  "N ∤ c/N" Kloosterman terms that derivation.md §5 names on the Petersson side.)
**Numerical evidence (checked):** at N = 211 (k = 2, 12,454 primes p ≤ 3N²) the exact twisted part has coefficient
1.004 ± 0.011 against S^R_loc (E3, χ² 10.8/12); at N = 139 (k = 4) S^R against 0: χ² 9.6/12; the per-prime variance of S_N
is the Sato–Tate size p/d (ratio 0.881, E5) and only 3.8 % of the random-sign variance of the t-sum: the cancellation that
Hypothesis T asks for is present, and it is spectral (Tr T_{p²} is a sum of d ≍ N eigenvalues bounded by Deligne), far
beyond what the arithmetic of the single discriminants shows. review60 reports 0.9994 ± 0.0027 over 22 levels (not
recomputed here).

## 6. Steps (e)–(f): polynomial weights and assembly (proved)
k > 2 enters only through U_{k−2}: in S^{CN} through G_m (bounded variation, ‖G_m‖ ≪_k √(Y+1); Lemma K applies), in S^R
through W_p (H4: ∫W_p = 0 for k ≥ 4 since ∫_0^π sin((k−1)θ)sin θ dθ = (π/2)[k=2]; the decay bound holds for every k).
Normalisation: d = (k−1)(N+1)/12 + O(k) for k ∈ K (template step (e), F3), so ρ_N = 1 + O(1/N); all four terms carry the
same factor ρ_N, which tends to 1. For even k ∉ K the new space is S_k(N) minus two copies of S_k(1): the extra term
−√p d^{−1}Σ_{g∈S_k(1)}λ_g(p²) needs the prime number theorem for Sym²g in windows of length hN² at x ≍ N² (not done).

## 7. What remains conditional

| statement | status | hypothesis |
|---|---|---|
| (a) trace identity for Tr T_{p²}, (b) Lemma B | **proved** (formula quoted, read: Kaplan–Petrow ESTFKL) | none |
| Lemma L (local means = ν), Lemma M (−2/(N²−1); E′_ℓ = 1) | **proved** | none |
| Theorem A: class-number + hyperbolic + identity parts, windows h ≥ N^{−1+δ}, error (log N)^{−A} | **proved** | fact 2047: real χ of conductor ≤ 4(Y+1)²N²; all χ mod q ≤ 16N^{2δ/5} |
| Proposition H: head of R = 24[k=2]y^{3/2}, windows h ≥ N^{−1/2+δ}, error N^{−δ/25} | **proved** | fact 2047: all χ mod q ≤ 16N^{1+2δ/5}; Weil bound, Bessel facts (quoted) |
| Hypothesis T: tail of R → 0 (moduli N^{1/16} < n < N^{4+ε}) | **NOT closed** | — |
| Theorem C: full statement ⇔ Hypothesis T | **proved** | fact 2047 |
| even k ∉ K; harmonic weights; windows h < N^{−1/2+δ} for the head; single primes | not done / not claimed | — |
| quoted, not opened | Davenport ch. 20 (via T3, T6), Cox Cor. 7.28 (via T2), Weil/IK 11.23, DLMF §18.17 and §10.17, Heath-Brown 1995, IK 5.15 | checked numerically where a check exists (F1–F3) |

## 8. Numerical checks (`sym2_check.py`; output `sym2_check_out.txt`)
B1–B2 own Hurwitz table = PARI; A1–A6 trace identity, decomposition, reduction, Gram route, dimensions; C1–C7 local means
(Lemma L, N-adic, twisted −2/(N²−1), E′_ℓ = 1), ν and M_k against derivation.md; F1–F4 Gegenbauer integral, Weil bound,
N-adic Fourier coefficients and their p-average; D1 the class-number key lemma at N = 10007; E0–E5 the one-level test of the
full statement at N = 211 (k = 2) and N = 139 (k = 4) for every prime p ≤ 3N², via the Gram route (independent of the
decomposition): χ² on 12 windows of width ¼: model 10.7 (zero 22.8, no hyperbolic term 204, ν ≡ 1 99.4, no R main term
8465) for k = 2; model 9.3 (zero 27.6, no hyperbolic 44.0, ν ≡ 1 29.4) for k = 4. At one level the test separates the model
from zero only at χ² ≈ 23–28: the per-prime noise is the Sato–Tate size √(p/d). The trace code covers k = 2, 4 only, so
k = 6–14 are proved but not checked. Final line of `sym2_check_out.txt`: **51 checks, 0 failed** (run time 284 s).

## 9. Sources
- Kaplan–Petrow, arXiv:1510.03980: **read** the statement of Theorem ESTFKL and the definitions of U_{k−2}(t, q), h_w, H
  (downloaded source, `dl_kp/`).
- Template mprime_theorem.md; derivation.md; review60/review_159_161.md: **read**.
- Kundu–Müller arXiv:2507.00738 §2, Burrin et al. arXiv:2603.25564 trace set-up, Lowry-Duda arXiv:2506.01640 §4–5: **read**
  (local copies): they use T_n∘W_N (Skoruppa–Zagier) or level 1; none evaluates the unweighted Σ_f λ_f(p²) at one prime level.
- Knightly–Li, Davenport, Iwaniec–Kowalski, Cox, DLMF, Heath-Brown (quadratic large sieve, Acta Arith. 72 (1995)): **not opened**.
