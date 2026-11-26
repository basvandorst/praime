# The murmuration density of one prime level with natural weights: a theorem given fact 2047

Folder: `scratchpad/mprime_thm/` (nothing in the repository was changed). Check script: `mprime_check.py`, output
`mprime_check_out.txt`. Labels: **proved** (written out here), **quoted** (exact reference; "read" = I opened the text,
"not opened" = cited from memory or from an earlier note), **checked** (a numbered assertion of the script),
**NOT closed** (a gap).

## 0. Verdict

**THEOREM under one named hypothesis: GRH for Dirichlet L-functions** (fact 2047), used for exactly two families:
(i) the real primitive characters χ_{Δ₀} of conductor ≤ 4(Y+1)N² attached to the discriminants −D with D = 4pN − r²N²,
and (ii) all characters to moduli ≤ 16N^{4δ/5} (primes in progressions). No GRH for L(s, f) or L(s, Sym²f) is used.
No average over the level: N is one fixed prime; the average is over primes p in a window of length ≥ N^{1/2+δ}.
All even k (k = 2 has the extra term −12y).

**The proof does not go through the Petersson formula.** It uses the Skoruppa–Zagier trace formula for T_p∘W_N at the
fixed prime level N (class numbers H(4pN − r²N²)) and the mean of L(1, χ_{−D}) over the primes p, local factor by local
factor. On the Petersson side, step (a) (harmonic weights with Hecke twists, uniform in ℓ ≤ N^η) is proved given
fact 2047, but step (b) (harmonic → natural weights) **stops** as a self-standing estimate: the available device
(Iwaniec–Luo–Sarnak §2, Lemma 2.12) needs GRH for L(s, Sym²f) (and, for prime-supported sequences, GRH for L(s, f)).
Once the theorem is proved by the trace route, step (b) follows on window averages (Corollary 3): the two routes agree.

## 1. Statements

**Setting.** N ≥ 5 prime, k ≥ 2 even, H_k(N) the newforms of S_k(Γ₀(N)) (all of S_k(Γ₀(N)) when S_k(SL₂(ℤ)) = 0,
i.e. k ∈ K := {2, 4, 6, 8, 10, 14}), a_f(1) = 1, λ_f(n) = a_f(n)n^{−(k−1)/2}, ε_f the root number
(ε_f = −i^k λ_f(N)√N: ILS (118) for k = 2, read; for k = 4 checked jointly with the trace formula, F2). Put
E_N(p) := |H_k(N)|^{−1} Σ_f ε_f √p λ_f(p) (for k = 2: (1/dim)Σ_f ε_f a_f(p)). For y₀ ≥ 0, h > 0, J := (y₀N, (y₀+h)N] and
⟨X⟩_J := Σ_{p∈J, p≠N} log p·X(p) / Σ_{p∈J, p≠N} log p. Let u_ℓ := 1/((ℓ−1)²(ℓ+1)) and

    ν′(r) := Π_{ℓ∤r}(1 − u_ℓ)·Π_{ℓ|r} ℓ²/(ℓ² − 1)      (= Σ_{δ|r, δ sqfree} c(δ) of densities.py; check C1),
    M_k(y) := (1/(k−1))·[2π√y + (−1)^{k/2−1}(12/π) Σ_{1≤r<2√y} ν′(r) U_{k−2}(r/(2√y)) √(4y − r²)] − 12y·[k = 2].

M_2 = M_prime of the task. M_k is continuous on [0, ∞) (each r-term vanishes at y = r²/4).

**Theorem 2 (natural weights, one prime level).** Assume GRH for Dirichlet L-functions. Fix k, Y ≥ 1, δ ∈ (0, ½), A ≥ 1.
Uniformly for N prime, 0 ≤ y₀ ≤ Y and N^{−1/2+δ} ≤ h ≤ 1,

    ⟨E_N⟩_J = (1/h) ∫_{y₀}^{y₀+h} M_k(y) dy + O_{k,Y,δ,A}((log N)^{−A}).

Hence (i) for h = h(N) → 0 with h ≥ N^{−1/2+δ}, ⟨E_N⟩_J → M_k(y₀) uniformly in y₀ ∈ [0, Y]; (ii) for a fixed
Φ ∈ C_c([0, ∞)), Σ_p log p·Φ(p/N)E_N(p)/Σ_p log p·Φ(p/N) → ∫ΦM_k/∫Φ (approximate Φ uniformly by step functions);
(iii) the same with weight 1 instead of log p (partial summation). Nothing is claimed for single primes p.

**Theorem 1 (step (a): harmonic weights with a Hecke twist).** k ∈ K, Σ^h_f the Petersson-weighted sum
(ω_f = Γ(k−1)/((4π)^{k−1}⟨f,f⟩)). Assume GRH for L(s, χ), χ mod m, m ≤ N^θ. For 0 < η < θ < ¼, 1 ≤ ℓ ≤ N^η,
0 ≤ y₀ ≤ Y, N^{−1/2+2θ} ≤ h ≤ 1 and g of bounded variation (‖g‖ := sup|g| + Var g on [y₀, y₀+h]):

    ⟨g(p/N)·√N Σ^h ε_f λ_f(ℓ)λ_f(p)⟩_J = (1/h)∫_{y₀}^{y₀+h} g(v)·2π Σ_{m≥1} (μ(m)c_m(ℓ)/(m φ(m))) J_{k−1}(4π√(ℓv)/m) dv
                                         + O_{k,Y}(‖g‖(N^{(η−θ)/2} + N^{−θ/2}) log³N).

With g(v) = √v: ⟨Σ^h ε_f λ_f(ℓ)√p λ_f(p)⟩_J → (1/h)∫M^h_ℓ, M^h_ℓ(v) = 2π√v Σ_m μ(m)c_m(ℓ)/(mφ(m)) J_{k−1}(4π√(ℓv)/m).
Unconditional variant: for ℓ ≤ (log N)^B and fixed h, the same with error o(1) (Siegel–Walfisz in place of GRH).

**Corollary 3 (step (b) on window averages).** Given fact 2047, k ∈ K, X := N^{η/2} with η, θ, h as in Theorem 1:
Y′_f := Σ_{(t,N)=1} λ_f(t²)/t, Y′_f(X) := Σ_{t≤X, (t,N)=1} λ_f(t²)/t, then
⟨|H_k(N)|^{−1}Σ_f ε_f √p λ_f(p)·(1 − Y′_f(X)/Y′_f)⟩_J = O_A((log N)^{−A}).

**Proposition 4 (step (d): Bessel form = closed form; unconditional).** For y > 0 and every even k,
M_k(y) = 2π Σ_{d sqfree} Σ_{s≥1, (s,d)=1} A((s/d)²)·(s/d)√y·J_{k−1}(4πs√y/d), A((s/d)²) = (6/π²)s^{−2}Π_{p|d} p²/((p−1)²(p+1)),
the double series converging absolutely.

## 2. Proof of Theorem 2 (the trace route)

**T1 (quoted; checked).** For N squarefree and a prime p ∤ N,
Σ_f ε_f√pλ_f(p) = ½H(4pN) + (−1)^{k/2−1}Σ_{1≤r<2√(p/N)} U_{k−2}(r√N/(2√p)) H(4pN − r²N²) − δ_{k,2}(p + 1),
H = Hurwitz class number. Source: Skoruppa–Zagier 1988, §2 (5), (7) (not opened), in the form of "On Murmurations",
arXiv:2609.29631, Lemma `lem:skoruppa-zagier-trace-formula` (read; it notes that Zubrilina 2310.07681 §2 prints U_{k−2}
outside the r-sum). For N prime the oldforms from level 1 do not contribute (W_N has trace 0 on each pair f(z), f(Nz)).
r²N = 4p is impossible. **Checked:** F1 (k = 2, ε_f = a_f(N): 918 pairs (N, p), levels 11–163) and F2 (k = 4,
ε_f = −a_f(N)/N: 376 pairs), exact rational equality with PARI trace forms.

**T2 (quoted classical; checked).** For D > 0, D ≡ 0, 3 (4), write −D = Δ₀f₀² (Δ₀ fundamental, χ := (Δ₀/·)). The Hurwitz sum
H(D) = Σ_{g|f₀} h′(Δ₀g²), the class number of orders h′(Δ₀g²) = h′(Δ₀)·g·Π_{ℓ|g}(1 − χ(ℓ)/ℓ) (Cox, *Primes of the form
x² + ny²*, Cor. 7.28; not opened) and Dirichlet's formula h′(Δ₀) = √|Δ₀|·L(1, χ)/π (h′ = 2h/w) give

    H(D) = (√D/π)·L*(D),  L*(D) := L(1, χ)·Π_{ℓ^a ∥ f₀} Q_ℓ,  Q_ℓ := (1 − χ(ℓ)/ℓ)Σ_{j<a}ℓ^{−j} + ℓ^{−a} ∈ [1 − 1/ℓ, 1 + 2/ℓ].

**Checked:** E1 (brute-force reduced forms = PARI qfbhclassno, D ≤ 3000), E2 (H(D) = h′(Δ₀)f₀ΠQ_ℓ exactly, D ≤ 3000).
Local form (proved): F_ℓ(D) := Q_ℓ/(1 − χ(ℓ)/ℓ) = Σ_{j<a}ℓ^{−j} + ℓ^{−a}/(1 − χ(ℓ)/ℓ) = Σ_{j≥0} c_{ℓ,j}(D)ℓ^{−j} with
c_{ℓ,j} = 1 (j < a), χ(ℓ)^{j−a} (j ≥ a), so |c_{ℓ,j}| ≤ 1; c_{ℓ,j}(D) depends only on D mod ℓ^{2j+2} (ℓ odd) and mod 2^{2j+6}
(ℓ = 2): if ℓ^{2j+2} | D (2^{2j+4} | D) then a > j and c = 1; otherwise v_ℓ(D) ≤ 2j+1 (≤ 2j+3), and v_ℓ(D), a and χ(ℓ) are
read off D mod ℓ^{v+1} (mod 2^{v+3}). F_N(4pN − r²N²) = 1 since N ∥ D for p ≠ N.

**T3 (proved, given GRH for L(s, χ_{Δ₀})).** Quoted consequence of GRH: θ(x, χ) := Σ_{ℓ≤x}χ(ℓ)log ℓ ≪ √x log²(qx) for χ ≠ χ₀
mod q (Davenport, *Multiplicative Number Theory*, ch. 20, (13)–(14) for ψ; not opened; the same citation as the refereed
family_proof, Theorem B). By partial summation Σ_ℓ χ(ℓ)ℓ^{−s} converges locally uniformly on Re s > ½, so it equals
log L(s, χ) − Σ_{ℓ, j≥2} χ(ℓ)^j/(jℓ^{js}) there (identity on Re s > 1; L ≠ 0 on Re s > ½), and Σ_{ℓ>z}χ(ℓ)/ℓ ≪ z^{−1/2}log²(qz).
Hence L(1, χ) = Π_{ℓ≤z}(1 − χ(ℓ)/ℓ)^{−1}·exp(O(z^{−1/2}log²(qz) + 1/z)) (the j ≥ 2 terms give 1/z). With ω(f₀) ≤ log D and Q_ℓ = 1 + O(1/ℓ),
Π_{ℓ|f₀, ℓ>z}Q_ℓ = exp(O(log D/z)). So, with P_z(D) := Π_{ℓ≤z}F_ℓ(D) ≪ (log z)³, uniformly for D ≤ 4(Y+1)N²:

    L*(D) = P_z(D) + O((log z)³·[z^{−1/2}log²(Nz) + log N/z]).

**T4 (proved).** P_z(D) = Σ_{P⁺(n)≤z} β_D(n)/n, β_D(n) := Π_{ℓ^j∥n} c_{ℓ,j}(D), |β_D(n)| ≤ 1, and β_D(n) depends only on
D mod Q(n) := Π_{ℓ^j∥n} ℓ^{2j+2}·16^{[2|n]} ≤ 16n⁴. Rankin with σ = 1/log z: Σ_{n>L, P⁺(n)≤z} 1/n ≤ L^{−σ}Π_{ℓ≤z}(1 − ℓ^{σ−1})^{−1}
≪ (log z)·exp(−log L/log z) (since ℓ^σ ≤ 1 + eσ log ℓ).

**T5 (local means; proved; checked E3).** Let E_ℓ(r) := ∫_{ℤ_ℓ^×} F_ℓ(4pN − r²N²) dp (Haar probability on the units).
- ℓ = N: E_N = 1 (T2).
- ℓ ∉ {2, N}, ℓ | r (r = 0 included): −D ≡ −4pN is a unit, χ(ℓ) = (−pN/ℓ) is ±1 with probability ½ each, a = 0:
  E_ℓ = ½[(1 − 1/ℓ)^{−1} + (1 + 1/ℓ)^{−1}] = ℓ²/(ℓ² − 1).
- ℓ ∉ {2, N}, ℓ ∤ r: u := r²N − 4p is Haar-uniform on ℤ_ℓ ∖ (r²N + ℓℤ_ℓ), with density ℓ/(ℓ−1), and −D = Nu. Over all of ℤ_ℓ
  (x := 1/ℓ; P(v_ℓ = 2a) = (1−x)x^{2a} with χ = ±1 equally, P(v_ℓ = 2a+1) = (1−x)x^{2a+1} with χ = 0) the mean of F_ℓ is
  Σ_a (1−x)x^{2a}[S_a + x^a/(1−x²)] + (1−x)x^{2a+1}[S_a + x^a] = 1/(1 − x²), S_a := (1 − x^a)/(1 − x). On the removed disk −D is a
  unit square, F_ℓ = (1 − 1/ℓ)^{−1}, measure 1/ℓ. So E_ℓ = (ℓ/(ℓ−1))[ℓ²/(ℓ²−1) − 1/(ℓ−1)] = ℓ(ℓ²−ℓ−1)/((ℓ−1)²(ℓ+1)) = 1 − u_ℓ (C2).
- ℓ = 2, r odd: −D ≡ 1 − 4 ≡ 5 (mod 8), F₂ = (1 + ½)^{−1} = 2/3. ℓ = 2, r = 2r′: −D = 4Δ′; r′ even: Δ′ uniform on ℤ₂^×, F₂ = 1
  (Δ′ ≡ 3 (4)), 2 (Δ′ ≡ 1 (8)), 4/3 (Δ′ ≡ 5 (8)): mean 4/3. r′ odd: −D = 8w, w uniform on ℤ₂; F₂(8·odd) = 1 and F₂(4Δ₁) = 1 + F₂(Δ₁)/2
  give the mean ½ + ½(1 + I/2) with I := mean of F₂(4w′) = 7/6 + I/8, I = 4/3; mean 4/3.
- Hence E(r) := Π_{ℓ≠N}E_ℓ(r) = ν′(r)·(1 + O(N^{−2})) for 1 ≤ r < N, and E(0) = (4/3)Π_{ℓ odd, ≠N} ℓ²/(ℓ²−1) = ζ(2)(1 − N^{−2}).
**Checked E3:** 208 cases (ℓ = 2, 3, 5, 7; N = 10007, 100003, 13 and N = ℓ; r = 0…12), exact enumeration of p mod ℓ^J with a
rigorous enclosure of the undetermined residues (width ≤ 1.5·10⁻⁵).

**T6 (primes in progressions; quoted GRH consequence).** θ(x; q, b) = x/φ(q) + O(√x log²(qx)) for (b, q) = 1 (Davenport ch. 20;
not opened). For G of bounded variation on [0, Y+1], ‖G‖ := sup|G| + Var G, partial summation gives
Σ_{p∈J, p≡b (q)} log p·G(p/N) = φ(q)^{−1}∫_J G(u/N)du + O(‖G‖√N log²(qN)). For n ≤ L the residue p mod Q(n) fixes D_r(p) mod Q(n), so
Σ_{p∈J} log p·G(p/N)β_{D_r(p)}(n) = e_r(n)∫_J G(u/N)du + O(‖G‖Q(n)√N log²N + L log N·sup|G|), with
e_r(n) := φ(Q)^{−1}Σ_{b∈(ℤ/Q)^×}β_{D_r(b)}(n) = Π_{ℓ^j∥n} E_ℓ[c_{ℓ,j}] (CRT), and Σ_{P⁺(n)≤z} e_r(n)/n = Π_{ℓ≤z}E_ℓ(r).

**Key lemma (proved).** For 0 ≤ r < 2√(Y+1), uniformly in J and G:
Σ_{p∈J} log p·G(p/N)·L*(D_r(p)) = E(r)∫_J G(u/N)du + O_A(‖G‖ hN (log N)^{−A}), D_r(p) := 4pN − r²N² (D₀ = 4pN).
*Proof.* T3 for each p (Σ_{p∈J}log p ≪ hN); T4: keep n ≤ L; T6 for each n ≤ L (total error ≪ ‖G‖L⁴√N log²N); the main terms
sum to Π_{ℓ≤z}E_ℓ(r) + O(Rankin) = E(r)(1 + O(1/z)) + O(Rankin) (E_ℓ = 1 + O(ℓ^{−2}) for ℓ > z > 2√(Y+1)). Choose L := N^{δ/5},
log z := log L/((A+1) log log N): Rankin ≪ (log N)^{−A}; the T3 error is ≪ exp(−c log N/(A log log N)); the T6 error is
≤ hN·N^{−δ/5}log²N since hN ≥ N^{1/2+δ}. ∎

**Assembly (proved).** Step (e): |H_k(N)| = (k−1)N/12·(1 + O_k(1/N)). With y_p := p/N, √(4pN) = 2N√y_p and √D_r = N√(4y_p − r²),
T1 and T2 give
E_N(p) = (1 + O(1/N))·{(12/((k−1)π))[√y_p L*(D₀) + (−1)^{k/2−1}Σ_r U_{k−2}(r/(2√y_p))√(4y_p − r²) L*(D_r)] − δ_{k,2}·12(y_p + 1/N)}.
Apply the key lemma with G₀(v) = √v and G_r(v) = 2√v·sin((k−1)θ_r(v)) for v > r²/4 (cos θ_r = r/(2√v); this equals
U_{k−2}(r/(2√v))√(4v − r²)), G_r = 0 otherwise: continuous, ‖G_r‖ ≪_k √(Y+1). Σ_{p∈J}log p·p/N = ∫_J u/N du + O(√N log²N).
With E(0)·12/π = 2π and E(r) = ν′(r)(1 + O(N^{−2})), the window mean is (1/h)∫M_k + O((log N)^{−A}). ∎ (Theorem 2)

## 3. Step (a): proof of Theorem 1 (proved, given GRH for characters mod m ≤ N^θ)

1. *Petersson with root number (quoted + proved).* For (n, N) = 1: ε_fλ_f(n) = −i^k√N λ_f(Nn), and Petersson (ILS (125), read;
   Iwaniec–Kowalski Prop. 14.5, not opened) gives √N Σ^h ε_fλ_f(n) = −2π Σ_{m≥1} S(Nn, 1; Nm) m^{−1} J_{k−1}(4π√(n/N)/m)
   (family_proof Lemma 2, refereed). S(Nn, 1; Nm) = −S(n, N̄; m) for (m, N) = 1 (ILS (157); checked A1) and S(Nn, 1; Nm) = 0 for
   N | m (proved: x = x₀ + N^{j−1}t kills the sum; checked A2). For p ∈ J with p > ℓ: λ_f(ℓ)λ_f(p) = λ_f(ℓp), so
   √N Σ^h ε_fλ_f(ℓ)λ_f(p) = 2π Σ_{(m,N)=1} S(ℓp, N̄; m) m^{−1} J_{k−1}(4π√(ℓp/N)/m). (Primes p ≤ ℓ cost O(ℓ log N) = o(hN).)
2. *Tail m > M₀ := N^θ (proved).* Weil |S(a, b; m)| ≤ τ(m)(a, b, m)^{1/2}m^{1/2} (quoted, IK Cor. 11.12, not opened; checked A4;
   here (ℓp, N̄, m) = 1) and |J_{k−1}(x)| ≤ x/2 (DLMF 10.14.1, 10.14.4, not opened; checked D0): ≪ √(ℓ(Y+1)) M₀^{−1/2}log M₀ per p.
3. *Head m ≤ M₀ (proved).* p > M₀ for all but O(M₀) primes, so K(p) := S(ℓp, N̄; m) is a function on (ℤ/m)^×, |K| ≤ τ(m)√m.
   T6 with G_m(u) := g(u/N)J_{k−1}(4π√(ℓu/N)/m) (‖G_m‖ ≪ ‖g‖(1 + √(ℓ(Y+1))/m), |J′| ≤ 1): Σ_{p∈J}log p K(p)G_m(p) =
   (Σ_a K(a)/φ(m))∫_J G_m + O(φ(m)τ(m)√m·‖G_m‖√N log²N), and Σ_{a∈(ℤ/m)^×}S(ℓa, b; m) = μ(m)c_m(ℓ) for (b, m) = 1 (family_proof
   Lemma 3(ii); checked A3). Summed over m ≤ M₀ with 2π/m: error ≪ ‖g‖√N M₀^{3/2}log³N.
4. *Collect.* Divide by Σ_{p∈J}log p = hN(1 + o(1)); extending the main m-sum to all m costs ≪ hN√ℓ τ(ℓ)M₀^{−2+ε}. Error
   ≪ ‖g‖(√ℓ M₀^{−1/2} + M₀^{3/2}N^{−1/2}h^{−1})log³N ≤ ‖g‖(N^{(η−θ)/2} + N^{−θ/2})log³N. ∎ **The uniformity in ℓ ≤ N^η that the
   gate (3.4) and family_proof (GAP 4) left open is now proved** (for this prime-side statement).
   Unconditional variant: M₀ = (log N)^{2B+1}, ℓ ≤ (log N)^B and Siegel–Walfisz (error N exp(−c√log N) per class): o(1).

## 4. Step (b): harmonic → natural weights

1. *Exact conversion (quoted + proved).* ω_f = 2π²/((k−1)N L(1, Sym²f)) (ILS Lemma 2.5 (2.36) and (124)–(127) for k = 2, read;
   family_proof Lemma 1, checked there to 10^{−11}) and Σ_t λ_f(t²)t^{−s} = L(s, Sym²f)/ζ^{(N)}(2s) with λ_f(N²) = 1/N give
   L(1, Sym²f) = ζ(2)Y′_f, hence **ω_f = 12/((k−1)N Y′_f)** and E_N(p) = ((k−1)N/(12|H_k(N)|))·Σ^h ε_f√pλ_f(p)Y′_f.
   Y′_f is a conditionally convergent series (ILS remark after Prop. 2.11).
2. *What the literature proves (read: ILS Prop. 2.11, (2.62)–(2.65), Lemma 2.12, Prop. 2.13).* ILS split Σ_m m^{−1}Δ_{k,M}(m², n)
   at m ≤ Y and bound the complement with r_f(KM; Y) := Σ_{m>Y, (m,KM)=1} λ_f(m²)/m ≪ Y^{−1/2}(kMY)^ε, which they derive
   "by the Riemann hypothesis for L(s, sym²(f))"; Lemma 2.12 also needs hypothesis (2.65) on the coefficient sequence, which for
   prime-supported sequences ILS take from GRH for L(s, f) ((2.66)–(2.67)). **Neither is a Dirichlet L-function: not covered by
   fact 2047, not used.** Kowalski–Michel's unconditional zero-density device (Duke 1999): **not opened.**
3. **Named gap (E_b), NOT closed by a direct argument:** ⟨|H_k(N)|^{−1}Σ_f ε_f√pλ_f(p)·r_f(X)/Y′_f⟩_J → 0 for X = N^{η/2},
   r_f(X) := Y′_f − Y′_f(X): the Sym² tail must be negligible against the window sum of ε_f a_f(p) on average over f.
4. *Corollary 3 closes (E_b) a posteriori (proved, given fact 2047).* By 1., |H|^{−1}Σ_f ε_f√pλ_f(p)Y′_f(X)/Y′_f =
   (1 + O(1/N))Σ_{t≤X,(t,N)=1} t^{−1}·Σ^h ε_fλ_f(t²)√pλ_f(p) (Hecke: p ∤ t). Theorem 1 at ℓ = t² ≤ N^η (uniform), summed with 1/t:
   (1/h)∫Σ_{t≤X} t^{−1}M^h_{t²} + O(N^{−c}log⁴N). The double series Σ_{t,m} converges absolutely (|c_m(t²)| ≤ (m, t²), |J| ≤ min(x/2,
   C_k x^{−1/2})), with tail Σ_{t>X} t^{−1}|M^h_{t²}(v)| ≪ X^{−1/2+ε} uniformly on [0, Y+1]; by (c) and (d) its full value is M_k.
   Subtracting Theorem 2 gives Corollary 3. So **the harmonic Hecke-twisted densities determine the natural density, as the gate
   claimed, and the truncation at t ≤ N^{η/2} is legitimate — but this is deduced from Theorem 2, not proved independently.**

## 5. Step (c): Hecke relation and the main-term algebra (proved; checked B1–B3)

λ_f(t²)λ_f(p) = λ_f(t²p) for p ∤ tN (multiplicativity). The natural Bessel series is
2π√y Σ_{t,m} μ(m)c_m(t²)/(t m φ(m)) J_{k−1}(4πt√y/m), (tm, N) = 1. Group by t/m = s/d in lowest terms, t = sg, m = dg. μ(dg) ≠ 0
forces d, g squarefree and coprime; for p | g, p | s²g² so c_p = p − 1; for p | d, p ∤ s²g² so c_p = −1; hence
**μ(dg)c_{dg}(s²g²)/φ(dg) = μ(g)/φ(d)** (B1, 2430 exact cases) and the term is (s d g²)^{−1}μ(g)/φ(d)·J_{k−1}(4πs√y/d). Summing g:
Σ_{(g,dN)=1} μ(g)/g² = (6/π²)Π_{p|d}(1 − p^{−2})^{−1}(1 − N^{−2})^{−1}. Writing the line as 2π A(ρ)√(ρy)J_{k−1}(4π√(ρy)), √ρ = s/d:
**A((s/d)²) = (6/π²)s^{−2}Π_{p|d} p²/((p−1)(p²−1)) = (6/π²)s^{−2}Π_{p|d} p²/((p−1)²(p+1))**, ×(1 − N^{−2})^{−1}; lines with N | s are
absent at level N and carry O(N^{−1/2}Y^{1/4}) in total. B2: the raw sum Σ_{t/m=s/d} μ(m)c_m(t²)/(φ(m)t²) (g ≤ 3000) matches for
13 pairs to 3·10⁻⁶; B3: d = 4 gives 0.

## 6. Step (d): Poisson in s (proved; checked D1–D3) = Proposition 4

With w(d) := Π_{p|d}p²/((p−1)²(p+1)), the Bessel form is (12/π)√y Σ_d (w(d)/d) Σ_{(s,d)=1} J_{k−1}(4πs√y/d)/s. Möbius:
Σ_{(s,d)=1} J(cs)/s = Σ_{e|d}(μ(e)/e)F_k(ce), F_k(x) := Σ_{s≥1} J_{k−1}(xs)/s. Poisson for g(u) := J_{k−1}(x|u|)/|u| (continuous,
|g| ≪ |u|^{−3/2}; ĝ continuous with compact support; Stein–Weiss Ch. VII, Cor. 2.6, not opened): ĝ(ξ) = 2∫₀^∞ J_ν(xu)cos(2πξu)du/u
= (2/ν)cos(ν arcsin(2π|ξ|/x)) for 2π|ξ| < x and 0 beyond (ν = k−1 odd: cos(νπ/2) = 0; Gradshteyn–Ryzhik 6.693.2, not opened),
and cos(ν arcsin t) = (−1)^{k/2−1}√(1−t²)U_{k−2}(t). So
F_k(x) = (1/(k−1))[1 + 2(−1)^{k/2−1}Σ_{1≤n<x/2π} √(1 − t_n²)U_{k−2}(t_n)] − (x/4)[k = 2], t_n = 2πn/x (D1: k = 2, 4, to 2·10⁻⁸).
Put x = 4π√y e/d, δ := d/e (squarefree), r := nδ. The n = 0 terms: Σ_d (w(d)/d)Σ_{e|d}μ(e)/e = Σ_d w(d)φ(d)/d² = Π_p p²/(p²−1) = ζ(2),
giving 2π√y/(k−1). The x/4 terms (k = 2): Σ_d (w(d)/d²)Σ_{e|d}μ(e) = 1, giving −12y. The n ≥ 1 terms: for fixed δ | r,
Σ_{(e,δ)=1} (w(eδ)/(eδ))(μ(e)/e) = (w(δ)/δ)Π_{p∤δ}(1 − u_p) = c(δ) (w(p)/p² = u_p), and Σ_{δ|r} c(δ) = ν′(r) (C1, C2). All
rearrangements are of absolutely convergent series (|F_k(x)| ≪ 1 + x; Σ_d w(d)τ(d)/d < ∞). D2: Möbius + Poisson per d, 24 cases
to 1·10⁻⁸; D3: the Bessel form (d ≤ 2·10⁴, s-sums exact) equals M_2 at y = 0.3, 1.3, 3.7, 7.9 to 3.8·10⁻⁴ (d-truncation). ∎

## 7. Step (e): bookkeeping (proved / quoted; checked F3)

S_k(SL₂(ℤ)) = 0 for k ∈ K, so S_k(Γ₀(N)) is new for N prime (for k = 2: S₂(Γ₀(N)) = S₂^new(N), there is nothing at level 1).
dim S₂(Γ₀(N)) = g₀(N) = (N+1)/12 − ν₂/4 − ν₃/3, ν₂ = 1 + (−1/N), ν₃ = 1 + (−3/N) (genus formula; Diamond–Shurman Thm 3.1.1, not
opened; F3: exact for all primes 5 ≤ N < 2000, |12 dim − N| ≤ 13). In general |H_k(N)| = (k−1)(N+1)/12 + O(k) (minus 2 dim S_k(1)),
so 12|H_k(N)|/((k−1)N) = 1 + O_k(1/N).

## 8. What remains conditional

| statement | status | hypothesis |
|---|---|---|
| Theorem 2: natural weights, one prime level, windows ≥ N^{1/2+δ}, all even k, error (log N)^{−A} | **proved** | GRH for Dirichlet L-functions only: χ_{Δ₀} real primitive, q ≤ 4(Y+1)N²; all χ mod q ≤ 16N^{4δ/5}. = fact 2047 |
| Theorem 1: harmonic, Hecke twist ℓ ≤ N^η, k ∈ K | **proved** | GRH for χ mod m ≤ N^θ (fact 2047); unconditional for ℓ ≤ (log N)^B, fixed h |
| Corollary 3: step (b) on window averages | **proved** | fact 2047 (through Theorems 1 and 2) |
| Proposition 4: Bessel form = closed form; step (c) algebra | **proved** | none |
| Step (b) as a direct estimate (E_b) | **NOT closed** | ILS Lemma 2.12 needs GRH for L(s, Sym²f) (+ GRH for L(s, f) for (2.65)): not available |
| Theorem 2 without GRH | **NOT closed** | would need L(1, χ_{−D(p)}) on average over primes p: Heath-Brown's quadratic large sieve over u = r²N − 4p plus a Bombieri–Vinogradov/BDH step (sketch only) |
| Windows shorter than N^{1/2}, single primes p | not claimed | — |
| Quoted, not opened | T1 (Skoruppa–Zagier; read in 2609.29631), Davenport ch. 20, Cox Cor. 7.28, Weil bound, Poisson/Stein–Weiss, GR 6.693.2, DLMF 10.14 | each checked numerically (F1–F2, E1–E2, A4, D0–D1) where a check exists |

## 9. Numerical checks (`mprime_check.py`, output in `mprime_check_out.txt`)

A1–A4 Kloosterman identities and the Weil bound (brute force); B1–B3 main-term algebra; C1–C3 Euler products (ν′, ζ(2), local
identities exact for p < 1000); D0–D3 Bessel bounds, Poisson (k = 2, 4), Möbius regrouping, Bessel form = closed form;
E1–E3 Hurwitz numbers, the product formula, local means E_ℓ(r); F1–F3 trace formula (k = 2, 4) and dimensions;
G1–G3 **one prime level, no average over N** (N = 10007, p ≤ 12N, h = ½: χ² 6.1 on 24 windows for M_2 against 31.9 for
Zubrilina's squarefree density; N = 100003, p ≤ 4N, h = ¼: χ² 2.1 on 16 against 194.5; k = 4: 3.4 vs 108.8 and 5.1 vs 1252.5;
means of πH(D_r)/√D_r over primes: z-scores |z| ≤ 0.64 against ζ(2), ν′(1), ν′(2), ν′(3)); the χ² are below the number of
windows because the residues of p mod small q are equidistributed better than random sampling (T6), so the per-prime scatter
overstates the error; H1–H2 the blind round-160 file (rms 0.058 against this note's M_2; M_2 = the registered model to 0.012).
Final line of `mprime_check_out.txt`: **28 checks, 0 failed** (run time about 2.5 min, PARI stack 300 MB).

## 10. Sources

- Zubrilina, arXiv:2310.07681 (v2 source): **read** Thm 1 (mainthm), §2 (trace formula setup), start of §3 (averages over N).
- "On Murmurations", arXiv:2609.29631: **read** abstract, §1 discussion of the misplaced U_{k−2}, Lemma (Skoruppa–Zagier).
- Skoruppa–Zagier, Invent. Math. 94 (1988), §2 (5), (7); Yamauchi 1973: **not opened**.
- Iwaniec–Luo–Sarnak, Publ. IHES 91 (2000), OCR of the published text: **read** §2 (Lemma 2.5, 2.7, Props. 2.8–2.11, (2.62)–(2.67),
  Lemma 2.12, Prop. 2.13); (118), (124)–(127), (157) as read for family_proof.
- Cowan arXiv:2408.12723 (abstract and Thm 1.4 statement, grep this session), Miller–Montague arXiv:0911.1830: read in the gate,
  not re-read. Cardi arXiv:2607.08982, arXiv:2609.00777: **abstract only** (squarefree levels averaged over N).
- Davenport ch. 20; Iwaniec–Kowalski Prop. 14.5, Cor. 11.12; Cox Cor. 7.28; Diamond–Shurman Thm 3.1.1; Stein–Weiss VII.2.6;
  Gradshteyn–Ryzhik 6.693.2; DLMF 10.14; Kowalski–Michel (Duke 1999): **not opened**.

**Prior art.** Every murmuration theorem with natural weights seen (Zubrilina; 2609.29631; Kundu–Müller as quoted there; Cardi
2607.08982; 2609.00777) averages over squarefree levels N in a short interval; the prime-level results seen (Cowan Thm 1.4 from
Miller–Montague) have harmonic weights. A natural-weight density for **one fixed prime level**, averaged over p only, was **not
found in print**; Theorem 2's method is Zubrilina's (trace formula + class-number averages) with the roles of N and p exchanged
(characters of the linear form r²N − 4p in the prime p instead of characters in N), which is why GRH for quadratic characters and
primes in progressions replace her Burgess/Pólya–Vinogradov step.


## Referee report (26 November 2026)

An independent referee (own code in [referee/](referee/)) re-derived the heart of the proof before reading it. **Verdict:
Theorem 2 is proved, given GRH for Dirichlet L-functions; no failing step, no hidden hypothesis.** Theorem 1 and Corollary 3
are sound (Corollary 3 is a subtraction, not circular, and not an independent proof of step (b)). Proposition 4 was not
recomputed.
- **Local means**, by exact enumeration with other code (p mod ℓ^J, N ≡ 1 and 7 mod 8): 2/3 and 4/3 at ℓ = 2; 15/16 and 9/8
  at ℓ = 3; 95/96 and 25/24 at ℓ = 5; F_N = 1 exactly. Their product is ν′(r).
- **The cut of the Euler product**: z = N^(δ/(5(A+1)·log log N)), expansion over z-smooth n, only n ≤ N^(δ/5) kept (moduli
  ≤ 16N^(4δ/5)), the tail by Rankin's trick for each discriminant separately: checked. With log z ≍ √(log N) both errors
  become exp(−c√(log N)), better than the stated (log N)^(−A).
- **The trace formula** T1 checked exactly, with root numbers from the functional equation, for every prime p ≤ 4N at
  k = 2 (N = 167, 173, 197), k = 4 (29, 31), k = 6 (17), k = 8 (19), k = 14 (11) and k = 12 (13, 17: old forms
  contribute 0).
- **Which GRH**: real primitive characters of conductor ≤ 4(Y+1)N² and all characters to moduli ≤ 16N^(4δ/5); no L(s, f),
  no L(s, Sym² f), no lower bound for L(1, χ).
- **One level, not used by the author** (N = 30011, all primes p ≤ 4N, 16 windows): χ² 3.2 against M_2 and 70.5 against the
  density of squarefree levels; weight 4: 1.5 and 388.8; at N = 30029: 4.0 and 69.5, 1.8 and 410.2. The errors are the
  scatter of single values, not calibrated: the gap between the two densities is what counts.
- **Slips to repair, none affecting the theorem**: line 77: the range of Q_ℓ is 1 ≤ Q_ℓ ≤ (ℓ+1)/(ℓ−1), not [1 − 1/ℓ,
  1 + 2/ℓ] (used only for ℓ > z); line 139: the primes p ≤ ℓ cost O(τ(ℓ)·ℓ·√N·log N); lines 221–225: the reason given for
  the low χ² is unproved; lines 90 and 109: bounds true but loose.
- **Prior art**: Zubrilina (arXiv 2310.07681) and its correction (2609.29631) average over squarefree levels; Cowan
  (2408.12723) has a fixed prime level under GRH with harmonic weights. A fixed prime level with natural weights was not
  found. The method is Zubrilina's with the roles of the level and the prime exchanged.
