# The murmuration density of one prime level, without GRH

Companion to `writeup/code/study/lowzero/universal/murmur/theorem/mprime_theorem.md` (Theorem 2, steps T1–T6, the
assembly; refereed). Check script: `uncond_check.py` (output `uncond_check_out.txt`: **46 checks, 0 failed**, 33 s).
Labels: **proved** (written out here), **quoted** (reference; *read* = I opened the text, *read (secondary)* = I read the
statement as restated in a later paper, *not opened*), **checked** (assertion of the script), **NOT closed**.

## 0. Verdict

**UNCONDITIONAL THEOREM** for windows of length hN with (log N)^(−B₀) ≤ h ≤ 1 (in particular h fixed), error
O((log N)^(−A)) for every A, constants **ineffective** (only through Siegel–Walfisz, step U5(ii)). The GRH steps T3–T6 and
the key lemma of the conditional proof are replaced by: the class number as a finite sum of L(1, Kronecker) values (U1); a
truncation of each L-series at M = N^(7/8) by **Burgess's bound**, for every p (U2); a factorisation of the terms into a
real character of the shifted prime 4p − r²N times a function of p mod a small modulus (U4); Siegel–Walfisz for the
conductors q ≤ (log N)^(B₁) and the **multiplicative large sieve** for (log N)^(B₁) < q ≤ N^(7/8) (U5). Because the
truncation is at N^(7/8) < hN·N^(−1/9), **range (iv) of the plan (conductors comparable to the length) never occurs**: no
shifted-prime bound of Friedlander–Gong–Shparlinski type is needed (not used, not opened).
**Correction to the plan, step (i):** for r ≥ 1 the main term is NOT carried by the principal n alone. The mean over primes p
of χ_q(4p − r²N) is â(ψ₀) = Π_{ℓ|q}(−(−r²N/ℓ)/(ℓ−1)), of modulus 1/φ(q), not 0 (p is a unit, so 4p − r²N misses the
class −r²N mod ℓ), and these terms are part of
the local densities: the principal-n sum is ≈ 1.639 for every r against E(1) = 0.6151, E(2) = 1.2303 (checks B3, B4).
Shorter windows (h < (log N)^(−B₀)) are **NOT closed**: the only obstacle is the small-modulus step U5(ii).

## 1. Statement

Setting, notation (E_N(p), ⟨·⟩_J with J = (y₀N, (y₀+h)N], log-weights, p ≠ N, ν′(r), M_k, H_k(N)) as in
mprime_theorem.md §1.

**Theorem U.** Fix k ≥ 2 even, Y ≥ 1, A ≥ 1, B₀ ≥ 0. Uniformly for N prime, 0 ≤ y₀ ≤ Y and (log N)^(−B₀) ≤ h ≤ 1,

    ⟨E_N⟩_J = (1/h) ∫_{y₀}^{y₀+h} M_k(y) dy + O_{k,Y,A,B₀}((log N)^(−A)).

No unproved hypothesis. The implied constant is ineffective. Consequences (i)–(iii) of Theorem 2 hold with "h ≥ N^(−1/2+δ)"
replaced by "h ≥ (log N)^(−B₀)": for h = h(N) → 0 no faster than a power of log N, ⟨E_N⟩_J → M_k(y₀) uniformly on [0, Y];
for fixed Φ ∈ C_c([0, ∞)) the Φ-weighted means converge to ∫ΦM_k/∫Φ; the same with weight 1 instead of log p.

## 2. Proof

Unchanged from mprime_theorem.md (none of them uses GRH): T1 (Skoruppa–Zagier trace formula; quoted, checked F1–F2 there),
T2 (H(D) = (√D/π)·L*(D), local factors F_ℓ), T5 (the local means E_ℓ(r); proved, refereed; rechecked here B1), step (e)
(dimensions), the Assembly. The assembly used GRH twice for the prime number theorem in J (Σ_{p∈J} log p and
Σ_{p∈J} log p·p/N, error √N log²N); replace by the PNT with the de la Vallée Poussin error N·exp(−c√log N) (quoted, not
opened), which is O(hN(log N)^(−A)) since h ≥ (log N)^(−B₀). Everything else in the assembly only needs:

**Key lemma U (proved, no hypothesis).** For 0 ≤ r < 2√(Y+1), (log N)^(−B₀) ≤ h ≤ 1 and G of bounded variation on [0, Y+1]
(‖G‖ := sup|G| + Var G) with G(v) = 0 for v ≤ r²/4 when r ≥ 1,

    Σ_{p∈J} log p · G(p/N) · L*(D_r(p)) = E(r) ∫_J G(u/N) du + O_{Y,A,B₀}(‖G‖ hN (log N)^(−A)),

D_r(p) := 4pN − r²N² = N·u, u := 4p − r²N, E(r) := Π_{ℓ≠N} E_ℓ(r) = ν′(r)(1 + O(N^(−2))), E(0) = ζ(2)(1 − N^(−2)).

Since p ≠ N, N ∥ D. Parameters: M := N^(7/8), L := m₀ := (log N)^(A+3), Q₀ := (log N)^(B₁), B₁ = B₁(A) fixed in U5.

**U1 (quoted classical; checked A1, A2).** For D > 0, D ≡ 0, 3 (4) let 𝓜(D) := {m ≥ 1 : m² | D, −D/m² ≡ 0, 1 (4)}. The
Hurwitz number counts all classes of forms of discriminant −D, so H(D) = Σ_{m∈𝓜(D)} h′(−D/m²), and Dirichlet's class number
formula for orders gives h′(Δ) = √|Δ|·L(1, ξ_Δ)/π for every discriminant Δ < 0, ξ_Δ := (Δ/·) (Kronecker), a character mod |Δ|
induced by χ_{Δ₀}: ξ_Δ(n) = χ_{Δ₀}(n)·1[(n, f) = 1] for Δ = Δ₀f² (Cox, *Primes of the form x² + ny²*, Thm 7.24, Cor. 7.28:
not opened). Hence **L*(D) = Σ_{m∈𝓜(D)} m^(−1) L(1, ξ_{−D/m²})** (a finite sum). Its local factor at ℓ (ℓ^a ∥ f₀) is
Σ_{i<a} ℓ^(−i) + ℓ^(−a)(1 − χ(ℓ)/ℓ)^(−1) = F_ℓ(D) of T2 (proved: −D/ℓ^(2i) is a discriminant iff i ≤ a, and ξ_{−D/ℓ^(2i)}(ℓ) = 0
for i < a, = χ(ℓ) for i = a). **A1:** all 750 D ≤ 1500, L-values exact from L(1, χ) = −q^(−1)Σ_a χ(a)ψ(a/q): max error
1.4·10^(−14); **A2:** PARI's H(D) = brute-force reduced forms for D ≤ 500.

**U2 (truncation for every p; proved from Burgess).** *Quoted:* for χ primitive mod q > 1 and every ε > 0,
|Σ_{a<n≤a+H} χ(n)| ≤ C(ε)·H^(1/2) q^(3/16+ε) (Burgess, r = 2, any modulus; r = 3 also any modulus, r ≥ 4 for cube-free q).
Burgess, Proc. LMS 13 (1963): not opened; statement **read (secondary)** in Treviño, arXiv:2511.17778, §1. Partial summation
(proved): |Σ_{n>X} χ(n)/n| ≤ 2C(ε) X^(−1/2) q^(3/16+ε). For Δ = Δ₀f², Möbius gives
Σ_{n>M} ξ_Δ(n)/n = Σ_{d|f} μ(d)χ(d)d^(−1) Σ_{n′>M/d} χ(n′)/n′. For d ≤ N^(1/32) use Burgess with q = |Δ₀| ≤ 4(Y+1)N²,
ε = 1/256: ≪_Y N^(−7/16+1/64+3/8+1/128) = N^(−5/128); for d > N^(1/32) use |Σ_{n′>X}χ(n′)/n′| ≤ |L(1, χ)| + log X + 1 ≪ log N
(L(1, χ) ≤ log q + 2, proved) and Σ_{d|f, d>N^(1/32)} 1/d ≤ τ(D)N^(−1/32) = N^(−1/32+o(1)). So for **every** p with u > 0 and
every m ∈ 𝓜(D): L(1, ξ_{Δ_m}) = Σ_{n≤M} ξ_{Δ_m}(n)/n + O_Y(N^(−1/32+o(1))), Δ_m := −D/m². Summed over m ≤ m₀ with 1/m and over
p with log p·|G|: O(‖G‖hN·N^(−1/32+o(1))). (Pólya–Vinogradov would need M = N^(1+ε) > hN, which puts conductors beyond
the large sieve: the plan's range (iv). Burgess removes it.) **C1–C3:** at N = 10007 the observed errors fall like M^(−1/2)
or faster, stay below the Burgess shape with constant 1 (max ratio 0.09), and the mean over all 1033 p ∈ (N, 2N] of the
series truncated at N^(7/8) matches E(r) (|z| ≤ 0.50, r = 0, 1, 2).

**U3 (large m, on average over p; proved).** |L(1, ξ_Δ)| ≤ |L(1, χ_{Δ₀})|·Π_{ℓ|f}(1 + 1/ℓ) ≪ log N·log log N, and
#{p ∈ J : m² | 4p − r²N} ≤ 4hN/m² + 1. So the terms m > m₀ contribute ≪ ‖G‖ log²N log log N·(hN/m₀² + log N) = O(‖G‖hN(log N)^(−A)).

**U4 (the terms as characters of the shifted prime; proved; checked A3–A6).** Put T(m, n; p) := 1[m ∈ 𝓜(D)]·ξ_{Δ_m}(n), so
L*(D) = Σ_m m^(−1) Σ_n T(m, n; p)/n.
(a) T(m, n; ·) is periodic mod 2nm² (proved: (Δ/n) depends on Δ mod 4n for Δ ≡ 0, 1 (4), from (Δ/n′) mod n′ for odd n′ and
(Δ/2) mod 8; m² | u and Δ_m mod 4n depend on u mod 4nm², i.e. on p mod nm²). A3: 9360 tests.
(b) r ≥ 1. Let q := Π{ℓ : v_ℓ(n) odd, ℓ ∤ 2mrN} (odd squarefree), q♮ := Π_{ℓ|q} ℓ^(v_ℓ(n)) = q·s², n″ := n/q♮, (n″, q) = 1.
Kronecker symbols are completely multiplicative in the lower entry, (Δ/ℓ^(odd)) = (Δ/ℓ), and (Δ_m/q) = (−N/q)(u/q) as (m, q) = 1:
**T(m, n; p) = W(p)·χ_q(4p − r²N)**, χ_q := (·/q) primitive real mod q, W(p) := (−N/q)T(m, n″; p) periodic mod K := 2n″m²,
(K, q) = 1, |W| ≤ 1. (If N | n then T = 0.)
(c) r = 0: q := Π{ℓ : v_ℓ(n) odd, ℓ ∤ 2N}; then χ_q(4p) = χ_q(p).
(d) For p ∤ q, χ_q(4p − c) = Σ_{ψ mod q} â(ψ)ψ(p) (c := r²N), â multiplicative over ℓ | q (CRT). For ℓ ∤ c: â_ℓ(ψ₀) =
−(−c/ℓ)/(ℓ − 1); for ψ ≠ ψ₀, (ℓ − 1)â_ℓ(ψ) is a Jacobi sum of the Legendre symbol and ψ̄, of modulus √ℓ, or 1 when ψ is the
Legendre symbol (quoted classical, |J(χ₁, χ₂)| = √ℓ when χ₁χ₂ ≠ ψ₀). So **|â(ψ₀)| = 1/φ(q), |â(ψ)| ≤ √q/φ(q)**. For ℓ | c
(r = 0): â_ℓ(ψ) = 1[ψ = (·/ℓ)]. Checked A4 (expansion), A5 (15 moduli up to 231, r = 1, 2, 4, 5, 7), A6.

**U5 (the n-sum in three pieces).** Per p, drop the n ≤ M with 2^(v₂(n)) > L or n₂ > L (n_odd = n₁n₂², n₁ squarefree): since
|T| ≤ 1 this costs ≪ log m₀·log N/L per p (proved). For the remaining n (the set 𝒩), n″ ≤ L³·m₀·2√(Y+1) and K ≤ 2n″m₀²
are ≤ (log N)^(C₃) with C₃ = C₃(A).
- **(ii) q ≤ Q₀ (proved, from Siegel–Walfisz).** *Quoted:* ψ(x; P, b) = x/φ(P) + O_B(x exp(−c_B√log x)) for (b, P) = 1,
  P ≤ (log x)^B, c_B ineffective (Walfisz 1936 via Siegel 1935; statement **read (secondary)**, Wikipedia
  "Siegel–Walfisz theorem"; Davenport ch. 22 not opened). Here n ≤ Q₀L²n″, so T(m, n; ·) has period P = 2nm² ≤ (log N)^(B₂).
  Partial summation (for x ≤ √N the trivial bound): Σ_{p∈J} log p·G(p/N)T(m, n; p) = e(m, n)∫_J G(u/N)du +
  O(P‖G‖N exp(−c√log N)), e(m, n) := mean of T(m, n; ·) over (ℤ/P)^×. There are (log N)^(O(1)) such (m, n); total error
  O(‖G‖hN(log N)^(−A)) since h ≥ (log N)^(−B₀).
- **(iii) Q₀ < q ≤ M (proved, from the large sieve).** *Quoted:* Σ_{d≤R} (d/φ(d)) Σ*_{ψ mod d} |Σ_{n∈I} a_n ψ(n)|² ≤
  (|I| + R²) Σ|a_n|² for any interval I (Montgomery–Vaughan 1973, Selberg; **read** as (1.1) in Friedlander–Iwaniec,
  arXiv:2409.19634). Fix (m, n″) (outer weight 1/(mn″)) and a dyadic block q ∈ (Q, 2Q]; the s-sum costs ≤ ζ(2). With
  S_b(ψ) := Σ_{p∈J, p≡b (K)} log p·G(p/N)ψ(p) and r ≥ 1:
  Σ_{q∼Q} q^(−1)|Σ_p log p G W χ_q(4p − c)| ≤ Σ_{b∈(ℤ/K)^×} [ Σ_{q∼Q} |S_b(ψ₀)|/(qφ(q)) + Σ_{1<d≤2Q} w_d Σ*_{ψ mod d}|S_b(ψ)| ]
  + O(‖G‖ Q^(1/2) φ(K) log²N), where ψ ≠ ψ₀ mod q is replaced by its primitive ψ* mod d | q (the error term: p | q) and
  w_d := Σ_{q∼Q, d|q} q^(−1/2)/φ(q) ≤ c₀Q^(−1/2)/φ(d) (proved: q = dq′, φ(dq′) ≥ φ(d)φ(q′), Σ_{q′≤X} q′/φ(q′) ≪ X).
  Cauchy–Schwarz over (d, ψ) and the large sieve with R = 2Q, |I| = hN:
  Σ_d w_d Σ*|S_b| ≤ c₀ (Q^(−1) Σ_{d≤2Q} 1/d)^(1/2)·((hN + 4Q²) Σ_{p≡b} log²p G²)^(1/2), and Cauchy–Schwarz over b gives a factor
  φ(K)^(1/2). Block total: ≪ ‖G‖ [hN/Q + φ(K)^(1/2) log^(3/2)N (hN Q^(−1/2) + (hN)^(1/2) Q^(1/2))]. For r = 0, â = 1[ψ = χ_q]
  and Cauchy–Schwarz over q (Σ_{q∼Q} q^(−2) ≤ Q^(−1)) gives the same bound. Summing the dyadic blocks Q₀ ≤ Q ≤ M and the outer
  (m, n″) (≪ (log N)^(C₃+2)):
  ≪ ‖G‖hN (log N)^(C₃+4) [Q₀^(−1/2) + (M/(hN))^(1/2)], and (M/(hN))^(1/2) ≤ N^(−1/16)(log N)^(B₀/2). Choose
  **B₁ := 2A + 2C₃ + 8**: the piece is O(‖G‖hN(log N)^(−A)). Checked on real data at N = 10007: D2 (large sieve, lhs/rhs 0.094),
  D3 (the whole chain L1 ≤ … ≤ L6 of this step, Q = 16, 32, 64, r = 1; Cauchy–Schwarz + large sieve for r = 0), D4 (block
  sums up to N^(7/8) below θ(J)(Q^(−1/2) + (Q/x)^(1/2)), r = 0, 1, 2).
- **(iv) q > M:** empty (n ≤ M). Empirical only (D5): |S_q| ≤ 3.1·(θ log θ)^(1/2) for 300 conductors near N and N^(1.1).

**U6 (the main term; proved; checked B1–B5).** T(m, n; p) = Π_ℓ T_ℓ(v_ℓ(m), v_ℓ(n); u) with T_ℓ depending on u ℓ-adically
(proved: (ab²/ℓ^j) = (a/ℓ^j) for ℓ ∤ b, also for ℓ = 2), so by CRT e(m, n) = Π_ℓ e_ℓ(v_ℓ(m), v_ℓ(n)). Pointwise
Σ_{i,j} T_ℓ(i, j; u)ℓ^(−i−j) = F_ℓ(D) (U1), so by dominated convergence Σ_{i,j} e_ℓ(i, j)ℓ^(−i−j) = ∫ F_ℓ dp = E_ℓ(r) (T5).
For ℓ ∤ 2rN: e_ℓ(0, j) = 1, −1/(ℓ−1), 1 − 1/(ℓ−1) (j = 0, odd, even ≥ 2) and |e_ℓ(i, ·)| ≤ ℓ^(1−2i)/(ℓ−1), so
Σ_{i,j}|e_ℓ|ℓ^(−i−j) = 1 + O(ℓ^(−2)): **Σ_{m,n} e(m, n)/(mn) = E(r), absolutely convergent**. The (m, n) left out in U2–U5
(m > m₀, n ∉ 𝒩, q > Q₀, using |e| ≤ |â(ψ₀)| = 1/φ(q)) contribute O((log N)^(−A)). For r ≥ 1, e_ℓ(0, odd) = −1/(ℓ−1) ≠ 0: the
non-principal n carry part of the main term. B1: local sums = closed forms of T5 (ℓ ≤ 13, r = 0..6, ≤ 2.5·10^(−6) with a
rigorous tail); B3: Σ_{m≤40, n≤2·10⁴} e(m, n)/(mn) = E(r) within 0.006 for r = 0..6 (E(0) = ζ(2)(1 − N^(−2)),
E(r) = ν′(r)(1 − u_N)^(−1)); B4: the principal-only sums 1.638–1.640 differ from E(r) by 0.16–1.02; B5: for r = 0 they agree.

**Proof of Key lemma U.** U1 (exact) + U2 (per p) + U3 + the per-p cut to 𝒩 + U5(ii) with U6 + U5(iii); each error is
O(‖G‖hN(log N)^(−A)). ∎ With the assembly of mprime_theorem.md (G₀(v) = √v, G_r(v) = 2√v sin((k−1)θ_r(v)), ‖G_r‖ ≪_k √(Y+1),
12E(0)/π = 2π(1 − N^(−2)), E(r) = ν′(r)(1 + O(N^(−2)))) this proves Theorem U. ∎

## 3. Inputs

| input | step | status |
|---|---|---|
| Skoruppa–Zagier trace formula at prime level, old forms trace 0 | assembly | quoted in mprime_theorem.md (read in arXiv:2609.29631); checked F1–F2 there; unconditional |
| H(D) = Σ_{m∈𝓜(D)} h′(−D/m²), h′(Δ) = √|Δ| L(1, ξ_Δ)/π | U1 | quoted classical (Cox, not opened); checked A1–A2 |
| Burgess r = 2, any primitive χ | U2 | quoted; read (secondary): Treviño arXiv:2511.17778 §1; Heath-Brown arXiv:1203.5219 abstract only |
| L(1, χ) ≪ log q; Möbius for imprimitive ξ_Δ; partial summation | U2, U3 | proved |
| Kronecker periodicity and multiplicativity | U4 | proved; A3 |
| Jacobi sums |J| = √ℓ (or 1) | U4 | quoted classical; A4–A6 |
| multiplicative large sieve, constant |I| + R² | U5(iii) | quoted; read: Friedlander–Iwaniec arXiv:2409.19634 (1.1); D2 |
| Siegel–Walfisz, moduli ≤ (log x)^B | U5(ii) | quoted; read (secondary): Wikipedia; **ineffective** |
| PNT with error x exp(−c√log x) | assembly | quoted, not opened |
| local means E_ℓ(r) | U6 | proved (T5, refereed); B1 |
| Friedlander–Gong–Shparlinski, Karatsuba, Rakhmonov | — | not used, not opened |

## 4. Remarks

- **Effectivity.** Burgess (C(ε) effective), the large sieve (explicit) and the PNT are effective; the only ineffective input
  is Siegel–Walfisz in U5(ii), i.e. a possible Siegel zero of a real character of conductor ≤ (log N)^(B₂). An effective
  version (Page's theorem, exceptional modulus ≤ (log N)^(2−ε)) is **not worked out**: the large-sieve losses (log N)^(C₃+4)
  force Q₀ far above (log N)².
- **Shorter windows: NOT closed.** U2–U5(iii) hold for h ≥ N^(−1/8+δ) (M = N^(7/8)); with Burgess r = 3 (any modulus;
  M = N^(2/3+δ)) for h ≥ N^(−1/3+2δ); with Burgess for all r on real characters (the odd part of Δ₀ is squarefree, the factor
  8 removed by splitting n mod 8; a sketch, not written out) for h ≥ N^(−1/2+2δ), the range of the GRH theorem. What is
  missing is U5(ii) in short windows: Σ_{p∈J, p≡b (P)} log p = hN/φ(P)·(1 + O((log N)^(−A))) for P ≤ (log N)^B and
  hN ≤ N^(1−η) — primes in progressions to small moduli in short intervals (Huxley-type results; not opened).
- **Rate.** (log N)^(−A) for every A, the same as under GRH; a power saving would need one in U5(ii).
- **Against the plan.** (i) is wrong for r ≥ 1 (U6, B4); (ii), (iii) are as proposed, with the composite and non-squarefree n
  handled by the split n = q♮·n″ (U4(b)) and the coprimality and 2-adic conditions put in W mod K; (iv) is avoided by
  truncating with Burgess instead of Pólya–Vinogradov (at N = 10007 the Pólya–Vinogradov tail bound exceeds 1 for every
  M ≤ N^(1.1), check C1).
- **Script history.** The first run had two empirical thresholds that failed (max error·√M ≤ 3: observed 5.28; error < 0.05 at
  M = N^(7/8): observed 0.083 for the worst of 120 samples); they were replaced by the scale-free tests C1 (decay) and C3 (mean
  over all p). The proof uses no numerical constant from these tests.
- **Prior art:** as in mprime_theorem.md §10; no new search was made.


## Referee report (26 November 2026)

An independent referee (own code in [referee_uncond/](referee_uncond/)), asked to assume the proof wrong and to attack the
two places most likely to break (weights of long period in the Siegel–Walfisz range; the bilinear structure of the large
sieve). **Verdict: proved unconditionally, given the quoted classical inputs (Burgess's bound with r = 2 for any modulus,
Siegel–Walfisz, the multiplicative large sieve, the prime number theorem with de la Vallée Poussin's error, the class-number
formula for orders) and the refereed steps T1, T2, T5 of mprime_theorem.md; no failing step, no hidden hypothesis.**
- The moduli in the Siegel–Walfisz range stay a power of log N because of the cut made for each p before any average
  (line 93: every n whose 2-part or odd square part exceeds (log N)^(A+3) is dropped, at a pointwise cost log M/L); the
  remaining periods are polylogarithmic.
- In the large sieve the sequence a_p does not depend on the conductor q (W depends on q only through a sign); the factor
  φ(K)^(1/2) comes from Cauchy–Schwarz over the classes mod K and is polylogarithmic; imprimitive ψ and ψ₀ are accounted for.
- Checked with other code: the class-number formula for orders on 300 discriminants (2·10⁻¹⁴); the factorisation of U4 in
  23,705 random cases; the main-term double series against ζ(2), ν′(1), ν′(2), ν′(3) (1.64493, 0.61488, 1.22888, 0.73798
  against 1.64493, 0.61513, 1.23027, 0.73816; the gap is the tail in m); at the level 20011, all 1,942 primes in (N, 2N]:
  the means by conductor block follow the local prediction (q = 1: +0.8201 against +0.8193; 1 < q ≤ 32: −0.2080 against
  −0.2075; every block above 32 at most 0.003).
- Slips, none affecting the theorem: lines 68–70: the truncation error is N^(−1/16+o(1)), not N^(−1/32+o(1)) (Burgess
  applies to every divisor; the split is unnecessary); lines 67–68: N divides the fundamental discriminant, so the
  conductor is never small (worth stating); line 81: the period is already n·m²; line 107: the term "p | q" is empty for
  r ≥ 1.
- Not opened by the referee: Burgess 1963, the large-sieve source, Walfisz, Cox. No new prior-art search; averages of
  L(1, χ) over shifted primes by the large sieve are classical, and for prime-power levels an unconditional statement of
  this kind is in print (arXiv 2603.25564).
