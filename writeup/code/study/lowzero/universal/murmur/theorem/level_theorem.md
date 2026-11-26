# Murmurations at one fixed squarefree level with natural weights: a theorem given fact 2047

Folder: `scratchpad/level_thm/` (nothing in the repository was changed). Check script `level_check.py`, output
`level_check_out.txt` (**17 checks, 0 failed**, 140 s under ulimit -v 1200000, PARI stack 300 MB). Blind test: `registration.txt`,
`predictions_blind.json`, `blind_test.py`. Labels: **proved** (written out here), **quoted** ("read" = I opened the text,
"not opened" = cited from memory or an earlier note), **checked** (assertion of the script), **NOT closed** (a gap).
Steps T1–T6, the key lemma and the assembly are those of `mprime_theorem.md` (Theorem 2; corrected copy
`mprime_theorem_corrected.md`), cited by label.

## 0. Verdict

**THEOREM 5 under GRH for Dirichlet L-functions only** (fact 2047), the same two families as Theorem 2: real primitive
characters χ_{Δ₀} of conductor ≤ 4(Y+1)N², and all characters to moduli ≤ 16N^{4δ/5}. It holds for **every squarefree
N ≥ 2, uniformly**, with a density M^(N) that depends on N only through the set of primes dividing N. The large prime factors
of N change it by O(Σ_{q|N, q large} 1/q); the small ones change it by O(1) factors (table in §3). Averaged over squarefree
N with the dimension as weight, M^(N) gives **exactly** Zubrilina's density: her constants α, β, γ, ν(r) are the averages of
the per-level local factors (Proposition 6). This explains the footnote of arXiv 2310.07681: families of levels weight the
local factors at the primes q | N differently.

## 1. Statements

**Setting.** N ≥ 2 squarefree, k ≥ 2 even, H_k^new(N) the newforms of S_k(Γ₀(N)), ε_f the root number,
ε_f = (−1)^{k/2}w_N(f), w_N = Π_{q|N}w_q, a_f(q) = −q^{k/2−1}w_q(f) (Atkin–Lehner; checked C5 against lfunrootres).
E_N(p) := |H_k^new(N)|^{−1}Σ_f ε_f√p λ_f(p); J := (y₀N, (y₀+h)N]; ⟨X⟩_J := Σ_{p∈J, p∤N} log p·X(p)/Σ_{p∈J, p∤N} log p.
For every prime ℓ (ℓ = 2 included) and r ≥ 0 let u_ℓ := 1/((ℓ−1)²(ℓ+1)) and

    e_ℓ(r) := ℓ²/(ℓ²−1) if ℓ | r (r = 0 included),   e_ℓ(r) := 1 − u_ℓ if ℓ ∤ r        (e₂ = 4/3 or 2/3),
    ν_N(r) := Π_{ℓ∤N} e_ℓ(r) = ν′(r)·Π_{q|N} e_q(r)^{−1}   (r ≥ 1; ν′ = ν_1 of Theorem 2),
    κ(N) := N/φ(N) = Π_{q|N} q/(q−1),     π(N) := Π_{q|N}(1 + 1/q) = κ(N)Π_{q|N}(1 − q^{−2}),
    M_k^(N)(y) := (1/(k−1))·[2π π(N)√y + (−1)^{k/2−1}(12/π)κ(N) Σ_{1≤r<2√y} ν_N(r)U_{k−2}(r/(2√y))√(4y − r²)] − 12κ(N)y·[k = 2].

**Theorem 5 (natural weights, one fixed squarefree level).** Assume GRH for Dirichlet L-functions. Fix k, Y ≥ 1,
δ ∈ (0, ½), A ≥ 1. Uniformly for squarefree N ≥ 2, 0 ≤ y₀ ≤ Y and N^{−1/2+δ} ≤ h ≤ 1,

    ⟨E_N⟩_J = (1/h)∫_{y₀}^{y₀+h} M_k^(N)(y) dy + O_{k,Y,δ,A}((log N)^{−A}).

Consequences (i)–(iii) of Theorem 2 carry over (h → 0, smooth weights Φ, weight 1). Nothing is claimed for single p.

**How M^(N) depends on N (proved, from the definition).** M^(N) is a function of the set P(N) of primes dividing N. For
q ∤ r, e_q(r)^{−1} = 1 + O(q^{−3}); only r < 2√(Y+1) occur. Let S ⊂ P(N) contain every prime factor ≤ 2√(Y+1) and put
ε := Σ_{q|N, q∉S}1/q ≤ 1. Then κ(N), π(N) and ν_N(r) are κ(S), π(S), ν_S(r) times 1 + O(ε), so M^(N) = M^(S) + O_Y(κ(S)ε) on
[0, Y], M^(S) being the same formula for the level Π_{q∈S}q. In particular N = mq with m fixed and q → ∞ gives M^(N) → M^(m). Completely general N is allowed: the theorem
holds with the exact M^(N); 1 ≤ κ(N) ≪ log log N and C₀ := Π_ℓ(1 − u_ℓ) ≤ ν_N(r) ≤ ζ(2), so M^(N) = O_Y(log log N).

**(a) Prime level (proved; checked D1).** For N = q prime, ν_q(r) = ν′(r)/(1 − u_q) (r < q), κ = 1 + 1/(q−1),
π(q) = 1 + 1/q: M^(q) = M_k + O_Y(1/q), and Theorem 5 is Theorem 2 (max |M^(q) − M_k| = 3.3·10⁻⁶ at q = 1000003, y ≤ 6).

**Proposition 6 ((b): the average over squarefree levels; proved; checked A1, A2).** For fixed k and y, as X → ∞,

    Σ^□_{N≤X} φ(N)·M_k^(N)(y) / Σ^□_{N≤X} φ(N)  →  M_zub,k(y) := (β/(k−1))√y + (α(−1)^{k/2−1}/(k−1))Σ_r ν(r)U_{k−2}(r/(2√y))√(4y−r²) − γy[k=2],

with Zubrilina's α, β, γ, ν(r) (arXiv 2310.07681 Thm 1, as restated with U_{k−2} inside the sum in arXiv 2609.29631 Thm 1.1;
both read). The same holds over N ∈ [X, X+Y] for Y ≥ X^{1/2+ε}. Since |H_k^new(N)| = (k−1)φ(N)/12 + O(2^{ω(N)}), the weight
φ(N) is the natural (dimension) weight of her theorem.

*Proof.* φ(N)M^(N)(y) = (N/(k−1))[2π√y Π_{q|N}(1−q^{−2}) + (−1)^{k/2−1}(12/π)Σ_r ν′(r)Π_{q|N}e_q(r)^{−1}U√(4y−r²)] − 12Ny[k=2].
*Mean-value lemma (proved).* Let g(N) = Π_{q|N}c(q) with c bounded and c(ℓ) − 1 = O(1/ℓ). Then
Σ_n μ²(n)g(n)n^{−s} = (ζ(s)/ζ(2s))G(s), G(s) = Π_ℓ(1 + (c(ℓ)−1)ℓ^{−s}/(1+ℓ^{−s})), so μ²g = μ² ∗ b with
b(ℓ^j) = (c(ℓ)−1)(−1)^{j−1}, and Σ_d|b(d)|d^{−1/2} = Π_ℓ(1 + O(ℓ^{−3/2})) < ∞. With Σ_{m≤x}μ²(m) = 6x/π² + O(√x):
Σ^□_{N≤X}g(N) = (6X/π²)G(1) + O(√X), G(1) = Π_ℓ(1 + (c(ℓ)−1)/(ℓ+1)); over (X, X+Y] the main term is (6Y/π²)G(1) and the error
is still O(√X), hence Y ≥ X^{1/2+ε}; the weight N follows by partial summation. The c needed are 1 − ℓ^{−1} (for φ(N)/N),
1 − ℓ^{−2} and e_ℓ(r)^{−1} (c − 1 = O(ℓ^{−3}) for ℓ ∤ r, bounded for the finitely many ℓ | r): all admissible.
*Euler factors (exact; A1 for every ℓ < 1000).* (ℓe_ℓ(0)+1)/(ℓ+1) = 1 + ℓ/((ℓ+1)²(ℓ−1))
(her A), (ℓe_ℓ(r)+1)/(ℓ+1) = (ℓ⁴−2ℓ²−ℓ+1)/(ℓ²−1)² for ℓ ∤ r (her B) and its ratio for ℓ | r is 1 + ℓ²/(ℓ⁴−2ℓ²−ℓ+1) (her ν);
with the mean of φ(N)/N, (ℓ+1−ℓ^{−1})/(ℓ+1), the ratios give β/2π = Π(ℓ³+ℓ²−1)/(ℓ(ℓ²+ℓ−1)), γ/12 = Πℓ(ℓ+1)/(ℓ²+ℓ−1) and
α/2π = Π(ℓ⁴−2ℓ²−ℓ+1)/(ℓ⁴−2ℓ²+ℓ) (using (12/π)/(2π)·ζ(2) = 1). ∎ **Brute force A2:** over squarefree N ≤ 3·10⁵ the weighted mean
of M^(N)(y) is 0.7829, 1.0077, −2.3039, 0.0041, 0.1453 at y = 0.1, 0.5, 1, 2, 3.3; M_zub: 0.7829, 1.0077, −2.3039, 0.0040, 0.1453.

**Reading.** Zubrilina averages over N for one prime P; Theorem 5 averages over p for one N. Both evaluate the same local
means: for ℓ ∤ N the local mean over p does not depend on N (S5 below), for ℓ | N it is 1; among squarefree N a prime ℓ divides
N with probability 1/(ℓ+1). Her density is the dimension-weighted mixture; the prime-level density (Theorem 2) is the member
with no small prime factor. (Logically, Proposition 6 is an identity between densities; combined with the uniformity of
Theorem 5 it says that the double average over (N, p) agrees with her theorem. It does not reprove her theorem.)

## 2. Proof of Theorem 5

**S1 = T1 for squarefree N (quoted; read; checked C1–C4).** For N squarefree and p ∤ N,
Σ_{f∈H_k^new(N)} ε_f√pλ_f(p) = ½H(4pN) + (−1)^{k/2−1}Σ_{1≤r<2√(p/N)} U_{k−2}(r√N/(2√p))H(4pN − r²N²) − δ_{k,2}(p + 1).
Source: Skoruppa–Zagier 1988 §2 (5), (7) (not opened); stated in this form in arXiv 2609.29631, Lemma
`lem:skoruppa-zagier-trace-formula` (read), and in Zubrilina §2 (read; U_{k−2} misplaced there). **Old forms (proved):**
SZ (7) is the trace on all of S_k(Γ₀(N)); for a newform g of level M | N, M < N, the old space V_g = span{g(dz): d | N/M}
is stable under T_p (scalar a_g(p)), and W_N maps g(dz) to a nonzero multiple of g(((N/M)/d)z) (Atkin–Lehner 1970, Lemma 26; not opened). N/M is
squarefree and > 1, so d ≠ (N/M)/d: W_N permutes 2^{ω(N/M)} lines without a fixed line and Tr(T_pW_N | V_g) = 0. (Zubrilina §2,
read: old forms from level M contribute only when N/M is a square.) **Checked:** C1 k = 2, 18 levels (6, 10, 22 have
S₂^new = 0: the formula gives 0), all p ≤ 300 with p = 2 for odd N, 1074 pairs; C2 k = 4, 243 pairs; C3 k = 6, 69 pairs —
exact rational equality with PARI trace forms of the new space, via ε_f = (−1)^{k/2+ω(N)}a_f(N)/N^{k/2−1}; C4: Tr(T_pW_N) on
the full and on the new space agree (80 traces). The term −δ_{k,2}(p+1) is the same for every squarefree N (C1).

**S2 = T2** (class-number factorisation H(D) = (√D/π)L*(D), L* = L(1, χ)ΠQ_ℓ): unchanged (Q_ℓ range as repaired).

**S3 = T3, S4 = T4**: unchanged; they hold for every D ≤ 4(Y+1)N² and do not involve N otherwise. GRH is used for χ_{Δ₀},
Δ₀ the fundamental part of −D_r(p), D_r(p) := 4pN − r²N² ≤ 4(Y+1)N²: conductor ≤ 4(Y+1)N², as before.

**S5 (local means at one level; proved; checked B1).** E_ℓ^(N)(r) := ∫_{ℤ_ℓ^×}F_ℓ(4pN − r²N²)dp.
- ℓ ∤ N (ℓ odd or ℓ = 2): **E_ℓ^(N)(r) = e_ℓ(r)**, by the computation of T5, which uses only that N is an ℓ-adic unit:
  for ℓ | r, −D ≡ −4pN is a unit and χ(ℓ) = (−pN/ℓ) is ±1 with probability ½; for ℓ ∤ r, −D = N(r²N − 4p) with r²N − 4p
  Haar-uniform off the disc r²N + ℓℤ_ℓ, and multiplying by the unit N maps Haar measure to Haar measure and the removed disc to
  a disc of unit squares N²r² + ℓℤ_ℓ; for ℓ = 2, N odd, the cases of T5 (−D ≡ 5 mod 8 for r odd; Δ′ = r′²N² − pN for r = 2r′)
  use only pN odd.
- ℓ | N, ℓ odd: D = N(4p − r²N) with 4p − r²N ≡ 4p ≢ 0 (mod ℓ), so ℓ ∥ D, ℓ ∥ Δ₀, ℓ ∤ f₀, χ(ℓ) = 0: F_ℓ = 1, **E_ℓ^(N) = 1**.
- ℓ = 2 | N, N = 2N′: −D = 4N′(r²N′ − 2p). r odd: m := N′(r²N′ − 2p) is odd, m ≡ 1 − 2pN′ ≡ 3 (mod 4); so −D = 4m,
  Δ₀ = 4m₀ (m = m₀g², g odd), 2 ∤ f₀, χ(2) = 0. r = 2r′ (r = 0 included): −D = 8N′(2r′²N′ − p) = 8·odd, Δ₀ = 8m₀, 2 ∤ f₀,
  χ(2) = 0. Both: F₂ = 1, **E₂^(N) = 1** (Zubrilina §2 makes the same observation for even N; read).
Hence E^(N)(r) := Π_ℓE_ℓ^(N)(r) = ν_N(r) exactly (r ≥ 1) and E^(N)(0) = ζ(2)Π_{q|N}(1 − q^{−2}). For ℓ > z both e_ℓ(r) and
the factor 1 at ℓ | N are 1 + O(ℓ^{−2}), as the key lemma needs. **B1:** 676 cases (ℓ = 2, 3, 5, 7; 13 levels: odd N in every
class mod 8, N even, 3 | N, 5 | N, 7 | N, N = 30030; r = 0…12), exact enumeration of p mod ℓ^J with a rigorous enclosure
(width ≤ 3.8·10⁻⁶). **E1 (data):** the mean over primes p ≤ 4N of πH(D_r)/√D_r on 5 levels per type near 10⁴ equals
E^(N)(0), ν_N(1), ν_N(2), ν_N(3) with |z| ≤ 0.77 (e.g. N = 2q: r = 1: 0.9224 against ν_N(1) = 0.9227; ν′(1) = 0.6151).

**S6 = T6** (primes in progressions; quoted GRH consequence): unchanged. D_r(p) is linear in p, so p mod Q(n) fixes
D_r(p) mod Q(n); the CRT factorisation e_r(n) = Π E_ℓ^(N)[c_{ℓ,j}] holds with the per-level local means of S5.
The primes p | N (at most log₂N of them) are left out, at cost O(log²N) ≪ hN.

**Key lemma (proved, uniform in squarefree N).** Σ_{p∈J, p∤N}log p·G(p/N)L*(D_r(p)) = E^(N)(r)∫_J G(u/N)du + O_A(‖G‖hN(log N)^{−A}).
The proof of Theorem 2 applies word for word: the error terms of T3, T4 (Rankin), T6 depend on N only through log N; the main
terms sum to Π_{ℓ≤z}E_ℓ^(N)(r) = E^(N)(r)(1 + O(1/z)) by S5.

**S7 Assembly (proved; dimension quoted, checked C6).** |H_k^new(N)| = (k−1)φ(N)/12 + O_k(2^{ω(N)}) (Martin, J. Number Theory 2005,
Thm 1, not opened; used in this form, O(N^ε), by arXiv 2609.29631 (read). C6: |dim − (k−1)φ(N)/12| ≤ 2^{ω(N)−1} for all
squarefree N ≤ 3000, k = 2, 4). Hence 12|H_k^new(N)|/((k−1)φ(N)) = 1 + O(N^{−1+ε}). With y_p = p/N, ½H(4pN) = (N√y_p/π)L*(D₀)
and H(D_r) = (N√(4y_p − r²)/π)L*(D_r), S1 gives
E_N(p) = (1 + O(N^{−1+ε}))·{(12κ(N)/((k−1)π))[√y_p L*(D₀) + (−1)^{k/2−1}Σ_r U_{k−2}(r/(2√y_p))√(4y_p − r²)L*(D_r)] − δ_{k,2}·12κ(N)(y_p + 1/N)}.
The key lemma with G₀, G_r of Theorem 2 and (12/π)ζ(2) = 2π, κ(N)Π(1 − q^{−2}) = π(N) give (1/h)∫M_k^(N); the error is
O(κ(N)(log N)^{−A}) = O((log N)^{−A+1}); rename A. ∎

**Which GRH.** Exactly that of Theorem 2: L(s, χ_{Δ₀}) for real primitive χ_{Δ₀} of conductor ≤ 4(Y+1)N² (T3), and L(s, χ) for
all χ mod q ≤ 16N^{4δ/5} (T6). No L(s, f), no L(s, Sym²f), no lower bound for L(1, χ). Covered by fact 2047.

## 3. (c) Five families, q → ∞ (k = 2; checked D1, D2; ν′(r) = C₀Π_{ℓ|r}(ℓ²/(ℓ²−1))/(1 − u_ℓ), C₀ = 0.615133)

    N = q :  M = 2π√y − 12y + (12/π)Σ_r ν′(r)√(4y−r²)
    N = 2q:  M = 3π√y − 24y + (24/π)Σ_r ν′(r)c₂(r)√(4y−r²),        c₂ = 3/2 (r odd), 3/4 (r even)
    N = 3q:  M = (8π/3)√y − 18y + (18/π)Σ_r ν′(r)c₃(r)√(4y−r²),    c₃ = 16/15 (3 ∤ r), 8/9 (3 | r)
    N = 6q:  M = 4π√y − 36y + (36/π)Σ_r ν′(r)c₂(r)c₃(r)√(4y−r²)
    N = 5q:  M = (12π/5)√y − 15y + (15/π)Σ_r ν′(r)c₅(r)√(4y−r²),   c₅ = 96/95 (5 ∤ r), 24/25 (5 | r)

| y | 0.1 | 0.2 | 0.5 | 1 | 2 |
|---|---|---|---|---|---|
| q (= Theorem 2) | 0.7869 | 0.4099 | 0.7925 | −1.6471 | 0.5008 |
| 2q | 0.5804 | −0.5851 | 1.7132 | −2.3662 | −1.9239 |
| 3q | 0.8492 | 0.1466 | 0.6833 | −3.1109 | 0.8318 |
| 6q | 0.3738 | −1.5801 | 2.1640 | −3.8991 | −1.8326 |
| 5q | 0.8843 | 0.3719 | 0.7994 | −2.3195 | 0.3872 |
| Zubrilina (average) | 0.7829 | 0.1093 | 1.0077 | −2.3039 | 0.0040 |

Weight 4 (same y): q 0.6623, 0.9366, 0.6977, 2.0944, 0.8652; 2q 0.9935, 1.4050, −0.1282, 3.1416, 2.8519; 6q 1.3246, 1.8733,
−0.7975, 4.1888, 3.3783 (`level_check_out.txt`). **How different:** a factor 2 in N changes the density by O(1) (y = 0.2:
0.41 for q, −0.59 for 2q, −1.58 for 6q; y = 2: 0.50 against −1.92); 3 | N by up to 1.5 (y = 1); 5 | N by up to 0.67 (y = 1).
The −12κy term (κ = 1, 2, 1.5, 3, 1.25) and the r = 1 line (κν_N(1) = 0.615, 1.845, 0.984, 2.953, 0.777) carry the differences.

## 4. Atkin–Lehner signs separately; non-squarefree N (remarks; NOT closed)

- For Q ∥ N, Q ≠ N: Σ_f w_Q(f)a_f(p) = Tr(T_pW_Q | S_k^new(N)). Two things change: (i) old forms contribute: on V_g (level M)
  W_Q acts by the scalar w_Q(g) on each line when Q | M, and permutes lines without fixed points otherwise, so the new trace
  needs an inversion over Q | M | N with weights 2^{ω(N/M)}; (ii) the Yamauchi/Skoruppa–Zagier trace of T_pW_Q has class numbers
  of discriminants r²Q² − 4pQ with a local factor counting solutions modulo N/Q (not opened). The natural variable is p/Q, not
  p/N. Not worked out; no claim.
- N not squarefree: old forms with N/M a square contribute (Zubrilina §2, read), q² | D_r(p) for q² | N so the local factor at
  q is a new computation, and SZ's general formula (not opened) has more terms. Not treated.

## 5. Numerical validation (discovery data, seen-type, declared as such)

Levels (5 per type, p ≤ 4N, all p ∤ N, log p weights, 16 windows of width ¼ in y ∈ (0, 4]): q = 10007, 10009, 10037, 10039,
10061; 2q with q = 5003, 5009, 5011, 5021, 5023; 3q with q = 3331, 3343, 3347, 3359, 3361; 6q with q = 1669, 1693, 1697, 1699, 1709.
Data from S1 with PARI qfbhclassno (exact). Errors: the scatter BETWEEN levels (σ_{T,w} = c_T s_w: c_T the rms of
data − M^(N) over the type's 80 window means, s_w a pooled shape over the 20 levels). It is 0.53 times the single-value standard
error (median over windows, N = 10006), confirming the referee. c_T = 0.212, 0.257, 0.285, 0.344 (q, 2q, 3q, 6q); c_T/κ =
0.21, 0.13, 0.19, 0.11 (even levels are quieter: F₂ ≡ 1 removes the largest local fluctuation).

The single levels named in the task (χ² on 16 windows with the between-level σ; rms over windows):

| N | type | M^(N) (own) | prime-level M_2 | Zubrilina | 2q-density | 3q-density | 6q-density |
|---|---|---|---|---|---|---|---|
| 10007 | q | **8.0** (0.157) | 8.1 | 70.3 (0.374) | 840 | 283 | 1756 |
| 10006 = 2·5003 | 2q | **14.7** (0.239) | 631 (1.160) | 423 (0.967) | — | 970 | 368 |
| 9993 = 3·3331 | 3q | **21.1** (0.348) | 209 (0.724) | 182 (0.723) | 865 | — | 1079 |
| 10014 = 6·1669 | 6q | **13.2** (0.307) | 694 (1.656) | 475 (1.350) | 221 | 663 | — |

All 20 levels (in-sample, 80 terms per type): own 57.2, 88.8, 65.3, 84.6; the best wrong density per type: Zubrilina 364
(type q), 6q-density 1898 (2q), Zubrilina 688 (3q), 2q-density 1177 (6q). The second counting code (`blind_test.py`, exact
rationals) reproduces all 320 window means to 5·10⁻¹⁵ (E4).

## 6. Blind test (frozen, not run)

40 levels nobody computed (types q, 2q, 3q, 6q; the first 10 primes q above 30030/m; N from 30034 to 30522), weight 2, all
p ≤ 4N, 16 windows. Predictions per level and window (own M^(N) with exact dimension; prime-level, Zubrilina and the three other
type densities as comparators) and σ_{T,w} from §5, unscaled for the larger N (conservative): `predictions_blind.json`,
sha256 2f67dc3a2f7fcf86327029c9a0db964cc0bbd21fb11a2fa0d6c6cd27f568a929. Decisions and commands: `registration.txt`.
Expected Δχ² from the predictions alone (Σ(own − alternative)²/σ²): smallest 680 (Zubrilina, type q), against a margin of 100.

## 7. Status

| step | status |
|---|---|
| Theorem 5 (every squarefree N, uniform, natural weights, error (log N)^{−A}) | **proved** given fact 2047 (GRH for Dirichlet L-functions only) |
| S1 trace formula on the new space, squarefree N | quoted (SZ not opened; 2609.29631 and Zubrilina §2 read); old-form vanishing proved (Atkin–Lehner lemma not opened); checked C1–C5 (1386 pairs) |
| S2–S4, S6, key lemma | as in Theorem 2 (refereed), uniformity in N proved here |
| S5 local factors at q \| N (ℓ = 2 included) and independence of N for ℓ ∤ N | **proved**; checked B1, E1 |
| S7 dimension of the new space | quoted (Martin Thm 1, not opened); checked C6 |
| (a) prime level = Theorem 2 | **proved**, checked D1 |
| (b) Proposition 6: dimension-weighted average = Zubrilina's density | **proved** (exact Euler factors + elementary mean value); checked A1, A2 |
| Atkin–Lehner signs w_Q, Q ≠ N; non-squarefree N | **NOT closed** (remarks only) |
| without GRH | **NOT closed** (as for Theorem 2) |

**Prior art.** Natural-weight murmuration theorems seen (Zubrilina; 2609.29631; Cardi 2607.08982 and 2609.00777, abstract only)
average over squarefree levels; Cowan 2408.12723 (fixed prime level, harmonic weights; read in an earlier gate). A density for
one fixed composite squarefree level, and the identification of Zubrilina's constants as averages of per-level factors, were
**not found in print**.


## Referee report (26 November 2026)

An independent referee (own code in [referee_level/](referee_level/)). **Verdict: Theorem 5 is proved, given GRH for
Dirichlet L-functions; no failing step.** Quoted and not proved here: the trace formula (checked exactly at 1,644 triples
(N, k, p): 18 levels with up to three prime factors, weights 2 to 12, levels with old forms included) and the dimension of
the new space up to O(2^ω). Old forms contribute exactly 0 (W_N permutes the lines g(dz) without a fixed line). Local
factors recomputed by enumeration in 396 cases (ℓ = 2, 3, 5, 7; eleven levels): 1 at every prime dividing N, the generic
mean otherwise; for even N the valuation of D at 2 is 2 or 3 and the factor is still 1, as the note says. The error is
uniform in N. Part (b): Zubrilina's constants from her source equal the φ-weighted averages to 7 digits. Own levels not
used by the author (2·10007, 3·6673, 6·3343, 5·4003): χ² 17.3, 21.5, 17.7, 14.1 for the level's own density against 435
to 1,623 for the prime-level density (72 for 5q) and 91 to 1,063 for the average over all levels.
**Prior art missed by the note:** arXiv 2603.25564, *Murmurations in the depth aspect*: one level at a time, N = ℓ^(2a+1)
with ℓ fixed and a → ∞, natural weights, mean over primes p with p/N in a fixed interval, by the Eichler–Selberg trace
formula, UNCONDITIONALLY with error (log N)^(−C) (a Bombieri–Vinogradov variant for square moduli and a mean-value bound
for real characters). So the method "one level, average over the primes" is in print for prime-power levels; the densities
M^(N) of squarefree levels and their average (part (b)) are not in it. It also suggests how to remove GRH here.
**Slips to repair:** line 33: the statement needs dim > 0 (N = 6, 10, 22 have no newforms of weight 2); line 41: "ε ≤ 1"
is false in general (write e^(O(ε))); lines 17–18: "explains the footnote" should be "suggests" (her footnote is about
non-squarefree levels); line 119: the deviation bound 2^(ω−1) fails beyond N = 3000 (N = 701,569: 25/3 against 8); only
O(2^ω) is used.
