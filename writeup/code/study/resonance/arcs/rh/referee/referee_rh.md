# Referee report on rh_proof.md (Hypothesis R_h: Gaussian primes in a residue class in sectors of width x^(−θ), θ > 1/4)

Folder: this folder. The repository was not changed. Labels: **proved** (checked line by line here),
**quoted (checked)** (statement read in the source), **quoted (secondary)** (known only from a source that cites it),
**computed** (my own scripts), **heuristic**.

## 1. Verdict

**PROVED UP TO one named gap.** Theorem 1 of rh_proof.md, and with it R_h of arcs_theorem.md §2 in the exact form used
there (norm ≤ x, not dyadic; sharp sector; fixed anchor φ, any orientation; uniform in φ and in ω ∈ [x^(−θ), 2π]),
follows correctly from the two stated hypotheses:

- **P_κ** (pointwise bound, κ < 1/4). Available from sources I read: Michel–Venkatesh, arXiv:0903.3591 (= Publ. IHÉS
  111), Theorem 1.1 and the theorem of §5.1 ("L(χ, 1/2) ≪_F C(χ)^{1/4−δ}", C the full analytic conductor,
  δ absolute, implied constant depending on F only) — **quoted (checked)**. Better, with an explicit exponent:
  Wu, *Burgess-like subconvexity for GL₁* (arXiv:1604.08551, Compositio 155 (2019)), Theorem 1.1:
  |L(1/2, χ)| ≪_{F,ε} C(χ)^{1/4−(1−2ϑ)/16+ε}, C(χ) = C(χ_∞)C(χ_fin), ϑ any Ramanujan–Petersson exponent — **quoted
  (checked)**; and Söhne (J. Number Theory 66 (1997) 225–250): ζ(1/2+it, λ) ≪_{K,ε} (g N𝔣)^{1/6+ε} + (N𝔣₀)^{1/2+ε} +
  (N𝔣/𝔣₀)^{1/4+ε} for any 𝔣₀ | 𝔣 — **quoted (secondary)**, verbatim in Wu's introduction (with g = C_∞(χ|·|^{it})) and
  in Heath-Brown's zbMATH review (Zbl 0885.11062). For fixed modulus this is κ = 1/6 + ε.
- **G_q** (Vinogradov–Korobov zero-free region with polylog bounds for L^{±1}, uniform in the infinity type, modulus q).
  **Not read in a primary source: this is the gap** (§2, item 1).

The note's arithmetic, its exponent A(κ) = 2 + 8κ, the corrected mean-value Lemma 5, and the counterexamples to the
published mean-value form and to Duke's conjecture are all correct (§3, §4).

Consequences for the exponent (proved, given G_q): θ < 1/(2+8κ); MV gives θ < 1/(4 − 8δ) (δ not explicit);
Wu with ϑ = 7/64 gives **θ < 64/231 = 0.2771** (ϑ = 7/64 over any number field is Blomer–Brumley, not opened; with
Kim–Shahidi ϑ = 1/9: 0.2769); Söhne gives **θ < 3/10** (Ricci's exponent, now for every fixed modulus). So the
note's Doubt 2 ("θ − 1/4 not quantified") can be removed.

**Title sentence of rh_proof.md** ("proof, with one quoted input not checked against its source"): accurate.

**Status of the Pythagorean-ray theorem (Theorem B of arcs_theorem.md):** *conditional on a citation* — on Coleman's
growth bound and zero-free region for Hecke L-functions with Grössencharacters mod h (Mathematika 37 (1990), Theorems
1–2), whose modulus is visible in two secondary sources but which I could not open. It is **no longer conditional on
GRH** (Theorem B′ is still the only version with a power saving), and it is **not yet unconditional** in the project's
sense (every input read).

## 2. Errors and gaps, ranked

1. **G_q for q ≠ 1 (Lemma 9, §3(c), §6 "σ ≤ 1 − η/2") — the one real gap; narrowed here.** It is used only in Lemma 9
   (bounds |L^{±1}| ≤ (log V)^B on Re s = 1 − η, V ≤ 3x, uniform in k), on the line where it is stated. New evidence:
   - zbMATH review of Coleman 1990 (Zbl 0721.11050, J. Kaczorowski; text obtained through the zbMATH API): "the author
     estimates (theorem 1) the order of magnitude of the Hecke L-functions with Größencharacters mod 𝔣 near the
     vertical line σ = 1. From this he derives an estimate of the zero-free region. In the special case of the Riemann
     zeta-function the results reduce to the classic theorems of H.-E. Richert and I. M. Vinogradov." — the growth bound
     (the missing part of G_q) **is stated for modulus 𝔣** — quoted (secondary: review only).
   - Sun–Zhang, arXiv:2008.09677 (LaTeX read, §"Zero-free region"), quote Coleman's Theorem 2 verbatim: "Let μλ^m mod 𝔣
     be a Hecke character … V = max(‖m‖, t). For a constant A > 0 depending only on K and sufficiently large V,
     L(σ+it, μλ^m) ≠ 0 for σ ≥ 1 − A/max(log N𝔣, (log V)^{2/3}(log log V)^{1/3}) apart from a possible exceptional zero
     when m = 0 and μ is a real character mod 𝔣" — quoted (secondary), modulus visible. Duke p. 222 (Urbialis, c₂ > 5/7)
     agrees (read).
   - Derivation (classical; not written in the note): a Richert-type bound L(σ+it, ξ) ≪ V^{c(1−σ)^{3/2}}(log V)^{2/3}
     in 1 − δ₀ < σ < 1 (Coleman Thm 1) plus the Titchmarsh §3.10–3.11 argument (Borel–Carathéodory / Landau lemma in
     discs of radius ≍ θ(V) = (log log V/log V)^{2/3}, with 1/|L(1+θ+it)| ≤ ζ_K(1+θ)) gives 1/L ≪ φ/θ ≍
     (log V)^{2/3}(log log V)^{1/3} for σ ≥ 1 − c(log V)^{−2/3}(log log V)^{−1/3}, uniformly in ξ — that is G_q with any
     a ∈ (2/3, 1). A zero-free region alone with convexity growth does **not** suffice (Borel–Carathéodory then gives
     only |1/L| ≤ exp(C(log V)^{1−a})), as the note says (lines 105–109); I re-derived this.
   - Exceptional zero: for fixed q there are finitely many real χ with k = 0; a possible exceptional zero is a fixed
     β < 1 and is excluded by taking c_q small. Harmless.
   - **What would close it:** (a) read Coleman 1990 Thm 1 and check that it is uniform in the infinity type k (V
     includes |k|) with an admissible dependence on 𝔣; or (b) switch the proof to the explicit-formula + zero-density
     route of Duke §3.1 (Lemmas 3.1.2–3.1.4, read) and Sun–Zhang, which needs **only the zero-free region**, plus a
     density estimate Σ_{|m|≤K} N(σ, T₀, χλ^m) ≪ K^{(2+8κ)(1−σ)}(log K)^B at height T₀ = (log x)^3. Sketch (heuristic,
     standard Montgomery zero-detection with mollifier length X = K and Y = K^{1+4κ}): type-I zeros are large values of
     polynomials of length in [K, Y], counted by Lemma 5 as ≪ Y^{2(1−σ)} = K^{(2+8κ)(1−σ)}; type-II zeros need
     |M_X(1/2+it)| ≥ Y^{σ−1/2}K^{−2κ}, counted as ≪ K^{1+4κ}Y^{1−2σ} = K^{(2+8κ)(1−σ)}; both match A = 2 + 8κ.
     Route (b) still rests on a secondary quotation of the zero-free region.

2. **The exponent is understated (Corollary, §3(e), Doubt 2) — not an error.** See §1: Söhne (secondary, two sources)
   gives θ < 3/10, Wu (read) gives θ < 64/231. Note that the note's Theorem 1 needs κ ≥ 1/8 only for case (L4); κ = 1/6
   is fine.

3. **Lemma 5 discussion: the T > M half of the published lemma is also false as printed (minor; the note says only the
   T ≤ M half fails).** Stucky's Lemma 4.2 asserts ∫_{−T}^{T}|Σ c(𝔞)λ^m(𝔞)N𝔞^{−it}|² dt = (T + O(N))Σ|c|². Taking
   c = 1 on the r ideals of one norm n₀, the m = 0 term alone is 2T r², against (MT + NM)·r: for M = 1 and T → ∞ the
   ratio tends to 2r, unbounded (r = 2^k for n₀ a product of k split primes) — proved (one line). The factor is
   N^{o(1)}, and for divisor-bounded coefficients Σ_n|b_n|² ≤ Σ_n r(n)Σ_{N𝔞=n}|c|² ≪ (log)^B, so the application is
   unaffected. Also: the false statement is **Stucky's Lemma 4.2 as printed**; whether "Ricci, Theorem C" in the thesis
   carries an extra hypothesis is unknown (thesis not opened). The note should say so.

4. **Small inaccuracies (no effect):** Step 1: Err must run over supp ψ₊ (the arc widened by Δ_a), not the arc.
   Step 2: the classes not prime to q contribute exactly 0 (χ(z) = 0), not O_q(log x). Table row (c), q = 1: Stucky's
   Lemmas 2.2–2.3 are themselves citations of Coleman ("follows from Theorems 1 and 2", "follows closely Titchmarsh
   3.11"), so "quoted (checked)" means "checked in Stucky", not in Coleman. The arXiv number 0904.2429 in my brief is
   Blomer–Harcos (verified via the arXiv API); the note's 0903.3591 for Michel–Venkatesh is correct.

No step fails. Every step of Lemmas 1–9, Proposition 7.1–7.4, the Main Lemma and the reduction in §5 is written out;
nothing is deferred "as in Stucky" (the class bookkeeping of §6 follows Stucky §4 but is written in full).

## 3. Re-derivations (done before reading the proofs)

**(a) Exponent from a pointwise bound.** For the family {χλ^m: |m| ≤ K} at bounded height the large sieve has
"family size" K (mean value (K + N)Σ|c|², up to the content weight) while the conductor is ≍ K², so a pointwise bound
of strength κ is K^{2κ} in family units. Montgomery's argument (ζ(1/2+it) ≪ t^c ⇒ A = 2 + 4c) with c ↦ 2κ gives
A = 2 + 8κ, and sector width x^(−θ) needs K^A ≤ x^{1−}, i.e. θ < 1/(2 + 8κ). Weyl (κ = 1/6) gives Ricci's 3/10,
convexity (κ = 1/4) gives exactly 1/4. Agrees with the note — proved (algebra of (L1)–(L4) checked line by line,
including the identity 2β(1−4κ)/(β−4κ) at the (L2) boundary) and computed (`exponent_scan.py`: sup = 2 + 8κ to 4
decimals for κ = 1/8, 1/6, 0.2, 0.2012, 0.24, 0.249, 1/4; maximisers at β = 1 + 4κ and β = 1/2 + 2κ; without the
exclusion the sup is unbounded, 1/(1−σ)).

**(b) Smooth vs sharp norm cut-off.** R_h has a sharp cut-off N z ≤ x. The note's Step 1 does not sandwich with
smooth majorants: it bounds |V − 1_{u≤x}| ≤ exp(−H²log²(u/x)/2) pointwise and counts **all lattice points** in each
polar shell (Lemma 2) with Λ ≤ log x. The edge (relative width 1/H) costs ≪ (log x)·ωx/H + √x log x, which is
o(ωx) because H = (log x)² (the log lost to Λ ≤ log x is paid by the second power of log) and ω ≥ x^(−θ), θ < 1/2.
So the edge is handled by a trivial lattice-point bound, not by the theorem and not by a sieve — and that is enough.
Heights: Ṽ(1/2+it) ≪ x^{1/2}e^{−t²/2H²}, so |t| ≤ T₁ = (log x)³ = H log x loses e^{−(log x)²/2} — proved. The
angular edge (width Δ_a = ωx^(−ε)) is a sandwich ψ₋ ≤ 1 ≤ ψ₊ whose main terms differ by O(x^(−ε)) relatively —
proved. The theorem is for N z ≤ x (full sector), uniform in φ (only |ψ̂(k)| enters) and in ω ∈ [x^(−θ), 2π]; all
error terms are uniform — proved.

**Where each input enters (checked):** P_κ only in Lemma 8, exactly on Re s = 1/2 after moving the Mellin contour from
Re s = 1 to Re s = 0 (no Phragmén–Lindelöf interpolation is needed). G_q only in Lemma 9, with V ≤ 3x, uniform in
|k| ≤ x. Both signs of m: Lemma 5 sums |m| ≤ M; the blocks K ≤ |k| < 2K contain both signs. The fixed character ξ_χ
is put into the coefficients (|ξ_χ| ≤ 1); Lemma 5 holds for arbitrary coefficients, so this is legitimate. J ≥ 10 is
fixed by (E2) (mild dependence on θ); g ≤ J/Δ + 1 is bounded. Short factors (x^{Δ/J} < N ≤ X), Möbius or not, are
handled by Lemma 5 with powering (R ≪ (KN)^{2(1−σ)}) and need G_q only for σ ≤ 1 − η/2; factors ≤ x^{Δ/J} are
bounded trivially (Z ≤ x^{2Δ}); long factors (N > X) are smooth 1/log factors and are covered by (L1)–(L4) for every
β > 0.

## 4. Verified, and to what precision

| item | how | result |
|---|---|---|
| Lemma 5 (all coefficients) | `mv_test.py`: exact LHS as c*Kc with K_ab = D_M(4Δθ)·2sin(TΔlogN)/ΔlogN (checked against direct m-sum + 400-point Gauss–Legendre in t: rel. diff 7.6e−15); worst case over **all** coefficient vectors = λ_max(W^{−1/2}KW^{−1/2}) | ≤ 0.1196 in all 24 configurations: N ∈ {500, 2500}, (M,T) ∈ {(1,1),(10,1),(300,1),(3000,1),(10⁵,1),(10⁵,0.1),(1,10³),(30,30),(10³,10³),(10,10⁴),(10³,10),(3,0.3)} — covers T ≪ M, T ≈ M, T ≫ M, N ≪ M, N ≫ M² — computed |
| 11 coefficient families | rational (k); content ≥ √N/6; one ray d(2+i); angle window 2/√N; primes; random signs; random phases; d(𝔞); μ; 1; log N | all LHS/Lemma-5 bound ≤ 0.119 — computed |
| published form (MT + N min(M,T))Σ\|c\|², T ≤ M | same script; `growth_test.py` with T = 1, M = 10N, c = 1 on (k), k² ∈ (N/2, N] | ratio = 1.05√N for N = 10³…10⁷ (3323 at N = 10⁷): **fails, growth ≍ √N** — computed; agrees with the note's 75, 110 (I get 110.1 at N = 2500, M = 10⁵) |
| same form for the coefficients that occur | `growth_test.py`: 1, μ, d, 1∗μ, 1∗1 on dyadic (N/2, N], N = 1000…6000, M = 10⁵, T = 1 | ratios 3.7–8.9, flat (μ, 1∗μ) or slowly growing (1, d): consistent with ≪ log-powers; the counterexample column grows 35 → 86 (×√N) — computed |
| Stucky's Theorem 1.1 still stands | Lemma 5 with c = b/√N𝔟, Σ\|b\|²cont/N𝔟 ≪ (log)^B (Lemma 6, proof checked) gives (M + N^g)T(log)^B = (MT + N^g min(M,T))(log)^B for T ≤ M | **proved**: Theorem A of arcs_theorem.md is unaffected |
| Duke's conjecture (Acta Arith. 52, p. 205) | read the page: "it seems reasonable to conjecture that if the range for \|m\| is \|m\| ≤ M then the left-hand side of (i) should be ≪ Σ_{N(𝔞)≤N}(N(𝔞) + qM^{n−1}T)\|c(𝔞)\|²", (i) = Σ_{χ mod 𝔮}Σ_{\|m\|≤T}∫_{−T}^{T}\|Σ_{N(𝔞)≤N}c(𝔞)χλ^m(𝔞)N(𝔞)^{it}\|²dt, c(𝔞) arbitrary | hypotheses: arbitrary c on N(𝔞) ≤ N, any integral ideal 𝔮 (𝔮 = (1) allowed; for 𝔮 ≠ 1 restrict to k ≡ 1 mod 𝔮, the ratio is still ≍ √N). For c = 1 on (k), k ≤ √N, ΣN(𝔞)\|c\|² + MTΣ\|c\|² ≤ (N + MT)Σ\|c\|², so the failure factor ≥ 1.05√N (T = 1, M ≥ 10N); computed 111.9 at N = 2500, M = 10⁵. **The counterexample is airtight** for the conjecture as printed (no log or ε factor would save it). Not found as known in one search; no novelty claimed |
| Michel–Venkatesh | own download of arXiv:0903.3591 source (byte-identical to the author's copy): Thm 1.1 ("for π an automorphic representation of GL₁(𝔸) or GL₂(𝔸) … L(π,1/2) ≪_F C(π)^{1/4−δ}"), the t-aspect remark, Theorem 5.1 ("L(χ,1/2) ≪_F C(χ)^{1/4−δ}, where C(χ) denotes the analytic conductor"), §local analytic conductor (archimedean C = Π(2+\|μ_i\|)^{deg k}), conventions (≪_F: constant depends on F only) | P_κ with κ = 1/4 − δ, δ absolute, uniform in k and t — quoted (checked) |
| Wu GL₁ | own download of arXiv:1604.08551 source: Thm 1.1 and the definition C(π) = C(π_∞)C(π_fin) | explicit κ = 1/4 − (1−2ϑ)/16 + ε — quoted (checked) |
| Söhne 1997 | Wu's introduction (verbatim formula with 𝔣, 𝔣₀) and zbMATH review Zbl 0885.11062 | κ = 1/6 + ε for fixed modulus — quoted (secondary) |
| Duke pp. 203–207, 218–223 | own download from matwbn.icm.edu.pl (aa5231.pdf), byte-identical to the author's copy; scanned pages read | Urbialis ZFR (p. 222) as the note says; Thm 3.1 only ≍ with δ < 1/(3n) = 1/6 for ℚ(i): does not contain R_h |

## 5. Not checked / could not open

- Coleman, Mathematika 37 (1990) 287–304: could not open (Wiley 403; sugaku.net 403; zbMATH web 403). Only the zbMATH
  review (via API) and Sun–Zhang's verbatim quotation of Theorem 2. Theorem 1's uniformity in the infinity type is not
  visible in the review.
- Söhne, J. Number Theory 66 (1997): could not open (ScienceDirect 403); secondary only.
- Koval'chik, Dokl. Akad. Nauk SSSR 219 (1974) 31–34: could not open (Math-Net.Ru full text "Unavailable"); zbMATH has no
  review; the longer version (Litovsk. Mat. Sb. 15 (1975) 133–151) is only listed in a bibliography (arXiv:2105.02492).
  **Whether Koval'chik already proves R_h is unknown.**
- Ricci's thesis (1976), Coleman Proc. LMS 1990, Coleman JNT 1992, Urbialis 1965, Harman's book Ch. 11, Blomer–Brumley:
  not opened.
- I did not re-run or audit rh_check.py beyond reading its output; my numbers agree with it where they overlap.

## 6. What to change in arcs_theorem.md

- §0 row "Theorem B": "proved assuming R_h" → "proved assuming G_h (Coleman 1990, Thms 1–2, Grössencharacters mod h;
  secondary), via rh/rh_proof.md Theorem 1; P_κ from Michel–Venkatesh (read) or Söhne (secondary, θ < 3/10)".
- §0 row "R_h for h > 1 unconditionally: GAP" → "proved up to G_h (one cited theorem not read in the original)";
  replace "needed > 1/4; published 5/24" by "this argument: any θ < 3/10 (Söhne, secondary) / 64/231 (Wu, read) /
  1/(4 − 8δ) (MV, read), given G_h".
- Title sentence: "(theorem under a hypothesis R_h that holds under GRH; unconditional GAP)" → "(theorem, conditional
  on a cited zero-free region with growth bounds for Hecke L-functions mod h, not read in the original)".
- §5 (e): Ricci's Weyl bound is not needed; any subconvex bound uniform in the infinity type suffices; Söhne's bound
  covers modulus 𝔣.
- §3 table: Stucky Lemma 4.2 is false as printed (T ≤ M by ≍ √N, T > M by a divisor factor); Theorem A stands with
  Lemma 5 + Lemma 6 of rh_proof.md in its place.
- Theorem B′ (GRH) stays as the version with an error term.

## 7. Files

- report: `referee_rh.md` (this file)
- `ideals.py` (ideals of ℤ[i]: angle, content, μ, d, primality), `mv_test.py` + `mv_test_out.txt` (Lemma 5, published
  form, Duke form, exact worst case), `growth_test.py` + `growth_test_out.txt` (growth in N), `exponent_scan.py` +
  `exponent_scan_out.txt` (A(κ))
- downloads (each in its own folder): `dl_wu/`, `dl_mv/`, `dl_duke/`, `dl_2008.09677/`, `dl_1506.09170/`,
  `dl_2105.02492/`, `dl_koval/` (failed)
