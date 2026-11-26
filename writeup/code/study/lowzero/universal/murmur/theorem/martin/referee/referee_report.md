# Referee report on martin_theorem.md (26 November 2026)

Brief: assume the proof wrong; re-derive the key lemma in the free prime q; check the use of Heath-Brown's
large sieve, the local Fubini lemma, the match with the printed conjecture, and uniformity. Own code in `scripts/`.

**Verdict: no failing step found.** Three small slips, all fixable without changing any result.

1. **Key lemma in q: holds, with slips.** The Jacobi-sum expansion was re-derived line by line: substituting
   x = (a/c)y gives the coefficient ψ̄(a/c)·χ_v(c)·J(χ_v ψ̄, χ_v)/(v − 1); the two special characters give
   −χ_v(−c)/(v − 1) for ψ₀ and −χ_v(a)/(v − 1) for χ_v; so |â(ψ)| ≤ √d/φ(d). The degenerate cases are complete
   (v | s, v | M, v | m and v = 2 go into n″; v | Q₀ gives T = 0; v = ℓ and v = q cannot occur since both exceed M_B).
   Burgess applies with conductor ≤ 4YβX² and gives X^(−1/16); the large sieve (length ≍ X, moduli ≤ X^(7/8)) and
   Siegel–Walfisz (moduli polylog) apply as in the one-level note, Steps 6.5–6.6.
   - Slip (a), Step 6.1: "#{q ∈ I : m² | D} ≤ 4|I|/m² + 1" is unproved as written, since |I| counts primes. Use the
     integer bound ≤ βX/(P m²) + 1; the tail is then ≪ X log²X loglog X / m₀², still admissible.
   - Slip (b): for s ≥ 1 the factor at ℓ is exactly (1 − 1/ℓ)^(−1), because Δ_s ≡ (sQ₀q)² mod ℓ: a relative error
     O(1/X), not O(ℓ^(−2)). Harmless (the Burgess tail absorbs it); the text should say so.
   - Slip (c): the definition of d, "v ∤ msP", is empty when s = 0; that case is treated separately anyway.
2. **Lemma T′: holds; Heath-Brown 1995 not read in the original** (download blocked). The quoted form matches the
   standard statement of Corollary 2. The family fits: the d_s are distinct negative fundamental discriminants
   (real primitive characters, conductor ≤ 4YPX); a_q is supported on primes, so n₁n₂ is a square only for
   n₁ = n₂; each d occurs for at most 6τ(4ℓQ₀) values of s. The saving X^(7/4+δ) against X²/log X is correct.
3. **Lemma F: holds.** No equidistribution in ℓ is needed: homogeneity of degree 2 and invariance of every local
   datum under unit squares give the q-mean at each fixed unit ℓ₀. Own enumeration (ℓ = 7, 11, 1201, 3001, 6007,
   9001; v = 2, 3, 5, 7; Q₀ = 1, 3; s ≤ 6): 322 cases, largest deviation 8.5·10⁻⁵.
4. **Statement against the printed conjecture: holds.** Read in the LaTeX source (arXiv v3). Σ′ restricts to
   levels coprime to ℓ, the denominator is not restricted, the normalisation is ℓ^(1−k/2) a_ℓ; continuity and
   x = 0 are fine; the δ-smoothed averages follow from uniformity in ℓ.
5. **Uniformity and effectivity: holds.** Only Siegel–Walfisz is ineffective. One unstated trivial case: if
   ℓ | P then A = 0 = M(0).

**Own numerical test** (P = 3, k = 2, q prime in [2000, 4000], X = 6000; 24 primes ℓ near x = 0.2, 0.5, 1.0, 1.5;
four eigenspaces): the trace-formula code matches PARI eigenforms with their Atkin–Lehner eigenvalues in 572
cases (levels 3q, q ≤ 41) and PARI trace forms at N = 11967 (7 of 7). χ² on 96 terms: the note's M^ε 37.2;
W_N term only 234; sign flips at 3, q or both 850–4831; W_q term only 1086; zero 1218. The terms that should
vanish have rms 0.341 at X = 1500 and 0.154 at X = 6000 (about X^(−0.57); the proof gives X^(−1/4)).

**Not read:** Heath-Brown 1995 (original); arXiv 2505.09637; Burgess, Siegel–Walfisz, Atkin–Lehner, Martin 2005,
Cox; Martin's proof of Prop. 5.5; Steps 1, 5, 7–9 of the one-level note; `martin_check.py` (only its output).
