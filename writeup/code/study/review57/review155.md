# Review of round 155 (linearised explicit formula for prime pairs) — adversary note

All code is mine (written before reading tz_counts.py / tz_test.py); the repo scripts were not imported.
Folder: this folder (the count files rc_full.npz, rc_disc.npz and the fit outputs are scratch files, not stored; rc_sieve.py rebuilds them).

## Verdict (one line)
Not broken for twins and cousins. The numbers come out the same with an independent sieve, an independent estimator, and on a second range.
The novelty sentence must be weakened. The structure "factor 2 × local density × mean of χ̄ over admissible classes" is a proved
theorem for the Goldbach average in residue classes (Bhowmik–Halupczok–Matsumoto–Suzuki 2019). The twin claim for one gap is its
parity-flipped analogue for an individual gap. That analogue is a genuine conjecture (major-arc dominance at order x^{1/2}), and
the data support it strongly. The triple (0,2,6) character reading is unstable across ranges.

## 1. Recount (computed)
* `rc_sieve.py`: a segmented, odd-only numpy sieve. Primality of p, p+2, p+4, p+6 and of 2p+1 is sieved directly on the p index.
  Weights are taken from the definition in twinzeros.md §3: ln p·ln(p+h), ln p·ln(2p+1), ln p·ln(p+2)·ln(p+6). Cells are (X_{k−1}, X_k].
* The whole holdout window was recounted, p ∈ (2,000,181,785; 19,995,038,405] (4,715 cells, 1.8·10¹⁰ in length). It also contains the
  sub-range (17,993,586,935; 19,995,038,405] (216 cells, 2.0·10⁹), done first.
  N2, N4, NSG, NT026, NT046: **identical in every cell** (totals 45,108,724; 45,102,596; 43,779,393; 4,258,577; 4,255,419).
  The float sums P2, P4, SG, T026, T046 agree to ≤ 6·10⁻¹⁵ relative per cell, and the class sums mod 3…28 to ≤ 5·10⁻¹⁴.
  hold.npz's `th`/`c*` are zero below 4·10⁹ by design (prime sums there come from disc.npz). Where they overlap (the sub-range), they agree.
* Second range (1e8, 2e9] on my own grid: π₂ is consistent with π₂(10⁸) = 440,312 to within the ≈100 pairs expected in the gap.

## 2. Refit with my own method (computed)
`rc_fit.py`: time-domain weighted least squares on per-cell increments, a random-walk noise model with white increments and variance ∝ (ln x)², and
regressors (θ − x)/√x plus **all 95 non-principal characters mod 420 free**, with 4 nuisance columns. The errors were estimated three ways:
OLS, a 40-block jackknife in ln x, and 22 circular shifts of the regressors (`rc_alt.py`).

| twins h=2, holdout | mine | registered |
|---|---|---|
| ζ | 2.7625 ± 0.117 (jk 0.121; shift-null sd 0.129) | 2.772 ± 0.119 (pred 2.6406) |
| χ₅ | 0.923 ± 0.104 (jk 0.106; shift sd 0.082) | — (pred 0.8802) |
| χ₋₃ | −0.067 ± 0.105 (jk 0.136; shift sd 0.109) | −0.034 ± 0.106 |
| χ₋₄ | +0.021 ± 0.106 (jk 0.094) | +0.038 ± 0.107 |
| loudness, q ∈ {5,7,15,21} | 1.074 ± 0.052 | 1.065 ± 0.045 (q ≤ 30) |
| loudness, q ∈ {35,105} (not in the registered fit) | 0.884 ± 0.135 | — |
| all 95 coefficients vs prediction | χ² 109.7 / 95 | — |
| 74 null coefficients | χ² 94.9 / 74 (p = 0.051) | 15.76 / 13 |

The registered errors are right; the shift-null sd is 0.129 against my 0.117 and their 0.119.
The largest null deviation is an odd character of conductor 105, j = (0,1,2,2), with z = +3.82; the look-elsewhere p over 74 nulls is ≈ 0.01. It does **not**
replicate on (1e8, 2e9] (z −0.22; `rc_one_*.out`), so I treat it as chance.
On (1e8, 2e9]: ζ 2.692 ± 0.087, χ₅ 0.814 ± 0.079, χ₋₃ +0.048 ± 0.080, loudness 1.054 ± 0.038, null χ² 89.5/74.
Cousins h=4 (holdout): ζ 2.685 ± 0.117, χ₅ **−0.851 ± 0.103** (pred −0.880), χ₂₁ +0.597 (pred +0.528), quartic 15 −0.925 (pred −0.880).
**The sign flip of χ₅ between h = 2 and h = 4 is confirmed (≈ 8σ each).**
The closed form β_h(χ) agrees with my own enumeration of admissible classes to 1.5·10⁻¹⁵ for all 95 characters mod 420 (h = 2, 4).
**Band split** (`rc_band.py`, twins holdout): ζ in [2,100) 1.3 ± 1.6; [100,300) 0.9 ± 1.4; [300,1000) 3.13 ± 0.45; [1000,2000) 2.76 ± 0.35;
[2000,4000) 2.78 ± 0.22; [4000,6434) 2.77 ± 0.18. **All the weight sits at γ ≥ 300**, scales x/6400 to x/300.
There the zeros are not resolved one by one. "Twins hear the zeros" means "twin counts track θ(x) − x at those scales". The individually resolved low zeros
carry nothing measurable.

## 3. Theorem or conjecture, reading by reading
* **Averages are theorems or identities.** Summed over all h ≤ x, Σ_h Σ_{n≤x−h} Λ(n)Λ(n+h) = (ψ(x)² − ΣΛ²)/2 exactly. This is
  consistent with the product rule once the members' ranges (0, y] and (h, h+y] are taken into account. A one-member (thinning) rule applied to all
  h is inconsistent with it: it gives ∫Δ instead of xΔ(x). The identity says nothing about a fixed h.
  The Goldbach average is a theorem under RH: Fujii 1991, Σ_{n≤N} G(n) = N²/2 − 2Σ N^{ρ+1}/(ρ(ρ+1)) + O(N^{4/3+ε}). That is the factor 2 on average.
  The Goldbach average in classes is a theorem with the full character structure (BHMS 2019, below).
* **The ζ coefficient 2S(h) for an individual gap is (d), a further hypothesis.** HL with error O(x^{1/2+ε}) says nothing about it. The
  Montgomery–Soundararajan route (CMP 2004, Thm 3: HL k-tuples with uniform O(N^{1/2+ε}) give moments for H ≤ N^{1/K}) would give,
  heuristically, the covariance identity Cov(R_h(I), ψ(I)) = 2S(h)·Var ψ(I) only for short intervals, H ≤ x^{1/3} or so. The data live at
  H ≥ x/6400 ≫ x^{1/2}, where only the x^{1/2}-order terms matter. Name of the hypothesis: the major-arc explicit formula for an individual pair
  count. It is equivalent to a local-density product law at scales ≥ x^{1/2}: the linear (χ₀,1)×(χ,ρ) cross terms of Pintz's (2.11), assumed
  for all zeros and not only exceptional ones, with minor arcs carrying no x^{1/2} lines.
* **The character coefficients μ(q)χ̄(h)/∏(p−2) and the odd cancellation:** these are bookkeeping of the same hypothesis, class by class. The bookkeeping
  is right (my enumeration). Under the hypothesis it is (c), HL in classes plus the explicit formula. Without it, (d).
* **The factor 2 is forced by every product-of-densities model, Cramér included.** The data test what follows. Conditional on the actual local prime density at scales
  x/6400 to x/300, the pair count responds as a product of the two members' densities. This holds even though the real fluctuations are rigid
  (sub-Poisson). Within the u-window I measured block-sum variance ratios of 0.5 to 0.9 for the twin residuals.
  Pseudo-data were built from the real prime regressors plus either the real twin residuals (shifted) or plain Poisson noise (`rc_alt.out`). The estimator is unbiased
  (claim planted → 2.62 ± 0.12, loudness 1.01). The real data are 23σ from independence (0.02), 12σ from one-member thinning (ζ S, χ₋₃ −S:
  real χ₋₃ is 11.7σ from −1.32), 10σ from a constant 4, 24σ from the original χ₋₃ = ±2S guess, and 23σ from "no L-lines" (loudness).
  Thinning and independence are consistent with every theorem for a fixed h. So the test has real content, but the content is the
  product law, not "zeros" as such.

## 4. Odd-character cancellation, second way (computed, `rc_class.py`)
Twins split by p mod q, holdout. Class-level form of the claim: E_b ≈ (Sφ(q)/A_q)[Δ(x;q,b) + Δ(x;q,b+2)].
* q = 5: own classes 1.830 ± 0.054 (pred 1.760), other classes 0.026 ± 0.054 (pred 0). χ² 17.7/12. The alternative "class follows the total"
  gives χ² 144–260 per class.
* q = 7: own 1.620 ± 0.039 (pred 1.584), other 0.017 ± 0.028, χ² 37.2/30. q = 105: own 4.272 ± 0.031 (pred 4.225), other 0.009 ± 0.008, χ² 34.7/30.
* q = 4: both classes have coefficient ≈ S on both Δ(4,1) and Δ(4,3) (1.34–1.43 ± 0.12; pred 1.320). So there is no χ₋₄, class by class.
* Odd characters mod 5 are **present in each class at the predicted size** and cancel only in the sum. b=1: Re 0.533 ± 0.046 / Im 0.504 (pred 0.440, 0.440).
  b=2: −0.425 / −0.338 (pred −0.440, −0.440). b=4: 0.039 / 0.020 (pred 0, 0). Summed: 0.147 / 0.186 ± 0.08 (pred 0), a 1.9σ and 2.4σ residual.
  The same residual appears in the registered fit (5.2: +1.55σ, −1.66σ); it is absent on (1e8, 2e9]. I treat it as chance.

## 5. Triples (computed, `rc_triple.py`)
My S₃ = 2.858249. Raw χ₋₃ on the holdout: (0,2,6) −2.87 ± 0.91, (0,4,6) +2.95 ± 0.91 (pred ∓2.858). The signs are right.
**The (0,2,6) joint loudness is 0.746 ± 0.087 (jk) on the holdout**; the registered value is 0.72 ± 0.09. On (1e8, 2e9] it is 1.102 ± 0.061. **The two disjoint ranges are 3.4σ apart.**
Conductor 21 is 0.18 ± 0.22 against 1.26 ± 0.17 (4.0σ apart). The mirror triple (0,4,6) is stable (0.98 ± 0.08 and 1.03 ± 0.06).
The jackknife does not inflate the errors. I have no explanation; with about 100 per-conductor readings the look-elsewhere p is ≈ 0.01. This sentence of the claim should be weakened.

## 6. Prior art
| source | how opened | verdict |
|---|---|---|
| Bhowmik, Halupczok, Matsumoto, Suzuki, "Goldbach representations in arithmetic progressions and zeros of Dirichlet L-functions", Mathematika 65 (2019) 57–97, arXiv:1704.06103v2 | relevant sections read (LaTeX source: intro, Thms 1–4) | **KNOWN IN SUBSTANCE (averaged Goldbach analogue, theorem).** Thm 4: Σ_{n≤x, n≡c (q)} G(n) = (𝔖_q(c)/2)x² − (2/φ(q)²) Σ_{a:(a(c−a),q)=1} Σ_χ χ̄(a) Σ_ρχ x^{ρχ+1}/(ρχ(ρχ+1)) + O(x^{2B*}(log qx)⁵). Thm 2: the same with (χ̄(a)+χ̄(b)) for summands ≡ a, b. Change of variables: Goldbach member c − a ↔ twin member a + h. Since a ↦ −a−h permutes the admissible set, mean χ̄(a+h) = χ(−1)·mean χ̄(a), so 2 → 1 + χ(−1). Summing over n adds one integration: x^{ρ+1}/(ρ(ρ+1)) ↔ x^ρ/ρ. The full S(h) of one gap corresponds to Granville's individual-n refinement J(n), which is conjectural. |
| Granville, "Refinements of Goldbach's conjecture, and the GRH", Funct. Approx. 37 (2007) 159–173, and Corrigendum 38 (2008) | could not open (Project Euclid serves a bot check); content as quoted in BHMS | KNOWN (as quoted). 1A: RH ⇔ Σ_{n≤x}(G − J) ≪ x^{3/2+o(1)}. 1B: GRH for the L-functions of odd squarefree divisors of q ⇔ the same over n ≡ 2 (q); this is the same odd-squarefree selection. Corrigendum (2): ΣG = x²/2 − 2Σ_{|γ|≤x} x^{ρ+1}/(ρ(ρ+1)) + … |
| Fujii 1991 I–III (Acta Arith. 58; Proc. Japan Acad. 67A) | could not open (J-STAGE HTTP 500; Project Euclid bot check); quoted in Goldston–Suriajaya | KNOWN: the n-average with −2Σ (ζ only) |
| Goldston–Suriajaya, arXiv:2110.14250 (Nagoya) | intro read (source) | KNOWN: Fujii, Bhowmik–Schlage-Puchta O(N log⁵N), Languasco–Zaccagnini O(N log³N), Ω(N log log N); no twins |
| Ikeda–Suriajaya, arXiv:2405.04315 | read via HTML | Goldbach over multiples of q: only ζ zeros; characters at O(N log³N), so no linear character terms for c ≡ 0. That is consistent with the class formula. |
| Suzuki, arXiv:1504.01967 (IJNT) | main theorem read (source) | Goldbach average in APs; no twins |
| Bhowmik–Ruzsa, arXiv:1711.06442 | source grep only | no twin or shift content |
| Pintz, arXiv:1804.05561 | Thm 1 (2.1)–(2.11), Main Lemma 1/1′ read (authors' text) | KNOWN IN SUBSTANCE: the linear terms are the (ρ₀=1,χ₀)×(ρ,χ) cross terms of (2.11), rigorous for exceptional zeros only. No parity or cancellation remark found (searched "χ(−1)", "odd", "cancel"). The MV Siegel-zero discussion (his §1) treats Goldbach with q \| m, where the odd linear term vanishes for the same reason. |
| Keating–Smith, J. Phys. A 52 (2019), arXiv:1903.07057 | relevant section read (source) | uses only α(h) (the HL main term); no remainder; NOT FOUND |
| Conrey–Keating, arXiv:1604.06124 | source grep only | no remainder terms found |
| Bogomolny–Keating 1995/96 | not opened (Keating–Smith summarises them as using α(h) only) | presumably NOT FOUND; unverified |
| Montgomery–Soundararajan, "Primes in short intervals", CMP 2004 (math/0409258) | Thm 3 and hypothesis (20) read | HL with uniform O(N^{1/2+ε}); nothing on x^{1/2} terms |
| Montgomery–Soundararajan, "Beyond pair correlation" (math/0003234) | could not read (glyph-encoded source) | unknown |
| Montgomery–Vaughan 1973, QJM 24 | could not open (paywalled; no author copy found) | unknown |
| Brent 1975, Math. Comp. 29 | abstract only (author PDF is a scan; no OCR here) | maxima and minima of the twin error; nothing more known |
| Brent, "Twin primes (seem to be) more random than primes" (slides 2014/2020) | read in full (text extracted) | says there "does not appear to be any analogous constraint" from zeros on twins. The claim refines this: the zero-driven part exists but is small (coherence 0.1–0.2). |
| Friedlander–Goldston 1995/96; Languasco–Perelli–Zaccagnini 2012/16/17; Languasco–Zaccagnini Cesàro averages; Heath-Brown 1983 (twins and Siegel zeros) | not opened | unknown. Heath-Brown is relevant to the odd cancellation: an odd real Siegel character would have no linear term for twins. |

**Individual gap h with coefficient 2S(h) and characters μ(q)χ̄(h)/∏(p − 2), with odd χ cancelling:** not found in print in the sources
opened. The verdict is provisional while MV 1973, Fujii III, Friedlander–Goldston, LPZ, Heath-Brown 1983 and Granville (the paper itself) are unopened.

## 7. Which sentences stand
* "Σ Λ(n)Λ(n+h) − S(h)x ≈ 2S(h)[(ψ − x) + Σ_χ …]": **stands as an empirically supported heuristic conjecture** at x ∈ [10⁸, 2·10¹⁰], γ ∈ [300, 6400].
* "Twins carry the zeros of ζ with coefficient 2S(2)": stands. Reword it as "track θ(x) − x at scales x/6400 to x/300". The low zeros carry no weight.
* "χ₅ at 1/3": stands (0.92 ± 0.10 and 0.81 ± 0.08). The cousins have the opposite sign (−0.85 ± 0.10).
* "No line of the mod 3 or mod 4 characters; odd characters cancel between the two members": stands for the total count. Per class the odd lines are present at the
  predicted size, which is a stronger check that passes.
* SG: counts verified; fit not redone (not rechecked).
* Triples: χ₋₃ signs stand. The joint loudness of (0,2,6) is −3σ on the holdout and 3.4σ away from (1e8, 2e9]. **Weaken to "unresolved".**
* Novelty: **weaken.** The structure is the parity-flipped, individual-gap analogue of a proved averaged Goldbach formula (BHMS 2019 Thm 4;
  Fujii 1991). The measurement is what is new, not the formula.

Honest size: a strong blind confirmation (12–23σ against the natural alternatives) of the circle-method or local-density product law for a single prime gap.
The averaged version of that law, in its Goldbach form, is a theorem.
