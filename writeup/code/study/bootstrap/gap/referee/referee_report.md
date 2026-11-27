# Referee report: gap_cert.md (zero below 9.225, real type, GRH + Ramanujan at 2), 27 Nov 2026

Own code only (repo validator read afterwards). Scripts and outputs in this folder.

## Verdict per point

1. **Lemma 2 (explicit formula): holds, but change the citation.**
   - Read verbatim (arXiv PDF text): Bober et al., Lemma 4.1, eq. (4.4). It writes Γ'/Γ(½(½+it)+μ) for Γ_R(s+μ). The correct argument is ½(½+it+μ). As printed, it is not the normalisation used here.
   - Miller (Duke 2002, arXiv math/0112196), eq. (3), read verbatim: (Γ'/2Γ)(¼+η/2+ir/2), with the u-form. It equals A_κ after x = 2u and ∫(e^{-t}/(1−e^{-t}) − e^{-t}/t)dt = γ. Cite Miller (3) and Rudnick–Sarnak.
   - Miller (3) asks for g differentiable, but G has kinks, so the smoothing step in the proof is needed. It is correct.
   - End to end, computed (`ref_zeros.py`): on 10 L-functions, Σ_ρ Ĝ(γ) over PARI zeros up to 201, plus the tail, equals W(G) with |diff| ≤ 2.3e-10. The L-functions: Δ, weight 16 level 1, 11a, 37a (central zero), χ₋₃, χ₋₄, χ₅, χ mod 5 of order 4, χ mod 7 of order 3 (both non-self-dual), Sym² 11a (degree 3, N = 121).
2. **Lemma 3: holds.**
   - A missing root contributes 0 ≤ M_p. Checked on 11a and 37a at their bad primes, on Sym² 11a at 11 (1 root of 3) and on χ at its modulus.
   - At p = 2 only m ≤ 4 enters, since G(5 log 2) = G(2a) = 0 exactly.
   - The maximum principle is correct for Re of a real polynomial.
   - W/d ≤ B holds for all 10 L-functions. Every zero with Ĝ(γ) < 0 lies in (−T, T); the largest is Δ's 9.2224.
3. **Float-with-margin steps: hold.** A sum of 4096 positive terms has relative error ≤ 4.6e-13, well below the 1e-9 margin. Redone in arb with no float steps: `ref_arb.py`.
4. **sup over κ: holds.**
   - Method: Gauss–Legendre with certified nodes plus a Bernstein-ellipse error bound (≤ 1e-46), and mean-value κ-intervals.
   - Range: κ ∈ [−1/2, 15], plus a tail argument giving F′ < 0 for κ > 14.26.
   - Result: sup F ∈ [−2.1015656e-3, −2.1015562e-3]. The upper end includes 1e-8 of slack by design.
5. **Independent numbers: hold.**
   - c0 = −0.02761853730648152 (±1e-37).
   - Upper bounds on B (Ramanujan radius R, or the bound R = p from Lemma 0):

     | R at p = 2 | R at 3 ≤ p ≤ 31 | T = 9.225 | T = 9.23 |
     |---|---|---|---|
     | 1 | 1 | ≤ −2.0511e-4 (float −2.05255e-4) | ≤ −9.2826e-4 |
     | 1 | p | ≤ −1.8115e-4 | **> 0, fails, as claimed** (+4.06e-4) |
     | 1 (also R = 1 at 5) | p | — | ≤ −6.28e-4 |
     | 1 | √p | — | ≤ −7.06e-4 |
     | 1.0365 | 1 | ≤ −2.4e-7 | — |
6. **Sanity: holds.** See points 1 and 2.
7. **Hypotheses: hold.**
   - Entire: there is no pole term.
   - The sign is right: c0 log N ≤ 0.
   - κ = 0 and κ ∈ (−1/2, 0) are covered: the pole of Γ_R(s+κ) at s = −κ is left of the line, and Gauss' integral needs κ > −1/2.

## Findings that change text

- **Headline B is mislabelled.** [−2.0525545e-4, −2.0525538e-4] is B for Ramanujan at all p ≤ 31. Under the Theorem's hypotheses (Ramanujan at 2 only, |α| ≤ p elsewhere) B ≈ −1.813e-4. Fix "In short" item 5 and commit b9c60a38's message.
- **Citation:** see point 1.
- **Cosmetic:** the scan and branch-and-bound pick maxima by `float()` comparison. That can lose 1 ulp (about 1e-19), which does not matter.

## Label

The theorem is computer-assisted. Two independent ball-arithmetic codes check it. It is conditional on GRH and on |α_j(2)| ≤ 1. The explicit formula is cited and its normalisation is checked numerically.

I would sign it under these hypotheses:
- L(s) = Σ a(n) n^{-s}, absolutely convergent for Re s > 1, equal to ∏_p ∏_{j ≤ d} (1 − α_j(p) p^{-s})^{-1}, with at most d nonzero roots per prime;
- Λ = N^{s/2} ∏_{j ≤ d} Γ_R(s+κ_j) L, with κ_j real > −1/2 and N ≥ 1;
- Λ entire of order 1, with Λ(s) = ε \overline{Λ(1−s̄)};
- GRH;
- |α_j(2)| ≤ 1.

The conclusion is then a zero with |γ| < 9.225.

## Novelty

Not found in print.
- Read in full: Miller 2002 and Bober et al. 2015 (arXiv PDFs).
- Seen as abstracts only: Zhao, arXiv 2503.15832 and 2511.13898; Bernard, arXiv 1404.6429.
- Web queries: four, on the lowest zero with local data, Δ/9.2224, and Zenodo/GitHub.
- Not seen: Bober's preprint "Bounds for large gaps", and Iwaniec–Kowalski Thm 5.12.
