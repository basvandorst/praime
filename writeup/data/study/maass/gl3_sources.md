# GL(3) level 1 Maass forms: sources for a positivity test

Companion to `gl3_key.json` (1050 spherical SL(3,Z) forms, conventions, checks). Page numbers are those of the arXiv PDFs. Every item below was read in the TeX source and the PDF, unless it says otherwise.

## Sources
| Source | Status |
|---|---|
| Farmer, Koutsoliotas, Lemurell (FKL), *Maass forms on GL(3) and GL(4)*, arXiv 1212.4545v1 (IMRN 2014) | read in full |
| Farmer, Koutsoliotas, Lemurell, Roberts (FKLR), *The landscape of L-functions: degree 3 and conductor 1*, arXiv 2303.00079v2 | read in full except Sec. 7 |
| github.com/davidfarmer/LandscapeData, `Degree3/degree3conductor1r0r0r0.dat` and README (FKLR ref. [8]) | parsed, 1050 entries |
| LMFDB API, degree 3, conductor 1 | 101 records read; later pages returned a reCAPTCHA, so I stopped |
| Miller, arXiv math/0006058v1 (preprint of Crelle 533, 2001) | abstract and introduction read, full source searched; the published version was not read |
| Miller, arXiv math/0112196v1 (preprint of Duke 112, 2002) | Secs. 1, 2, 6 and 8 read |
| Bian, Bull. LMS 42 (2010) 827-842, doi 10.1112/blms/bdq038 | **not read**: the publisher returned 403 and the paper is not on arXiv. I saw only its Crossref metadata. |

## Conventions
- FE (FKL eq. (3.2), p.3; FKLR eqs. (2.3)-(2.5), p.3-4, and (3.1), p.7):
  Λ(s,F) = Γ_R(s+iλ1) Γ_R(s+iλ2) Γ_R(s−i(λ1+λ2)) L(s,F), with Γ_R(s) = π^{−s/2}Γ(s/2), and Λ(s,F) = conj Λ(1−conj s, F), with N = 1 and ε = 1 (FKL p.4).
  - Normalisation: λ1 ≥ λ2 ≥ 0 and λ3 = −λ1−λ2. Each point stands for the pair {F, F~}. The dual form has conj(a_p) and parameters −λ_j (FKL p.4; FKLR p.4 and p.7).
- Euler product (FKLR eq. (1.3), p.2): L(s) = ∏_p (1 − a_p p^{−s} + conj(a_p) p^{−2s} − p^{−3s})^{−1}. FKL p.13 gives a_{p²} = a_p² − conj(a_p).
  - Here a_p = A(p,1) = α1+α2+α3 with α1α2α3 = 1. The condition |α_i| = 1 is the unproved Ramanujan conjecture, and the search did not impose it (FKLR Sec. 2.3).
- Laplace eigenvalue = 1 + λ1² + λ1λ2 + λ2² (FKLR eq. (6.5), p.19). Miller's normalisation is the same: λ = (n³−n)/24 − Σμ_j²/2 (Duke p.25).
- AFE form (FKL eq. (4.2), p.10): Λ(s) = Q^s ∏ Γ(s/2 + iλ_j/2) L(s) with Q = π^{−3/2}.
  - The GitHub file stores x = λ/2. Its gamma template {−XX[1], XX[2], −XX[1]−XX[2]} is unbalanced; read literally, it fails the FE by a factor of 10^5 to 10^9 (computed).
- Rankin–Selberg. This is not in FKL or FKLR; the Γ-list below is derived here, and the poles are classical but not read in this run.
  - Λ(s,F×F~) = Γ_R(s)³ ∏_{d} Γ_R(s+id) Γ_R(s−id) L(s,F×F~), with d ∈ {λ1−λ2, 2λ1+λ2, λ1+2λ2}. These are the same three differences as in the FKLR Plancherel density (p.18).
  - Λ(s,F×F~) = [Γ_R(s)ζ(s)] · [Γ_R(s)² ∏_d Γ_R(s±id) L(s, Ad F)], with simple poles at s = 0 and 1.
  - Miller's positivity inequality for π⊗π~ is Duke eq. (13), p.32: ∫g(x)(e^{x/2}+e^{−x/2})dx + 2Re Σ_{j,k} l(μ_jk) + g(0) log D ≥ 0. Here l(η) = (1/2π)∫h(r) Γ_R'/Γ_R(1/2+η+ir) dr (eq. (3), p.5-6).

## Excluded regions
- **Miller, Duke Sec. 6** (read on p.26-29).
  - Criterion 6.1 (p.26): Σ_j s(r_j) < 0 implies that no cusp form has Im μ_j = r_j. The function s comes from the test function g_{1,1/2}, so no primes enter.
  - Prop. 6.1 (p.26): s(r) < 0 for |r| < 7.2, which gives Σr_j² > 77.76 for n = 3.
  - Thm 6.3 (p.27), with d = 174 (p.29): λ ≥ 87.625. For tempered forms this means λ1²+λ1λ2+λ2² ≥ 87.
- **FKL p.4** states "λ1²+λ1λ2+λ2² ≥ 80", and FKL Fig. 3.1 (p.5) attributes the black region to [Mil1] (Crelle). The arXiv version of [Mil1] has no such region.
- **FKL Fig. 3.1, digitised here.** Axes calibrated on the 124 plotted forms (rms 0.05 units).
  - The black region ends at r ≈ 9.3 within 15° of the axes (λ1²+λ1λ2+λ2² ≈ 94-109).
  - It grows to a spike along λ1 = λ2 that ends near (17.2, 17.2).
  - The Sec. 6 disk criterion cannot produce that spike. I did not find its source.
- **FKL Thm 3.1** (p.4; proof in Sec. 7, p.17-18): assuming |a_p| ≤ 3, the lighter region in Fig. 3.1 contains no forms.
- **The first forms and the boundary.** Neither FKL nor FKLR remarks on this. Digitised: the five forms with the smallest Laplace eigenvalue lie 0.13 to 0.52 units from the FKL (Ramanujan-conditional) region and 1.4 to 4.1 units from the black region.

## First forms
- Lexicographic order: FKL Table 3.1 (p.5) starts with (13.59658451, 4.76468206), then (14.14163558, 2.38038848).
- By Laplace eigenvalue, (14.14163558, 2.38038848) comes first (240.315), then (16.40312474, 0.17112189) and (13.59658451, 4.76468206).
- FKLR Fig. 4.1 (p.13) uses (14.141, 2.380) as its example.
- Completeness claims: FKL p.4 says the list should be complete for λ1 ≤ 20, and may be complete up to λ1 ≤ 24. FKLR p.12 expects to have found almost all forms only up to the first few Plancherel contours.
- **Zeta zeros:** neither paper relates the spectral parameters to zeros of ζ. The only remark about the first zero is the Sp(4) zero at 14.496 (FKL p.6). The difference 14.1416 − 14.1347 = 0.0069 is mine, and I treat it as a coincidence: with about 50 forms per unit of λ1, a gap that small is to be expected (heuristic).

## Checks done here (computed)
- **Own AFE code** (mpmath and numpy, two codes) for FKL Thm 4.1.
  - It reproduces L(s,χ_{−4}) to 6e-14.
  - It reproduces FKL eqs. (4.8) and (4.9) (p.12) and (7.1) (p.17) to every printed digit, and the Sec. 7 constants 0.1866 and −4.1854.
  - It does not reproduce the printed ranges ±1.33 and ±0.301; my own d_3(n) bound gives ±1.03 and ±0.27.
  - FKL's example "λ = (8.4i, …)" works only as Γ_R(s+8.4i), that is, as Γ(s/2+4.2i).
- **FE test of all 1050 entries** at s = 1/2+i with g = e^{−iβs}, β ∈ {0, 1/4, 1/2}.
  - 802 entries pass with the stored a_p, 242 pass only after conjugation, and 6 self-dual entries pass both ways.
  - The key stores the oriented coefficients.
  - A second check with the mpmath code agrees on the sampled forms.
  - With these orientations the data match FKL Table 3.1 to below 1e-8, and all 101 LMFDB records agree.
- **Accuracy and the Ramanujan bound.**
  - Accuracy of a_p falls along each list. For the self-dual entries |Im a_p| ≤ 2e-11 for p ≤ 13; it reaches 1e-4 at p = 23 and up to 0.3 for p ≥ 29.
  - All a_p with p ≤ 29 lie in the SU(3) deltoid, so they are consistent with Ramanujan. 67 entries have a_p outside the deltoid at p ≥ 31, 62 of them at p = 37.

## Missing
- No per-number accuracy is given in the source.
- Bian's paper was not read.
- I could not find the source of the diagonal spike in Miller's region.
- The non-spherical degree-3 conductor-1 types (r0r1r1, r-δ-c2κ) are in the same repository but were not parsed.
