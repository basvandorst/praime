# Referee report: certified Lieb-Oxford bounds 1.604335 (A) and 1.597856 (B)

Verdict: no failing step found. Both bounds hold as computer-assisted results, given LLS Theorem 1 and Lemma 10, and assuming Arb and sympy are correct. The claim to be the best rigorous bound and the prior-art section do not hold (finding 1).

## Findings
1. **Prior art.** Ye Li, "Improvement on the Lieb-Oxford Bound", Phys. Rev. A, published 25 Nov 2026, DOI 10.1103/mqsm-pwgn. The abstract (OpenAlex) gives C = 1.5342 inside a generalised LLS framework: a signed nu gives 1.5448, and an extra correction term gives 1.5342. The paper has a GitHub repo, pkupkupku123/improvement-lieb-oxford-bound (created 18 Oct 2026). Its verifier works in exact rational arithmetic and bounds three pieces: origin, middle and tail. It also proves pi <= 355/113. I saw the README, VALIDATION.md and the verifier only through fetch summaries; I did not read the paper and did not run or check the verifier. VALIDATION.md states that the four paper certificates were not fully re-verified when the repo was packaged. The May 2026 document "certified ... below 1.64" still does not resolve (DNS failure); I saw only search snippets.
2. **Float value for Theorem B.** The note gives J = 1.37702062, c = 1.59785372. Both are artefacts of the trapezoid rule on a linear grid. Integrating my g with the note's own scheme gives 1.37702068. The true value is J = 1.3770192 (to within 2e-9) and c = 1.5978532. The certified J is therefore 6.2e-6 above the true J (4.5e-6 relative), not 4.7e-6.
3. **Wrong kappa_mu in the text.** Note section 6 gives kappa_mu = 0.0319135871783. The code and my own computation give 0.031770351651336, at b = 0.5703437.
4. **Reproducibility.** Note section 10 points to files in scratch. out_*.txt, meas_K12.json and explore*.py are not in the repo, and make_combined.py hard-codes a /tmp path.

## Checks
- Point 1 holds, checked in the LaTeX source:
  - Hypotheses: Theorem 1 needs only radial probability measures mu and nu with D(mu,mu) finite. Lemma 10 adds support in the unit ball and mu({0}) = 0.
  - A single sphere is allowed for nu: Newton gives D(eta,eta) <= D(rho,rho).
  - Re-derived: Psi = Phi_mn + Phi_nm - Phi_nn, the factor 4 pi int f r^-5, and the scaling with t^3 = 2J/D.
  - Version: arXiv v3 is final. No v4 and no erratum found.
- Point 2 holds. My own Psi (Newton plus Archimedes, ball by exact shell quadrature):
  - It agrees with Monte Carlo of (19) to within 3 sigma, with exact (78) to 1.6e-15 at 20000 points, and with (43) in float arithmetic.
  - Theorem A: J = 1.6993340347, c = 1.6043348213.
  - Theorem B: J = 1.3770191667, c = 1.5978531596. D(mu,mu) recomputed as an exact rational matches.
  - Every segment integral is at or below its certified value. S1 agrees with the exact rational to 1e-13.
- Point 3 holds:
  - Cell-by-cell test: my g against chord + Haa+ d^2/2 at 75832 points (every ball interval) and 10104 points (every 4th sphere interval). No violations. The smallest margin is 1.9e-11 relative.
  - Dense sampling of Psi_aa reaches at most 0.992 x Haa+ (A) and 0.921 x Haa+ (B).
  - Support cuts, both tails and both S0 bounds: re-derived. The sphere window bound is two-sided, so pairs with negative coefficients are covered.
  - q_B(7/8,1) = 420633/20971520. My g equals Psi(a,a)/2 on [1/1024, 7/8] to 4e-15 and leaves the diagonal at a = 0.8891.
  - One configuration of (78), the partial overlap, re-derived symbolically by my own route.
- Point 4 holds: 1.6043348788155 <= 1.604335 and 1.5978555539888 <= 1.597856, at 50 digits.
- My instrumented copy of the script reproduces all 81 checks with the same numbers. The embedded measure equals meas_K48.json.

## Not checked
- The correctness of Arb and sympy.
- The sympy steps other than one configuration.
- The 3/4 of the sphere intervals that the cell test skipped.

## Label
Computer-assisted with one rigorous code. An independent float code agrees, and the cell-by-cell attack failed to break it.

## Statement I would sign
Given LLS 2022, Theorem 1 and Lemma 10, lo_certify.py proves c_LO <= 1.604335 (mu = nu = ball) and c_LO <= 1.597856 (the 48-sphere mu, nu = sigma_{2839/4096}), with f = g. The method's values are 1.6043348 and 1.5978532. These are not the smallest bounds claimed as rigorous: Phys. Rev. A (25 Nov 2026) claims 1.5342.
