# Referee report on theory/global/global.md (26 November 2026)

Brief: assume each proof wrong. Own checks (quadrature of the transform of F2 at 28 points, agreement 5e-14; Prop T on four twisted examples with a smooth bump of support 0.95, agreement 5e-11; on the odd integers the formula holds for support 0.45 and fails at 0.9 by the predicted -F(1/2)/2).

| result | verdict |
|---|---|
| Lemma G | no failing step found |
| Prop T / Theorem E | gaps but fixable |
| Theorem 1 | gaps but fixable |
| Theorem 2 | no failing step found |
| Theorem B / Corollary 5.1 | no failing step found |
| Corollary F | no failing step found |
| R7 and the list G1-G3 | gaps but fixable (wording) |

1. **Lemma G.** The transform G2 of F2 is >= 0 on |xi| >= 1, vanishes there exactly on Z, and (1 + xi^2)|G2| <= 0.283. The pairing uses dilated mollifications whose transforms are bounded by C/(1 + xi^2); dominated convergence applies since int d tau/(1 + x^2) is finite. "tau on Z" to "tau = alpha(comb - delta_0)" is sound (a 1-periodic distribution vanishing on an interval of length 2 is zero). In the remark after the lemma the two inequalities are swapped (triangle: alpha + beta >= 0; F2: -alpha - beta >= 0); the proof does not use the remark.
2. **Prop T and Theorem E.** Tate's local functional equation, the residues (r F(0)/2 at s = 1, zeta_P(0) F^(0) at s = 0), the bounds on Re s = 1 - c', Phragmen-Lindelof in the strip and the role of c' < 3 (the growth condition forces all trivial zeros) check out. To add: Theorem E must assume that zeta_P converges in some half-plane; Landau's theorem (a_n >= 0) with finite order then gives absolute convergence on Re s > 1, which the Euler product, the growth condition and int d tau/(1 + x^2) need.
3. **Theorem 1.** Holds. In part (b) Fubini fails for each half-kernel separately; it holds for the symmetrised kernel w/(w^2 + t'^2). Wording fix.
4. **Theorem 2.** Holds.
5. **Theorem B, Corollary 5.1, Corollary F.** Hold. The ratio is +-e^(p0 (s - 1/2)) R(s) with p0 <= 0 and R a finite Blaschke product, |R| <= 1 on Re s >= 1/2; so the growth condition holds with N = 0 and Theorem E applies.
6. **R7 and the open cases.** (ii) must count zeros of Xi_P with multiplicity. G1-G3 are not exhaustive: they miss "zeta_P meromorphic of finite order with poles in Re s < 1/2 away from -2k" (shape zeta(s)/(1 - 2^(-s))); G3 should read "Xi_P meromorphic of finite order with poles". G2 is narrower than what is proved: the proved complement is "Lambda_P(1-s)/Lambda_P(s) not polynomially bounded on any Re s >= c'"; "positive exponential type" rests on a heuristic factorisation. The LP experiment used a cap d <= 50 and a t-grid; its claim is conditional on both.

Integrality of the multiplicities is used nowhere; the claim that every proved statement holds for real weights is right. Positivity enters through C_P >= 0, tau >= 0 and Landau's theorem; a_0 = 1 with nothing else at 1 is used in Corollary 4.1 and Theorem E.

Not checked: R5, R6, the LP numbers, the Hermite figures of Corollary 4.3, Prop 5(iv) of proofs.txt, the literature (Ryavec, the equality case of Cohn-Elkies, Hamburger/Bochner-type converse theorems).
