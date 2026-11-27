# Referee report on exact3.md (27 Nov 2026)

**Verdict: no failing step.** Theorem 4 and the table hold. There is one gap in rigour (never triggered) and six overstatements.

- **A. No failing step** (re-derived). F = b(i d/du)φ decays super-exponentially, so |g_M − g_B| ≤ Ce^{−2M}e^{−2|u|} uniformly. This beats e^{−u/2}; the prime sum converges. A₁ goes through Ω = O(log) and H¹.
- **B. Sound; the sum over m is complete** (factor 2 of the tail right). Gap: `poly_min_lower`, branch `elif m > 0`, takes the endpoint minimum even when the ball |dc/m| straddles rad. Repair: use it only if lower(|dc/m|) ≥ rad, else Pc − dc²/(2m). Interior minima are caught: at degree 16, P_5 has its minimum at x = −0.99992, and both codes give 1.92e−23.
- **C. No failing step.** Constants re-derived.
- **D. No sign error.** maass.md §2, engine.py and exact3 agree. The χ₋₄ control gives +Σ h(γ).
- **E. Direction correct.**

**Recomputed** (computed; mpmath 60 digits; Fourier inversion of Ξ²B from mpmath.dirichlet):

- B = 1: Defect/g(0) = 0.04531890667, q* ≥ 2.74004536.
- Degree 8: 5.136229585e−7, q* ≥ 3 − 3.081736168e−6, carrier 13.
- Degree 22: 6.349200222e−23, q* ≥ 3 − 3.809520133e−22, carrier 41.
- All 12 table rows agree.
- Lemma 1: relative residual ≤ 1.2e−59, also for t², 1 − t²/10 (not ≥ 0) and (t² − 64.64)².
- Control ĝ = Ξ_{χ₋₄}²: residual 1.48132342e−4, equal to the sum over the 28 zeros below 70 to within 1e−30.

**Stated more strongly than proved**
1. The title "3 − ε" reads as "for every ε"; only ε = 3.81e−22 is shown.
2. The lower ends of the "certified" intervals leave out the tail.
3. §6 says "g(log 4) < 0 makes P_2 increasing", but at degree 0 g(log 4) = +0.0118 g(0).
4. §6 "binds" is heuristic, and the run whose log was not kept cannot be checked.
5. Lemma 1 is the classical explicit formula with Johns' kernel.
6. If no family is feasible (unconditionally open), the theorem is vacuous.

**Not checked:** arb, sympy, DLMF, Iwaniec–Kowalski, the search codes, the §6 dual paragraph.

**Files:** [referee/ref_exact3.py](referee/ref_exact3.py), [referee/ref_exact3_out.txt](referee/ref_exact3_out.txt) (106 checks, 0 failed).
