# Referee report on exact5.md (27 Nov 2026)

**Verdict: no failing step.** Theorem and table hold. One caveat is false.

- **A. No failing step.** Theta law at four new x; τ(χ₅) = √5; no n = 0 term. The Mellin integral from θ₀ equals Λ(s) at s = −1.5 + 0.7i, so it is entire.
- **B. No failing step.** Λ(0⁺) = 2 log φ. ψ_R crosses only the strip [1/2, 2]. A₀ in u-form equals the Fourier side to 4e−39.
- **C. No failing step.** Closed form, constant 4, d(N) and recurrence agree with Fourier inversion to 1e−50.
- **D. No failing step.** Two inherited sentences are false for κ = 0, not load-bearing: Ce^{−u/2}, not Ce^{−u}; φ decays at rate π/5, not exp(−e^{2|u|}).
- **E. Direction correct.** For a compact f = (1−(u/1.2)²)⁸, Q = 0.0259012896181, equal to the zero sum.

**Recomputed** (computed; mpmath 50 digits, Fourier inversion of Ξ²B, mpmath L agrees with arb):
- B = 1: Defect/g(0) = 9.512241947e−3, q ≥ 4.905776697.
- Degree 8: Defect/g(0) = 2.206347658e−7, q ≥ 5 − 2.206347171e−6, carrier 19.
- All ten rows agree.
- Lemma 1: relative residual ≤ 4e−48 for seven kernels.
- Negative controls: the κ = 1 kernel gives 0.429. The kernel Ξ_{χ₈}² gives 8.2024e−4, equal to the sum over the zeros of L(s, χ₅) to 4e−42.

**Stated too strongly**
1. **False:** "no feasible family is known". μ_p = δ₀ is feasible unconditionally at q ≥ e^{−Ω(0)} = 215.33, since Ω(t) ≥ Ω(0). The theorem is therefore not vacuous. The same holds for exact3 (9.305), which corrects point 6 of the earlier report.
2. GRH is sufficient for the upper bound, not shown to be necessary.
3. In short 1 omits the condition B ≥ 0.
4. The "stall": degree 18 gains 2e−7 relative.

**Not checked:** arb, sympy, DLMF, the search codes. validate_exact5.py was not run.

Code: [referee5/ref_exact5.py](referee5/ref_exact5.py), output [referee5/ref_exact5_out.txt](referee5/ref_exact5_out.txt) (114 checks, 0 failed).
