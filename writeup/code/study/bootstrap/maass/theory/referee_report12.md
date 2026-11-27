# Referee report on exact12.md (27 Nov 2026)

**Verdict: no failing step.** Theorem 6, the assembly and all ten rows hold. My own code reproduces every number I tested. One caveat is false, and three statements are too strong.

- **A. No failing step.**
  - I re-proved Lemma 5. k′ ≤ 1 because e^{−v} ≥ 1 − v.
  - B = 1 covers (0, 7.91467), because g = φ*φ > 0 makes A′(s) = ∫u e^{−su/2} g/(1−e^{−u}) du positive. Degree 18 covers [4, 12 − 3.69306e−8).
  - The cosine bound points the right way, since an odd truncation lies below cos. Fubini applies and the s-integration is right.
  - From my own moments the bound gives A′(12)/g(0) ≥ 0.160073, 0.123009 and 0.121747, and ε = 4.08533, 8.78443e−3 and 3.69306e−8. Each bound lies below the exact Δ_B, whose crossings are ε = 3.76519, 7.40594e−3 and 3.09043e−8.
- **B. No failing step.**
  - Λ(10) from the Dirichlet series equals the Mellin integral to 1e−25.
  - Ξ agrees with PARI's lfunlambda(η²⁴, 6+it) to 1e−45 relative at six t, including the first zero.
  - τ*τ is the Dirichlet convolution and the constant 8 is right: the Bessel series equals Fourier inversion at log 1, log 2 and log 3. The misreading Σ τ(N)²K₀ does not.
  - Lemma 1 has residual/g(0) ≤ 2e−47 for five kernels, (t² − 85)² and t² among them. The controls x₂ + 10⁻³ and τ(3) + 1 break it.
- **C. No failing step.**
  - αβ = 1 always (Hecke polynomial), so α^m + β^m = 2T_m(x_p) holds algebraically.
  - All three branch-and-bound branches give valid lower bounds. Dropping a₀ is harmless.
  - My minima: Defect/g(0) = 0.3583922917, 5.40296313e−4 and 2.248086463e−9. The argmin of P₂ is −0.268127, with Defect₂/g(0) = 1.18e−7.
- **D. Direction correct.** A feasible family gives A^(κ) ≥ 2Σ min P_p = A^(12) − 2 Defect.

**Stated too strongly**
1. **False:** "may be vacuous". The uniform measure dθ/π sets every coefficient to 0, and Re ψ(x+iy) ≥ ψ(x). So it is feasible unconditionally once ψ(κ/2) ≥ log 2π, i.e. κ ≥ KU. Hence 12 − 3.7e−8 ≤ κ* ≤ KU. This is the same point as report 5.
2. "The restriction matters" and "κ ≥ 4 is needed" are false. Re ψ′(1/2+iy) = π²/(2cosh²πy) > 0, and Re ψ′ is bounded and harmonic on x > 1/2. By Poisson, Re ψ′ > 0 for x ≥ 1/2. So Δ_B is monotone on [1, 12] for every B. (Proved here; not searched.)
3. Deligne point (i): Theorem 6 does not need Defect_p ≥ 0. The real use is the bound |2T_m(x_p)| ≤ 2 for the terms with p^m > 700. List that use.
4. Write 3.694e−8 or 3.7e−8, not 3.69e−8.
5. The header of validate_exact12.py is stale.

**Not checked:** the internals of arb, Lemma 3, the Pólya citation, or a rigorous quadrature error for my numbers (NG 20 against 24: 2e−31). validate_exact12.py was not run.

Code: [referee12/referee12_check.py](referee12/referee12_check.py), output [referee12/referee12_check_out.txt](referee12/referee12_check_out.txt) (NCHK checks, 0 failed).
