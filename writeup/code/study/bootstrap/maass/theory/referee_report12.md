# Referee report on exact12.md (27 Nov 2026)

**Verdict: no failing step.** Theorem 6 and all ten rows hold. One caveat is false.

- **A. No failing step.**
  - Lemma 5 re-proved. B = 1 (g = φ*φ > 0) covers (0, 7.91467), and degree 18 covers [4, 12 − 3.69306e−8).
  - The cosine bound points the right way. My own moments reproduce A′(12)/g(0) ≥ 0.160073, 0.123009, 0.121747 and ε = 4.08533, 8.78443e−3, 3.69306e−8. Each bound lies below the exact Δ_B, whose crossings are ε = 3.7652, 7.4059e−3, 3.0904e−8.
- **B. No failing step.**
  - Ξ equals PARI's lfunlambda(η²⁴, 6+it) to 1e−45 relative, including at the first zero.
  - τ*τ is the Dirichlet convolution and the constant 8 is right. The Bessel series matches Fourier inversion, while Σ τ(N)²K₀ does not.
  - Lemma 1: residual/g(0) ≤ 2e−47 for five kernels. Changing x₂ by 10⁻³ or τ(3) by 1 breaks it.
- **C. No failing step.**
  - αβ = 1 always, so the identity is algebraic.
  - The branch and bound is valid; dropping a₀ is harmless.
  - Defect/g(0) = 0.3583922917, 5.40296313e−4, 2.248086463e−9. P₂ has its argmin at −0.268127.
- **D. Direction correct.**

**Stated too strongly**
1. **False:** "may be vacuous". The uniform dθ/π kills every coefficient and Re ψ(x+iy) ≥ ψ(x). So it is feasible unconditionally for κ ≥ 13.5531, and 12 − 3.7e−8 ≤ κ* ≤ 13.5531. This is the same point as report 5.
2. "κ ≥ 4 is needed" is false. Re ψ′(1/2+iy) = π²/(2cosh²πy) > 0, and Re ψ′ is bounded and harmonic on x > 1/2, so Re ψ′ > 0 for x ≥ 1/2 (Poisson; not searched in print). Monotonicity holds on [1, 12] for every B.
3. Deligne point (i): Theorem 6 does not need Defect_p ≥ 0. The real use is |2T_m(x_p)| ≤ 2 in the tail p^m > 700.
4. Write 3.694e−8, not 3.69e−8.
5. The header of validate_exact12.py is stale.

**Not checked:** arb internals, Lemma 3, the Pólya citation, rigorous error bounds for my quadrature (20 against 24 nodes per panel: 2e−31). validate_exact12.py was not run.

Code: [referee12/referee12_check.py](referee12/referee12_check.py), output [referee12/referee12_check_out.txt](referee12/referee12_check_out.txt) (70 checks, 0 failed).
