# Round 165 review (own Legendre/x-space Weil form, exact Toeplitz-moment SDP, own barrier solver, float64)

**Findings**
1. **The root number is missing.** Odd forms have L(1/2,f) = 0, so under GRH Q_f ≥ (∫φ)². The true forms stay feasible. At X = 8 (odd) it gives t* < 0 on 11.5–15.75 except at 12.1730 (+7e-6) and 14.3585 (+4e-4). Example: 13.0 goes from +0.152 to −0.051. The feasible region [12.2, 20] shrinks to [16.0, 20].
2. **Proved:** Q_odd − Q_even = ½∫|φ̂|² Σ_ρ sech(π(t−ρ)) dt ≥ 0, with ρ = ±r (and ±3r for Sym³). So t*_odd ≥ t*_even, and the odd scan admits every even eigenvalue (your scans agree).
3. **L(f) alone gives the whole t* at X = 8.** This holds at r = 9, 9.53 and 10. Sym² alone is feasible, and Sym³ and Sym⁴ were never binding.
4. **D2 and D5 fail because of the r-grid.** The island is [9.533695256, 9.533696035] (N = 32). Grid points give −0.0022 (at 9.55) and −0.008 (at 13.80), both below −0.002. D5 misses for p = 7 (off by 0.053).
5. No bug found in engine.py or maass.py; one measure per prime, as 3d requires.

**Task 2** (X = 8; t* at N = 24 / 48; a_p at N = 48)

| eps, r | t* | a_2 a_3 a_5 a_7 |
|---|---|---|
| 1, 3 | −1.60994 / −1.60995 | −.66 1.20 .03 2.00 |
| 1, 6 | −.78683 / −.78684 | .10 −.67 −.16 −1.17 |
| 1, 9.0 | −.20387 / −.20388 | −1.10 −.52 −.42 −.02 |
| 1, 9.53369526 | 1.7e-9 / 0+ | −1.0683 −.4562 −.2907 −.7396 |
| 1, 10 | −.05718 / −.05756 | −.94 −.73 −.13 −1.00 |
| 1, 11 | −.10482 / −.10509 | −.43 −1.16 .06 −.55 |
| 1, 12.5 | .05425 / .05408 | .47 −1.06 .12 .80 |
| 0, 12 | −.59099 / −.59175 | .57 −1.67 −.21 .52 |
| 0, 13.77975135 | 2.2e-6 / 0+ | 1.5493 .2469 .7371 −.2594 |
| 0, 15 | −.53479 / −.53490 | 1.15 .07 −.44 .46 |

Your values reproduce within 7e-4; your −0.0003 at r₁ is an angle-grid artifact.

**Task 1.** Normalisation checked against ζ zeros (4e-14). λ_min with the true a_p (N = 24): form #1 1.2e-13 (L), 6.3e-5 (Sym²); form #3 3.2e-10, 0.031. With the wrong parity, #1 gives −0.94 but #3 gives +1.9e-4 (not rejected).

**3a.** Odd r = 13, 15, 18 at N = 96: .1520, .0853, .5199 (N = 24: .1521, .0854, .5229). No sign change.

**3b.** With eps = 0, t* lies in [−0.66, −0.54] on 9–12.5, with no peak at 9.53.

**Not checked:** a primary source for the Gamma factors (derived only), the holomorphic control, my own X ≥ 16 scan, and precision beyond float64.
