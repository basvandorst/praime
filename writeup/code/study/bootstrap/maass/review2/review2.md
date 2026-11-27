# Second review of maass.md

Method: own engine (weil2.py, Legendre basis), checked against the zeros of zeta, L(χ₋₃), L(χ₅) and L(Δ) to within 4e-13. Exact Toeplitz moments per prime, with no angle grid. Own barrier solver; a dual bound certifies every t*.

| | note | here |
|---|---|---|
| q*, χ odd, X = 8 / 14 / 24 | 2.9999977 / 2.9999999996 / failed | 2.9999981 / 3.0000000 / 3.0000000 |
| q*, χ even, X = 8 / 14 / 24 | 4.9998922 / 4.99999988 / 4.9999999994 | 4.9998949 / 4.9999999 / 5.0000000 |
| p and p² decoupled, X = 8 / 24 | "about 2.98, 4.985" | 2.999998, 4.998633 / 3.000000, 5.000000 |
| t* at κ = 11.75, 12, 14 | −0.0366, −0.00012, 0.165 | −0.03659, +2.2e-7, 0.16537 |
| t* at κ = 14, 17.75, 18, with the central zero | −0.470, −0.026, +1e-5 | −0.47049, −0.02602, +1.5e-4 |
| smallest weight at X = 8 / X = 16 | 12.00, 18.00 | 11.9999983, 17.99853 / 12.0000000, 17.9999986 |
| odd forms, central zero, X = 16, at 12.17, 14.36, 14.45, 16.0 | −.00062, −.00022, −.0098, +.0076 | −.00055, −.00017, −.0097, +.0075 |
| GL(3) (13.5966, 4.7647): t*, a₂ | +0.0016, 1.040−0.382i | 0.001731, 1.039−0.383i |

## Corrections

1. **In short 1, §3, §10.1.** The bound is exact only from X = 14 (odd) and X = 24 (even). At X = 8 a measure that is not the character does better by 6.3e-7 (odd) and 2.0e-5 (even) in log q; at X = 14 (even) by 2.0e-8. These gaps are certified and stable for N = 16 to 40. At X = 40 (even), p = 37 sits at −0.91, so "at every prime" is false.
2. **§3, "locates the loss".** Refuted. Decoupling every prime power leaves the values at 3 and 5. The 2.981 problem lets the positions vary freely and has a pole.
3. **In short 3, §8.** The "pieces" are pieces of {t* ≥ −0.002}. In fact t* < 0 at 9.54, 12.17, 14.36, 13.78, 17.69, 19.32 and 19.43. The true islands are much narrower than the 0.01 grid, so pieces narrower than the step cannot be excluded.
4. **In short 2.** "To 3 decimals" holds at X = 16 only. At X = 8 the optimum is not the true form: a₅ and a₇ are off by 0.007 and 0.03 at κ = 18. The match the note shows at support 8 came from a suboptimal LP point.
5. **§10.** The proof also needs GRH for Sym³, Sym⁴ and L(χ). Real angles are Deligne's theorem for the holomorphic forms. For the Maass forms they were checked: |a_p| ≤ 1.70.
6. **§2, minimax.** The dual optimum has rank ≥ 2 in each parity block, so the bound is attained by positive combinations of test functions, not by single ones.
7. **§2, "closed under conjugation".** This is false for L(F) on GL(3), and the formula does not need it.
8. **§4, levels.** Every prime gets the Steinberg option, not only the primes dividing the level. This is a valid relaxation, and the local factors are correct.
9. **§7, "on the edge".** t* = +0.0017 at the form: the form lies inside the feasible set and is not the optimum.

The Gamma factors match PARI's lfunsympow at weight 2. The functional-equation test for Sym²Δ gives −27 bits, against −2 for Γ_R(s).

Not checked: the full fine scans, X = 32, G2/G3, rounds 166 and 169. Novelty: not found in print. I searched for Chirre–Pereira–de Laat, Booker–Hiary–Keating, Chenevier and Chenevier–Taïbi, and saw snippets only.
