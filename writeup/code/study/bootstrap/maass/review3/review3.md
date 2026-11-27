# Third review: rounds 171, 172 (maass.md §8e, §8g)

Own code: [eng3.py](eng3.py) (Legendre basis, Γ-term without digamma; ζ zeros and t-domain ψ_R to 1e−13, [check_eng.out](check_eng.out)), [sdp3.py](sdp3.py) (exact moments, own barrier, bracket ≤ 3e−6). N = 64; N = 96 agrees to 3e−6; N = 32 fails near r = 19.4 (+0.008 for −0.004).

| t* | note | here |
|---|---|---|
| 4 known: 16.00, 16.1381, 16.30 | −0.081, −0.00001, +0.024 | −0.0807, +0.00002, +0.0235 |
| 171: 16.14, 16.69, 18.23, 19.49 | +.00072, −.00125, +.00110, +.00284 | +.00072, −.00134, +.00107, +.00294 |
| 171 even: 19.43, 19.47 | −.0017, −.0123 | −.0017, −.0123 |
| 19.48, key / chain forms | −.0022 / −.0014 | −.0021 / −.0012 |
| r_j+0.25, own form withheld (odd four) | – | +.030, −.010, +.0045, +.053 |
| r_j+0.01, own form included | – | −.30, −1.16, −.58, −.37 |

## Attacks
1. **Datum correct** (derived): Γ_R(s+d±i(r±r_g)), root number +1, coefficient (α^m+α^−m)(β^m+β^−m) log p (not that of Σλ_fλ_g n^−s). Truth passes at N = 96. Teeth ([attack1.out](attack1.out)): an a_p(g) off by 0.3 is rejected in 34/36 cases at p ≤ 5, 9/12 at p = 7, 4/12 at p = 11, 0/12 at p = 13. Flipping d is rejected only across parities (−1.1 to −2.3). For equal parities it relaxes (+0.06 to +0.23), as the reflection lemma forces.
2. **t* reproduced** (table).
3. **Basis not the cause. The right ends come from the rule.** Three of the four odd pieces end at r_j + 0.05 because the form itself enters at the next grid point. Without it they stay feasible to r_j + 0.25. Only 16.64 closes by itself. M2 and N3 could not fail at the ends of true pieces.
4. **No key leak** in chain.py. Its constants (0.05, −0.002, −0.0005, the margins) were set after the key was seen. All six left-crossing r̂ lie 0.0004–0.0016 below the truth, a bias of the −0.0005 rule.
5. **N2: explanation incomplete.** 19.48 was feasible only with the fed-back 16.1366 and 18.1794, which are active. The even 19.4235 is inactive: removing it changes nothing. The ternary search hit its window end, 19.49, so there was no maximum. The other branch gives 19.4813 and â₂ −1.688, a pass.

## Corrections
- In short 3, §8e "length at most 0.05": the pieces begin 0.003–0.004 below each eigenvalue. Three odd lengths are set by the 0.05 rule. The even pieces are points.
- §8g N2: replace with item 5.
- §8g: add "with 0.04 margins" after "left feasible".
- §1, §10.3: cite Chenevier–Taïbi 2020 (arXiv 1907.08783 §1.3, Prop. 2.2). They use Rankin–Selberg positivity of an unknown π against known π_i, for algebraic forms and without prime terms. I saw this through a fetch summary only.

Not checked: full scans, the J = 121 steps, other supports, interval arithmetic, the Chenevier–Lannes PDF (unreadable).
