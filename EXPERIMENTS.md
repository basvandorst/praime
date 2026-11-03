# Experiments on the logic of the prime and twin prime sequences

Branch `experiments`. Puzzle results live on `main`; this branch holds the research runs, each under two minutes on one core. Odlyzko's first 100,000 zeros (`zeros1`, 1.8 MB) are fetched with `curl -sL -o zeros1 http://www.dtc.umn.edu/~odlyzko/zeta_tables/zeros1` and kept outside the repo.

| Date | Question | Result | Code |
|---|---|---|---|
| 2026-11-03 | Is the twin coin the square of the prime coin, and are the two coins independent? | Twin-survivor density = V₂ · B₁² to 10⁻³; twin coin = ρ² · ∏_{q>y}(1 − 1/(q−1)²) to 5·10⁻⁴; independent at y = x^0.4 to 10⁻³ at 10⁹, 10¹², 10¹⁵; parity census = independent coins (fact 4027) | writeup/code/twins_coin.py |
| 2026-11-03 | Does the coin sequence compress? What periods hide in its failures? | zstd within 0.7 % of the entropy bound, same size as a shuffled coin; the failures are divisible by exactly the 499 primes in (y, 2y] | writeup/code/coin_logic.py |
| 2026-11-03 | Do past twin gaps predict the next one beyond the wheel? | p mod 5005 cuts the gap entropy 3.81 → 2.89 bits; the last two gaps only to 3.66; adding the last gap to p mod 385 gains 0.007 bits | writeup/code/coin_logic.py |
| 2026-11-03 | Do the zeros rebuild the primes? | 100,000 zeros return all 1229 primes below 10⁴, no false positive; Λ(n) to 0.004 below 100 (fact 2047) | writeup/code/zeros_to_primes.py |
| 2026-11-03 | Zero spacings GUE? Prime spectrum dips at the zeros? | GUE bin by bin (var 0.162 vs 0.178; Poisson 1); dips of exactly −ln X at the zeros, 0 between them | writeup/code/zeros_gue_spectrum.py |
| 2026-11-03 | Do twins carry the zeros once (independent coin) or twice (both members)? | 1.85 ± 0.20 above 10², 1.94 ± 0.29 above 10⁴ (circle method: 2; coin: 1, excluded at 4σ). The repeating per-zero pattern between 10⁸ and 10⁹ is shared logarithmic shot noise (fact 4028) | writeup/code/zeros_twin_spectrum.py |
| 2026-11-03 | Are the Dirichlet L-zeros mod 3, 5 visible in the twin spectrum? | Not resolved: expected 0.4 and 0.1 of the prime dip, 20 L-zeros each give ±0.3 | scratch (e7_lzeros.py) |
| 2026-11-03 | Cousins, sexy pairs, triplets, Sophie Germain pairs at the zeros | cousins 1.72 ± 0.17, sexy 1.65 ± 0.12 (predicted 2), triplets 2.04 ± 0.58 (predicted 3); Sophie Germain: fitted −0.907 − 0.799 cos(γ ln 2) against the predicted −1 − 0.707 cos(γ ln 2), corr −0.31 over 300 zeros, twins show none (0.017) | writeup/code/zeros_tuple_spectrum.py |
| 2026-11-03 | Fujii's formula for Goldbach: Σ_{N≤X} r(N) − X²/2 = −2 Σ_ρ X^(ρ+1)/(ρ(ρ+1))? | corr 0.993 and slope 1.016 over 400 values of X in [10⁴, 2·10⁶]; the factor 2 is exact here | writeup/code/zeros_goldbach_fujii.py |
