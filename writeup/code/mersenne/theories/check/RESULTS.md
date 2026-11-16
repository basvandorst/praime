# Independent check of items A-H

All scripts are in this directory and were written from scratch. None reads or imports anything else in `mexp/theory/`. Each ran on one core under `ulimit -v 1500000; timeout 900 nice`.

| Item | Verdict | Our value | Method (script) |
|---|---|---|---|
| A | AGREE | counts 13,517 / 13,615 / 13,617 / 13,610 (total 54,359); survivors 10,796 / 10,423 / 11,054 / 8,922; shares 0.7987 / 0.7656 / 0.8118 / 0.6555 | `pow(2, p, 2kp+1) == 1` for every k = 1..100, with no mod-8 or primality filter; primes from sympy.primerange ([a_c.py](a_c.py)) |
| B (relative, p = 10^6, K = 3000) | AGREE | 1.0488 / 1.0073 / 1.0645 / 0.8794 (absolute 0.6844 / 0.6574 / 0.6947 / 0.5739, mean 0.65258) | exact product over 2,880 residue combinations mod 5, 7, 11, 13; tail over r > 13 computed directly to 2·10^7 (0.9822987) and matched to C2/∏(3..13); r \| k correction (r-1)/(r-2) ([b_model.py](b_model.py)) |
| B (absolute, p = 2.6e5, K = 100) | CLOSE, not exact | 0.7903 / 0.7555 / 0.8051 / 0.6489 (claimed ≈ 0.795 / 0.761 / 0.809 / 0.657) | same model. The claimed numbers are what the model gives at p ≈ 4·10^5: 0.7954 / 0.7614 / 0.8097 / 0.6578 ([b_var.py](b_var.py), [b_var2.py](b_var2.py)) |
| C | AGREE | 11 / 20 / 13 / 5 (49 exponents > 3, list length 51) | direct count ([a_c.py](a_c.py)) |
| D | AGREE | byte 223: 14, byte 251: 14, byte 255: 0. First 251 at j = 1,875,108 (p = 45,002,591). Up to j = 416,666: 223 occurs 5 times, 251 never. HL: C7 = 323.83 for both patterns, prediction 13.218 (to 10^9) and 1.081 (to 10^7) | (1) our own segmented sieve over all integers, aligned to 24, full byte histogram ([de_sieve.py](de_sieve.py)); (2) separate enumeration: p = 24j-1 filtered by residues mod 5..997, then gmpy2.is_prime on the 7 members and on the excluded slot. Both give the same 14 + 14 positions ([d_tuples.py](d_tuples.py)) |
| E | AGREE (pairs with p+2 ≤ 10^9) | up 1,711,773, down 1,712,729, sum 3,424,502. All twin pairs below 10^9: 3,424,506. N(p² − 2 prime, 5 ≤ p ≤ 31,622) = 651 | full sieve to 10^9 (not reduced to 10^8) ([de_sieve.py](de_sieve.py), [d_tuples.py](d_tuples.py)) |
| F | AGREE | 21 / 105 / 585 / 3,845; with q = 1 mod 3: 852; model sum 3,828.68 | smallest-prime-factor sieve on q − 1, then pow(2, l, q) == 1 for each odd prime l \| q − 1 ([f_g.py](f_g.py)); cross-check to 10^5 with sympy.n_order gives 585 ([f_cross.py](f_cross.py)) |
| G | AGREE | 1,593; 529 with 6p+1 \| 2^p − 1; 529 with 6p+1 = x² + 27y²; the two sets are identical | pow, plus a brute-force search over y ([f_g.py](f_g.py)) |
| H | AGREE | mean Δ = −66.3685 / −65.8531 / −65.3337 / +0.4585; mean S = 55.5926 / 55.4568 / 55.6308 / 0 | exact: counted ordered prime pairs with p+q in [2·10^6, 4·10^6] with searchsorted, split by class (no FFT). 300 random n per class, counted directly: Δ −75.0±10.5, −60.4±7.4, −47.2±9.7, +2.1±7.8; S 59.1±1.5, 54.0±1.3, 55.9±1.5, 0 (all within 2 SE) ([h_goldbach.py](h_goldbach.py)) |

## Notes

- **B absolute.** The claimed "roughly" values do not come from p = 2.6·10^5. They match this model at p ≈ 4·10^5 to within 0.001. At p = 2.6·10^5 every class is lower, by 0.004 to 0.008. The relative values at p = 10^6 agree to 3 decimals.
- **E boundary.** 10^9 + 7 and 10^9 + 9 are both prime. They form an "up" pair (10^9 + 7 = 11 mod 12). If the upper limit really is about 10^9 + 20, up becomes 1,711,774 and the sum 3,424,503. The claimed numbers count pairs with p + 2 ≤ 10^9. Counted the same way, up_all + down_all + 1 for (3,5) = 3,424,506, which matches the known π₂(10^9).
- **D, the excluded slot.** The 8-tuple 0, 2, 6, 8, 12, 14, 18, 20 covers every residue mod 5, so it is inadmissible. So for p > 5, a byte of 223 or 251 is exactly an admissible prime septuplet with p = 23 mod 24, and byte 255 never occurs (the sieve also finds 0).
- **D, the 251 positions.** j = 1875108, 3901083, 5196503, 6443133, 9213348, 12532013, 13211433, 20634723, 21227903, 30384638, 30672338, 34798908, 35800048, 37955663.
- **D, the 223 positions.** j = 235, 3700, 35655, 273375, 314185, 795400, 3300350, 5169420, 6008930, 11462735, 17408920, 19300285, 29442515, 35325525.
