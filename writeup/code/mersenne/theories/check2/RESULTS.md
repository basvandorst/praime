# Independent check 2 (items A-E)

All scripts in this directory were written from scratch, without reading or importing the earlier theory scripts or their outputs.
Outputs: A_out.txt, B_out.txt, C_out.txt, D_out.txt, E_zeros_out.txt, E_race_out.txt.

## A. Wheel-24 bytes (a_bytes.py, a_hl.py), full range J = 41,666,666 (to 10^9)

Sieve: eight wheel arrays indexed by j, one per offset; each prime p from 5 to 31,623 strikes j = -o/24 mod p, starting at 24j+o >= p^2.
Check: the bits add up to 50,847,526 = pi(10^9) - 8. At J = 416,666 the bits add up to sympy's primepi(10^7) - 8 = 664,571.
HL: S(U) over all subsets. r = 2 and 3 give the factor 3^|U|. For r = 5..43, nu_r is the number of distinct o mod r in U. For r >= 47, nu = |U|; the product runs to 10^7, plus a tail term. Sanity check: S({-1,+1}) = 12*C2, to 3e-10.
E(T) by Moebius over supersets, with the sums of 1/ln(24j+9)^k over j.

| claim | ours | verdict |
|---|---|---|
| 27 impossible byte values, E = 0, never observed | 27, E exactly 0, observed 0 | AGREE |
| chi2 (229 values), J = 41,666,666 | 174.18 | AGREE |
| chi2 (224 values with E >= 5), J = 416,666 | 196.20 (all 229 possible: 215.08) | AGREE |
| entropy observed / predicted | 4.91877 / 4.91886 bits | AGREE |
| MI observed, plug-in | 0.03257 bits (41,666,665 pairs) | AGREE |
| MI predicted, 16-position HL | 0.03210 bits (L = 1/ln(24j+9) or 1/ln(24j+21), mixed over j; no negative cells) | AGREE |

Note: the predicted joint law evaluated at the single point j = J/2, without mixing over j, gives MI = 0.03047. The 0.03210 includes the mixture over j.

## B. Lane words, Mertens predictor (b_lanes.py)

prod_{r<=47}(1-1/r) = 0.138704; prod_{47<r<2e6}(1-1/(r-1)^2) = 0.99604630. Classes come from a sieve, with spot checks against sympy.

| range | words | log-loss | verdict |
|---|---|---|---|
| c = 1,000,002..1,999,998 | 166,667 | 0.8315 bits | AGREE |
| c = 9,000,000..9,999,996 | 166,667 | 0.8256 bits (same with the tail product to 2e7) | AGREE |

## C. Primes with prime order of 2 (c_order.py)

- A = C2 * sum m g / k^2. The sum factors as (1/2 + sum_{a>=2} 4^-a) * prod_{r odd} (1 + (r-1)/((r-2)(r^2-1))) = (7/12) prod_{r>=3} (1 + 1/((r-2)(r+1))).
  So A = (7/12) prod_{r>=3} r(r^2-r-1)/((r-1)^2(r+1)) = **0.538241**. The direct partial sum to k = 2e6 gives the same. AGREE.
- Share of 3|k: the Euler factor at r = 3 is 1 + 2/(1*8) = 5/4, and the terms with 3|k contribute 1/4 of it. The share is exactly **1/5** (proved; the numeric value is 0.200000). AGREE.
- N(x): Gauss-Legendre in u = ln t, 30 nodes per piece of length 0.5, vectorised over k; spot checks against mpmath.quad to 1e-8.
  N(10^6) = **3,841.6**, 3|k share **0.2202**; N(10^7) = **26,690.2**, 3|k share **0.2162**. AGREE.
- Direct count: for each prime q <= 10^7, test pow(2, l, q) == 1 for each odd prime l | q-1. Result: **26,560** primes, of which **5,727** have q = 1 mod 3. AGREE.
  Cross-check to 10^6 with sympy n_order: 3,845 = our count to 10^6.

## D. Goldbach mod 8, even n in [1,000,000, 2,000,000] (d_goldbach8.py)

Ordered pair counts by FFT convolution of the class indicator arrays mod 8, rounded. Spot checks by direct loops at 5 values of n. S(n) by shifted sieve sums, with spot checks.

| n mod 8 | count | mean Delta | mean S | verdict |
|---|---|---|---|---|
| 0 | 125,001 | -36.29 | 29.09 | AGREE |
| 2 | 125,000 | -37.77 | 29.00 | AGREE |
| 4 | 125,000 | -39.31 | 29.08 | AGREE |
| 6 | 125,000 | -25.75 | 29.14 | AGREE |

## E. Real characters mod P (e_zeros.py, e_race.py)

chi = (./P) = Kronecker (D/.), with D = P* = +-P (checked for n < 200). Zeros come from PARI 2.17.2 lfunzeros.
V is computed three ways:
1. mpmath: L'/L(1) = -ln P + sum chi(a) gamma_1(a/P) / sum chi(a) psi(a/P), from the Hurwitz expansion. No PARI.
2. PARI lfun(L,1,1)/lfun(L,1).
3. 2*sum_{0<gamma<400} 1/(1/4+gamma^2) over all PARI zeros, plus the tail integral of (1/pi) ln(Pt/2pi)/(1/4+t^2) from 400 to infinity. The zero counts match N_smooth(400) within about 2.

| P | first zero | V (route 1) | V (route 2) | V (route 3) | claim zero / V | verdict |
|---|---|---|---|---|---|---|
| 3 | 8.0397 | 0.11323 | 0.11323 | 0.11323 | 8.0397 / 0.113 | AGREE |
| 5 | 6.6485 | 0.15656 | 0.15656 | 0.15656 | 6.6485 / 0.157 | AGREE |
| 7 | 4.4757 | 0.25524 | 0.25524 | 0.25524 | 4.4757 / 0.255 | AGREE |
| 11 | 2.4772 | 0.50751 | 0.50751 | 0.50752 | 2.4772 / 0.508 | AGREE |
| 13 | 3.1193 | 0.39665 | 0.39665 | 0.39665 | 3.1193 / 0.397 | AGREE |
| 17 | 3.7281 | 0.38716 | 0.38716 | 0.38716 | 3.7281 / 0.387 | AGREE |
| 19 | 1.5161 | 1.05267 | 1.05267 | 1.05267 | 1.5161 / 1.053 | AGREE |
| 23 | 2.8713 | 0.56931 | 0.56931 | 0.56930 | 2.8713 / 0.569 | AGREE |
| 29 | 1.7938 | 0.87927 | 0.87927 | 0.87927 | 1.7938 / 0.879 | AGREE |
| 31 | 2.0350 | 0.80973 | 0.80973 | 0.80972 | 2.0350 / 0.810 | AGREE |
| 37 | 2.0572 | 0.81445 | 0.81445 | 0.81446 | 2.0572 / 0.814 | AGREE |
| 41 | 2.4753 | 0.68277 | 0.68277 | 0.68277 | 2.4753 / 0.683 | AGREE |
| 43 | 0.8364 | 2.53822 | 2.53822 | 2.53822 | 0.8364 / 2.538 | AGREE |
| 47 | 2.2174 | 0.84588 | 0.84588 | 0.84588 | 2.2174 / 0.846 | AGREE |

### Race, log-time share of x in [1000, 10^8] with the +1 class strictly ahead

The primes come from an odd-only sieve: pi(10^8) = 5,761,455. The counts are constant between consecutive primes, so each gap is weighted by ln(p_{k+1}/p_k), clipped to [1000, 10^8].

| P | chi(2) | incl. p = 2 (definition as stated) | odd primes only | claim |
|---|---|---|---|---|
| 3 | -1 | 0.00000 | 0.00000 | 0.000 |
| 5 | -1 | 0.00000 | 0.00000 | 0.000 |
| 7 | +1 | 0.00060 | 0.00035 | 0.000 |
| 11 | -1 | 0.04010 | 0.04760 | 0.048 |
| 13 | -1 | 0.02011 | 0.02307 | 0.023 |
| 17 | +1 | 0.00486 | 0.00450 | 0.004 |
| 19 | -1 | 0.11807 | 0.13696 | 0.137 |
| 23 | +1 | 0.06150 | 0.05342 | 0.053 |
| 29 | -1 | 0.11008 | 0.11973 | 0.120 |
| 31 | +1 | 0.10113 | 0.09071 | 0.091 |
| 37 | -1 | 0.07984 | 0.08415 | 0.084 |
| 41 | +1 | 0.05028 | 0.04801 | 0.048 |
| 43 | -1 | 0.20021 | 0.20413 | 0.204 |
| 47 | +1 | 0.10150 | 0.08817 | 0.088 |

**Verdict: AGREE only if p = 2 is left out of the counts.** The claimed race numbers are the odd-prime race: all 14 match to the third decimal (0.00450 against 0.004 for P = 17 is a rounding edge).
As stated ("primes p <= x, p != P"), p = 2 counts, because chi(2) = (2/P) = +-1. Then 11 of the 14 values differ (all P >= 7; at P = 7 only in the 4th decimal).
Including 2 shifts the walk by chi(2). It lowers the share when (2/P) = -1 (P = 3, 5 mod 8) and raises it when (2/P) = +1.
The write-up should say "odd primes p != P", or the numbers should be recomputed with p = 2 included.

Coordinator note on E3: the claims were for odd primes, and the write-up will say so. With that definition all 14 race values AGREE.

## F. Goldbach strands against the prime race, N = 4,000,000 (f_strands.py, F_out.txt)

Method:
- Ordered pair counts come from FFT convolution (length 2^23) of the class indicators mod 12 (1, 5, 7, 11) and mod 8 (1, 3, 5, 7). The rounding error is below 0.01. Two values of n per class were checked by a direct loop.
- pi_a(n) is a cumulative sum of the class indicator.
- S(n) and Q(n) are sums of shifted sieve arrays over the primes s >= s0, with s0 = 5 for mod 12 and 3 for mod 8.
- pred = S/Q * SIG. For n = 10 mod 12, S(n-2) replaces S(n) and Q(n) stays. Using Q(n-2) as well gives the same values to 4 decimals.
- 40 log bins with edges exp(linspace(ln(N/16), ln N, 41)); n >= lower edge, n < upper edge, Q(n) > 0. The last edge evaluates to 3999999.9999999986, so n = 4,000,000 itself is not in the last bin.

| class | corr (claim) | MAD pred (claim) | MAD -S (claim) | verdict |
|---|---|---|---|---|
| mod 12, n = 0 | 0.9930 (0.993) | 3.988 (3.99) | 19.411 (19.41) | AGREE |
| mod 12, n = 2 | 0.9930 (0.993) | 2.267 (2.27) | 11.065 (11.07) | AGREE |
| mod 12, n = 6 | 0.9893 (0.989) | 2.338 (2.34) | 10.493 (10.49) | AGREE |
| mod 12, n = 10 | 0.9897 (0.990) | 2.290 (2.29) | 10.865 (no claim) | AGREE |
| mod 8, n = 0 | 0.9859 (0.986) | 1.820 (1.82) | 8.920 (8.92) | AGREE |
| mod 8, n = 2 | 0.9871 (0.987) | 1.686 (1.69) | 8.354 (8.35) | AGREE |
| mod 8, n = 4 | 0.9879 (0.988) | 2.333 (2.33) | 12.231 (12.23) | AGREE |
| mod 8, n = 6 | 0.9890 (0.989) | 1.603 (1.60) | 7.019 (7.02) | AGREE |

## G. Mod-12 Goldbach strands, 60 log bins from 10^6 to 10^9 (g_sieve.py, g_bins.py, g_check.py)

Method:
- **Sieve.** A segmented odd-only sieve below 10^9 gives the sorted prime arrays of the classes 1, 5, 7 and 11 mod 12 (50,847,532 primes in all, which is pi(10^9) - 2) and a packed bit array.
- **Delta, exact.** For a class pair (a, b), the sum of #(a,b)(n) over a bin is G_ab(hi) - G_ab(lo), where G_ab(t) = #{(p, q): p = a, q = b mod 12, p + q < t} = sum over p of pi_b(t-1-p). This is one vectorised searchsorted per edge.
- **Checks on G.** G matches brute force at t = 20,000. Its symmetry G_ab = G_ba holds. A window of 15 n near 5.01e8 for n = 0 and n = 6 mod 12 matches direct per-n pair counting: sums -12,978 and -15,686, equal both ways.
- **Mean S(n), exact.** s^2 = 1 mod 12, so r = n - s^2 lies in class c - 1. For each prime s up to 31,622, the bin count is a count of class primes in [lo - s^2, hi - s^2). For n = 10 mod 12, S(n) = 0 identically (r = 9 mod 12), so S(n-2) uses the class 7.
- **pred, sampled.** pred = S*SIG/Q, from uniform random n of the class in each bin (seed 12345). Two runs: 2,000 and 20,000 n per bin and class. S(n) comes from bit lookups, pi_a(n) from searchsorted, Q(n) = #{primes 5 <= s, s^2 <= n}.

### Bin means of Delta (exact)

| bin start | n = 0 | n = 2 | n = 6 | n = 10 | claims (0 / 2 / 6) | verdict |
|---|---|---|---|---|---|---|
| 1,000,000 | -13.4526 | -27.2979 | -41.8708 | -14.2088 | -13.45 / -27.3 / -41.9 | AGREE (n = 10: no claim given) |
| 100,000,000 | -300.0555 | -260.4755 | -221.0979 | 39.5088 | -300.06 / -260.5 / -221.1 | AGREE |
| 501,187,234 | -751.0252 | -604.8795 | **-458.8903** | 146.2968 | -751.03 / -604.9 / **-537.6** | n = 0, 2 AGREE; **n = 6 DISAGREE** |
| 891,250,938 | -426.9025 | -275.1461 | -123.4930 | 151.7867 | -426.90 / -275.1 / -123.5 | AGREE |

The claimed -537.6 is the n = 6 mean of the previous bin, [446,683,592, 501,187,234), where we get -537.6015. Most likely the value was copied from the wrong row; the correct value for the bin starting at 501,187,234 is -458.89. Our sampled pred, -470.1 +- 1.6, is also near -459, not near -538.

### Mean of -S(n), n = 0 mod 12 (exact)

| bin start | ours | claim | verdict |
|---|---|---|---|
| 1,000,000 | -38.6158 | -38.62 | AGREE |
| 100,000,000 | -211.3965 | -211.40 | AGREE |
| 501,187,234 | -395.5298 | -395.53 | AGREE |
| 891,250,938 | -495.4729 | -495.47 | AGREE |

### pred bin means, sampled with 20,000 n per bin (mean +- standard error)

| bin start | n = 0 | n = 2 | n = 6 | n = 10 (S(n-2)) |
|---|---|---|---|---|
| 1,000,000 | -21.02 +- 0.08 | -32.66 +- 0.11 | -44.03 +- 0.15 | -11.58 +- 0.05 |
| 100,000,000 | -285.70 +- 0.95 | -250.14 +- 0.79 | -214.26 +- 0.70 | 35.66 +- 0.24 |
| 501,187,234 | -737.12 +- 2.39 | -601.88 +- 1.96 | -470.13 +- 1.62 | 132.97 +- 0.55 |
| 891,250,938 | -438.44 +- 1.50 | -304.02 +- 0.97 | -171.99 +- 1.04 | 133.05 +- 0.85 |

The 2,000-sample run, in the same order: -20.60 +- 0.25, -32.53 +- 0.34, -44.56 +- 0.49, -11.19 +- 0.17; -288.91 +- 3.07, -250.27 +- 2.55, -217.20 +- 2.26, 36.40 +- 0.77; -742.92 +- 7.55, -603.76 +- 6.11, -468.87 +- 5.14, 134.61 +- 1.78; -432.89 +- 4.71, -301.03 +- 3.03, -177.82 +- 3.40, 132.71 +- 2.64.

### Correlation over 60 bins: exact mean Delta against sampled mean pred

| class | 2,000 per bin | 20,000 per bin | claim | verdict |
|---|---|---|---|---|
| n = 0 | 0.9968 | 0.9972 | 0.997 | AGREE |
| n = 2 | 0.9977 | 0.9975 | 0.997 | AGREE within sampling noise (0.9975 is at the rounding edge) |
| n = 6 | 0.9966 | 0.9968 | 0.997 | AGREE |
| n = 10 (S(n-2)) | 0.9953 | 0.9954 | 0.995 | AGREE |

The large scratch arrays (G_bits.npy, G_p*.npy) were deleted after the run. The bin results are kept in G_results.npz (20,000-sample run) and G_results_2000.npz.
