# Two new prime2.html layouts: mersieve and gaptowers

N = 500,000 throughout (41,538 primes). All files are in this directory.

Tags show which script produced a number and which script rechecked it:

| tag | produced by | rechecked by |
|---|---|---|
| [A] | `ms_main.py` (numpy, vectorised over p, all k <= 1000) | `ms_check.py` (built-in `pow`, sympy primes, p outer / k inner): first-k table and full factor list identical |
| [B] | `ms_analyse.py` (numpy + gmpy2) on the [A] tables | `ms_analyse_check.py` (pure Python + `sympy.isprime`, its own C2, on the `ms_check.py` tables): every number identical |
| [C] | `ms_speed.c` timings (`ms_speed_run.py`, best of 3, gcc -O2, this machine) | every output table compared with `chk_first.txt` by `ms_speed_check.py` |
| [D] | `gt_main.py` (numpy) | `gt_check.py` (pure Python, sympy, `random`; histogram by explicit inclusion-exclusion) |
| [E] | `gt_five.py` (sympy) | count also in [D] |
| [F] | `ms_extra.py` (numpy, `ms_first.txt`, `ms_analyse.out`) | `ms_extra_check.py` (pure Python, sympy, `chk_first.txt`, `ms_analyse_check.out`) |
| [X] | `ms_explore.py`, `ms_explore2.py`: exploration only, not rechecked | - |

Outputs: `ms_analyse.out`, `ms_analyse_check.out`, `ms_check.out`, `gt_main.out`, `gt_check.out`, `ms_speed_run.out`, `ms_speed_check.out`. Data: `ms_rows.csv` (row counts k = 1..1000), `ms_first.txt`, `ms_pairs.txt`, `gt_gapcounts.csv`, `gt_gallagher_hist.txt`.

---

## Layout 1: mersieve

Definition used: for each prime p <= N, the smallest k in 1..K with q = 2kp+1 dividing M_p = 2^p - 1 and q < M_p; otherwise p is a survivor. The smallest k does not depend on K, so a row's count is the same for every K >= k.

### Two facts used throughout (proved)

- **Every divisor of M_p is 1 mod 2p** (p an odd prime). A prime r | M_p has ord_r(2) = p, so p | r - 1, and r is odd. Products of such r are also 1 mod 2p.
- **The smallest k always gives a prime q.** If q = 2kp+1 | M_p is composite, each prime factor r of q divides M_p, r = 2jp+1 with j < k, and r < q < M_p. So a page may skip composite q (for example q with a factor 3..13) without changing any row.

### Q1. Row counts [A, B]

| K | primes in rows | survivors | survivor fraction | non-empty rows |
|---:|---:|---:|---:|---:|
| 16 | 7,166 | 34,372 | 0.8275 | 12 |
| 64 | 9,805 | 31,733 | 0.7640 | 48 |
| 256 | 12,137 | 29,401 | 0.7078 | 192 |
| 1000 | 14,230 | 27,308 | 0.6574 | 694 |

All 31 Mersenne exponents <= N are survivors for every K (checked) [A, B]. For K = 16, 64 and 256 every row k that is not 2 mod 4 holds at least one prime (12 = 16 - 4, 48 = 64 - 16, 192 = 256 - 64).

Rows k = 1..64 (smallest k), rows not listed are empty [B]:

| k | count | k | count | k | count | k | count |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 2189 | 17 | 99 | 33 | 80 | 49 | 25 |
| 3 | 1309 | 19 | 79 | 35 | 64 | 51 | 48 |
| 4 | 838 | 20 | 199 | 36 | 142 | 52 | 43 |
| 5 | 541 | 21 | 133 | 37 | 37 | 53 | 25 |
| 7 | 264 | 23 | 65 | 39 | 89 | 55 | 40 |
| 8 | 412 | 24 | 230 | 40 | 88 | 56 | 60 |
| 9 | 367 | 25 | 79 | 41 | 27 | 57 | 54 |
| 11 | 172 | 27 | 99 | 43 | 34 | 59 | 20 |
| 12 | 489 | 28 | 98 | 44 | 62 | 60 | 119 |
| 13 | 124 | 29 | 58 | 45 | 78 | 61 | 19 |
| 15 | 292 | 31 | 51 | 47 | 26 | 63 | 47 |
| 16 | 169 | 32 | 78 | 48 | 102 | 64 | 42 |

Rows 65..256 hold 2,332 primes together, rows 257..1000 hold 2,093 [B]. Row 1 is 22.3% of the primes drawn at K = 64 [B]. Rows 12, 24, 36, 48 and 60 are fuller than their neighbours: k = 0 mod 4 admits primes p of both classes mod 4 (Q2), and when 3 | k the number q = 2kp+1 is never divisible by 3, so it is prime more often (see Q4).

### Q2. The mod-4 law [B]: zero exceptions

Checked on all 17,199 factors q = 2kp+1 <= 2000p+1 of all M_p (not only the smallest): q is never 3 or 5 mod 8 (10,352 are 7 mod 8, 6,847 are 1 mod 8); no factor has k = 2 mod 4; k = 1 mod 4 occurs only with p = 3 mod 4 (5,843 factors); k = 3 mod 4 only with p = 1 mod 4 (4,509 factors); k = 0 mod 4: 6,847 factors, both classes.

Proof. (1) q | 2^p - 1 with p odd gives 2 = 2^(p+1) = (2^((p+1)/2))^2 mod q, so 2 is a square mod q, so q = ±1 mod 8 (second supplement to quadratic reciprocity). (2) q = 1 + 2kp mod 8 is ±1 mod 8 exactly when kp = 0 or 3 mod 4; for odd p this kills k = 2 mod 4 (kp = 2), and for odd k it forces kp = 3, that is p = 3 mod 4 when k = 1 mod 4 and p = 1 mod 4 when k = 3 mod 4.

### Q3. Row 1 is Euler's rule [B]

Row 1 = {p > 3 : p = 3 mod 4 and 2p+1 prime}: both sets have 2,189 elements and are equal. Proof both ways: if p = 3 mod 4 and q = 2p+1 is prime, then q = 7 mod 8, 2 is a square mod q, and Euler's criterion gives 2^p = 2^((q-1)/2) = 1 mod q. Conversely, 2p+1 | M_p forces 2p+1 prime (its prime factors are 1 mod 2p, hence >= 2p+1) and p = 3 mod 4 (Q2 with k = 1). p = 3 is excluded because 7 = M_3 is not a proper factor. (Consequence: an exponent p = 3 mod 4 needs 2p+1 composite; true for all 51 known exponents [B].)

### Q4. How the survivor fraction depends on K and p

**Model (derivation).** Let q = 2kp+1.
1. q can divide M_p only if q = ±1 mod 8 (Q2).
2. If q is prime, its unit group is cyclic of order 2kp, and q | M_p means 2 lies in the subgroup of 2k-th powers (index 2k). Given that 2 is already a square (step 1), the chance is 1/k (heuristic: 2 behaves like a random square; for k = 1 it is certain, which is Euler's rule). Tested directly below.
3. Exact-q model: P(survive K) = product over k <= K with q prime, q = ±1 mod 8, q < M_p, of (1 - 1/k). Uses the true primality of each q.
4. Averaged model: replace [q prime] by its Hardy-Littlewood probability. For p ranging over primes, q is odd (factor 2), never divisible by an odd prime r | k (factor r/(r-1)), and divisible by an odd r not dividing k with chance 1/(r-1) instead of 1/r (factor r(r-2)/(r-1)^2). So P(q prime) ≈ 2·C2·prod_{odd r | k} (r-1)/(r-2) / ln q, with C2 = 0.66016186 the twin prime constant (computed [B] as a product over primes < 10^6).
5. Crude model: the weight [q = ±1 mod 8]·2·C2·prod(...) averages to 1/2 · 2 = 1 over k (k mod 4 admits half the k; the average of prod_{odd r|k}(r-1)/(r-2) is prod_r (1 + 1/(r(r-2))) = 1/C2). So the expected number of factors is about sum_{k<=K} 1/(k ln(2kp)) ≈ ln(ln(2Kp)/ln(2p)), and P(survive) ≈ exp(-that) = ln(2p)/ln(2Kp).

**Direct test of the 1/k rule** (admissible = q prime, q = ±1 mod 8, q < M_p; all p <= N) [B]:

| k | admissible q | q divides M_p | expected (sum of 1/k) | sd |
|---|---:|---:|---:|---:|
| 1 | 2,189 | 2,189 | 2,189.0 | 0 |
| 3 | 4,028 | 1,309 | 1,342.7 | 29.9 |
| 4 | 3,911 | 962 | 977.8 | 27.1 |
| 5 | 2,658 | 541 | 531.6 | 20.6 |
| 7 | 2,239 | 293 | 319.9 | 16.6 |
| 8 | 3,733 | 446 | 466.6 | 20.2 |
| 9 | 3,736 | 421 | 415.1 | 19.2 |
| 11-12 | 9,260 | 750 | 786.8 | 26.8 |
| 13-16 | 10,345 | 712 | 695.0 | 25.5 |
| 17-32 | 38,074 | 1,547 | 1,596.3 | 39.1 |
| 33-64 | 82,809 | 1,748 | 1,770.2 | 41.6 |
| 65-256 | 450,098 | 3,221 | 3,238.5 | 56.7 |
| 257-1000 | 1,623,173 | 3,060 | 2,988.5 | 54.6 |

(k = 2 mod 4 has no admissible q.) All 17,199 factors are admissible, as Q2 requires. Exploratory check on a bigger set [X, `ms_explore2.py`]: for 22 values of k up to 64 and all primes q < 4·10^7 in the class q = 1+2k mod 4k, ±1 mod 8, the fraction with 2 a 2k-th power is 1/k within |z| < 2.9 (3 of 44 tests have |z| > 2, about what chance gives).

**Survivors: observed vs the three models** [B] (fractions in brackets; z for the exact-q model, sd from the sum of P(1-P)):

| K | p band | primes | observed | exact-q model | z | averaged (HL) | crude ln(2p)/ln(2Kp) |
|---:|---|---:|---:|---:|---:|---:|---:|
| 16 | p < 10^4 | 1,229 | 946 (0.7697) | 931.5 (0.7579) | +1.23 | 938.5 (0.7636) | 928.8 (0.7558) |
| 16 | 10^4-10^5 | 8,363 | 6,708 (0.8021) | 6,783.1 (0.8111) | -2.65 | 6,790.6 (0.8120) | 6,727.3 (0.8044) |
| 16 | 10^5-5·10^5 | 31,946 | 26,718 (0.8363) | 26,558.1 (0.8313) | +3.02 | 26,636.5 (0.8338) | 26,397.4 (0.8263) |
| 16 | all | 41,538 | 34,372 (0.8275) | 34,272.7 (0.8251) | +1.62 | 34,365.7 (0.8273) | 34,053.5 (0.8198) |
| 64 | p < 10^4 | 1,229 | 838 (0.6819) | 837.1 (0.6811) | +0.07 | 843.1 (0.6860) | 828.5 (0.6742) |
| 64 | 10^4-10^5 | 8,363 | 6,149 (0.7353) | 6,205.3 (0.7420) | -1.63 | 6,218.1 (0.7435) | 6,128.4 (0.7328) |
| 64 | 10^5-5·10^5 | 31,946 | 24,746 (0.7746) | 24,540.4 (0.7682) | +3.16 | 24,621.4 (0.7707) | 24,288.5 (0.7603) |
| 64 | all | 41,538 | 31,733 (0.7640) | 31,582.8 (0.7603) | +2.00 | 31,682.6 (0.7627) | 31,245.5 (0.7522) |
| 256 | p < 10^4 | 1,229 | 757 (0.6159) | 757.0 (0.6160) | -0.00 | 763.0 (0.6208) | 748.0 (0.6086) |
| 256 | 10^4-10^5 | 8,363 | 5,641 (0.6745) | 5,705.5 (0.6822) | -1.68 | 5,718.3 (0.6838) | 5,627.5 (0.6729) |
| 256 | 10^5-5·10^5 | 31,946 | 23,003 (0.7201) | 22,747.7 (0.7121) | +3.51 | 22,830.1 (0.7146) | 22,491.9 (0.7041) |
| 256 | all | 41,538 | 29,401 (0.7078) | 29,210.2 (0.7032) | +2.28 | 29,311.3 (0.7056) | 28,867.4 (0.6950) |
| 1000 | p < 10^4 | 1,229 | 686 (0.5582) | 691.7 (0.5628) | -0.36 | 697.2 (0.5673) | 682.8 (0.5556) |
| 1000 | 10^4-10^5 | 8,363 | 5,203 (0.6221) | 5,282.5 (0.6316) | -1.96 | 5,294.9 (0.6331) | 5,209.2 (0.6229) |
| 1000 | 10^5-5·10^5 | 31,946 | 21,419 (0.6705) | 21,211.4 (0.6640) | +2.67 | 21,289.1 (0.6664) | 20,967.6 (0.6563) |
| 1000 | all | 41,538 | 27,308 (0.6574) | 27,185.5 (0.6545) | +1.38 | 27,281.2 (0.6568) | 26,859.6 (0.6466) |

Reading: over all p, all three models get the survivor fraction to within 1.3 percentage points (within 1.6 in each band); the crude formula is always a little low overall (by 0.8 to 1.3 points). The fraction rises with p (bigger p: the first possible factors 2p+1, 6p+1, ... are larger, so less likely prime) and falls slowly with K (6.4, 5.6 and 5.0 points from K = 16 to 64 to 256 to 1000). Not explained: the exact-q model leaves a surplus of survivors among p > 10^5 of about 3 sd (z = +3.0 to +3.5; the K values are nested, so this is one effect, not four). [X, `ms_explore.py`] shows it is a deficit of factors with 3 <= k <= 64 in that band (6,085 found, 6,328 expected, z = -3.3), not clustering of several factors on the same p; in the band 10^4-10^5 the deviation goes the other way (z = +1.4). Summed over all k <= 1000 the 1/k table gives 17,199 factors against 17,318.0 expected (z = -1.0), and for 3 <= k <= 64 it gives 8,729 against 8,902.0 (z = -2.0) [F]. I treat it as a fluctuation plus model error I could not pin down.

### Q5. Classes and the exponents [B]

| K | class | primes | survivor fraction | exact-q model |
|---:|---|---:|---:|---:|
| 64 | p = 1 mod 4 | 20,731 | 0.7907 | 0.7852 |
| 64 | p = 3 mod 4 | 20,806 | 0.7373 | 0.7355 |
| 1000 | p = 1 mod 4 | 20,731 | 0.6788 | 0.6764 |
| 1000 | p = 3 mod 4 | 20,806 | 0.6361 | 0.6326 |

Cause: p = 3 mod 4 can be hit at k = 1 (2p+1, certain when prime), p = 1 mod 4 first at k = 3 (6p+1, chance 1/3). Mod 8 adds nothing beyond mod 4 (K = 64: 0.7898, 0.7385, 0.7916, 0.7360 for p = 1, 3, 5, 7 mod 8; `ms_analyse.py`, rechecked in `ms_extra_check.py`).

The exponents: 17 of the 30 odd exponents <= N are 1 mod 4 (binomial P(X >= 17 | 30, 1/2) = 0.29); of all 50 odd known exponents 31 are 1 mod 4 (P(X >= 31 | 50, 1/2) = 0.06). The direction matches the class effect above, but it is not significant. Mod 8 counts of the 50: 19, 9, 12, 10 for 1, 3, 5, 7. Neighbours [F]: of the 42 primes just before or after the 21 exponents in [100, N] (none of them an exponent, no repeats), 32 are survivors at K = 64, against 32.1 expected at the overall rate 76.4% and 29.8 at the rate of each one's band. The other 10 sit in rows 1, 1, 3, 3, 12, 12, 15, 28, 36, 57 (primes 216071, 131, 601, 4441, 4243, 11197, 113, 2207, 44491, 2179): ordinary rows. Conclusion: in this picture "survivor" is the only thing that marks an exponent; there is no pattern in the rows around them.

### Q6. Build time [C]

Cost = number of modpows times ~19 squarings (p < 2^19; the "times 2" step is an add and a compare, no mulmod). Modpow counts for the straightforward loop (all k up to the first hit) [B, C]: K = 16: 587,028; 64: 2,162,388; 94: 3,104,158; 128: 4,150,518; 256: 7,975,502; 1000: 28,891,764 (C counts are a few hundred lower: the C loop stops early for tiny p, where q reaches M_p). Skipping q not ±1 mod 8 halves these; also skipping q with a factor 3..13 (safe, see the facts above) cuts them to about 19% (K = 64: 409,202; K = 256: 1,523,620).

**Exactness limit.** With plain JavaScript numbers, a·b is exact only below 2^53. Here a, b < q <= 2Kp+1, so plain `a*b % q` is exact for all p <= 500,000 only when K <= 94 (2^26.5 / 10^6 = 94.9). Confirmed [C]: at K = 256 the plain-double versions give wrong rows for 418 primes (fmod) and 443 primes (floor), first at p = 211,493 (fmod) and p = 203,383 (floor); all exact versions agree with the Python tables.

C timings, best of 3, seconds (filter 0 = every k, 1 = skip q not ±1 mod 8, 2 = also skip q with a factor 3..13); one other job was using the second core:

| K | mulmod | filter 0 | filter 1 | filter 2 |
|---:|---|---:|---:|---:|
| 64 | `fmod(a*b, q)` (what JS `%` on large doubles does) | 0.845 | 0.461 | 0.239 |
| 64 | `t = a*b; r = t - q*floor(t/q)`, one correction | 0.408 | 0.227 | 0.142 |
| 64 | uint64 `%` (reference) | 0.405 | 0.230 | 0.143 |
| 64 | Math.imul trick (below) | 0.529 | 0.288 | 0.164 |
| 94 | floor version | 0.586 | 0.326 | 0.211 |
| 256 | split into 16-bit halves, two fmods | 5.744 | 3.065 | 1.430 |
| 256 | uint64 `%` (reference) | 1.525 | 0.849 | 0.540 |
| 256 | Math.imul trick | 1.989 | 1.368 | 0.848 |
| 1000 | uint64 `%` / Math.imul trick | - | - | 2.243 / 2.368 |

Math.imul trick (exact for q < 2^30, so for every K <= 1073 at this N): `qh = Math.floor(a*b/q); r = (Math.imul(a,b) - Math.imul(qh,q)) | 0; if (r < 0) r += q; else if (r >= q) r -= q;`. The double quotient is off by at most 1, so the true remainder lies in [-q, 2q), which fits in 32 signed bits, and the low 32 bits of a·b - qh·q determine it.

The per-run output tables were deleted after `ms_speed_check.py` validated them (37 exact runs OK, the 6 plain-double runs at K >= 256 WRONG; see `ms_speed_check.out`); rebuild with `gcc -O2 -o ms_speed ms_speed.c -lm` and `ms_speed_run.py`.

Recommendation: default K = 64 with filter 1 or 2 and the floor mulmod (0.14 to 0.23 s in C). K = 256 needs the imul trick and takes 0.85 s in C even with filter 2, so it is likely over 1 s in a browser. Not measured: JavaScript itself (no browser or node timing was made here); V8 usually runs such double loops at 1 to 2 times the C time, but that is an estimate.

### Note text for mersieve (K = 64; every number from [A]/[B])

Row k holds the primes p for which 2kp+1 is the smallest proper factor of 2^p − 1 of that form, and the top row holds the survivors with no such factor up to k = 64: 9,805 of the 41,538 primes below 500,000 get a row and 31,733 (76.4%) survive, among them all 31 Mersenne exponents in range. Every factor of 2^p − 1 is 1 or 7 mod 8, because 2 = 2^(p+1) is a square modulo it; so rows 2, 6, 10, … stay empty, rows 1, 5, 9, … hold only primes 3 mod 4 and rows 3, 7, 11, … only primes 1 mod 4. Row 1 is Euler's rule: 2p+1 divides 2^p − 1 exactly when p > 3 is 3 mod 4 and 2p+1 is prime, which happens for 2,189 primes, the fullest row. A prime q = 2kp+1 that is 1 or 7 mod 8 divides 2^p − 1 with chance about 1/k, so the rows get sparser as k grows and the survivors fall only slowly as K grows (82.7% at K = 16, 70.8% at 256, 65.7% at 1000; the rule of thumb ln(2p)/ln(2Kp) is within 1.3 percentage points). Primes 1 mod 4 survive more often than primes 3 mod 4 (79.1% against 73.7%), because their first possible factor is 6p+1, which divides only one time in three, not 2p+1, which divides whenever it is prime.

---

## Layout 2: gaptowers

### Q1. Gap distribution up to 500,000 [D]

41,537 gaps between consecutive primes <= N (mean 12.04; ln N = 13.12). The last column, p = 499,979, is cut at N: it shows 22 numbers, its true gap to 500,009 is 30. With w = 1415 the 41,538 columns fill 30 bands.

| gap | count | gap | count | gap | count | gap | count |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 26 | 539 | 52 | 29 | 76 | 2 |
| 2 | 4565 | 28 | 609 | 54 | 62 | 78 | 4 |
| 4 | 4558 | 30 | 861 | 56 | 15 | 82 | 1 |
| 6 | 7474 | 32 | 237 | 58 | 21 | 84 | 1 |
| 8 | 3081 | 34 | 267 | 60 | 39 | 86 | 2 |
| 10 | 3842 | 36 | 331 | 62 | 8 | 90 | 1 |
| 12 | 4237 | 38 | 163 | 64 | 7 | 96 | 1 |
| 14 | 2281 | 40 | 195 | 66 | 16 | 100 | 1 |
| 16 | 1512 | 42 | 200 | 68 | 5 | 112 | 1 |
| 18 | 2512 | 44 | 70 | 70 | 11 | 114 | 1 |
| 20 | 1214 | 46 | 67 | 72 | 2 | | |
| 22 | 1101 | 48 | 67 | 74 | 1 | | |
| 24 | 1278 | 50 | 44 | | | | |

Jumping champion: 6 (7,474 times), then 2 (4,565) and 4 (4,558). Twin pairs (gap 2, from (3,5) on): 4,565. Multiples of 6 stand out (6, 12, 18, 24, 30 are each fuller than their neighbours), as the singular series 2·C2·prod_{r|g}(r-1)/(r-2) predicts.

Maximal gaps (gap: first prime) [D]: 1: 2, 2: 3, 4: 7, 6: 23, 8: 89, 14: 113, 18: 523, 20: 887, 22: 1129, 34: 1327, 36: 9551, 44: 15683, 52: 19609, 72: 31397, 86: 155921, 96: 360653, 112: 370261, 114: 492113.

### Q2. Gaps after the Mersenne exponents [D]

Gap after each of the 31 exponents <= N, as gap / ln p: 2: 1 (1.44), 3: 2 (1.82), 5: 2 (1.24), 7: 4 (2.06), 13: 4 (1.56), 17: 2 (0.71), 19: 4 (1.36), 31: 6 (1.75), 61: 6 (1.46), 89: 8 (1.78), 107: 2 (0.43), 127: 4 (0.83), 521: 2 (0.32), 607: 6 (0.94), 1279: 4 (0.56), 2203: 4 (0.52), 2281: 6 (0.78), 3217: 4 (0.50), 4253: 6 (0.72), 4423: 18 (2.14), 9689: 8 (0.87), 9941: 8 (0.87), 11213: 26 (2.79), 19937: 12 (1.21), 21701: 12 (1.20), 23209: 18 (1.79), 44497: 4 (0.37), 86243: 6 (0.53), 110503: 24 (2.07), 132049: 10 (0.85), 216091: 12 (0.98).

Test: S = sum of gap/ln p. Null: each exponent replaced by a uniformly random prime among the 101 primes around it by index (fewer at the start), so small exponents are compared with small primes.

| set | S | null mean | null sd | z | two-sided p (simulation) | mean percentile of gap in its window |
|---|---:|---:|---:|---:|---:|---:|
| all 31 exponents | 36.423 | 32.728 | 3.792 | +0.97 | 0.327 (2·10^5 draws; recheck 0.329 with 10^5 draws, other RNG) | 0.476 |
| the 21 exponents >= 100 | 21.249 | 21.474 | 3.294 | -0.07 | 0.990 (recheck 0.989) | 0.502 |

The gaps after the exponents are ordinary. The mean of gap/ln p over all primes 100..N is 1.0016, and the percentile column is from `gt_main.py` only.

### Q3. Gallagher: primes in windows of length ln n [D]

Windows: start at n = 100,000, length L = round(ln n), next window at n + L, last window ends at 499,994. That gives 31,848 windows (14,029 of length 12, 17,819 of length 13) holding 31,946 primes; mean 1.0031 per window (expected L/ln n averaged over windows: 1.0040, `gt_main.py` only), variance 0.6938.

| primes in window | observed | Poisson(1) | Hardy-Littlewood prediction |
|---:|---:|---:|---:|
| 0 | 9,495 | 11,716.2 | 9,473.7 |
| 1 | 14,180 | 11,716.2 | 14,172.4 |
| 2 | 6,835 | 5,858.1 | 6,877.9 |
| 3 | 1,258 | 1,952.7 | 1,248.9 |
| 4 | 78 | 488.2 | 74.0 |
| 5 | 2 | 97.6 | 1.0 |
| 6+ | 0 | 18.9 [F] | 0.0 |

Chi-square on bins 0, 1, 2, 3, >= 4 (4 df): 1,804.6 against Poisson(1), 0.72 against Hardy-Littlewood. Poisson with the observed mean is no better (1,800.8 on bins 0..4, >= 5; `gt_main.py` only).

Hardy-Littlewood prediction (derivation): average over the window start, the chance that the numbers at offsets T ⊂ {0, ..., h-1} are all prime is S(T)/L^|T| with S(T) = prod_p (1 - ν_p(T)/p)/(1 - 1/p)^|T| (ν_p = number of residues mod p that T covers; S = 0 if T covers all residues mod some p). Inclusion-exclusion gives P(exactly j primes) = sum_t (-1)^(t-j) C(t, j) A_t / L^t, A_t = sum of S(T) over t-subsets. Used per window with h = 12 or 13 and L = ln of the window's midpoint. A_2 for h = 12 is 48.41 = 10·S(2) + 8·S(4) + 6·S(6) + 4·S(8) + 2·S(10) with S(2) = S(4) = S(8) = 1.3203, S(6) = 2.6406, S(10) = 1.7604. The variance it predicts is 0.6917 (from the pair sum, and again from the full histogram), against 0.6938 observed and 1.0040 for Poisson. The recheck caught a bug in the first version of `gt_main.py` (one prime factor of S(T) skipped); the numbers above are after the fix, and both scripts agree.

Largest window counts [E]: two windows hold 5 primes, both prime quintuplets with offsets 0, 4, 6, 10, 12 (starting at 340,927 and 375,247). A window of 12 numbers can never hold 5 primes this high up (A_5 = 0 for h = 12: every 5-subset of 12 consecutive numbers covers all residues mod 2, 3 or 5).

### Note text for gaptowers (every number from [D]/[E])

Each column starts at a prime and stacks the numbers up to the next prime, so its height is the prime gap: up to 500,000 there are 41,537 complete gaps, and the most common height is 6 (7,474 columns), ahead of 2 (4,565 twin pairs) and 4 (4,558). Heights that are multiples of 6 stand out, and the record columns (maximal gaps) start at 2, 3, 7, 23, 89, 113, 523, 887, 1129, 1327, 9551, 15683, 19609, 31397, 155921, 360653, 370261 and 492113, where the tallest column of this picture, 114, starts. The columns of the 31 Mersenne exponents below 500,000 are ordinary: measured in units of ln p their heights add up to 36.4, against 32.7 ± 3.8 for random primes of the same size. Cut 100,000 to 500,000 into 31,848 windows of length ln n (12 or 13 numbers): Gallagher's theorem says the number of primes per window tends to a Poisson law with mean 1, but here 29.8% of the windows hold no prime (Poisson: 36.8%) and 44.5% hold exactly one (36.8%). The Hardy–Littlewood prime-pattern constants predict the actual counts 9,495 / 14,180 / 6,835 / 1,258 / 78 / 2 almost exactly (9,474 / 14,172 / 6,878 / 1,249 / 74 / 1); the two windows with 5 primes are prime quintuplets starting at 340,927 and 375,247.

---

## What I could not verify

- JavaScript speed: only C timings were made; the browser figure is an estimate (1 to 2 times C).
- The 3-sd survivor surplus for p > 10^5 in the exact-q model (Q4) is unexplained.
- `ms_explore.py` and `ms_explore2.py` are exploratory and have no second script.
- The p-values in Q2 (gaps) are by simulation (two independent simulations agree to 0.002).
