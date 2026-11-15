# Polyspiral and chebrace: the maths checked

All scripts and outputs are in this directory. Each number carries a tag "made by → rechecked by":

| tag | script (output) | what it does |
|---|---|---|
| M1 | `spiral_main.py` (`spiral_main.txt`, `.json`) | layout 1: numpy sieve; Bateman–Horn C by Legendre symbols |
| C1 | `spiral_check.py` (`spiral_check.txt`) | layout 1, independent: rays rebuilt from T(m), `sympy.isprime`, C by brute-force root counting mod every p < 10^4, reducibility by `sympy.factor_list` |
| F1 | `spiral_factors.py` (`spiral_factors.txt`) | local factors p ≤ 13 of the bright rays, brute force (their full products are the C of M1 = C1) |
| M2 | `race_main.py` (`race_main.txt`, `.json`) | layout 2: numpy, sampled on a geometric grid of 10^5 points |
| C2 | `race_check.py` (`race_check.txt`) | layout 2, independent: pure-Python sieve, squares by the local (CRT) criterion, prime-by-prime walk, exact log-scale integrals with `mpmath.li`, exact prime-power identity in Fractions |
| T2 | `race_theory.py` (`race_theory.txt`, `.json`) | zero variances V(χ) via L′/L(1, χ) (mpmath) and the Gaussian bias estimate |
| Z2 | `race_check_zeros.py` (`race_check_zeros.txt`) | rechecks T2: L(1), L′(1) by the Dirichlet series; V by summing the actual zeros (4 characters) |

Every comparison between a main script and its check agreed exactly (integers) or to the stated tolerance.

---

## Layout 1: polyspiral (N = 500,000)

### Set-up (derived; checked M1 → C1)

* Ring m holds P·m numbers; corner j of ring m is c_j(m) = T(m−1) + j·m = (P/2)m² + (j − P/2)m + 1 (C1 checks this for P = 3, 4, 6, 40, m < 50).
* The ray (P, j, s) is f(m) = P·m(m−1)/2 + j·m + 1 + s. It is integer-valued for every P, j, s, because m(m−1)/2 is an integer. Rays are counted from m = 1 with 2 ≤ f(m) ≤ N. There are 38 values of P, Σ P = 817 corners and 13 shifts, so 10,621 rays.
  *Caveat:* for m < |s| the number c_j(m) + s has not reached the straight part yet (it wraps round the previous corner). That is at most 6 points per ray, and they are included as the brief asked.
* j = P, s = 0 is the centred P-gonal numbers 1 + P·m(m+1)/2. For P = 6 they are 3m² + 3m + 1 = (m+1)³ − m³, so this ray carries exactly the **cuban primes**: 126 of them up to 500,000 (C1 counts the cuban primes directly; the ray count in M1 is 126).

### The Bateman–Horn constant: three corrections to the brief

1. **p = 2 needs m mod 4, not m mod 2.** For odd P, f(m) mod 2 is c, j + c, 1 + c, 1 + j + c for m ≡ 0, 1, 2, 3 (mod 4), with c = 1 + s. So exactly 2 of every 4 values are odd, and the 2-factor is exactly 1 for **every** ray with odd P (F1, C1). For even P = 2q, f = q·m² + (j − q)·m + c is always odd when c is odd and j is even, giving a 2-factor of 2. As a result, the best C over odd P is 2.861, and over even P it is 5.911 (M1 → C1).
2. **Reducible rays.** When disc(2f) = (2j − P)² − 8P(1 + s) is a square, f factors over Q and Bateman–Horn does not apply. The truncated product then only creeps towards 0 as the prime bound grows.
   * For s = −1 at every corner, 2f = m·(P(m−1) + 2j). **The line one step before each corner is composite from m = 3 on.** That is 817 rays.
   * 720 other rays have a square discriminant, which includes 0 (for example 4m² + 4m + 1 = (2m+1)² for P = 8).
   * In total 1,537 reducible rays hold 731 primes, at most 2 per ray (M1 → C1).

   They are left out of the test and the rankings below. (With them left in, the z-scores would have mean −0.243 and sd 0.734: M1 → C1.)
3. **C = 0 rays.** 1,713 rays are always even or always divisible by an odd p | P. They hold 115 primes in total, at most 1 per ray (M1 → C1).

How C is computed: for odd p | P, f is linear mod p, so ω ∈ {0, 1, p}. For p ∤ 2P, ω = 1 + (D/p). The product runs over p < 10^4. Extending it to p < 10^6 changes C by at most 0.49 % on all 63 reported rays (M1 only, not rechecked).

### Do the predictions hold in aggregate? (irreducible rays with C > 0: 7,864; M1 → C1)

| values counted | observed | predicted C·Σ1/ln f | obs/pred | z mean | z sd |
|---|---|---|---|---|---|
| all (f ≥ 2) | 209,337 | 212,855.8 | 0.9835 | −0.084 | 0.584 |
| f ≥ 10^3 | 191,999 | 192,781.0 | 0.9959 | −0.019 | 0.605 |
| f ≥ 10^4 | 165,761 | 166,215.8 | 0.9973 | −0.013 | 0.626 |

By size of C (all values / values ≥ 10^3): [0,1) 0.9833/0.9972, [1,2) 0.9829/0.9957, [2,3) 0.9824/0.9935, [3,4) 0.9925/1.0036, [4,10) 0.9856/0.9945 (M1 → C1).

* **Mean: yes.** The 1.6 % excess on all values comes from values below 1,000, where C/ln f overstates the chance of being prime. Above 1,000 the total is right to 0.4 %, with no trend in C.
* **Spread: no, sd ≈ 0.6 instead of 1.** A model where each value is prime independently, with probability min(C/ln f, 1), would give sd 0.923 (M1 → C1). So the counts are *sub-Poisson*. The likely reason (heuristic): the small primes hit a ray in an exactly periodic pattern, so their removals do not fluctuate like coin flips. This is the same effect as the Montgomery–Soundararajan variance for primes in short intervals. Use sqrt(pred) as an upper bound for the noise on the page, not as the noise itself.

### Brightest ray for each P (M1 → C1; "pts" = number of values on the ray)

| P | brightest by C: f(m) [j, s] | C | obs/pred/pts | most primes: [j, s] | C | obs/pred |
|---|---|---|---|---|---|---|
| 3 | (3m²+3m−2)/2 [3,−2] | 2.404 | 128/133.6/576 | same | | |
| 4 | 2m²−1 [2,−2] | 3.699 | 174/174.1/499 | same | | |
| 5 | (5m²+m−2)/2 [3,−2] | 2.246 | 98/97.4/447 | same | | |
| 6 | 3m²+3m+1 [6,0] (cuban) | 3.357 | 126/129.2/407 | same | | |
| 7 | (7m²+7m+2)/2 [7,0] | 2.300 | 74/81.9/377 | [2,−5] | 1.940 | 80/69.0 |
| 8 | 4m²+7 [4,6] | 3.942 | 129/131.4/353 | same | | |
| 9 | (9m²+9m+8)/2 [9,3] | 2.450 | 73/76.5/332 | [3,−5] | 2.404 | 83/75.0 |
| 10 | 5m²+5m+1 [10,0] | 4.427 | 130/131.3/315 | same | | |
| 11 | (11m²−9m+2)/2 [1,0] | 2.190 | 64/64.9/301 | same | | |
| 12 | 6m²+6m+1 [12,0] | 4.147 | 114/112.3/288 | same | | |
| 13 | (13m²−11m+2)/2 [1,0] | 2.573 | 65/70.4/277 | same | | |
| 14 | 7m²−3m+7 [4,6] | 3.990 | 106/100.8/267 | same | | |
| 15 | (15m²+15m+4)/2 [15,1] | 2.372 | 61/57.2/257 | same | | |
| 16 | 8m²−6m−1 [2,−2] | 4.718 | 109/110.0/249 | same | | |
| 17 | (17m²−9m+8)/2 [4,3] | 2.520 | 55/57.9/242 | [17,1] | 2.506 | 59/56.8 |
| 18 | 9m²+9m−1 [18,−2] | 4.837 | 110/106.6/235 | same | | |
| 19 | (19m²−11m+8)/2 [4,3] | 2.571 | 52/55.9/229 | same | | |
| 20 | 10m²+7 [10,6] | 3.769 | 79/79.1/223 | [4,0] | 3.656 | 80/78.1 |
| 21 | (21m²+15m+2)/2 [18,0] | 2.453 | 52/49.9/217 | same | | |
| 22 | 11m²+9m−1 [20,−2] | 3.935 | 82/78.2/212 | same | | |
| 23 | (23m²+19m−8)/2 [21,−5] | 2.387 | 51/46.6/208 | same | | |
| 24 | 12m²+12m−1 [24,−2] | 4.147 | 80/78.8/203 | same | | |
| 25 | (25m²−15m−8)/2 [5,−5] | 2.807 | 53/52.1/199 | same | | |
| 26 | 13m²−11m−1 [2,−2] | 4.452 | 81/81.1/195 | same | | |
| 27 | (27m²−21m−2)/2 [3,−2] | 2.706 | 47/52.1/192 | [15,−5] | 2.539 | 49/46.0 |
| 28 | 14m²−8m+1 [6,0] | 4.439 | 86/79.9/189 | same | | |
| 29 | (29m²+25m−8)/2 [27,−5] | 2.609 | 46/45.2/185 | [4,3] | 2.518 | 48/44.6 |
| 30 | 15m²+15m−1 [30,−2] | 3.786 | 62/64.3/182 | [12,4] | 3.490 | 65/59.8 |
| 31 | (31m²−13m+8)/2 [9,3] | 2.527 | 42/42.8/179 | [14,1] | 2.512 | 44/42.4 |
| 32 | 16m²+2m−1 [18,−2] | 4.718 | 82/78.1/176 | same | | |
| 33 | (33m²+33m+8)/2 [33,3] | 2.645 | 43/42.7/173 | same | | |
| 34 | 17m²−7m+7 [10,6] | 5.589 | 96/90.1/171 | same | | |
| 35 | (35m²−25m+4)/2 [5,1] | 2.419 | 41/39.0/169 | same | | |
| 36 | 18m²+18m+1 [36,0] | 4.535 | 70/70.2/166 | same | | |
| 37 | (37m²−33m−8)/2 [2,−5] | 2.861 | 46/43.5/163 | same | | |
| 38 | 19m²+19m−1 [38,−2] | 5.911 | 90/88.7/161 | same | | |
| 39 | (39m²+39m−2)/2 [39,−2] | 2.326 | 34/34.5/159 | [9,3] | 2.241 | 37/33.9 |
| 40 | 20m²−16m−5 [4,−6] | 4.492 | 69/65.6/157 | same | | |

### Top 15 rays over all P (M1 → C1)

By number of primes:

| # | P | f(m) [j,s] | pts | obs | pred | C |
|---|---|---|---|---|---|---|
| 1 | 4 | 2m²−1 [2,−2] | 499 | 174 | 174.1 | 3.699 |
| 2 | 4 | 2m²+2m−1 [4,−2] | 499 | 131 | 131.8 | 2.764 |
| 3 | 10 | 5m²+5m+1 [10,0] | 315 | 130 | 131.3 | 4.427 |
| 4 | 4 | 2m²+2m+1 [4,0] | 499 | 130 | 129.8 | 2.742 |
| 5 | 8 | 4m²+7 [4,6] | 353 | 129 | 131.4 | 3.942 |
| 6 | 3 | (3m²+3m−2)/2 [3,−2] | 576 | 128 | 133.6 | 2.404 |
| 7 | 6 | 3m²+3m+1 [6,0] | 407 | 126 | 129.2 | 3.357 |
| 8 | 6 | 3m²+m−1 [4,−2] | 408 | 126 | 126.0 | 3.224 |
| 9 | 8 | 4m²+2m−1 [6,−2] | 353 | 124 | 118.8 | 3.542 |
| 10 | 8 | 4m²+4m−1 [8,−2] | 353 | 120 | 123.4 | 3.699 |
| 11 | 8 | 4m²−2m−1 [2,−2] | 352 | 118 | 117.3 | 3.542 |
| 12 | 6 | 3m²−m−1 [2,−2] | 407 | 117 | 123.6 | 3.224 |
| 13 | 4 | 2m²+2m+7 [4,6] | 499 | 117 | 120.7 | 2.567 |
| 14 | 12 | 6m²+6m+1 [12,0] | 288 | 114 | 112.3 | 4.147 |
| 15 | 3 | (3m²+3m−4)/2 [3,−3] | 575 | 114 | 110.8 | 2.044 |

By C (density = obs/pts):

| # | P | f(m) [j,s] | C | obs/pts | density |
|---|---|---|---|---|---|
| 1 | 38 | 19m²+19m−1 [38,−2] | 5.911 | 90/161 | 0.559 |
| 2 | 34 | 17m²−7m+7 [10,6] | 5.589 | 96/171 | 0.561 |
| 3 | 34 | 17m²+7m+7 [24,6] | 5.589 | 87/171 | 0.509 |
| 4 | 34 | 17m²−15m−1 [2,−2] | 5.178 | 79/170 | 0.465 |
| 5 | 34 | 17m²+15m−1 [32,−2] | 5.178 | 77/171 | 0.450 |
| 6 | 18 | 9m²+9m−1 [18,−2] | 4.837 | 110/235 | 0.468 |
| 7–10 | 16, 16, 32, 32 | 8m²∓6m−1, 16m²∓2m−1 | 4.718 | 109, 109, 75, 82 of 249, 249, 176, 176 | 0.426–0.466 |
| 11 | 36 | 18m²+18m+1 [36,0] | 4.535 | 70/166 | 0.422 |
| 12–13 | 40 | 20m²±16m−5 [36,−6], [4,−6] | 4.492 | 55, 69 of 157 | 0.350, 0.439 |
| 14–15 | 26 | 13m²±11m−1 [24,−2], [2,−2] | 4.452 | 74, 81 of 195 | 0.379, 0.415 |

Mirror rays (j, s) and (P − j, s) are the same polynomial with m → −m, so they always have the same C.

### Answers to question 2

* **Brightest single line.**
  * By C and by share of primes: **P = 38**, 19m² + 19m − 1 (two places before the top corner): C = 5.911, 90 primes among 161 values.
  * A near tie is **P = 34**, 17m² − 7m + 7: C = 5.589, 96 primes among 171 values, the highest share at 0.561, with its mirror 17m² + 7m + 7.
  * By count: **P = 4**, 2m² − 1, with 174 primes among 499 values, because small P gives the longest rays (about √(2N/P) values).
  * Why 19m² + 19m − 1 is bright: its local factors are 2.000, 1.500, 1.250, 1.167, 1.100, 1.083 at p = 2, 3, 5, 7, 11, 13. So it is never even and never a multiple of 3, 5, 7, 11 or 13, and the product up to 13 is already 5.21 (F1). Here disc(2f) = 1748 = 4·19·23.
* **Famous polynomials: none, and none is possible.**
  * The Euler family needs disc = −163·k². For odd P, disc(2f) = (2j − P)² − 8P(1+s) is ≡ 1 (mod 8), while −163k² is ≡ 5 (mod 8) for odd k and even for even k.
  * For even P = 2q it would need q(1+s) = x² + x + 41 (always prime, ≥ 41) with q ≤ 20 and 1 + s ≤ 7, which is impossible (derived).
  * The 2m² − 199 family needs a discriminant with squarefree part 398. P = 4 rays have (j−2)² − 8(1+s) ≤ 44.
  * The computation agrees: 0 rays have squarefree discriminant part −163 or 398 (M1 → C1). Every ray has C ≤ 5.911, below C(m² + m + 41) = 6.656 and C(2m² − 199) = 7.296 (p < 10^4, M1 → C1).
  * The underlying reason: the constant term 1 + s is tiny (−5 … 7).
* **Recommendations:**
  * **Default P = 10.** Its corner ray 5m² + 5m + 1 (the centred decagonal numbers) has 130 primes among 315 values (41 %, C = 4.427). That is the most primes on one line for any P ≥ 5, and the rays are still long.
  * **P = 4** for the line with the most primes: 2m² − 1, 174 of 499. Two more P = 4 rays have 131 and 130.
  * **P = 6** for the story: the corner ray is exactly the cuban primes (126 of 407), with 3m² + m − 1 (126) and 3m² − m − 1 (117) beside the other corners.
  * **P = 38** (or 34) for the single brightest line (56 % primes).
  * Odd P look flat: no line beats C = 2.861.

---

## Layout 2: chebrace (to 10^7, and 5·10^5 for the page)

### Derivation (proved unless marked; identities checked C2)

1. The units mod P form a group G. R is the subgroup of squares. The squaring map has kernel of size r = #{b : b² ≡ 1}, so |R| = φ/r. Every square has exactly r square roots, every non-square has none. So c(a) = r on R and 0 on Nn, and Σ_{a∈R} c(a) = φ. (M2 asserts this for all P; C2 rebuilds R from the local criterion.)
2. Summing the brief's formula: π_R ≈ |R|·Li(x)/φ − ½·Li(√x) and π_N ≈ |Nn|·Li(x)/φ. Hence h = |R|π_N − |Nn|π_R ≈ ½·|Nn|·Li(√x). The Li(x) terms cancel exactly by the choice of weights.
3. Exact version. Let Π(x; a) = Σ_{p^k ≤ x, p^k ≡ a} 1/k (Riemann's prime-power count).
   * Because Nn is a union of cosets of R, p^k is a square exactly when k is even or p is a square. C2 finds 0 violations for all p^k ≤ 10^7 and all P.
   * Hence the identity, exact for all x:
     **h(x) = H_Π(x) + |Nn|·Σ_{k even} π*(x^{1/k})/k − Σ_{k odd ≥ 3} h(x^{1/k})/k**
     where H_Π = |R|·Π(x; Nn) − |Nn|·Π(x; R) and π* skips primes dividing P. C2 checks it exactly with Fractions: 0 failures for P ∈ {3, 4, 8, 15, 19, 24, 40} and x ∈ {10^3, 12345, 10^5, 999983}.
   * So h(x) = ½|Nn|·Π*(√x) + H_Π(x) + (tiny odd-power terms), and H_Π = −|R|·Σ_{χ real, χ ≠ χ₀} Π(x, χ).
   * Under GRH each Π(x, χ) is −Σ_ρ Li(x^ρ) plus lower-order terms, which average to 0 in ln x. Π*(y) ≈ Li(y) by the prime number theorem. This gives **E[h(x)] ≈ ½|Nn|·Li(√x)** (conditional on GRH for the "mean 0" part).
4. Normalised (Rubinstein–Sarnak): (ln x/√x)·h/|R| ≈ (r − 1) + Σ_χ Σ_γ x^{iγ}/(½ + iγ).
   * The bias r − 1 is the number of nonprincipal real characters.
   * Under GRH + LI the oscillating part has mean 0 and variance V(P) = Σ_χ V(χ*), with V(χ_d) = Σ_γ 1/(¼ + γ²) = ln(|d|/π) + ψ((1+a)/2) + 2L′(1, χ_d)/L(1, χ_d). The last formula comes from the Hadamard product (Davenport ch. 12), with a = 0 for d > 0 and a = 1 for d < 0.
   * Gaussian estimate: P(h < 0) ≈ Φ(−(r−1)/√V(P)).

### Question 1: table to 10^7 (M2 → C2)

"avg Li" is the average of h/(½|Nn|Li(√x)) over ln x in [10^3, 10^7]. M2 samples 10^5 grid points and C2 integrates exactly; they differ by at most 3.8·10^-5. "avg Π*" uses ½|Nn|Π*(√x) instead. "#p<0" counts the primes p ≤ 10^7 after which h < 0. "log dens" is the share of ln x in [10^3, 10^7] with h < 0.

| P | \|R\| | \|Nn\| | r | avg Li | avg Π* | #p<0 (share %) | log dens | first h<0 | first >1000 | last ≤10^7 | h(10^7) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 3 | 1 | 1 | 2 | 0.959 | 0.996 | 0 | 0 | — | — | — | 190 |
| 4 | 1 | 1 | 2 | 0.952 | 0.998 | 239 (0.036) | 0.0006 | 26861 | 26861 | 633797 | 218 |
| 5 | 2 | 2 | 2 | 1.031 | 1.063 | 0 | 0 | — | — | — | 612 |
| 6 | 1 | 1 | 2 | 0.914 | 1.008 | 0 | 0 | — | — | — | 189 |
| 7 | 3 | 3 | 2 | 0.974 | 0.997 | 13 (0.002) | 0.0004 | 2 | 48673 | 49757 | 534 |
| 8 | 1 | 3 | 4 | 0.978 | 1.027 | 0 | 0 | — | — | — | 674 |
| 9 | 3 | 3 | 2 | 0.959 | 0.996 | 0 | 0 | — | — | — | 570 |
| 10 | 2 | 2 | 2 | 0.986 | 1.083 | 0 | 0 | — | — | — | 610 |
| 11 | 5 | 5 | 2 | 0.975 | 0.992 | 19882 (2.99) | 0.029 | 5 | 7789 | 2532319 | 1700 |
| 12 | 1 | 3 | 4 | 0.935 | 1.034 | 0 | 0 | — | — | — | 533 |
| 13 | 6 | 6 | 2 | 1.006 | 1.025 | 18319 (2.76) | 0.015 | 2083 | 2083 | 6696259 | 1572 |
| 14 | 3 | 3 | 2 | 1.019 | 1.113 | 3 (0.0005) | 0.00005 | 48731 | 48731 | 49739 | 537 |
| 15 | 2 | 6 | 4 | 0.929 | 1.004 | 0 | 0 | — | — | — | 1050 |
| 16 | 2 | 6 | 4 | 0.978 | 1.027 | 0 | 0 | — | — | — | 1348 |
| 17 | 8 | 8 | 2 | 1.065 | 1.084 | 1052 (0.16) | 0.0013 | 2 | 17863 | 2374529 | 592 |
| 18 | 3 | 3 | 2 | 0.914 | 1.008 | 0 | 0 | — | — | — | 567 |
| 19 | 9 | 9 | 2 | 0.868 | 0.882 | 216584 (32.6) | 0.148 | 11 | 1061 | 8365061 | 1872 |
| 20 | 2 | 6 | 4 | 0.943 | 1.033 | 0 | 0 | — | — | — | 1634 |
| 21 | 3 | 9 | 4 | 0.952 | 1.024 | 0 | 0 | — | — | — | 1611 |
| 22 | 5 | 5 | 2 | 0.930 | 1.000 | 21203 (3.19) | 0.038 | 3 | 7649 | 2533507 | 1695 |
| 23 | 11 | 11 | 2 | 0.970 | 0.984 | 38192 (5.75) | 0.069 | 2 | 3209 | 2709677 | 3894 |
| 24 | 1 | 7 | 8 | 0.948 | 1.052 | 0 | 0 | — | — | — | 1481 |
| 25 | 10 | 10 | 2 | 1.031 | 1.063 | 0 | 0 | — | — | — | 3060 |
| 26 | 6 | 6 | 2 | 0.961 | 1.040 | 19126 (2.88) | 0.018 | 3 | 2063 | 6696281 | 1566 |
| 27 | 9 | 9 | 2 | 0.959 | 0.996 | 0 | 0 | — | — | — | 1710 |
| 28 | 3 | 9 | 4 | 0.977 | 1.066 | 0 | 0 | — | — | — | 1599 |
| 29 | 14 | 14 | 2 | 1.116 | 1.135 | 20596 (3.10) | 0.084 | 719 | 7621 | 9182039 | 952 |
| 30 | 2 | 6 | 4 | 0.914 | 1.058 | 0 | 0 | — | — | — | 1048 |
| 31 | 15 | 15 | 2 | 0.952 | 0.965 | 78180 (11.8) | 0.106 | 2 | 1861 | 4505687 | 6000 |
| 32 | 4 | 12 | 4 | 0.978 | 1.027 | 0 | 0 | — | — | — | 2696 |
| 33 | 5 | 15 | 4 | 0.980 | 1.048 | 0 | 0 | — | — | — | 3845 |
| 34 | 8 | 8 | 2 | 1.110 | 1.210 | 862 (0.13) | 0.0010 | 435949 | 435949 | 2374523 | 600 |
| 35 | 6 | 18 | 4 | 0.946 | 1.008 | 2 (0.0003) | 0 | 11 | — | 751 | 5766 |
| 36 | 3 | 9 | 4 | 0.935 | 1.034 | 0 | 0 | — | — | — | 1599 |
| 37 | 18 | 18 | 2 | 1.049 | 1.061 | 48337 (7.27) | 0.068 | 11 | 4831 | 3131321 | 8820 |
| 38 | 9 | 9 | 2 | 0.823 | 0.883 | 220496 (33.2) | 0.171 | 7 | 1051 | 8365099 | 1863 |
| 39 | 6 | 18 | 4 | 0.987 | 1.053 | 561 (0.084) | 0.0009 | 873913 | 873913 | 1006441 | 3798 |
| 40 | 2 | 14 | 8 | 0.948 | 1.033 | 0 | 0 | — | — | — | 4258 |

**Is the average ≈ 1? Yes, to leading order.**

* Over the 38 moduli the average of avg Li is 0.970 (median 0.960, range 0.823 at P = 38 to 1.116 at P = 29) (M2; per-P values rechecked by C2).
* With Π*(√x) the average is 1.029 (median 1.027), and Π*(√x)/Li(√x) itself averages 0.944 (M2).
* The scatter between moduli comes from the zeros (P = 19, 38, 29, 34) and from two lower-order effects:
  * the primes dividing P drop out;
  * an even real character adds about |R|·ln ln x, from the trivial zero of L(s, χ) at s = 0 (heuristic).
* Data agree with that heuristic. The 14 moduli without an even real character average 0.940 (Li) / 0.987 (Π*). The 24 with one average 0.987 / 1.054 (C2, from per-P values that M2 and C2 both computed).

**Small-x artefacts.**

* For odd P, the races mod P and mod 2P are the same except for the prime 2. h_P − h_2P is constant for every prime ≤ 10^7: +|R| when 2 is a non-square mod P, −|Nn| when it is a square (C2, P = 3 … 19 odd).
* So "first h < 0 = 2" for P = 7, 17, 23, 31 only says that 2 is a square mod P. Use the "first > 1000" column instead.

### Same at N = 5·10^5 (M2 → C2)

| P | avg Li | avg Π* | #p<0 | first h<0 | last ≤ 5·10^5 |
|---|---|---|---|---|---|
| 3, 5, 6, 8, 9, 10, 12, 15, 16, 18, 20, 21, 24, 25, 27, 28, 30, 32, 33, 36, 39, 40 | see below | | 0 | — | — |
| 4 | 0.937 | 1.001 | 1 | 26861 | 26861 |
| 7 | 0.959 | 0.990 | 13 | 2 (48673 after 1000) | 49757 |
| 11 | 0.922 | 0.944 | 1364 | 5 (7789) | 247259 |
| 13 | 1.052 | 1.078 | 101 | 2083 | 150221 |
| 14 | 1.022 | 1.152 | 3 | 48731 | 49739 |
| 17 | 1.061 | 1.087 | 187 | 2 (17863) | 440431 |
| 19 | 0.910 | 0.928 | 3436 | 11 (1061) | 208223 |
| 22 | 0.859 | 0.952 | 1668 | 3 (7649) | 247279 |
| 23 | 0.905 | 0.923 | 5271 | 2 (3209) | 411937 |
| 26 | 0.990 | 1.099 | 152 | 3 (2063) | 150223 |
| 29 | 1.081 | 1.106 | 11115 | 719 (7621) | 459607 |
| 31 | 0.927 | 0.945 | 2606 | 2 (1861) | 149899 |
| 34 | 1.124 | 1.262 | 152 | 435949 | 440399 |
| 35 | 0.976 | 1.062 | 2 | 11 (none after 1000) | 751 |
| 37 | 1.070 | 1.086 | 1818 | 11 (4831) | 195887 |
| 38 | 0.847 | 0.928 | 3602 | 7 (1051) | 208231 |

avg Li for the never-negative moduli at 5·10^5: 3: 0.945, 5: 1.052, 6: 0.882, 8: 0.977, 9: 0.945, 10: 0.989, 12: 0.908, 15: 0.917, 16: 0.977, 18: 0.882, 20: 0.930, 21: 0.932, 24: 0.940, 25: 1.052, 27: 0.945, 28: 0.973, 30: 0.897, 32: 0.977, 33: 0.960, 36: 0.908, 39: 1.000, 40: 0.899. The mean over all 38 moduli is 0.961 (M2 → C2). **22 of the 38 moduli never go below 0 up to 5·10^5; 21 never do up to 10^7** (P = 39 first drops at 873,913) (M2 → C2).

### Question 2: which P are lopsided? (T2; V(χ) rechecked Z2)

Characters and zero variances for the real primitive characters (Z2 agrees to 9·10^-11 by the Dirichlet series):

| d | −3 | −4 | 5 | −7 | −8 | 8 | −11 | 12 | 13 | −15 | 17 | −19 | −20 | 21 | −23 | −24 | 24 | 28 | 29 | −31 | 33 | −35 | 37 | −39 | −40 | 40 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| V | .113 | .156 | .157 | .255 | .316 | .235 | .508 | .330 | .397 | .459 | .387 | **1.053** | .639 | .614 | .569 | .797 | .552 | .573 | .879 | .810 | .566 | 1.110 | .814 | .823 | 1.404 | .697 |

Summing the actual zeros up to height 60, plus the density tail, gives V = 0.1129 / 0.1554 / 0.1564 / 1.0525 for d = −3 / −4 / 5 / −19, against 0.1132 / 0.1556 / 0.1566 / 1.0527 (Z2). The first zeros are at 8.0397, 6.0209, 6.6485 and **1.5161** (Z2).

z = (r − 1)/√V(P), ranked: 24: 4.43; 12, 36: 3.88; 40: 3.69; 8, 16, 32: 3.57; 15, 30: 3.51; 20: 3.08; 21: 3.03; 28: 3.02; 3, 6, 9, 18, 27: 2.97; 33: 2.75; 39: 2.60; 4: 2.54; 5, 10, 25: 2.53; 35: 2.43; 7, 14: 1.98; 17, 34: 1.61; 13, 26: 1.59; 11, 22: 1.40; 23: 1.33; 31, 37: 1.11; 29: 1.07; 19, 38: 0.97 (T2).

* **The divisors of 24 do not stand out as a group.**
  * P = 24 is the most lopsided (z = 4.43, Gaussian chance of h < 0 ≈ 5·10^-6) and 12 is second (3.88). Both never go negative to 10^7.
  * But 3 and 6 (z = 2.97) and 4 (2.54) are mid-table, and P = 4 does go negative (239 primes). Nine moduli that do not divide 24 score higher than 3 and 6 (36, 40, 16, 32, 15, 30, 20, 21, 28). Sixteen moduli that do not divide 24 never go negative to 10^7: 5, 9, 10, 15, 16, 18, 20, 21, 25, 27, 28, 30, 32, 33, 36, 40.
  * R = {1} only maximises r for a given φ. What matters is r − 1 (the number of real characters) against V(P), and V grows with the conductors.
* The ranking explains the data. The Spearman correlation between z and the observed log density of h < 0 is 0.868 (p = 1.7·10^-12). Every never-negative modulus has z ≥ 2.53. The two that go negative closest to that line are P = 39 (2.60) and P = 4 (2.54), so near z ≈ 2.5 the outcome to 10^7 is luck (T2 with M2 data).
* The Gaussian estimate is rough. It gives 0.0056 (P = 4) and 0.0015 (P = 3), while Rubinstein–Sarnak's exact densities are 0.0041 and 0.0010 (literature, quoted in FACTS.md 1602/1606; not recomputed).
* **The least lopsided are P = 19 and 38.** h < 0 after 32.6 % and 33.2 % of the primes ≤ 10^7. The cause is that L(s, χ₋₁₉) has a zero at height γ = 1.5161. This one pair contributes 0.785 of V = 1.053 (Z2). Its oscillation has amplitude 2/|ρ| ≈ 1.25 against a bias of 1, with period 2π/γ ≈ 4.1 in ln x (derived from Z2's zero).

### Question 3 (M2 → C2)

* **P = 4:** the first prime with π(p; 4, 1) > π(p; 4, 3) is **26,861**, which confirms Leech 1957. Up to 10^7 the squares lead after 239 primes, the last of them 633,797. Up to 5·10^5 they lead only after 26,861 itself.
* **P = 3:** **no sign change up to 10^7** (h ≥ 0 at every prime). π(10^7; 3, 2) − π(10^7; 3, 1) = 190. The first lead, at 608,981,813,029 (Bays–Hudson 1978), is far out of range and was not checked here.

---

## Text for the page's explanatory notes (every number above)

**Polyspiral.**

* Each ring is a P-sided polygon, so the numbers a fixed step from a corner lie on a straight line and are the values of a quadratic, f(m) = (P/2)m² + (j − P/2)m + 1 + s. The corner line itself is the centred P-gonal numbers, and for P = 6 these are the cuban primes (m+1)³ − m³.
* The Bateman–Horn conjecture says how bright each line should be: a constant C, which measures how well f(m) avoids small prime factors, times Σ 1/ln f(m). On the 7,864 lines where it applies it predicts 192,781 primes above 1,000, and 191,999 are there.
* The brightest line, 19m² + 19m − 1 for P = 38, is never even and never a multiple of 3, 5, 7, 11 or 13 (C ≈ 5.9): 90 of its 161 numbers are prime.
* The line one step before each corner is always dark, because there f(m) = m × (something) factors.
* For odd P exactly two in every four numbers on every line are even, so no line gets brighter than C ≈ 2.9.

**Chebyshev race.**

* For each modulus P the walk steps up at every prime that is not a square mod P and down at every prime that is. The step sizes are chosen so that it would hover around zero if primes were spread evenly over the classes.
* Instead it drifts upward, by about ½·|Nn|·Li(√x). Counting prime powers spreads evenly over the classes, but every prime squared lands in a square class, so the primes themselves fall behind there. This is Chebyshev's bias, explained by Rubinstein and Sarnak in 1994.
* Up to 500,000 the walk never drops below zero for 22 of the 38 moduli. For P = 4 the only prime after which the squares lead is 26,861.
* How one-sided a race is depends on how many real characters P has, against the noise from the zeros of their L-functions. P = 24 is the most one-sided. P = 19 and 38 are below zero after a third of the primes up to 10^7, because L(s, χ₋₁₉) has an unusually low zero at height 1.52.

## Not verified here

* Bateman–Horn, the Rubinstein–Sarnak mean-zero statement and the V(χ) = Σ 1/(¼ + γ²) identification all rest on conjectures (BH; GRH + LI). They are tested only empirically, in the ranges above.
* The exact Rubinstein–Sarnak densities: only the Gaussian approximation was computed.
* The mod-3 sign change at 608,981,813,029.
* The ln ln x term from even characters and the sub-Poisson explanation are heuristics.
* C with the product to p < 10^6 (0.49 % at most) comes from M1 alone.
