# quadform and eisenstein layouts: the numbers behind the pictures

N = 500,000, π(N) = 41,538 (numpy sieve in `qf_count.py`; `sympy.primerange` in `qf_check.py`, which asserts the same count).
Every number below is printed by a script in this directory. Each result was rechecked by a second script that uses a different method; the pairs are listed here and named again at each table.
Run with `(ulimit -v 1000000; timeout 900 nice /home/bas/praime/.venv/bin/python <script>)`; the order is `qf_count.py` → `qf_check.py` → `make_table.py`, then the rest in any order. `make_table.py` and `make_tables2.py` only format the JSON files.

| script (method A) | what | recheck (method B) | result of the recheck |
|---|---|---|---|
| `qf_count.py` | Q1–2: direct search over (a, b) with the page's smallest-b rule; h(−4P) by counting reduced forms; the residue-class test mod 4P | `qf_check.py`: Cornacchia for each prime (root from `sympy.sqrt_mod`), h from PARI `qfbclassno`, class group from PARI `quadclassunit`, idoneal test by genus count (Cox §3) | 0 mismatches over 100 P (`qf_check.out`) |
| `qf_extra.py` | Q3: power residues, Tonelli–Shanks, hand-written root counting mod p, compared with direct search | `qf_extra_check.py`: representation by PARI `qfbsolve`; criterion from the class polynomial H₋₄ₚ (PARI `polclass`, Cox Thm 9.2); root counts by `sympy.polys.galoistools`; same-field test by PARI `nfisisom` | all agree (`qf_extra_check.out`) |
| `eis.py` | Q4: direct search, a ≥ b ≥ 0, smallest b | `eis_check.py`: PARI `qfbsolve` for x² + xy + y², orbit under the 12 automorphisms | count, set and every drawn position agree |
| `lines.py` | Q5: lines through the lattice points of the drawn region, Bateman–Horn with root counts by brute force mod p | `lines_check.py`: prime positions from PARI `qfbsolve`, a plain loop along each line, root counts from the Legendre-symbol formula | 0 mismatches on 64 lines, top-5 rankings agree |
| `picture_stats.py` | Q6: density of the dots, background constant, row/column/diagonal contrasts for every P | `picture_check.py`: lattice points by a closed sum, counts from Cornacchia, h from PARI, Bateman–Horn constants by brute force | 0 mismatches over 100 P |
| `subsets.py` | Q6: how the pictures nest (P = 16, 64 in P = 1; P = 12, 27 in P = 3); P = 58 | `subsets_check.py`: PARI `qfbsolve`, sympy | all true |
| `bias.py` | why every count falls below 1/(2h) | `bias_check.py` (from the Cornacchia counts and the genus count) | 0 mismatches |

## Conventions

* **Drawn region R**: for P = 1, 0 ≤ b ≤ a and a² + b² ≤ N, an octant, because the smallest b always has b ≤ a. For P ≥ 2, a, b ≥ 0 and a² + Pb² ≤ N. For eisenstein, a ≥ b ≥ 0 and a² + ab + b² ≤ N, a 30° wedge of the hexagonal lattice.
* **Primes and positions**: a prime that is drawn has exactly one lattice point in R, so its position does not depend on the smallest-b rule. Only composites have several points. `lines.py` asserts this, and `eis_check.py` and `lines_check.py` confirm the positions independently.
* **Primes dividing 4P**: these are 2 and the prime divisors of P, and they are excluded from the density comparison and from the residue-class test. 2 is drawn only for P = 1 and 2 (1 + 1 and 0 + 2). An odd prime q dividing P is drawn only when q = P (a = 0, b = 1), which happens for the 24 odd primes P ≤ 100. The table column "share among p not dividing 4P" divides the drawn primes p ∤ 4P by the number of primes p ∤ 4P up to N.

## 1–2. Primes drawn, class numbers, congruences (all P = 1..100)

Source: `qf_count.py` → `qf_count.json`. Rechecked by `qf_check.py`, with 100 of 100 rows agreeing (`make_table.py`).

{{SUMMARY12}}

{{table.md}}

### Why every count is a little below 1/(2h)

In all 100 cases the fraction (excluding p | 4P) is below 1/(2h). This is a Chebyshev-type bias. The primes of the principal form are exactly the primes whose Frobenius in Gal(L/Q) is the identity, where L is the ring class field and G is generalized dihedral of order 2h. Counting prime powers with weight 1/k gives no bias. Removing the prime squares p² ≤ x takes away the most from the identity class, because every element of order ≤ 2 squares to it. The mean deficit is then about ½·(I − 1)/(2h)·π(√N), where I = #Cl[2] + h = #{g : g² = 1}. This is the Rubinstein–Sarnak heuristic carried over to Chebotarev sets (Ng 2000; Bellaïche 2016). It is a heuristic, and the value at a single x fluctuates.

{{BIAS}}

## 3. The classical criteria for non-idoneal P (every prime ≤ 500,000)

Method A is `qf_extra.py`. Method B is `qf_extra_check.py`, whose representations come from PARI `qfbsolve` and are independent of method A.

{{q3_table.md}}

Notes:
* **P = 11 and 23 need "p ≠ P"**. The only disagreements anywhere are p = 11 for P = 11 and p = 23 for P = 23. 11 = 0² + 11·1² and 23 = 0² + 23·1² are drawn (a = 0), but p divides the discriminant of the cubic (−44 and −23), so the cubic has a double root there instead of three roots. The verified statements are:
  * p = x² + 23y² ⟺ p ≠ 23 and x³ − x − 1 has three roots mod p.
  * p = x² + 11y² ⟺ p ≠ 11 and x³ − x² + x + 1 has three roots mod p. An equivalent form is (−11/p) = 1 and the cubic has a root, which also agrees on all 41,536 primes p ≠ 2, 11.
* **Why these cubics**: x³ − x² + x + 1 has discriminant −44 and is exactly PARI's `polredabs` of the class polynomial H₋₄₄ (the minimal polynomial of j(√−11)). For P = 23, `polredabs(H₋₉₂)` = x³ − x² + 1, which `nfisisom` shows defines the same field as x³ − x − 1 (h(−92) = h(−23) = 3, so the ring class field of Z[√−23] is the Hilbert class field of Q(√−23)). In the same way, `polredabs(H₋₁₀₈)` = x³ − 2, `polredabs(H₋₂₅₆)` = x⁴ − 2, and `polredabs(H₋₅₆)` defines the same field as (x² + 1)² − 8. This is how Gauss's cubic and biquadratic rules for 2 and Cox's x² + 14y² rule fit Cox Thm 9.2 (p = x² + ny² ⟺ (−n/p) = 1 and the class polynomial has a root mod p, for p ∤ n·disc).
* The rules for P = 27, 64 and 14 hold for every prime, including 2, 3 and 7, where both sides are false.
* **Chebotarev for x³ − x − 1** (p ≠ 23): {{CHEB}}. The densities 1/3, 1/2, 1/6 are those of S₃.

## 4. eisenstein

Source: `eis.py`. Rechecked by `eis_check.py`.

{{EIS}}

## 5. Bright straight lines (P = 1, 2, 3 and eisenstein)

Source: `lines.py` → `lines.json`. Every line below was rechecked by `lines_check.py` (column "rechecked").

**Method.** A line is a set of lattice points of R with db·a − da·b = key, for a primitive step (da, db). Along the line, n = f(t) is a quadratic. The Bateman–Horn prediction is C(f)·Σ 1/ln f(t), summed over the line's lattice points in R with f ≥ 2. Here C(f) = ∏_{p ≤ 10,000} (1 − ω(p)/p)/(1 − 1/p), and ω(p) is the number of roots of f mod p.

**Contrast.** What the eye sees is the density of primes on a line compared with the area around it. A lattice point of R with norm n is prime with probability about B/ln n, where B = (share of primes drawn)/(lattice points of R per unit of n):
* P = 1: B = (1/2)/(π/8) = 4/π.
* P ≥ 2: B = (1/(2h))·4√P/π.
* eisenstein: B = (1/2)/(π/(6√3)) = 3√3/π.

The empirical B (primes in R divided by Σ 1/ln n over R) agrees to within 1%. A line with C/B = 3 is lit three times as densely as its surroundings.

{{lines_table.md}}

**What a viewer sees, and why**:
* **Parity sets the stripes.** For P odd (1, 3), n ≡ a + b (mod 2), so only the points with a + b odd can be prime, a checkerboard. Each 45° diagonal a − b = c has a single parity. Odd diagonals carry every prime and even ones none, which gives the striped octant of P = 1. For P even, n ≡ a (mod 2), so every even column is empty and P = 2 shows vertical stripes. For eisenstein, n is even only when a and b are both even, and n ≡ (a − b)² (mod 3). The lines a − b ≡ 0 (mod 3) are therefore dark, and the lines a − b = c with 3 ∤ c gain a factor 3/2, which is why a − b = 7 and 91 are bright.
* **Offsets divisible by split primes make a line brighter.** Take P = 1 and the line a − b = c, where n = b² + (b + c)². If q ≡ 1 (mod 4) divides c, the line meets the dark multiples of q once per period instead of twice, a factor (q − 1)/(q − 2). If q ≡ 3 (mod 4) divides c, it meets them once instead of never, a factor (q − 1)/q. The brightest P = 1 lines therefore have offsets 5, 25, 65 = 5·13 and 85 = 5·17. The same holds for P = 2 (columns a = 9, 33 = 3·11, 51 = 3·17, 57 = 3·19, 99, where 3, 11, 17 and 19 split), P = 3 (a ± b = 91, 273, 637: 7 and 13 split) and eisenstein (b = 7, 21, 49, 13 and a − b = 91).
* **Landau's line does not stand out.** b = 1 in P = 1 (a² + 1) holds 89 primes against a prediction of 93.6, but its C/B is only 1.08, so it is hardly brighter than its surroundings. The row b = 1 of P = 2 (a² + 2) has C/B = 0.79 and is darker than the background.
* **Steeper lines (step 3) reach C/B ≈ 3.3–4.0**, for example P = 1 a − 3b = ±221 and P = 3 a ± 3b = ±91, 637. Their points are √10 or more apart, so they read as dotted rather than solid lines.

## 6. Which P to show

Per-P picture numbers come from `picture_stats.py` → `picture_stats.json` and are rechecked by `picture_check.py`. The nesting facts come from `subsets.py`, rechecked by `subsets_check.py`.

{{PIC}}

### Recommendation

**Default P = 1. I agree with your choice.**
* It ties for the most dots: 20,732, that is 2 and every p ≡ 1 (mod 4), Fermat's theorem.
* It has the cleanest structure: the checkerboard turns into 45° diagonals. Odd diagonals carry every prime, and the brightest, a − b = 5, holds 169 primes (Bateman–Horn predicts 170.2). a − b = 65 is lit 3.13 times as densely as its surroundings.
* The octant (b ≤ a) is a thin wedge. P = 4 shows the same 20,731 odd primes over a full quadrant, with vertical stripes because a must be odd. That is worth one sentence in the note.

**Values for the note**:
1. **P = 5**: the first picture in which a quadratic-residue condition is not enough. (−5/p) = 1 for p ≡ 1, 3, 7, 9 (mod 20), but only p ≡ 1, 9 (mod 20) are drawn: 10,345 dots, a fraction of 0.2490 against 1/(2h) = 1/4. The others (p ≡ 3, 7 mod 20) belong to 2a² + 2ab + 3b²; every prime with (−5/p) = 1 lies in exactly one of the two forms (`subsets.py`, `subsets_check.py`). A congruence still decides it, because 5 is idoneal.
2. **P = 27**: the P = 3 picture thinned to every third row. Its 6,898 dots are exactly the P = 3 primes with 3 | b (33.27% of them). They are also exactly the primes p ≡ 1 (mod 3) for which 2 is a cubic residue (Gauss), with 0 exceptions among 41,538 primes. 18 of the 36 classes mod 108 are mixed, so no congruence can draw this picture.
3. **P = 64**: Gauss's biquadratic rule. Its 5,139 dots are the P = 1 primes whose even coordinate is divisible by 8 (24.8% of p ≡ 1 mod 4; in the P = 1 octant, 2,549 lie on rows b ≡ 0 and 2,590 on columns a ≡ 0 mod 8). 32 of the 128 classes mod 256 are mixed.
4. **P = 58**: the densest picture of all 100. It has 10,380 dots on 51,952 lattice points, so 20.0% of the points are lit, against 10.5% for P = 1 and 15.8% for the runner-up P = 37. The reason is that h(−232) = 2 is unusually small for P = 58, which gives the largest background constant: B = 2√58/(2π) = 2.42. P = 58 is idoneal (28 classes mod 232). The other form class, 2a² + 29b², contains the prime-generating polynomial 2x² + 29, which is prime for x = 0..28.

Alternatives:
* **P = 23**: the cleanest cubic rule, x³ − x − 1 with three roots, from the Hilbert class field of Q(√−23).
* **P = 11**: the smallest non-idoneal P (h = 3, 10 of 20 classes mixed), with rule x³ − x² + x + 1.

The most line-like pictures are P = 2 (vertical stripes, columns a = 33 and 99 lit 3.52 times as densely as the background) and P = 14 (column a = 15 lit 4.38 times as densely, the largest contrast among the rows, columns and diagonals with index ≤ 20 of any P).

## Not verified / caveats

* **Primality and ranges.** All primality tests are below 500,000, so they are exact (sieve, sympy, PARI). Every statement in Q1–Q4 is checked only up to 500,000. The equivalences themselves are theorems (Gauss, Euler, Cox Thm 9.2). Our runs confirm them on this range; they do not prove them.
* **The bias explanation in §1–2 is a heuristic.** The observed deficits (4 to 101 per P) scatter widely around the predicted mean (31.5 to 55.1), with correlation 0.36. Only the sign (100 of 100) and the total (4,693 against 4,126) fit well.
* **Choices in the line scan.** "Brightest" lines are taken among directions with steps up to 3 and lines with at least 100 lattice points (for the contrast lists). The per-P contrast maxima in §6 cover only rows, columns and diagonals with index ≤ 20.
* **Bateman–Horn.** The constants use primes up to 10,000. The product left out is 1 + O(1/(10⁴ ln 10⁴)) per line.
* **Citations from memory.** The literature references (Cox's section and theorem numbers, Ng 2000, Bellaïche 2016, Rubinstein–Sarnak) were not checked against the texts. The mathematics they are cited for was checked numerically as stated.
* **What "striking" means here.** It is judged only by the numbers above (dot density, line contrast, nesting, whether a congruence decides). I did not render the pictures.
