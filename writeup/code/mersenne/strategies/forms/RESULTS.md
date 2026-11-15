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

* **Class numbers.** h(−4P) from counting reduced forms equals PARI `qfbclassno` for all 100 P (`qf_check.py`). h = 1 for P = 1, 2, 3, 4, 7; the largest is h = 12 (P = 89).
* **Densities.** Excluding p | 4P, every fraction is within 2.81% (relative) of 1/(2h). The largest gap is at P = 41, the smallest (0.039%) at P = 10. All 100 fractions are *below* 1/(2h); see the bias note below.
* **Congruence test.** Every coprime class mod 4P holds at least 198 primes ≤ N, and no class is empty, so the test has teeth. The represented primes are a union of classes mod 4P for exactly 36 values of P. These are exactly Euler's idoneal numbers ≤ 100 from your list: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 18, 21, 22, 24, 25, 28, 30, 33, 37, 40, 42, 45, 48, 57, 58, 60, 70, 72, 78, 85, 88, 93. For every non-idoneal P, between 12.5% and 50.0% of the classes are mixed.
* **Four idoneal tests agree for all 100 P**: (i) the congruence test; (ii) h(−4P) equals the number of genera 2^(μ−1) (genus count as in Cox §3, `qf_check.py`); (iii) the class group has exponent ≤ 2 (PARI `quadclassunit`); (iv) every reduced form is ambiguous (b = 0, b = a or a = c; `qf_count.py`). This is genus theory: a congruence mod 4P fixes only the genus, so it decides exactly when each genus is a single class (D. Cox, *Primes of the form x² + ny²*, §3, genus theory and Euler's convenient numbers; section numbers in this file are cited from memory).
* **Rules for some idoneal P** (classes mod 4P, from `qf_count.json`): P = 1: p ≡ 1 (mod 4); P = 2: p ≡ 1, 3 (mod 8); P = 3: p ≡ 1, 7 (mod 12); P = 5: p ≡ 1, 9 (mod 20); P = 7: p ≡ 1, 9, 11, 15, 23, 25 (mod 28); P = 58: 28 classes mod 232.

| P | h(-4P) | primes drawn | / 41,538 | share among p not dividing 4P | 1/(2h) | congruence-decided | idoneal |
|---:|---:|---:|---:|---:|---:|:---:|:---:|
| 1 | 1 | 20,732 | 0.4991 | 0.4991 | 0.5000 | yes | yes |
| 2 | 1 | 20,753 | 0.4996 | 0.4996 | 0.5000 | yes | yes |
| 3 | 1 | 20,734 | 0.4992 | 0.4992 | 0.5000 | yes | yes |
| 4 | 1 | 20,731 | 0.4991 | 0.4991 | 0.5000 | yes | yes |
| 5 | 2 | 10,345 | 0.2490 | 0.2490 | 0.2500 | yes | yes |
| 6 | 2 | 10,321 | 0.2485 | 0.2485 | 0.2500 | yes | yes |
| 7 | 1 | 20,706 | 0.4985 | 0.4985 | 0.5000 | yes | yes |
| 8 | 2 | 10,334 | 0.2488 | 0.2488 | 0.2500 | yes | yes |
| 9 | 2 | 10,329 | 0.2487 | 0.2487 | 0.2500 | yes | yes |
| 10 | 2 | 10,380 | 0.2499 | 0.2499 | 0.2500 | yes | yes |
| 11 | 3 | 6,858 | 0.1651 | 0.1651 | 0.1667 | no | no |
| 12 | 2 | 10,329 | 0.2487 | 0.2487 | 0.2500 | yes | yes |
| 13 | 2 | 10,307 | 0.2481 | 0.2481 | 0.2500 | yes | yes |
| 14 | 4 | 5,151 | 0.1240 | 0.1240 | 0.1250 | no | no |
| 15 | 2 | 10,327 | 0.2486 | 0.2486 | 0.2500 | yes | yes |
| 16 | 2 | 10,334 | 0.2488 | 0.2488 | 0.2500 | yes | yes |
| 17 | 4 | 5,151 | 0.1240 | 0.1240 | 0.1250 | no | no |
| 18 | 2 | 10,357 | 0.2493 | 0.2493 | 0.2500 | yes | yes |
| 19 | 3 | 6,881 | 0.1657 | 0.1656 | 0.1667 | no | no |
| 20 | 4 | 5,167 | 0.1244 | 0.1244 | 0.1250 | no | no |
| 21 | 4 | 5,135 | 0.1236 | 0.1236 | 0.1250 | yes | yes |
| 22 | 2 | 10,342 | 0.2490 | 0.2490 | 0.2500 | yes | yes |
| 23 | 3 | 6,911 | 0.1664 | 0.1664 | 0.1667 | no | no |
| 24 | 4 | 5,137 | 0.1237 | 0.1237 | 0.1250 | yes | yes |
| 25 | 2 | 10,344 | 0.2490 | 0.2490 | 0.2500 | yes | yes |
| 26 | 6 | 3,427 | 0.0825 | 0.0825 | 0.0833 | no | no |
| 27 | 3 | 6,898 | 0.1661 | 0.1661 | 0.1667 | no | no |
| 28 | 2 | 10,320 | 0.2484 | 0.2485 | 0.2500 | yes | yes |
| 29 | 6 | 3,423 | 0.0824 | 0.0824 | 0.0833 | no | no |
| 30 | 4 | 5,094 | 0.1226 | 0.1226 | 0.1250 | yes | yes |
| 31 | 3 | 6,888 | 0.1658 | 0.1658 | 0.1667 | no | no |
| 32 | 4 | 5,142 | 0.1238 | 0.1238 | 0.1250 | no | no |
| 33 | 4 | 5,104 | 0.1229 | 0.1229 | 0.1250 | yes | yes |
| 34 | 4 | 5,140 | 0.1237 | 0.1237 | 0.1250 | no | no |
| 35 | 6 | 3,422 | 0.0824 | 0.0824 | 0.0833 | no | no |
| 36 | 4 | 5,142 | 0.1238 | 0.1238 | 0.1250 | no | no |
| 37 | 2 | 10,284 | 0.2476 | 0.2476 | 0.2500 | yes | yes |
| 38 | 6 | 3,428 | 0.0825 | 0.0825 | 0.0833 | no | no |
| 39 | 4 | 5,154 | 0.1241 | 0.1241 | 0.1250 | no | no |
| 40 | 4 | 5,146 | 0.1239 | 0.1239 | 0.1250 | yes | yes |
| 41 | 8 | 2,524 | 0.0608 | 0.0607 | 0.0625 | no | no |
| 42 | 4 | 5,167 | 0.1244 | 0.1244 | 0.1250 | yes | yes |
| 43 | 3 | 6,916 | 0.1665 | 0.1665 | 0.1667 | no | no |
| 44 | 6 | 3,414 | 0.0822 | 0.0822 | 0.0833 | no | no |
| 45 | 4 | 5,133 | 0.1236 | 0.1236 | 0.1250 | yes | yes |
| 46 | 4 | 5,155 | 0.1241 | 0.1241 | 0.1250 | no | no |
| 47 | 5 | 4,090 | 0.0985 | 0.0984 | 0.1000 | no | no |
| 48 | 4 | 5,137 | 0.1237 | 0.1237 | 0.1250 | yes | yes |
| 49 | 4 | 5,145 | 0.1239 | 0.1239 | 0.1250 | no | no |
| 50 | 6 | 3,447 | 0.0830 | 0.0830 | 0.0833 | no | no |
| 51 | 6 | 3,403 | 0.0819 | 0.0819 | 0.0833 | no | no |
| 52 | 4 | 5,146 | 0.1239 | 0.1239 | 0.1250 | no | no |
| 53 | 6 | 3,437 | 0.0827 | 0.0827 | 0.0833 | no | no |
| 54 | 6 | 3,416 | 0.0822 | 0.0822 | 0.0833 | no | no |
| 55 | 4 | 5,147 | 0.1239 | 0.1239 | 0.1250 | no | no |
| 56 | 8 | 2,553 | 0.0615 | 0.0615 | 0.0625 | no | no |
| 57 | 4 | 5,128 | 0.1235 | 0.1235 | 0.1250 | yes | yes |
| 58 | 2 | 10,380 | 0.2499 | 0.2499 | 0.2500 | yes | yes |
| 59 | 9 | 2,270 | 0.0546 | 0.0546 | 0.0556 | no | no |
| 60 | 4 | 5,133 | 0.1236 | 0.1236 | 0.1250 | yes | yes |
| 61 | 6 | 3,416 | 0.0822 | 0.0822 | 0.0833 | no | no |
| 62 | 8 | 2,556 | 0.0615 | 0.0615 | 0.0625 | no | no |
| 63 | 4 | 5,160 | 0.1242 | 0.1242 | 0.1250 | no | no |
| 64 | 4 | 5,139 | 0.1237 | 0.1237 | 0.1250 | no | no |
| 65 | 8 | 2,543 | 0.0612 | 0.0612 | 0.0625 | no | no |
| 66 | 8 | 2,552 | 0.0614 | 0.0614 | 0.0625 | no | no |
| 67 | 3 | 6,840 | 0.1647 | 0.1647 | 0.1667 | no | no |
| 68 | 8 | 2,561 | 0.0617 | 0.0617 | 0.0625 | no | no |
| 69 | 8 | 2,531 | 0.0609 | 0.0609 | 0.0625 | no | no |
| 70 | 4 | 5,142 | 0.1238 | 0.1238 | 0.1250 | yes | yes |
| 71 | 7 | 2,914 | 0.0702 | 0.0701 | 0.0714 | no | no |
| 72 | 4 | 5,137 | 0.1237 | 0.1237 | 0.1250 | yes | yes |
| 73 | 4 | 5,133 | 0.1236 | 0.1236 | 0.1250 | no | no |
| 74 | 10 | 2,049 | 0.0493 | 0.0493 | 0.0500 | no | no |
| 75 | 6 | 3,426 | 0.0825 | 0.0825 | 0.0833 | no | no |
| 76 | 6 | 3,417 | 0.0823 | 0.0823 | 0.0833 | no | no |
| 77 | 8 | 2,531 | 0.0609 | 0.0609 | 0.0625 | no | no |
| 78 | 4 | 5,106 | 0.1229 | 0.1229 | 0.1250 | yes | yes |
| 79 | 5 | 4,105 | 0.0988 | 0.0988 | 0.1000 | no | no |
| 80 | 8 | 2,554 | 0.0615 | 0.0615 | 0.0625 | no | no |
| 81 | 6 | 3,436 | 0.0827 | 0.0827 | 0.0833 | no | no |
| 82 | 4 | 5,129 | 0.1235 | 0.1235 | 0.1250 | no | no |
| 83 | 9 | 2,258 | 0.0544 | 0.0543 | 0.0556 | no | no |
| 84 | 8 | 2,552 | 0.0614 | 0.0614 | 0.0625 | no | no |
| 85 | 4 | 5,120 | 0.1233 | 0.1233 | 0.1250 | yes | yes |
| 86 | 10 | 2,021 | 0.0487 | 0.0487 | 0.0500 | no | no |
| 87 | 6 | 3,421 | 0.0824 | 0.0824 | 0.0833 | no | no |
| 88 | 4 | 5,146 | 0.1239 | 0.1239 | 0.1250 | yes | yes |
| 89 | 12 | 1,687 | 0.0406 | 0.0406 | 0.0417 | no | no |
| 90 | 8 | 2,553 | 0.0615 | 0.0615 | 0.0625 | no | no |
| 91 | 6 | 3,405 | 0.0820 | 0.0820 | 0.0833 | no | no |
| 92 | 6 | 3,439 | 0.0828 | 0.0828 | 0.0833 | no | no |
| 93 | 4 | 5,124 | 0.1234 | 0.1234 | 0.1250 | yes | yes |
| 94 | 8 | 2,561 | 0.0617 | 0.0617 | 0.0625 | no | no |
| 95 | 8 | 2,554 | 0.0615 | 0.0615 | 0.0625 | no | no |
| 96 | 8 | 2,559 | 0.0616 | 0.0616 | 0.0625 | no | no |
| 97 | 4 | 5,138 | 0.1237 | 0.1237 | 0.1250 | no | no |
| 98 | 8 | 2,555 | 0.0615 | 0.0615 | 0.0625 | no | no |
| 99 | 6 | 3,398 | 0.0818 | 0.0818 | 0.0833 | no | no |
| 100 | 4 | 5,166 | 0.1244 | 0.1244 | 0.1250 | no | no |

### Why every count is a little below 1/(2h)

In all 100 cases the fraction (excluding p | 4P) is below 1/(2h). This is a Chebyshev-type bias. The primes of the principal form are exactly the primes whose Frobenius in Gal(L/Q) is the identity, where L is the ring class field and G is generalized dihedral of order 2h. Counting prime powers with weight 1/k gives no bias. Removing the prime squares p² ≤ x takes away the most from the identity class, because every element of order ≤ 2 squares to it. The mean deficit is then about ½·(I − 1)/(2h)·π(√N), where I = #Cl[2] + h = #{g : g² = 1}. This is the Rubinstein–Sarnak heuristic carried over to Chebotarev sets (Ng 2000; Bellaïche 2016). It is a heuristic, and the value at a single x fluctuates.

`bias.py` (rechecked by `bias_check.py`) uses π(√N) = 126. All 100 of 100 P have a positive deficit (expected count n_excl/(2h) minus drawn). The observed deficits run from 4.0 to 101.0 primes (median 45.7). The predicted mean runs from 31.5 to 55.1 (median 39.4). Summed over all P, 4693 were observed against 4126 predicted, with correlation 0.36 across P. So the sign and the size are explained, but the scatter from P to P is not.

## 3. The classical criteria for non-idoneal P (every prime ≤ 500,000)

Method A is `qf_extra.py`. Method B is `qf_extra_check.py`, whose representations come from PARI `qfbsolve` and are independent of method A.

| P | h(-4P) | criterion (method A, qf_extra.py) | primes checked | agree | exceptions | method B (qf_extra_check.py): class polynomial H_{-4P} | method B: simple polynomial |
|---:|---:|---|---:|---:|---|---|---|
| 27 | 3 | p = 1 (mod 3) and 2^((p-1)/3) = 1 (mod p) | 41,538 | 41,538 | 0 | 41,525 of 41,525 agree (skipped p dividing 4P*disc H: [2, 3, 5, 11, 17, 23, 29, 41, 47, 59, 71, 89, 101]) | p = 1 mod 3 and x^3 - 2 has a root: 41,538 of 41,538 agree; same field as polredabs(H) = x^3 - 2: True |
| 64 | 4 | p = 1 (mod 4) and 2^((p-1)/4) = 1 (mod p) | 41,538 | 41,538 | 0 | 41,517 of 41,517 agree (skipped p dividing 4P*disc H: [2, 3, 7, 11, 19, 23, 31, 43, 47, 59, 71, 83, 103, 107, 131, 139, 191, 211, 227, 239, 251]) | p = 1 mod 4 and x^4 - 2 has a root: 41,538 of 41,538 agree; same field as polredabs(H) = x^4 - 2: True |
| 14 | 4 | (-14/p) = 1 and (x^2+1)^2 = 8 (mod p) solvable | 41,536 | 41,536 | 0 (set aside: p=2: drawn=no, criterion=no; p=7: drawn=no, criterion=no) | 41,527 of 41,527 agree (skipped p dividing 4P*disc H: [2, 7, 11, 17, 29, 31, 37, 41, 43, 47, 53]) | (-14/p) = 1 and (x^2+1)^2 - 8 has a root: 41,538 of 41,538 agree; same field as polredabs(H) = x^4 - 2*x^3 + x^2 - 2*x + 1: True |
| 11 | 3 | x^3 - x^2 + x + 1 has three roots mod p | 41,536 | 41,536 | 0 (set aside: p=2: drawn=no, criterion=no; p=11: drawn=yes, criterion=no) | 41,530 of 41,530 agree (skipped p dividing 4P*disc H: [2, 7, 11, 13, 17, 19, 29, 41]) | x^3 - x^2 + x + 1 has three roots: 41,537 of 41,538 agree, exceptions [11]; same field as polredabs(H) = x^3 - x^2 + x + 1: True |
| 23 | 3 | x^3 - x - 1 has three roots mod p | 41,536 | 41,536 | 0 (set aside: p=2: drawn=no, criterion=no; p=23: drawn=yes, criterion=no) | 41,526 of 41,526 agree (skipped p dividing 4P*disc H: [2, 5, 7, 11, 17, 19, 23, 53, 61, 79, 83, 89]) | x^3 - x - 1 has three roots: 41,537 of 41,538 agree, exceptions [23]; same field as polredabs(H) = x^3 - x^2 + 1: True |

Further consistency rows of qf_extra.py:

| statement | primes checked | agree | disagree |
|---|---:|---:|---:|
| P=11 (-11/p)=1 and x^3-x^2+x+1 has a root | 41,536 | 41,536 | 0 |
| P=27 x^3-2 three roots | 41,536 | 41,536 | 0 |
| P=64 (p=1 mod 4) and x^4-2 has a root | 41,537 | 41,537 | 0 |
| P=14 x^4+2x^2-7 has a root and (-14/p)=1 | 41,536 | 41,536 | 0 |
| x^3-x-1 root-count distribution (p != 23) | {'0': 13809, '1': 20818, '3': 6910} | | |

Notes:
* **P = 11 and 23 need "p ≠ P"**. The only disagreements anywhere are p = 11 for P = 11 and p = 23 for P = 23. 11 = 0² + 11·1² and 23 = 0² + 23·1² are drawn (a = 0), but p divides the discriminant of the cubic (−44 and −23), so the cubic has a double root there instead of three roots. The verified statements are:
  * p = x² + 23y² ⟺ p ≠ 23 and x³ − x − 1 has three roots mod p.
  * p = x² + 11y² ⟺ p ≠ 11 and x³ − x² + x + 1 has three roots mod p. An equivalent form is (−11/p) = 1 and the cubic has a root, which also agrees on all 41,536 primes p ≠ 2, 11.
* **Why these cubics**: x³ − x² + x + 1 has discriminant −44 and is exactly PARI's `polredabs` of the class polynomial H₋₄₄ (the minimal polynomial of j(√−11)). For P = 23, `polredabs(H₋₉₂)` = x³ − x² + 1, which `nfisisom` shows defines the same field as x³ − x − 1 (h(−92) = h(−23) = 3, so the ring class field of Z[√−23] is the Hilbert class field of Q(√−23)). In the same way, `polredabs(H₋₁₀₈)` = x³ − 2, `polredabs(H₋₂₅₆)` = x⁴ − 2, and `polredabs(H₋₅₆)` defines the same field as (x² + 1)² − 8. This is how Gauss's cubic and biquadratic rules for 2 and Cox's x² + 14y² rule fit Cox Thm 9.2 (p = x² + ny² ⟺ (−n/p) = 1 and the class polynomial has a root mod p, for p ∤ n·disc).
* The rules for P = 27, 64 and 14 hold for every prime, including 2, 3 and 7, where both sides are false.
* **Chebotarev for x³ − x − 1** (p ≠ 23): 0 roots for 13,809 primes (33.2%), 1 root for 20,818 primes (50.1%), 3 roots for 6,910 primes (16.6%). The densities 1/3, 1/2, 1/6 are those of S₃.

## 4. eisenstein

Source: `eis.py`. Rechecked by `eis_check.py`.

* **Primes drawn: 20,734**, a fraction 0.49916 of π(N) (0.49915 without 3). Density 1/2, since h(−3) = 1.
* **The set is exactly {3} ∪ {p ≡ 1 (mod 3)}**: True in `eis.py`, True in `eis_check.py`. There are 20,733 primes ≡ 1 and 20,804 primes ≡ 2 (mod 3), so the same Chebyshev bias shows here.
* **Positions.** 3 is drawn at (a, b) = (1, 1), on the edge a = b of the wedge; its orbit under the 12 automorphisms has 6 points. Every other drawn prime has an orbit of 12 points with exactly one in the wedge (orbit sizes {'6': 1, '12': 20733}). Every drawn prime's position in `eis.py` matches `eis_check.py` (0 mismatches).
* 93,120 of the numbers 1..N are drawn at all (Löschian numbers).

## 5. Bright straight lines (P = 1, 2, 3 and eisenstein)

Source: `lines.py` → `lines.json`. Every line below was rechecked by `lines_check.py` (column "rechecked").

**Method.** A line is a set of lattice points of R with db·a − da·b = key, for a primitive step (da, db). Along the line, n = f(t) is a quadratic. The Bateman–Horn prediction is C(f)·Σ 1/ln f(t), summed over the line's lattice points in R with f ≥ 2. Here C(f) = ∏_{p ≤ 10,000} (1 − ω(p)/p)/(1 − 1/p), and ω(p) is the number of roots of f mod p.

**Contrast.** What the eye sees is the density of primes on a line compared with the area around it. A lattice point of R with norm n is prime with probability about B/ln n, where B = (share of primes drawn)/(lattice points of R per unit of n):
* P = 1: B = (1/2)/(π/8) = 4/π.
* P ≥ 2: B = (1/(2h))·4√P/π.
* eisenstein: B = (1/2)/(π/(6√3)) = 3√3/π.

The empirical B (primes in R divided by Σ 1/ln n over R) agrees to within 1%. A line with C/B = 3 is lit three times as densely as its surroundings.

#### quadform P = 1  (background B = 1.2732, empirical 1.2644; 20,732 primes on 196,955 lattice points)

*5 lines with the most primes (b, a, a - b or a + b fixed)*

| line | n along the line | points | primes | Bateman-Horn | obs/pred | C | C/B | rechecked |
|---|---|---:|---:|---:|---:|---:|---:|:---:|
| a - b = 5 | 2b^2 + 10b + 25 | 498 | 169 | 170.2 | 0.993 | 3.6561 | 2.87 | yes |
| a - b = 65 | 2b^2 + 130b + 4225 | 467 | 164 | 163.3 | 1.004 | 3.9884 | 3.13 | yes |
| a - b = 85 | 2b^2 + 170b + 7225 | 456 | 156 | 154.3 | 1.011 | 3.8998 | 3.06 | yes |
| a - b = 25 | 2b^2 + 50b + 625 | 488 | 154 | 161.2 | 0.955 | 3.6561 | 2.87 | yes |
| a - b = 35 | 2b^2 + 70b + 1225 | 483 | 137 | 135.5 | 1.011 | 3.1338 | 2.46 | yes |

*5 brightest relative to the background (same directions, >= 100 points)*

| line | n along the line | points | primes | Bateman-Horn | obs/pred | C | C/B | rechecked |
|---|---|---:|---:|---:|---:|---:|---:|:---:|
| a - b = 65 | 2b^2 + 130b + 4225 | 467 | 164 | 163.3 | 1.004 | 3.9884 | 3.13 | yes |
| a - b = 325 | 2b^2 + 650b + 105625 | 311 | 99 | 100.1 | 0.989 | 3.9884 | 3.13 | yes |
| a + b = 845 | 2b^2 - 1690b + 714025 | 267 | 78 | 82.5 | 0.945 | 3.9884 | 3.13 | yes |
| a + b = 325 | 2b^2 - 650b + 105625 | 163 | 58 | 58.4 | 0.993 | 3.9884 | 3.13 | yes |
| a - b = 85 | 2b^2 + 170b + 7225 | 456 | 156 | 154.3 | 1.011 | 3.8998 | 3.06 | yes |

*brightest in any direction (steps up to 3), >= 100 points*

| line | n along the line | points | primes | Bateman-Horn | obs/pred | C | C/B | rechecked |
|---|---|---:|---:|---:|---:|---:|---:|:---:|
| a - 3b = -221 | 10b^2 - 1326b + 48841 | 178 | 57 | 63.6 | 0.896 | 4.2543 | 3.34 | yes |
| a - 3b = 221 | 10b^2 + 1326b + 48841 | 157 | 56 | 55.0 | 1.017 | 4.2543 | 3.34 | yes |
| a - 3b = -377 | 10b^2 - 2262b + 142129 | 145 | 49 | 48.9 | 1.002 | 4.1361 | 3.25 | yes |
| a - 3b = 377 | 10b^2 + 2262b + 142129 | 108 | 36 | 35.6 | 1.011 | 4.1361 | 3.25 | yes |
| a + 3b = 481 | 10b^2 - 2886b + 231361 | 121 | 45 | 43.7 | 1.030 | 4.1024 | 3.22 | yes |

*reference lines*

| line | n along the line | points | primes | Bateman-Horn | obs/pred | C | C/B | rechecked |
|---|---|---:|---:|---:|---:|---:|---:|:---:|
| b = 1 | a^2 + 1 | 707 | 89 | 93.6 | 0.950 | 1.3710 | 1.08 | yes |
| a - b = 1 | 2b^2 + 2b + 1 | 500 | 130 | 129.8 | 1.002 | 2.7420 | 2.15 | yes |
| b = 2 | a^2 + 4 | 706 | 92 | 91.3 | 1.007 | 1.3710 | 1.08 | yes |

#### quadform P = 2  (background B = 0.9003, empirical 0.8960; 20,753 primes on 278,269 lattice points)

*5 lines with the most primes (b, a, a - b or a + b fixed)*

| line | n along the line | points | primes | Bateman-Horn | obs/pred | C | C/B | rechecked |
|---|---|---:|---:|---:|---:|---:|---:|:---:|
| a = 33 | 2b^2 + 1089 | 500 | 140 | 144.2 | 0.971 | 3.1669 | 3.52 | yes |
| a = 51 | 2b^2 + 2601 | 499 | 138 | 136.5 | 1.011 | 3.0402 | 3.38 | yes |
| a = 99 | 2b^2 + 9801 | 496 | 138 | 137.9 | 1.001 | 3.1669 | 3.52 | yes |
| a = 57 | 2b^2 + 3249 | 499 | 137 | 135.0 | 1.015 | 3.0178 | 3.35 | yes |
| a = 9 | 2b^2 + 81 | 500 | 131 | 133.1 | 0.984 | 2.8502 | 3.17 | yes |

*5 brightest relative to the background (same directions, >= 100 points)*

| line | n along the line | points | primes | Bateman-Horn | obs/pred | C | C/B | rechecked |
|---|---|---:|---:|---:|---:|---:|---:|:---:|
| a = 561 | 2b^2 + 314721 | 305 | 77 | 80.3 | 0.959 | 3.3780 | 3.75 | yes |
| a = 627 | 2b^2 + 393129 | 232 | 57 | 60.0 | 0.950 | 3.3532 | 3.72 | yes |
| a = 33 | 2b^2 + 1089 | 500 | 140 | 144.2 | 0.971 | 3.1669 | 3.52 | yes |
| a = 99 | 2b^2 + 9801 | 496 | 138 | 137.9 | 1.001 | 3.1669 | 3.52 | yes |
| a = 297 | 2b^2 + 88209 | 454 | 113 | 118.3 | 0.955 | 3.1669 | 3.52 | yes |

*brightest in any direction (steps up to 3), >= 100 points*

| line | n along the line | points | primes | Bateman-Horn | obs/pred | C | C/B | rechecked |
|---|---|---:|---:|---:|---:|---:|---:|:---:|
| 3a + 2b = 969 | 22t^2 - 1292t + 104329 | 162 | 47 | 48.4 | 0.971 | 3.5767 | 3.97 | yes |
| 3a - 2b = 969 | 22t^2 + 1292t + 104329 | 108 | 27 | 31.3 | 0.862 | 3.5767 | 3.97 | yes |
| 3a + 2b = 2091 | 22t^2 - 2788t + 485809 | 132 | 29 | 35.3 | 0.822 | 3.4646 | 3.85 | yes |
| 3a + 2b = 2193 | 22t^2 - 2924t + 534361 | 106 | 23 | 28.1 | 0.817 | 3.4604 | 3.84 | yes |
| a = 561 | 2b^2 + 314721 | 305 | 77 | 80.3 | 0.959 | 3.3780 | 3.75 | yes |

*reference lines*

| line | n along the line | points | primes | Bateman-Horn | obs/pred | C | C/B | rechecked |
|---|---|---:|---:|---:|---:|---:|---:|:---:|
| b = 1 | a^2 + 2 | 708 | 54 | 49.2 | 1.097 | 0.7125 | 0.79 | yes |
| a = 1 | 2b^2 + 1 | 500 | 69 | 68.1 | 1.012 | 1.4251 | 1.58 | yes |
| a - b = 1 | 3b^2 + 2b + 1 | 408 | 52 | 55.0 | 0.946 | 1.4251 | 1.58 | yes |

#### quadform P = 3  (background B = 1.1027, empirical 1.0959; 20,734 primes on 227,298 lattice points)

*5 lines with the most primes (b, a, a - b or a + b fixed)*

| line | n along the line | points | primes | Bateman-Horn | obs/pred | C | C/B | rechecked |
|---|---|---:|---:|---:|---:|---:|---:|:---:|
| a + b = 427 | 4b^2 - 854b + 182329 | 409 | 95 | 91.2 | 1.042 | 2.7308 | 2.48 | yes |
| b = 7 | a^2 + 147 | 708 | 93 | 88.3 | 1.053 | 1.3427 | 1.22 | yes |
| a + b = 399 | 4b^2 - 798b + 159201 | 400 | 92 | 93.6 | 0.983 | 2.8433 | 2.58 | yes |
| a + b = 637 | 4b^2 - 1274b + 405769 | 381 | 92 | 87.4 | 1.053 | 2.9294 | 2.66 | yes |
| a + b = 413 | 4b^2 - 826b + 170569 | 409 | 91 | 88.5 | 1.029 | 2.6398 | 2.39 | yes |

*5 brightest relative to the background (same directions, >= 100 points)*

| line | n along the line | points | primes | Bateman-Horn | obs/pred | C | C/B | rechecked |
|---|---|---:|---:|---:|---:|---:|---:|:---:|
| a + b = 637 | 4b^2 - 1274b + 405769 | 381 | 92 | 87.4 | 1.053 | 2.9294 | 2.66 | yes |
| a - b = 91 | 4b^2 + 182b + 8281 | 329 | 81 | 83.8 | 0.967 | 2.9294 | 2.66 | yes |
| a - b = -91 | 4b^2 - 182b + 8281 | 284 | 69 | 69.9 | 0.987 | 2.9294 | 2.66 | yes |
| a + b = 273 | 4b^2 - 546b + 74529 | 274 | 68 | 70.5 | 0.965 | 2.9294 | 2.66 | yes |
| a - b = 273 | 4b^2 + 546b + 74529 | 265 | 63 | 63.6 | 0.991 | 2.9294 | 2.66 | yes |

*brightest in any direction (steps up to 3), >= 100 points*

| line | n along the line | points | primes | Bateman-Horn | obs/pred | C | C/B | rechecked |
|---|---|---:|---:|---:|---:|---:|---:|:---:|
| a - 3b = -91 | 12b^2 - 546b + 8281 | 196 | 77 | 76.4 | 1.008 | 4.3941 | 3.99 | yes |
| a + 3b = 637 | 12b^2 - 3822b + 405769 | 213 | 71 | 78.1 | 0.909 | 4.3941 | 3.99 | yes |
| a + 3b = 1183 | 12b^2 - 7098b + 1399489 | 211 | 67 | 72.0 | 0.931 | 4.3941 | 3.99 | yes |
| a - 3b = 91 | 12b^2 + 546b + 8281 | 181 | 65 | 68.7 | 0.947 | 4.3941 | 3.99 | yes |
| a - 3b = -637 | 12b^2 - 3822b + 405769 | 129 | 48 | 45.4 | 1.056 | 4.3941 | 3.99 | yes |

*reference lines*

| line | n along the line | points | primes | Bateman-Horn | obs/pred | C | C/B | rechecked |
|---|---|---:|---:|---:|---:|---:|---:|:---:|
| b = 1 | a^2 + 3 | 708 | 84 | 76.4 | 1.099 | 1.1189 | 1.01 | yes |
| a = 1 | 3b^2 + 1 | 409 | 63 | 65.3 | 0.964 | 1.6783 | 1.52 | yes |
| a - b = 1 | 4b^2 + 2b + 1 | 354 | 70 | 74.8 | 0.936 | 2.2378 | 2.03 | yes |
| a = 2 | 3b^2 + 4 | 409 | 66 | 66.1 | 0.998 | 1.6783 | 1.52 | yes |

#### eisenstein  (background B = 1.6540, empirical 1.6415; 20,734 primes on 151,715 lattice points)

*5 lines with the most primes (b, a, a - b or a + b fixed)*

| line | n along the line | points | primes | Bateman-Horn | obs/pred | C | C/B | rechecked |
|---|---|---:|---:|---:|---:|---:|---:|:---:|
| b = 7 | a^2 + 7a + 49 | 697 | 170 | 171.6 | 0.991 | 2.6853 | 1.62 | yes |
| b = 21 | a^2 + 21a + 441 | 676 | 163 | 161.9 | 1.007 | 2.6853 | 1.62 | yes |
| b = 49 | a^2 + 49a + 2401 | 633 | 155 | 147.0 | 1.054 | 2.6853 | 1.62 | yes |
| a - b = 7 | 3b^2 + 21b + 49 | 405 | 150 | 151.3 | 0.991 | 4.0280 | 2.44 | yes |
| b = 13 | a^2 + 13a + 169 | 688 | 146 | 151.8 | 0.962 | 2.4412 | 1.48 | yes |

*5 brightest relative to the background (same directions, >= 100 points)*

| line | n along the line | points | primes | Bateman-Horn | obs/pred | C | C/B | rechecked |
|---|---|---:|---:|---:|---:|---:|---:|:---:|
| a - b = 91 | 3b^2 + 273b + 8281 | 362 | 132 | 137.2 | 0.962 | 4.3941 | 2.66 | yes |
| a - b = 133 | 3b^2 + 399b + 17689 | 340 | 120 | 122.9 | 0.976 | 4.2649 | 2.58 | yes |
| a - b = 217 | 3b^2 + 651b + 47089 | 295 | 92 | 101.5 | 0.906 | 4.1669 | 2.52 | yes |
| a - b = 259 | 3b^2 + 777b + 67081 | 272 | 94 | 92.1 | 1.021 | 4.1430 | 2.50 | yes |
| a - b = 301 | 3b^2 + 903b + 90601 | 249 | 88 | 83.2 | 1.058 | 4.1262 | 2.49 | yes |

*brightest in any direction (steps up to 3), >= 100 points*

| line | n along the line | points | primes | Bateman-Horn | obs/pred | C | C/B | rechecked |
|---|---|---:|---:|---:|---:|---:|---:|:---:|
| a - b = 91 | 3b^2 + 273b + 8281 | 362 | 132 | 137.2 | 0.962 | 4.3941 | 2.66 | yes |
| a + 2b = 637 | 3b^2 - 1911b + 405769 | 213 | 76 | 75.7 | 1.004 | 4.3941 | 2.66 | yes |
| 2a + b = 1183 | 3t^2 + 3t + 349873 | 197 | 64 | 67.3 | 0.951 | 4.3941 | 2.66 | yes |
| 2a + b = 637 | 3t^2 + 3t + 101443 | 106 | 33 | 40.1 | 0.824 | 4.3941 | 2.66 | yes |
| a - b = 133 | 3b^2 + 399b + 17689 | 340 | 120 | 122.9 | 0.976 | 4.2649 | 2.58 | yes |

*reference lines*

| line | n along the line | points | primes | Bateman-Horn | obs/pred | C | C/B | rechecked |
|---|---|---:|---:|---:|---:|---:|---:|:---:|
| b = 1 | a^2 + a + 1 | 706 | 140 | 150.8 | 0.928 | 2.2378 | 1.35 | yes |
| a - b = 1 | 3b^2 + 3b + 1 | 408 | 126 | 129.2 | 0.976 | 3.3566 | 2.03 | yes |

Over the 64 distinct lines above: 5,905 primes observed, 5,991.5 predicted by Bateman-Horn (ratio 0.986); per line obs/pred between 0.817 and 1.099.

**What a viewer sees, and why**:
* **Parity sets the stripes.** For P odd (1, 3), n ≡ a + b (mod 2), so only the points with a + b odd can be prime, a checkerboard. Each 45° diagonal a − b = c has a single parity. Odd diagonals carry every prime and even ones none, which gives the striped octant of P = 1. For P even, n ≡ a (mod 2), so every even column is empty and P = 2 shows vertical stripes. For eisenstein, n is even only when a and b are both even, and n ≡ (a − b)² (mod 3). The lines a − b ≡ 0 (mod 3) are therefore dark, and the lines a − b = c with 3 ∤ c gain a factor 3/2, which is why a − b = 7 and 91 are bright.
* **Offsets divisible by split primes make a line brighter.** Take P = 1 and the line a − b = c, where n = b² + (b + c)². If q ≡ 1 (mod 4) divides c, the line meets the dark multiples of q once per period instead of twice, a factor (q − 1)/(q − 2). If q ≡ 3 (mod 4) divides c, it meets them once instead of never, a factor (q − 1)/q. The brightest P = 1 lines therefore have offsets 5, 25, 65 = 5·13 and 85 = 5·17. The same holds for P = 2 (columns a = 9, 33 = 3·11, 51 = 3·17, 57 = 3·19, 99, where 3, 11, 17 and 19 split), P = 3 (a ± b = 91, 273, 637: 7 and 13 split) and eisenstein (b = 7, 21, 49, 13 and a − b = 91).
* **Landau's line does not stand out.** b = 1 in P = 1 (a² + 1) holds 89 primes against a prediction of 93.6, but its C/B is only 1.08, so it is hardly brighter than its surroundings. The row b = 1 of P = 2 (a² + 2) has C/B = 0.79 and is darker than the background.
* **Steeper lines (step 3) reach C/B ≈ 3.3–4.0**, for example P = 1 a − 3b = ±221 and P = 3 a ± 3b = ±91, 637. Their points are √10 or more apart, so they read as dotted rather than solid lines.

## 6. Which P to show

Per-P picture numbers come from `picture_stats.py` → `picture_stats.json` and are rechecked by `picture_check.py`. The nesting facts come from `subsets.py`, rechecked by `subsets_check.py`.

| P | h | dots | lattice points of R | dots/points | B | row b=1: C/B | best row b ≤ 20 (C/B) | worst row (C/B) | mixed classes mod 4P |
|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| 1 | 1 | 20,732 | 196,955 | 0.1053 | 1.273 | 1.08 | b = 20 (1.44) | b = 3 (0.72) | 0 of 2 |
| 2 | 1 | 20,753 | 278,269 | 0.0746 | 0.900 | 0.79 | b = 18 (1.58) | b = 5 (0.63) | 0 of 4 |
| 3 | 1 | 20,734 | 227,298 | 0.0912 | 1.103 | 1.01 | b = 14 (1.22) | b = 5 (0.81) | 0 of 4 |
| 4 | 1 | 20,731 | 196,884 | 0.1053 | 1.273 | 1.08 | b = 20 (1.44) | b = 3 (0.72) | 0 of 8 |
| 5 | 2 | 10,345 | 176,139 | 0.0587 | 0.712 | 0.74 | b = 18 (1.48) | b = 11 (0.67) | 0 of 8 |
| 7 | 1 | 20,706 | 148,920 | 0.1390 | 1.684 | 1.17 | b = 11 (1.30) | b = 15 (0.62) | 0 of 12 |
| 11 | 3 | 6,858 | 118,869 | 0.0577 | 0.704 | 0.72 | b = 15 (1.93) | b = 7 (0.62) | 10 of 20 |
| 14 | 4 | 5,151 | 105,383 | 0.0489 | 0.596 | 0.70 | b = 15 (1.88) | b = 11 (0.64) | 6 of 24 |
| 16 | 2 | 10,334 | 98,613 | 0.1048 | 1.273 | 1.08 | b = 20 (1.44) | b = 3 (0.72) | 0 of 32 |
| 23 | 3 | 6,911 | 82,319 | 0.0840 | 1.018 | 0.80 | b = 18 (1.61) | b = 5 (0.64) | 22 of 44 |
| 27 | 3 | 6,898 | 76,010 | 0.0908 | 1.103 | 1.01 | b = 14 (1.22) | b = 5 (0.81) | 18 of 36 |
| 37 | 2 | 10,284 | 64,985 | 0.1583 | 1.936 | 1.21 | b = 19 (1.28) | b = 15 (0.64) | 0 of 72 |
| 41 | 8 | 2,524 | 61,745 | 0.0409 | 0.510 | 0.68 | b = 15 (1.83) | b = 13 (0.63) | 20 of 80 |
| 58 | 2 | 10,380 | 51,952 | 0.1998 | 2.424 | 1.23 | b = 16 (1.23) | b = 15 (0.65) | 0 of 112 |
| 64 | 4 | 5,139 | 49,493 | 0.1038 | 1.273 | 1.08 | b = 20 (1.44) | b = 3 (0.72) | 32 of 128 |
| 67 | 3 | 6,840 | 48,381 | 0.1414 | 1.737 | 1.20 | b = 17 (1.28) | b = 15 (0.64) | 66 of 132 |
| 89 | 12 | 1,687 | 41,995 | 0.0402 | 0.500 | 0.69 | b = 15 (1.83) | b = 11 (0.62) | 44 of 176 |

* **Densest pictures** (dots per lattice point of R): P = 58 (0.200), P = 37 (0.158), P = 67 (0.141), P = 7 (0.139), P = 28 (0.138), P = 25 (0.131). **Sparsest**: P = 59 (0.044), P = 41 (0.041), P = 89 (0.040). `picture_check.py` confirms P = 58 as the densest.
* **Strongest single lines** (rows, columns and diagonals of index ≤ 20 with ≥ 100 points, over all P): column a = 15 for P = 14 (C/B = 4.38), diagonal a − b = 5 for P = 41 (4.38), row b = 15 for P = 71 (1.95). In this range, rows (fixed b) never exceed 1.95.
* **Large P.** The page draws (a, b√P), so the picture becomes a stack of horizontal rows √P apart (P = 58: 93 rows). Row brightness varies by about a factor 2 to 3 (best and worst rows above), following the primes that divide b.
* **How the pictures nest** (`subsets.py`, rechecked by `subsets_check.py`; every flag True in both):
  * P = 4 has the same 20,731 odd primes as P = 1 (all p ≡ 1 mod 4).
  * P = 16 has the 10,334 primes whose even coordinate in the P = 1 picture is divisible by 4 (⟺ p ≡ 1 mod 8: True).
  * P = 64 has the 5,139 with the even coordinate divisible by 8 (True).
  * P = 12 has the P = 3 primes with b even (10,329).
  * P = 27 has the P = 3 primes with 3 | b: 6,898, a share of 0.3327 of the 20,733 primes ≡ 1 (mod 3).
* **P = 5.** (−5/p) = 1 for p ≡ 1, 3, 7, 9 (mod 20); x² + 5y² takes 1, 9, 2x² + 2xy + 3y² takes 3, 7; each such prime is in exactly one form: True (rechecked: True).
* **P = 58.** 2x² + 29 is prime for x = 0..28 (29 values; checked: True). x² + 58y² takes 28 classes mod 232, and 2x² + 29y² takes 28 others (disjoint: True).

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
