# Praime
## Beating prime puzzle records with LLM reasoning 
This repository is an experiment in how far LLM can go with hard math: puzzles, logical reasoning, and problems nobody has cracked yet. The test case is Carlos Rivera's Prime Puzzles site, which has collected 1200+ prime challenges since 1999, many with records that stood for decades. 

## Results 30-10-2026
**39** beaten records, **26** answers to open questions (a few of them partial or conditional), **9** results on the conjectures pages (2 disproved, 4 proved, 1 explained, 2 open questions answered) and **1** mistake found on the pages (puzzle 1209).

#### Beaten records
| # | Puzzle | Current record | Our result |
  |---|---|---|---|
  | [65](https://www.primepuzzles.net/puzzles/puzz_065.htm) | Equal power sums from consecutive primes | degree 3 | [degree 5](results/65.md) |
  | [87](https://www.primepuzzles.net/puzzles/puzz_087.htm) | Prime magic squares A+B+1 / A+B+C | 1263 / 3327 | [1041 / 789](results/87.md) |
  | [92](https://www.primepuzzles.net/puzzles/puzz_092.htm) | Pile of prime spheres (n = 5) | total 2,675,763 | [total 16,561](results/92.md) |
  | [99](https://www.primepuzzles.net/puzzles/puzz_099.htm) | Prime with p² = R(R(p)²) | 50 digits | [3001 digits](results/99.md) |
  | [134](https://www.primepuzzles.net/puzzles/puzz_134.htm) | Carrousel prime, one active digit | 50 digits | [201 digits](results/134.md) |
  | [227](https://www.primepuzzles.net/puzzles/puzz_227.htm) | floor(tʲ) prime for j = 1, 2, … | 100 | [239](results/227.md) |
  | [249](https://www.primepuzzles.net/puzzles/puzz_249.htm) | Pandigital prime magic square, sum | 1530758629047 | [961547036829](results/249.md) |
  | [307](https://www.primepuzzles.net/puzzles/puzz_307.htm) | Largest Smith number (Costello method) | 7.8M digits | [41M digits](results/307.md) |
  | [324](https://www.primepuzzles.net/puzzles/puzz_324.htm) | Largest self-descriptive prime | 70 digits | [74 digits](results/324.md) |
  | [340](https://www.primepuzzles.net/puzzles/puzz_340.htm) | Prime p contained in p² | 1339 digits | [3925 digits](results/340.md) |
  | [342](https://www.primepuzzles.net/puzzles/puzz_342.htm) | Die of squares with prime corner sums | 1330 | [658](results/342.md) |
  | [356](https://www.primepuzzles.net/puzzles/puzz_356.htm) | Prime sets with distinct pair sums | 751 / 1409 / 2273 | [661 / 1297 / 2213](results/356.md) |
  | [359](https://www.primepuzzles.net/puzzles/puzz_359.htm) | First N primes in a circle, no shared digits | 108 | [202 (maximum)](results/359.md) |
  | [373](https://www.primepuzzles.net/puzzles/puzz_373.htm) | Largest self prime | 2²⁴⁰³⁶⁵⁸³ − 1 | [2⁸²⁵⁸⁹⁹³³ − 1](results/373.md) |
  | [398](https://www.primepuzzles.net/puzzles/puzz_398.htm) | Prime no digit insertion keeps prime | 1000 digits | [1057 digits](results/398.md) |
  | [419](https://www.primepuzzles.net/puzzles/puzz_419.htm) | σ(N) as a power of its prime sum | P = 31 | [P = 47](results/419.md) |
  | [438](https://www.primepuzzles.net/puzzles/puzz_438.htm) | Same as 398 | 200 digits | [1057 digits](results/438.md) |
  | [450](https://www.primepuzzles.net/puzzles/puzz_450.htm) | Consecutive-prime sum with most prime factors | 195,005 | [1,290,006](results/450.md) |
  | [476](https://www.primepuzzles.net/puzzles/puzz_476.htm) | Odd Collatz terms all prime | 19 | [30](results/476.md) |
  | [522](https://www.primepuzzles.net/puzzles/puzz_522.htm) | Three runs of consecutive primes | k = 14 | [k = 15](results/522.md) |
  | [553](https://www.primepuzzles.net/puzzles/puzz_553.htm) | Product of primes = sum p..q | 319 | [333](results/553.md) |
  | [597](https://www.primepuzzles.net/puzzles/puzz_597.htm) | Prime/semiprime values of a quadratic | 9567 | [9743](results/597.md) |
  | [607](https://www.primepuzzles.net/puzzles/puzz_607.htm) | Largest prime power missing a digit | 7022 digits | [13.2M digits](results/607.md) |
  | [793](https://www.primepuzzles.net/puzzles/puzz_793.htm) | Every digit pair adjacent (K = 6, 8, 10) | Vantieghem's values | [smaller for all three](results/793.md) |
  | [808](https://www.primepuzzles.net/puzzles/puzz_808.htm) | Palprime palindromic in many bases | 82 | [227](results/808.md) |
  | [830](https://www.primepuzzles.net/puzzles/puzz_830.htm) | Minimal N (K = 2012 / 2013) | 4.86·10¹⁴ / 3.87·10¹¹ | [3.58·10¹³ / 1.74·10¹¹](results/830.md) |
  | [905](https://www.primepuzzles.net/puzzles/puzz_905.htm) | Primes in a row from a recursion | 19 | [28](results/905.md) |
  | [908](https://www.primepuzzles.net/puzzles/puzz_908.htm) | Descending prime remainders | 14 | [34](results/908.md) |
  | [957](https://www.primepuzzles.net/puzzles/puzz_957.htm) | Pandigital expression primes | 489–753 digits | [1685–2989 digits](results/957.md) |
  | [1001](https://www.primepuzzles.net/puzzles/puzz_1001.htm) | Anti-divisor ratio Q/P | 4.781 | [5.520](results/1001.md) |
  | [1028](https://www.primepuzzles.net/puzzles/puzz_1028.htm) | Σ(A·n mod B) prime generator | 41 | [49](results/1028.md) |
  | [1097](https://www.primepuzzles.net/puzzles/puzz_1097.htm) | Modular prime generator | 65 | [130](results/1097.md) |
  | [1099](https://www.primepuzzles.net/puzzles/puzz_1099.htm) | Floor-function prime generator | 120 | [570](results/1099.md) |
  | [1135](https://www.primepuzzles.net/puzzles/puzz_1135.htm) | Primes with no alternating digit parity | 23 | [30](results/1135.md) |
  | [1190](https://www.primepuzzles.net/puzzles/puzz_1190.htm) | Largest 10-digit Friedman prime | 9999999001 | [9999999967](results/1190.md) |
  | [1209](https://www.primepuzzles.net/puzzles/puzz_1209.htm) | Unlucky primes in a row | 45 (invalid, really 26) | [34](results/1209.md) |
  | [1217](https://www.primepuzzles.net/puzzles/puzz_1217.htm) | f(p(n+1), p(n)) prime for successive n | 30,800 | [103,518](results/1217.md) |
  | [1230](https://www.primepuzzles.net/puzzles/puzz_1230.htm) | Primes with even-length period of 1/p | 29 | [30](results/1230.md) |
  | [1251](https://www.primepuzzles.net/puzzles/puzz_1251.htm) | Cycles of "add b, take largest prime factor" | 72 | [477](results/1251.md) |

  ### (First) answers to open questions
| # | Open question | Status on the page | Our answer |
|---|---|---|---|
| [75 Q2](https://www.primepuzzles.net/puzzles/puzz_075.htm) | Smallest palprime with digit sum 2000 | unanswered | [223 digits (PRP)](results/75.md) |
| [77](https://www.primepuzzles.net/puzzles/puzz_077.htm) | Prime pine with 10 rows | open since 2012 | [first one, 55 primes](results/77.md) |
| [81 c)](https://www.primepuzzles.net/puzzles/puzz_081.htm) | Least Sophie-Germain magic squares | unanswered | [magic sums 145060413 / 143168367](results/81.md) |
| [87 Q3](https://www.primepuzzles.net/puzzles/puzz_087.htm) | 4×4 versions | solutions never published | [constants 548 / 570](results/87.md) |
| [92](https://www.primepuzzles.net/puzzles/puzz_092.htm) | Piles for n = 6, 7 | none listed | [first piles](results/92.md) |
| [106](https://www.primepuzzles.net/puzzles/puzz_106.htm) | Prime magic squares summing to 2000 | no solutions | [least k = 4, most k = 26](results/106.md) |
| [134 Q4](https://www.primepuzzles.net/puzzles/puzz_134.htm) | 100-digit carrousel prime, one inactive digit | unanswered | [first example](results/134.md) |
| [170 Q10](https://www.primepuzzles.net/puzzles/puzz_170.htm) | Prime n with n+1 a base-2 pseudoprime | unanswered | [impossible (proof)](results/170.md) |
| [172](https://www.primepuzzles.net/puzzles/puzz_172.htm) | Prime rules for K = 9, 10 | table stops at K = 8 | [L = 113 and L = 144](results/172.md) |
| [179 Q1](https://www.primepuzzles.net/puzzles/puzz_179.htm) | Prime alphametic | unanswered | [TROIS + ONZE = SEPT + SEPT](results/179.md) |
| [189](https://www.primepuzzles.net/puzzles/puzz_189.htm) | Open path cases up to n = 120 | several open | [all settled, threshold 57](results/189.md) |
| [299 Q2](https://www.primepuzzles.net/puzzles/puzz_299.htm) | Smallest increasing run of 9 | only "≥ 7·10⁸" | [35014031359](results/299.md) |
| [325 Q4](https://www.primepuzzles.net/puzzles/puzz_325.htm) | Zeisel number with 10,000+ digits | no contribution | [10,008 digits](results/325.md) |
| [435](https://www.primepuzzles.net/puzzles/puzz_435.htm) | Solve φ(m) + σ(m) = 5m | congruences only (López PDF) | [m = 60493590969525342720](results/435.md) |
| [458](https://www.primepuzzles.net/puzzles/puzz_458.htm) | Nested prime sets of σ(n), n, φ(n) | unanswered | [explicit n](results/458.md) |
| [490](https://www.primepuzzles.net/puzzles/puzz_490.htm) | Other odd n with 2ⁿ − 2ᵏ + 1 a square | unanswered | [none (proof)](results/490.md) |
| [493](https://www.primepuzzles.net/puzzles/puzz_493.htm) | Is pᵖ = xa² + yb² always solvable? | unanswered | [yes (proof)](results/493.md) |
| [516 Q1](https://www.primepuzzles.net/puzzles/puzz_516.htm) | Other n = P + π(P), P = digit product | no contributions | [59 is the only one (proof)](results/516.md) |
| [631](https://www.primepuzzles.net/puzzles/puzz_631.htm) | Three prime gaps all arbitrarily large? | no answers | [yes (Shiu's theorem)](results/631.md) |
| [639](https://www.primepuzzles.net/puzzles/puzz_639.htm) | Composite solutions for every m | no contributions | [yes, under the prime-quadruple conjecture](results/639.md) |
| [677](https://www.primepuzzles.net/puzzles/puzz_677.htm) | Primes in a row with square sums | no solution | [tried values impossible, range narrowed](results/677.md) |
| [890 Q1](https://www.primepuzzles.net/puzzles/puzz_890.htm) | Magic square from Pythagorean triples | no contributions | [impossible (proof)](results/890.md) |
| [900](https://www.primepuzzles.net/puzzles/puzz_900.htm) | Why the numerator of ∏(p²+1)/(p−1) ends in 25 | "a conjecture" (Vantieghem) | [proved; same question as conjecture 80](results/conj-80.md) |
| [1080](https://www.primepuzzles.net/puzzles/puzz_1080.htm) | Third term a(3) | none known | [a(3) ≤ 1.9·10³² (exactly 3 solutions)](results/1080.md) |
| [1145](https://www.primepuzzles.net/puzzles/puzz_1145.htm) | Covering set behind the listed k | no contributions | [found; listed k is its minimum (bounded search)](results/1145.md) |
| [1204](https://www.primepuzzles.net/puzzles/puzz_1204.htm) | Explain the 51.6% | zero answers | [Hardy–Littlewood explanation](results/1204.md) |


## Results 01-11-2026
**1** beaten record and **4** more answers to open questions (two of them partial).

#### Beaten records
| # | Puzzle | Current record | Our result |
|---|---|---|---|
| [882](https://www.primepuzzles.net/puzzles/puzz_882.htm) | A(Z)ᵢB prime for i = 0..k | k = 9 (Taneja) | [k = 10: 95773717(6)ᵢ1](results/882.md) |

#### Answers to open questions
| # | Open question | Status on the page | Our answer |
|---|---|---|---|
| [489](https://www.primepuzzles.net/puzzles/puzz_489.htm) | Other q with σ(σ(φ(q+1))) = q+1 | searched to 10⁶; Firoozbakht's family open beyond p = 1278 | [all solutions with q+1 = 2^a·3^b found (proof); family ends at p = 5; none new ≤ 10⁹](results/489.md) |
| [542](https://www.primepuzzles.net/puzzles/puzz_542.htm) | Other n with σ(n) = q^k (q odd prime, k > 1); primes p with σ(p^m) = q^n | no contributions | [only 81 and 400 up to 10¹⁸; Q2 needs n free of 2 and 3](results/542.md) |
| [870](https://www.primepuzzles.net/puzzles/puzz_870.htm) | A 4×4 prime nested magic square, grade 2 | "I couldn't find a single one" (Rivera) | [three: types R, L and LR, inner primes < 1000](results/870.md) |
| [1238](https://www.primepuzzles.net/puzzles/puzz_1238.htm) | Rules for the reversed period of 1/p; relation to the forward rules | no contributions | [works iff p \| 10^N + 10^J − 1; reverse lags = −forward lags (mod T); rules for 154 primes < 1000, none possible for 11 (proof)](results/1238.md) |

#### Conjectures
| # | Conjecture | Status on the page | Our result |
|---|---|---|---|
| [21](https://www.primepuzzles.net/conjectures/conj_021.htm) | Least N needing R consecutive primes (± sums), R = 11, 12 | table stops at R = 10 (Nash, 2000) | [R = 11: 113313, R = 12: 72598584](results/conj-21.md) |
| [39](https://www.primepuzzles.net/conjectures/conj_039.htm) | Is there n with d(n) = 9 (Jacobsthal gap)? | open (Firoozbakht) | [yes, first at 36465 = 3·5·11·13·17](results/conj-39.md) |
| [55](https://www.primepuzzles.net/conjectures/conj_055.htm) | Sand's variant: every n > 5 is p² + q + d, d ≤ 2 | "checked until 10⁶" | [false: 305, and infinitely many (mod 6)](results/conj-55.md) |
| [76](https://www.primepuzzles.net/conjectures/conj_076.htm) | Six inequalities for p_n^(1/n) | only via Firoozbakht's conjecture (to 4·10¹⁸) | [all six proved unconditionally](results/conj-76.md) |
| [80](https://www.primepuzzles.net/conjectures/conj_080.htm) | Why the numerator of ∏(p²+1)/(p−1) ends in 25 | "certainly true for large n", checked to n = 500,000 | [proved for all n (explicit bound + sieve to 8·10⁹)](results/conj-80.md) |
| [99](https://www.primepuzzles.net/conjectures/conj_099.htm) | Strong semiprimes / primes → ln 3 | no contributions | [proved (prime number theorem); data explained](results/conj-99.md) |
| [111](https://www.primepuzzles.net/conjectures/conj_111.htm) | Shielded primes have density 1/e² | data only, no argument | [explained: Flory's dimer theorem (1939)](results/conj-111.md) |
| [112](https://www.primepuzzles.net/conjectures/conj_112.htm) | Lucas–Lehmer type test for 2^k + 3 | no contributions | [proved; no pseudoprime for n ≤ 10⁴](results/conj-112.md) |
| [115](https://www.primepuzzles.net/conjectures/conj_115.htm) | Triangle test: passes iff prime or 9 | tested to 3000 | [false: the prime 3917 fails](results/conj-115.md) |


## Facts about primes
[FACTS.md](FACTS.md) collects 200 facts about primes, each with its formula and how we know it: a theorem, our own computation (one pass over all primes below 1.4·10¹¹), or a source in this archive. A few of them:
- Every prime p > 3 is 6k ± 1 and has p² ≡ 1 (mod 24). The gap after it is ≡ 2 (mod 6) only when p ≡ 5 (mod 6).
- Below 1.4·10¹¹ the most common gap is 6 and the largest is 468 (after 127976334671). Every gap after p ≥ 11 is below 1.32·√p, so Legendre's and Andrica's conjectures hold there.
- Consecutive primes avoid repeating their last digit: among the first 10⁸ primes, a prime ending in 1 is followed by another one ending in 1 only 18.5 % of the time.
- Chebyshev's bias: π(x; 4, 1) > π(x; 4, 3) first happens at x = 26861, and below 1.4·10¹¹ at only 0.44 % of the primes.
- Euler–Rodríguez chains: k consecutive primes of the form A + x(x + 1). The longest known has 14 (2003); a chain of 15 is still open.
- The longest run of consecutive primes on a quadratic below 1.4·10¹¹ is 27797667649 − x(x + 1), x = 0 … 11. Chains and their mirror images are equally common, as Hardy–Littlewood predicts.
- n² + n + 41 is provably the best polynomial of its kind (40 primes; Rabinowitsch–Heegner). Exact formulas for all primes exist (Willans, Mills, Fridman's 2.920050977…), but each of them encodes the primes it produces.
- Open: twin primes, Goldbach, Legendre, Cramér, primes n² + 1, a CPAP of 11 primes, and a formula for p_n that is not circular.

## Working notes
- [CLAUDE.md](CLAUDE.md): contains the instructions for the LLM: where the data is, how to crack a puzzle, the rules.
- [IDEAS.txt](IDEAS.txt): contains tricks to try next and open puzzles that look reachable.
- [FINDINGS.txt](FINDINGS.txt): contains the results that surprised us.
- [FACTS.md](FACTS.md): a long list of facts about primes (residues, counts, gaps, runs of consecutive primes, formulas, open problems), each with its formula and how we checked it.
- [writeup/](writeup/INDEX.md): contains our attempts on open puzzles, and an index of all puzzles and conjectures with their status.

## Models

The work started on **Claude Fable 5.1** at medium effort, which found about half of the results. When the credits ran out, **Claude Opus 5.5** at high effort did the rest.

Both got a short prompt and the tips in [CLAUDE.md](CLAUDE.md). When a model got stuck, two things helped. Telling it to think harder often made it find a new way. Asking why other solutions worked let it use the same idea on the new puzzle.

## Personal
I wanted to see what happens when you use AI on some of the hardest small problems in number theory. Primes look random, but they're not. So The Prime Puzzles turned out to be the perfect playground. Every record showed a bit more of the rules primes secretly follow, and how far clever reasoning via a LLM gets you before you need raw computing power. What's next? the prime distribution itself :)

Found any mistakes or have nice ideas on how to improve the reasoning? Open a PR or send a mail basvandorst [at] gmail [dot] com
