# Praime
This repository is an experiment in how far AI can go with hard math: puzzles, logical reasoning, and problems nobody has cracked yet. The test case is Carlos Rivera's Prime Puzzles site, which has collected hundreds of prime challenges since 1999, many with records that stood for decades. No racks of computers grinding through billions of candidates, just careful reasoning about why each old record worked and a trick that made the hard part happen by itself. 

**The result: 39 records broken and 25 open questions answered.**

## results
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
| [1080](https://www.primepuzzles.net/puzzles/puzz_1080.htm) | Third term a(3) | none known | [1.9·10³²](results/1080.md) |
| [1145](https://www.primepuzzles.net/puzzles/puzz_1145.htm) | Covering set behind the listed k | no contributions | [found; listed k is its minimum](results/1145.md) |
| [1204](https://www.primepuzzles.net/puzzles/puzz_1204.htm) | Explain the 51.6% | zero answers | [Hardy–Littlewood explanation](results/1204.md) |


## LLM Instructions
### Puzzle reasoning
- **Ask for reasoning before code.** "Explain why the old record works, then propose a construction" gets far better ideas than "find a better one".
- **Forbid brute force.** Say explicitly: computation is only for checking a construction. Otherwise the model drifts into giant searches.
- **Demand several angles.** "Try at least three different approaches and write down why each one failed" stops it from giving up after the first dead end.
- **Give it the toolbox.** List tricks that worked before (fix a remainder, reuse record primes, stack building blocks). It applies them to new puzzles surprisingly well.
- **Split the work.** Separate agents for page ranges (e.g. puzzles 1–66, 67–133, …) scout in parallel. One agent per hard puzzle digs deeper.
- **Never trust a claim without a re-check.** Re-run every claimed beat in a separate script yourself. Several "records" didn't survive.
- **Make it read the whole page.** Models often read only the top or a truncated page. The fix: decode it properly (these pages are Windows-1252) and read to the last line. An audit like this withdrew six false claims.
- **Watch for overclaiming.** "Proved", "minimal" and "first" need evidence. Ask "what exactly was checked, and what is assumed?"
- **Set hard limits.** Per-command CPU caps, no background jobs, low priority. LLM agents will happily start a 5-hour run.
- **Keep a running log.** One results file and one submission file, updated after every check. That makes the final write-up easy.
### Code writing
- **Read the full puzzle page.** Download it, read it to the bottom, and note the best answer and who found it. Later contributions often beat the one at the top.
- **Factor the record numbers.** Run them through `sympy.factorint` and check their remainders mod 3, 4, 5, 8 and 40. Patterns show up fast.
- **Turn the condition into a remainder.** Ask whether "even period", "unlucky" or "no digit 6" depends only on p mod m. If it does, build a formula that is always ≡ that value mod m.
- **Search inside that class only.** Fix the remainder in the coefficients, then search just those polynomials. Far fewer candidates, much higher hit rate.
- **Plug in known record objects.** Before computing anything, try the Prime Pages records, PrimeGrid's AP-27 and the Dress–Landreau polynomials. They often win outright.
- **Sieve before testing.** Throw out candidates with small factors first using a `numpy` boolean sieve, then run `gmpy2.is_prime` only on the survivors.
- **Vectorise exhaustive checks.** When proving "this is the minimum" or "this can't exist", use `numpy` arrays instead of Python loops. Whole ranges get checked in seconds.
- **Keep runs short.** Cap every job (`timeout`, `nice`) so one bad idea can't eat a day of CPU.
- **Verify with a second, separate script.** Recompute the result from scratch with `sympy.isprime` and plain integer arithmetic before calling it a record.
- **Mark probable primes "PRP".** Huge numbers that pass the tests but aren't formally proven get the label, as the site expects.


## Personal
I wanted to see what happens when you pair AI and point them at some of the hardest little problems in number theory. Primes look random, but they're not. So The Prime Puzzles turned out to be the perfect playground. Every record showed a bit more of the rules primes secretly follow, and how far clever reasoning via a LLM gets you before you need raw computing power. What's next? the prime distribution itself :)
