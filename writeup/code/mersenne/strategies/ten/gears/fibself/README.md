# Fibonacci fools itself (prime39.html, FACTS 6121)

The Fibonacci wheel of the prime gears page says n is prime when n divides F(n − ε), ε = (n/5). This folder recounts its mistakes (Fibonacci pseudoprimes) and shows why Fibonacci numbers are among them: on n = F(k) the test is the same as k | F(k) − ε, so every composite F(p) and F(2p) gets through (Duparc 1955, E. Lehmer 1964, Parberry 1970).

| file | what it does |
|---|---|
| `fibpsp.py N` | every Fibonacci pseudoprime below N (numpy, fast doubling), by n mod 5, and the pair with the base-2 Fermat test → `fibpsp_1e+08.txt` (1,460) |
| `fibfib.py` | the direct test on F(k) for k ≤ 1500, and the rule "k prime or twice a prime" with its three exceptions 572, 646, 754 |
| `indices.py` → `indices.json` | all k ≤ 10⁵ with k \| F(k) − ε: 9,589 primes, 5,131 doubles of primes, 59 others (terms of OEIS A279072) |
| `ref39.py` → `ref39.json` | reference for the page: the grid k ≤ 1500 by gmpy2 on F(k) itself, the mistakes below 10⁷, the examples (F(4181), 874 digits) |
| `gen39.py` | `prime39.src.html` + the json files → `prime39.html`; recomputes the recount table of section 5 (machine 6, the ring-5 list in `../tour/`) and asserts it |
| `compare39.py page39.json` → `compare39.txt` | the page's own results (run in headless Chrome with `t39.js`) against `ref39.json`: 0 mismatches |

The golden-wheel row on prime3.html now says that its 192 mistakes below 10⁷ are for n ≡ 2, 3 (mod 5), where the Fibonacci wheel makes the same 192 (`../gen_panel3.py`).
