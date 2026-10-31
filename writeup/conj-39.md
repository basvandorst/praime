# Conjecture 39: Yes, d(n) = 9 occurs, first at n = 36465

> **In short:** On the page of Ingala's conjecture, d(n) is the largest gap between consecutive integers coprime to n (Jacobsthal's function). Farideh Firoozbakht proposed a table in which d(n) depends only on how many primes divide n and whether 2 or 3 is among them. The table never produces 9, so she asked whether d(n) = 9 is possible. It is. The smallest n is **36465 = 3·5·11·13·17**: the integers 829 and 838 are coprime to it, and the eight numbers in between are not. Her table cannot be right in general, because d(n) depends on which primes divide n. For example, d(3·5·7·11) = 7 but d(3·7·11·13) = 6.

**Conjecture:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_039.htm)

**Status:** answered (Firoozbakht's question) · smallest n (exhaustive over n ≤ 40000) · her table disproved

## The question

S_n is the set of integers coprime to n, and d(n) is the largest difference between consecutive elements of S_n. Ingala's conjecture (d(n) < 2p for the largest prime factor p) was refuted on the page by Knjzek, Calligy and Andersen. Firoozbakht then wrote down a table of d(n) by the "shape" of n. With n = p^k·q^l·… and 3 < p < q < …, her table reads:

```
p                           1        2^a·p·q                   6
p^k (k>1)                   2        3^a·p·q·r                 7
p·q                         3        p·q·r·s·t                 7
2^a·p                       4        2^a·p·q·r                 8
p·q·r                       4        2^a·3^b·p·q               10
3^a·p·q                     5        2^a·p·q·r·s               10
p·q·r·s                     5        3^a·p·q·r·s               11
```

She asked: "Does there exist a number n such that d(n) = 9?"

## What the page had vs ours

| | On the page | Ours |
|---|---|---|
| Is d(n) = 9 possible? | open (Firoozbakht) | **yes; smallest n = 36465 = 3·5·11·13·17** |
| The shape table | proposed as "seems true" | **false**: d(3·7·11·13) = 6, not 7; d(3·5·11·13·17) = 9, not 11 |

## How we found it

d(n) only depends on the set of primes dividing n, and the gaps repeat with period n. So for each squarefree n we scanned one period and took the largest gap. The first n with each value of d:

```
d     2   3   4   5    6   7     8    10   11     12    14    22     9
n     2   15  6   105  30  1155  770  210  15015  6006  2310  30030  36465
```

Every value from 2 to 12 shows up early, except 9. It first appears at 36465.

## The numbers

```
n = 36465 = 3 · 5 · 11 · 13 · 17
829                       coprime to n (prime)
830 = 2·5·83     831 = 3·277    832 = 2⁶·13    833 = 7²·17
834 = 2·3·139    835 = 5·167    836 = 2²·11·19  837 = 3³·31
838 = 2·419               coprime to n
```

So 829 and 838 are consecutive elements of S_n, 9 apart, and no gap in S_n is larger.

Why the table fails: d(n) is decided by how the residue classes 0 mod p can be lined up to cover a run of integers, and that depends on the actual primes. With 3·5·7·11, the multiples of 3, 5, 7 and 11 can cover 6 consecutive integers (d = 7). With 3·7·11·13 they can only cover 5 (d = 6). Both have the same shape in her table.

Two smaller slips on the same page, both read off the page's own table of d(p#), whose values agree with OEIS A048670:
- Dan Lima's "d(p#(k)) < 2p(k−1) for k = 9, 12" is reversed. From the page's own table, d(23#) = 40 > 38 = 2·19 and d(37#) = 66 > 62 = 2·31. So his stronger conjecture d(p#(k)) ≤ 2p(k−1) already fails at 23#.
- n = 36465 also violates the table's claimed value 11 for the shape 3·p·q·r·s.

## Validation

This script computes d(n) by scanning one period, finds the first n with each value up to 22 among n ≤ 40000 (squarefree n suffice, since d(n) only depends on the primes dividing n), and checks the claims above. It runs in a few seconds.

```python
from math import gcd
import numpy as np
from sympy import factorint

def d(n):                                    # largest gap between consecutive integers coprime to n
    hit = np.zeros(n + 2, dtype=bool)        # one period 1..n+1 is enough
    for p in factorint(n):
        hit[::p] = True
    xs = np.nonzero(~hit[1:])[0] + 1
    return int(np.diff(xs).max()) if len(xs) > 1 else 1

first = {}
for n in range(2, 40001):
    if max(factorint(n).values()) == 1:      # squarefree n suffice: d depends only on the primes
        first.setdefault(d(n), n)
assert first[9] == 36465 and all(first[k] < 36465 for k in (2, 3, 4, 5, 6, 7, 8, 10, 11, 12))
assert 36465 == 3 * 5 * 11 * 13 * 17
assert gcd(829, 36465) == gcd(838, 36465) == 1 and all(gcd(x, 36465) > 1 for x in range(830, 838))
assert d(3 * 5 * 7 * 11) == 7 and d(3 * 7 * 11 * 13) == 6      # same "shape", different d
print("d(n) = 9 first at n = 36465 = 3*5*11*13*17")
```
