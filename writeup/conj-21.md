# Conjecture 21: The next two entries of the R(N) table

> **In short:** R(N), the fewest primes in a signed sum of consecutive primes equal to N, had its table of first occurrences up to R = 10. The method: (1) walk all sign patterns of every block of up to 12 primes below 10⁸ in Gray-code order; (2) mark every sum; (3) confirm with a second program. This gives R = 11 first at N = 113313 and R = 12 at 72598584.

**Conjecture:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_021.htm)

**Status:** answered (Q5: two new table entries, exhaustive search to 10⁸) · Q4.2: no constant bound (heuristic, not a proof)

## The question

Write N as a sum of consecutive primes with signs, none of them larger than N, for example 14 = 13 − 11 + 7 + 5. R(N) is the smallest number of primes in such an expression. Question 5 asks to continue the table of the least N for each value of R. Question 4.2 asks whether R(N) stays below a constant Kr ("my personal bet is that … it is not greater than 20").

```
R(N) | least N | minimal expression
1    | 2       | 2
2    | 8       | 5+3
3    | 6       | 5+3−2
4    | 14      | 13−11+7+5
5    | 21      | 17+13−11+7−5
6    | 32      | 19+17−13+11−7+5
7    | 111     | 101+97−89+83−79−73+71
8    | 1024    | +239+241−251+257+263+269−271+277          (Nash)
9    | 3153    | −997+1009−1013−1019+1021+1031+1033+1039+1049   (Nash)
10   | 56924   | +14153−14159−14173+14177+14197+14207−14221+14243+14249+14251   (Nash)
11   | ?
12   | ?
```

## What the page had vs ours

| | On the page | Ours |
|---|---|---|
| R = 11 | ? (Nash searched to 65,521) | **113313** = +16103 −16111 +16127 −16139 +16141 +16183 +16187 +16189 +16193 +16217 +16223 |
| R = 12 | ? | **72598584** = +9074783 +9074789 +9074803 +9074827 +9074839 +9074843 +9074851 +9074861 +9074867 −9074873 +9074881 −9074887 |
| R = 13 | | none below 10⁸ |
| Q4.2: a constant bound Kr? | Rivera bets Kr ≤ 20, Nash expects Kr exists | **no constant bound**: R grows like 2·log₂(ln N) (heuristic) |

## How we found it

A block of r consecutive primes has 2^r sign patterns. Walking through them in Gray-code order changes one sign at a time, so each new sum costs a single addition. For every block of up to 12 primes below 10⁸, the program marks every sum that lands in [largest prime of the block, 10⁸]. The least unmarked N for each r gives the table. That is about 5·10¹⁰ additions, and it takes 52 seconds in C.

The search reproduces Nash's entries for R = 1 … 10 exactly. Two independent programs agree:
- A numpy version, which multiplies blocks by a matrix of sign patterns, matches every count per R up to 10⁷. Run on its own up to N = 72598584, it finds that every N from 5 to 72598583 needs at most 11 primes and that 72598584 needs 12. So minimality rests on two separately written exhaustive programs.
- A targeted search over all blocks near N/s, for each sign balance s, confirms that 113313 needs exactly 11 primes (44 expressions) and 72598584 needs exactly 12 (48 expressions).

## The math

**Parity.** A block of odd primes with signs has the parity of its length. Blocks that contain 2 have all their primes small, so for N > 238 every expression uses only odd primes. So **R(N) ≡ N (mod 2)**. That is why the table alternates: 111, 3153 and 113313 are odd, while 1024, 56924 and 72598584 are even.

**Why R has no constant bound (Q4.2).** Take a sign pattern with k plus signs and r − k minus signs, and let s = 2k − r. A block near x sums to about s·x, so only blocks near N/s can hit N. As the block slides one prime to the right, its sum moves by about s·ln N. Sums of r odd primes all have the parity of r, so each block hits N with chance about 2/(s·ln N). Summing over patterns:

```
expected number of r-prime expressions of N  ≈  Σ_s C(r, (r+s)/2) · 2/(s·ln N)  ≈  c · 2^r / ln N
```

The chance that N has no expression with r primes is then about exp(−c·2^r/ln N). Among the numbers up to X, the first that needs more than r primes appears when this chance is about 1/X:

```
c · 2^r / ln X ≈ ln X     →     2^r ≈ (ln X)² / c     →     record R(X) ≈ 2·log₂(ln X) + constant
```

This grows without bound, so no constant Kr exists. The data fit it. For the record values (R, X) = (10, 56924), (11, 113313), (12, 72598584), the constant c = (ln X)²/2^r comes out at 0.12, 0.066 and 0.080. With c ≈ 0.08, the model puts R = 13 near 10¹¹, R = 14 near 10¹⁵ and R = 20 near 10¹²⁶. So Rivera's "Kr ≤ 20" holds for every N anyone can compute, but not for all N.

Counts up to 10⁸ (how many N have each R):

```
R      1        2        3        4         5         6         7         8         9        10      11    12
count  5761455  3001132  6395506  11613020  16907897  23472989  18494661  11611888  2437649  300957  2835  9
```

## Validation

This script checks the two new expressions, confirms with a targeted search that neither number has a shorter one, and rebuilds the table for R = 1 … 11 up to 2·10⁵ with numpy (about 20 seconds). Minimality of 72598584 needs the full search to 10⁸: [code/conj-21.c](code/conj-21.c) (`./conj21 100000000 12`, about a minute).

```python
import itertools
import numpy as np
from sympy import isprime, nextprime, primerange

def consecutive(ps):
    return all(isprime(p) for p in ps) and all(nextprime(ps[i]) == ps[i+1] for i in range(len(ps) - 1))

new = {113313:   [+16103, -16111, +16127, -16139, +16141, +16183, +16187, +16189, +16193, +16217, +16223],
       72598584: [+9074783, +9074789, +9074803, +9074827, +9074839, +9074843, +9074851, +9074861,
                  +9074867, -9074873, +9074881, -9074887]}
for N, e in new.items():
    assert sum(e) == N and consecutive([abs(t) for t in e]) and max(map(abs, e)) <= N

def shortest(N, rmax):
    # least r with a signed block of r consecutive primes (all <= N) summing to N
    for r in range(1, rmax + 1):
        for s in range(2 - r % 2, r + 1, 2):              # s = (#plus) - (#minus), same parity as r
            k = (r + s) // 2
            ps = list(primerange(max(2, N // s - 500 * r), min(N, N // s + 500 * r) + 1))
            for i in range(len(ps) - r + 1):
                B = ps[i:i + r]
                for plus in itertools.combinations(range(r), k):
                    if sum(B[j] if j in plus else -B[j] for j in range(r)) == N:
                        return r
    return None
assert shortest(113313, 11) == 11 and shortest(72598584, 12) == 12

# The table R = 1..11 (least N), exhaustive up to 2*10^5.
M = 200000
P = np.array(list(primerange(2, M + 1)), dtype=np.int64)
best = np.full(M + 1, 99)
for r in range(1, 12):
    S = np.array(list(itertools.product((1, -1), repeat=r))).T
    n = len(P) - r + 1
    for lo in range(0, n, 2000):                          # chunks keep memory small
        B = np.stack([P[lo + j: min(n, lo + 2000) + j] for j in range(r)], axis=1)
        sums = B @ S
        v = sums[(sums >= B[:, -1:]) & (sums <= M)]
        best[v] = np.minimum(best[v], r)
least = [int(np.argmax(best == r)) for r in range(1, 12)]
assert least == [2, 8, 6, 14, 21, 32, 111, 1024, 3153, 56924, 113313]
print("table verified: R = 11 first at 113313; 72598584 needs exactly 12 primes")
```
