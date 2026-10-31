# Conjecture 111: Why the shielded primes have density 1/e²

> **In short:** It was conjectured that the density of shielded primes, survivors of a sieve that removes consecutive pairs in order of gap size, tends to 1/e². The explanation: (1) the sieve is Flory's random sequential adsorption of dimers; (2) for random gaps a site stays empty with probability exactly e⁻²; (3) the primes approach it up to 10¹⁵. This explains 1/e², without a proof for the primes.

**Conjecture:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_111.htm)

**Status:** explained (Q3) · exact theorem for random gaps (Flory), checked numerically for primes up to 10¹⁵; not a proof for the primes

## The question

Take the first k primes. Go through the gap sizes 1, 2, 4, 6, … in increasing order. For each size, scan left to right and remove every pair of *originally consecutive* primes with that gap, provided both are still present. Primes whose neighbours were removed earlier are "shielded" and survive forever. The density of survivors is 0.1358 for k = 10⁷ and 0.13556 for k ≈ 10¹⁰ (Gennady Gsev), slowly decreasing. The conjecture is that it tends to 1/e² = 0.135335…; Q3 asks for a proof. Emmanuel Vantieghem asked for "some reasonable arguments" for the value e⁻².

## What the page had vs ours

| | On the page | Ours |
|---|---|---|
| Q1, Q2 | counts to k = 10¹¹ (Gsev, Cavegn) | first six counts reproduced exactly |
| Q3: why 1/e² | no argument | **Flory's theorem**: for random gaps the density is exactly e⁻² (proof below) |
| Is the limit really e⁻²? | cumulative density 0.13552 at 10¹¹, still falling | local density near 10¹² and 10¹⁵ equals e⁻² within noise |

## How we found it

Each prime is a site on a line, and each gap is an edge between neighbouring sites. The sieve takes edges in order of length. It removes an edge's two sites only if both are still there, and a removed site can never be used again. That is a *greedy matching*: every site is used at most once, and shorter edges go first.

If the edge lengths were random numbers with no ties, "in order of length" would just mean "in random order". Dimers dropped in random order on a line, sticking only where both sites are free, is one of the oldest models in statistical physics: Flory's 1939 model of reactions between neighbouring groups on a polymer chain. Its answer is famous: when nothing more fits, a fraction exactly e⁻² of the sites is left empty.

## The math

**Flory's theorem.** Give every edge an independent arrival time, uniform in [0, 1], and accept an edge if both its sites are still free when it arrives. Then every site ends up empty with probability e⁻².

*Proof.* First take a half-line of sites 1, 2, 3, … Let q(t) be the chance that site 1 is still free at time t. Site 1 can only be taken by edge (1, 2). That edge, arriving at time s < t, is accepted exactly when site 2 is still free at time s in the half-line 2, 3, …, which has the same law. So

```
q(t) = 1 − ∫₀ᵗ q(s) ds     ⇒     q(t) = e^(−t)
```

Now take a site 0 with edges to −1 and to 1, arriving at times a and b. Say a < b. The edge to −1 is refused only if site −1 is already taken (chance 1 − e^(−a)). Then the edge to 1 is refused only if site 1 is already taken (chance 1 − e^(−b)). The two half-lines are independent, and the case b < a is the same. So

```
P(site 0 stays empty) = E[(1 − e^(−a))(1 − e^(−b))] = (∫₀¹ (1 − e^(−t)) dt)² = (1/e)² = e^(−2)   ∎
```

**From random times to prime gaps.** The sieve uses the gaps as arrival times. Two things differ from Flory's model:

- **Ties.** Equal neighbouring gaps are resolved left to right. An isolated tie does not change the density on average, because swapping the order of the two edges gives the mirror image of an equally likely configuration. Only clusters of ties matter. Numerically, resolving ties right to left instead changes the survivors at k = 10⁷ from 1,358,010 to 1,358,422, only 0.03%.
- **Correlations.** Consecutive prime gaps are not independent. For example, a gap of 2 or 4 is never followed by a gap of the same size (apart from 3, 5, 7), because one of p, p + g, p + 2g is divisible by 3. In the Cramér/Gallagher picture these effects fade as the gaps grow like ln p.

So the limit should be exactly e⁻², with corrections that shrink as the primes grow. That is what the numbers show.

## The numbers

Cumulative counts (first k primes), as on the page:

```
k        survivors   density
10²      14          0.14
10³      134         0.134
10⁴      1374        0.1374
10⁵      13494       0.13494
10⁶      135530      0.135530
10⁷      1358010     0.1358010
```

Random edge weights instead of prime gaps (pure Flory model, k = 10⁷, six runs): mean density 0.135303, run-to-run spread ±0.000065. So e⁻² = 0.135335 is confirmed, and ±0.00007 is the noise level for windows of 10⁷ sites.

The cumulative density hides the trend, because the small primes dominate it. The *local* density, measured on 10⁷ consecutive primes starting at X, is much closer to e⁻²:

```
X        local density   excess over e⁻²     same for a random mod-6 model
10⁸      0.1359234       +0.000588            −0.000336
10¹⁰     0.1355724       +0.000237
10¹²     0.1353920       +0.000057            −0.000140
10¹⁵     0.1353560       +0.000021
```

From 10¹² on, the excess is below the noise level. A random-prime model (n ≡ ±1 mod 6 taken with probability 3/ln n) approaches e⁻² from *below* while the real primes approach it from *above*. The leftover deviation therefore comes from fine gap statistics, and it vanishes as the primes grow.

**What remains open.** A proof for the actual primes would need the fine statistics of consecutive prime gaps (a uniform version of the Hardy–Littlewood conjectures), which is far beyond current methods. What we give is the exact theorem behind the constant, and evidence that the primes follow it.

## Validation

This script reproduces the counts for k = 10² … 10⁶, checks Flory's integral, and simulates the random model on 4·10⁶ sites. It runs in a few seconds. The C program for the windows near 10⁸ … 10¹⁵ is [code/conj-111.c](code/conj-111.c).

```python
import math, random
import numpy as np
from sympy import prime

def survivors(weights):
    # greedy matching on a path: take edges in order (weight, position), keep an edge if both ends are free
    order = np.lexsort((np.arange(len(weights)), weights))
    dead = bytearray(len(weights) + 1)
    for i in order:
        if not dead[i] and not dead[i + 1]:
            dead[i] = dead[i + 1] = 1
    return dead.count(0)

# 1. The page's counts, from the first k primes.
N = prime(10**6)
sieve = np.ones(N + 1, dtype=bool); sieve[:2] = False
for i in range(2, int(N**0.5) + 1):
    if sieve[i]: sieve[i*i::i] = False
P = np.nonzero(sieve)[0]
for k, expected in [(100, 14), (1000, 134), (10**4, 1374), (10**5, 13494), (10**6, 135530)]:
    assert survivors(np.diff(P[:k])) == expected, k

# 2. Flory's constant: (integral of 1 - e^-t over [0,1])^2 = e^-2.
integral = sum(1 - math.exp(-(j + 0.5) / 10**5) for j in range(10**5)) / 10**5
assert abs(integral**2 - math.exp(-2)) < 1e-9

# 3. Random arrival times: the empty fraction is e^-2 (noise about 1e-4 for 4*10^6 sites).
random.seed(1)
dens = survivors(np.array([random.random() for _ in range(4*10**6 - 1)])) / (4*10**6)
assert abs(dens - math.exp(-2)) < 5e-4
print("counts reproduced; random model density %.5f vs e^-2 = %.5f" % (dens, math.exp(-2)))
```
