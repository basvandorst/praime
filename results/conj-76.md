# Conjecture 76: All six inequalities for p_n^(1/n) are true

> **In short:** Mahdi Meisami conjectured six inequalities for a_n = p_n^(1/n) (the n-th root of the n-th prime), such as a_n · a_m > a_(nm) and a_n + a_m > 1 + a_(n+m). On the page they were only known to follow from Firoozbakht's unproved conjecture, which has been checked up to 4·10¹⁸. We prove all six unconditionally. They all follow from one lemma, p_(a+b) ≤ p_a · p_b, which itself needs only Bertrand's postulate, a standard explicit bound on p_n, and a check of 148 small cases.

**Conjecture:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_076.htm)

**Status:** proved (all six statements, for all n, m ≥ 1)

## The question

With p_n the n-th prime and a_n = p_n^(1/n), prove or refute, for all natural n and m:

```
Conj1:   a_n · a_m > a_(nm)            Conj3:   a_n · a_m > a_(n+m)
Conj2:   a_n + a_m > a_(nm)            Conj4:   a_n + a_m > a_(n+m)
Stronger-Conj2:  a_n + a_m > 1 + a_(nm)      Stronger-Conj4:  a_n + a_m > 1 + a_(n+m)
```

## What the page had vs ours

| | On the page | Ours |
|---|---|---|
| All six | follow from Firoozbakht's conjecture (a_n decreasing), so true for n·m < 4·10¹⁸ (Kourbatov) | **proved for all n, m**, without Firoozbakht's conjecture |

## How we found it

Firoozbakht's conjecture says a_n is decreasing, which makes all six obvious. But a weaker fact is enough: a_(n+m) and a_(nm) are never larger than the larger of a_n and a_m. Since every a_k is above 1, multiplying by the other factor, or adding the other term minus 1, can only increase it. So everything comes down to comparing p_(n+m) with p_n · p_m, a much cruder statement than any gap conjecture.

## The math

**Basic facts.** a_k > 1 because p_k ≥ 2. Bertrand's postulate gives p_(k+1) < 2p_k, so p_k ≤ 2^k and a_k ≤ 2, with equality only for k = 1.

**Lemma.** p_(a+b) ≤ p_a · p_b for all a, b ≥ 1.

*Proof.* If a = 1, this is Bertrand: p_(b+1) < 2p_b. If 2 ≤ a ≤ b, then p_a ≥ 3 and

```
p_(a+b) ≤ p_(2b) < 3·p_b ≤ p_a · p_b.
```

It remains to show p_(2b) < 3p_b for b ≥ 2. For b ≥ 150 use Rosser's bound n·ln n < p_n (all n ≥ 1) and the Rosser–Schoenfeld bound p_n < n(ln n + ln ln n) (n ≥ 6):

```
p_(2b) < 2b(ln 2b + ln ln 2b)  ≤  3b·ln b  < 3·p_b
```

The middle step is equivalent to f(b) = ln b − 2 ln 2 − 2 ln ln 2b ≥ 0. Here f(150) = 0.142, and f is increasing for b ≥ 4 (its derivative is (1 − 2/ln 2b)/b). For 2 ≤ b < 150 the inequality p_(2b) < 3p_b is checked directly. ∎

**Consequence 1: a_(n+m) ≤ max(a_n, a_m).** By the lemma,

```
a_(n+m)^(n+m) = p_(n+m) ≤ p_n · p_m = a_n^n · a_m^m ≤ max(a_n, a_m)^(n+m).
```

**Consequence 2: a_(nm) ≤ min(a_n, a_m).** By induction on n, using the lemma with a = (n−1)m and b = m, we get p_(nm) ≤ p_(nm−m) · p_m ≤ … ≤ p_m^n. So a_(nm)^(nm) ≤ a_m^(nm), i.e. a_(nm) ≤ a_m. By symmetry a_(nm) ≤ a_n as well.

**The six statements.** Let M = max(a_n, a_m). Since both a_n and a_m exceed 1:

```
a_n · a_m > M           (multiply M by the other factor, which is > 1)
a_n + a_m − 1 > M       (add the other term minus 1, which is > 0)
a_n + a_m > 2 ≥ a_k     for every k
```

- Conj1 and Conj3: a_n · a_m > M ≥ a_(nm), a_(n+m) by the two consequences.
- Stronger-Conj2 and Stronger-Conj4: a_n + a_m > 1 + M ≥ 1 + a_(nm), 1 + a_(n+m).
- Conj2 and Conj4 follow from the stronger forms, or directly from a_n + a_m > 2 ≥ a_k.

All inequalities are strict. ∎

The proof gives more than was asked: a_(nm) ≤ min(a_n, a_m) and a_(n+m) ≤ max(a_n, a_m) for all n, m. Firoozbakht's conjecture would put a_(n+m) below the minimum too; that stronger statement remains open.

## Verify it yourself

This script checks the finite part of the proof (p_(2b) < 3p_b for 2 ≤ b < 150) and that f is positive from b = 150 on. As sanity checks it also tests the lemma for a, b < 600, the explicit bounds on the whole range used, and all six inequalities for n, m < 400. It runs in about ten seconds.

```python
import math
from sympy import primerange

P = [0] + list(primerange(2, 2 * 10**7))          # P[k] = k-th prime
K = len(P) - 1

# The finite part of the proof, plus a sanity check far beyond it.
assert all(P[2*b] < 3 * P[b] for b in range(2, 150))
assert all(P[2*b] < 3 * P[b] for b in range(2, K // 2))

# The explicit bounds used for b >= 150 (Rosser 1939; Rosser-Schoenfeld 1962), on the range we have.
assert all(n * math.log(n) < P[n] for n in range(1, K + 1))
assert all(P[n] < n * (math.log(n) + math.log(math.log(n))) for n in range(6, K + 1))
f = lambda b: math.log(b) - 2 * math.log(2) - 2 * math.log(math.log(2 * b))
assert f(150) > 0 and all(f(b + 1) > f(b) for b in range(4, 10**5))

# The lemma p_(a+b) <= p_a * p_b, and the six inequalities.
assert all(P[x + y] <= P[x] * P[y] for x in range(1, 600) for y in range(1, 600))
a = [0] + [P[k] ** (1.0 / k) for k in range(1, 400 * 400)]
for n in range(1, 400):
    for m in range(1, 400):
        assert a[n] * a[m] > a[n*m] and a[n] + a[m] > 1 + a[n*m]
        assert a[n] * a[m] > a[n+m] and a[n] + a[m] > 1 + a[n+m]
print("lemma and all six inequalities verified on the tested ranges; finite part of the proof checked")
```
