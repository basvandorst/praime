# Problem 79: 223 numbers with three prime factors whose reciprocals sum to 1

> **In short:** The record was 688 numbers with three prime factors (squares like 12 = 2²·3 allowed) whose reciprocals sum to 1. The construction: (1) the sum is an integer exactly when, for every prime p, p²·(the sum of the terms divisible by p) is 0 mod p², and cubes p³ can never occur; (2) start from the smallest such numbers and repair the primes from the largest down; (3) close 2, 3, 5, 7 exactly. This gives 223 terms; at least 189 are needed.

**Puzzle:** [primepuzzles.net](https://www.primepuzzles.net/problems/prob_079.htm)

**Status:** beaten · exact (rational arithmetic)

## The question

Write 1 as a sum of reciprocals of distinct numbers that each have exactly three prime factors, counted with multiplicity (12 = 2·2·3 counts). The first solution on the page had 914 terms.

**Q.** Can you get a solution with fewer than 914 terms?

## Record vs ours

| | Terms | Largest |
|---|---|---|
| Record (Dec 2020) | 688 (266 not squarefree) | — |
| **Ours** | **223** (91 not squarefree) | 977309 |
| Lower bound | 189 | — |

## How we found it

This is the method of [problem 76](prob-76.md) and [puzzle 1020](1020.md), with one more layer. Once squares are allowed, the condition at a prime p lives modulo p², not modulo p. It splits into two steps: first the terms p²·a fix the residue mod p, then the terms p·a·b fix the next digit.

## The math

**Step 1 (proved). The criterion.** For a prime p, let S_p be the sum of 1/t over the terms divisible by p; every other term has no p in its denominator. A term p²q contributes q⁻¹ to p²·S_p, and a term p·m with p ∤ m contributes p·m⁻¹, both taken mod p². So the total is an integer iff

```
R_p := Σ_(p²q ∈ S) q⁻¹ + p · Σ_(pm ∈ S, p ∤ m) m⁻¹ ≡ 0 (mod p²)    for every prime p
```

A cube p³ would add 1/p³, which nothing else can cancel (it is the only term with p³), so cubes never occur.

**Step 2 (proved). Lower bound 189.** The reciprocals of the 189 smallest numbers with three prime factors, cubes excluded (12, 18, 20, 28, 30, …, 782), are the first to reach 1.

**Step 3 (construction).**

1. *Bulk:* all such numbers t ≤ X₀ = 850.
2. *Repair* every prime v > 7 from the largest down. If R_v ≢ 0 (mod v), toggle a subset of the terms v²·a (a < v prime) to fix the residue mod v. Then fix the residue mod v² with one or two terms v·a·b (a ≤ b < v). The best option by Δ(terms) − X₀·Δ(sum) is taken.
3. *Close* 2, 3, 5, 7: their residues R₂ mod 4, R₃ mod 9, R₅ mod 25 and R₇ mod 49 form one of 44100 states. The moves are the 16 terms built only from 2, 3, 5, 7, and gadgets: sets of terms v·m (m a product of two of 2, 3, 5, 7) whose contributions cancel mod v. A dynamic programme finds the cheapest combination; its cost also counts the change in sum, so that it lands on 1 rather than 0.

With seed 7 the result has 223 terms on 35 primes, 91 of them not squarefree, and the exact sum is 1. Over the settings tried (X₀ = 700 … 1000, 10 seeds each) the best results were 223, 226, 227, 228 and 229.

**What this does not settle.** The minimum lies between 189 and 223.

## Validation

[code/prob79_verify.py](code/prob79_verify.py) reads the 223 numbers from [data/prob79_solution.txt](data/prob79_solution.txt). It checks that they are distinct, each with exactly three prime factors and no cubes, that their reciprocals sum to exactly 1, that R_p ≡ 0 (mod p²) for every prime, and the lower bound 189. [code/prob79_build.py](code/prob79_build.py) rebuilds the set in a few seconds. Run after `pip install sympy numpy`:

```
python3 prob79_verify.py
python3 prob79_build.py 850 7 1.0 out.txt
```
