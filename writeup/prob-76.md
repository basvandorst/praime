# Problem 76: 711 odd semiprimes whose reciprocals sum to 1

> **In short:** The record was 2587 distinct odd semiprimes with reciprocals summing to 1. The construction: (1) the sum is an integer exactly when, for every prime p, the partners q of p satisfy Σ q⁻¹ ≡ 0 (mod p); (2) start from the smallest semiprimes and repair this condition prime by prime from the largest down; (3) close the last primes 3…13 exactly. This gives 711 terms; at least 551 are needed.

**Puzzle:** [primepuzzles.net](https://www.primepuzzles.net/problems/prob_076.htm)

**Status:** beaten · exact (rational arithmetic, no primality question)

## The question

Find ODD semiprimes whose reciprocals sum up to 1, or prove that this is impossible. "Send your best result."

## Record vs ours

| | Terms | Largest term |
|---|---|---|
| Current record (E. Vantieghem, Jan 2021) | 2587 | — |
| **Ours** | **711** | 265189 |
| Lower bound (proved) | 551 | — |

Earlier results on the page: 2775 terms (Nov 2020), then 2767 (Dec 2020).

## How we found it

The page's method was a greedy loop that adds and removes terms until the sum is exactly 1, which "needs more luck than intelligence". We looked at the denominators instead. View each semiprime pq as an edge between the primes p and q. The prime p leaves the denominator of the sum exactly when the reciprocals of its partners cancel mod p. These are local conditions, one per prime, and they can be fixed one prime at a time. When all hold, the sum is an integer, and a sum near 1 is then exactly 1. The criterion and the fact that squares never occur are known from [puzzle 897](https://www.primepuzzles.net/puzzles/puzz_897.htm), where they are used for semiprimes of any parity (sum 1 with 47 terms, sums 2 and 3) and for numbers with three prime factors. What is new here is the odd case of problem 76 and the construction: top-down repair with upward gadgets and an exact closing step.

## The math

Let S be a finite set of odd semiprimes and Σ = Σ_(s∈S) 1/s.

**Step 1 (proved; also on puzzle 897). No square p² can be used.** Write the terms that contain p as 1/p² (if p² ∈ S) and (1/p)·A with A = Σ_(pq∈S, q≠p) 1/q, whose denominator is coprime to p. Their sum is (p·A + 1)/p², and p·A + 1 ≢ 0 (mod p), so Σ would keep p² in its denominator. Hence every solution uses only squarefree semiprimes pq with p < q.

**Step 2 (proved; also on puzzle 897). The local criterion.** For a prime p, call the q with pq ∈ S the partners of p. The terms containing p add up to (1/p)·Σ_q 1/q, and all other terms have no p in the denominator. So p divides no denominator of Σ iff

```
Σ_(q partner of p) q⁻¹ ≡ 0 (mod p)        ("p is balanced")
```

If every prime is balanced, Σ is an integer, and 0 < Σ < 2 then means Σ = 1. Consequences: every prime of a solution has at least 2 partners (a single q⁻¹ is never 0 mod p), and the largest prime has at least 3 (with two smaller partners q₁ + q₂ ≡ 0 (mod p) and q₁ + q₂ < 2p would force q₁ + q₂ = p, which is even).

**Step 3 (proved). Lower bound 551.** The reciprocals of the 551 smallest odd squarefree semiprimes (15, 21, 33, …, 2827) are the first to reach 1, so no solution has fewer terms.

**Step 4 (construction).** Primes are vertices and terms are edges; F = {3, 5, 7, 11, 13}.

1. *Bulk:* all odd squarefree semiprimes pq ≤ X₀ = 3950.
2. *Repair, from the largest prime v > 13 down:* if v is not balanced, toggle edges from v to smaller primes. For v ≤ 47 all subsets of these edges are tried; for larger v, combinations of one to three edges. An "upward gadget" is also allowed: an edge (v, u) to a larger prime u, together with edges (u, f), f ∈ F, chosen so that u stays balanced. The option with the smallest Δ(edges) − 0.8·X₀·Δ(sum) is taken. Smaller primes are repaired later, so each step leaves every larger prime balanced.
3. *Closing:* only the primes of F can still be unbalanced. Their residues form one of 3·5·7·11·13 = 15015 states. The available moves are the toggles of the 10 edges inside F, and gadgets (v, Q): for a prime v and Q ⊆ F with Σ_(q∈Q) ±q⁻¹ ≡ 0 (mod v), toggle all edges (v, q), q ∈ Q, which keeps v balanced. Example: 71 divides the numerator of 1/3 + 1/5 + 1/7 = 71/105. A dynamic programme over the 15015 states, taking at most one move per v, finds the combination with the fewest edges that balances all of F.

After step 3 every prime is balanced, so Σ is an integer, and it lies between 0 and 2, so Σ = 1. With seed 42 the result has 711 terms on 106 primes, every prime has at least 3 partners, and the largest term is 265189 = 359·739. The run is deterministic for a given (X₀, seed, factor 0.8). Over 9 parameter settings × 60 seeds the best results were 711, 711, 713 terms.

**What this does not settle.** The minimum lies between 551 and 711. Reaching the lower bound would need almost all of the 551 smallest semiprimes, but many of them have a large prime with only one or two partners, which step 2 forbids. So the true minimum is probably well above 551.

## Validation

[code/prob76_verify.py](code/prob76_verify.py) reads the 711 numbers from [data/prob76_solution.txt](data/prob76_solution.txt). It checks that they are distinct odd squarefree semiprimes, that their reciprocals sum to exactly 1 (exact fractions), that every prime is balanced with at least 2 partners, and the lower bound 551. [code/prob76_build.py](code/prob76_build.py) rebuilds the set from scratch in about 2 seconds. Run after `pip install sympy`:

```
python3 prob76_verify.py
python3 prob76_build.py 3950 42 0.8 out.txt
```
