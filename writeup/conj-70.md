# Conjecture 70: Every prime p > 7 has two consecutive primitive roots (it is a theorem)

> **In short:** Emmanuel Vantieghem conjectured that for every prime p > 7 there is an n with n and n + 1 both primitive roots mod p. The page knew it only for p > 10¹⁹ (Szalay) and up to the 20-millionth prime (computation). It is in fact a theorem. In 2015 Cohen, Oliveira e Silva and Trudgian proved that every finite field with q > 169 elements (q odd) has even *three* consecutive primitive elements. The primes 7 < p < 170 are checked directly. So the conjecture holds for every prime p > 7.

**Conjecture:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_070.htm)

**Status:** proved (by a published theorem for p > 169, plus a direct check of the 35 primes 7 < p < 170)

## The question

Q1: prove that every prime p > 7 has two consecutive primitive roots, or find a counterexample. (For p = 7 the primitive roots are 3 and 5, so 7 is a genuine exception.)

## What the page had vs ours

| | On the page | Ours |
|---|---|---|
| p > 10¹⁹ | true (Szalay 1975) | |
| p < 20,000,000th prime | checked by computer (Vantieghem) | |
| The gap in between | "I do not believe that this can be done in reasonable time" | **closed**: a 2015 theorem covers all p > 169 |

## How we found it

The question "are there consecutive primitive roots?" is a well-studied problem about primitive elements in finite fields. S. D. Cohen worked on it from 1985 on ([Proc. Amer. Math. Soc. 93 (1985), 189–197](https://www.ams.org/proc/1985-093-02/S0002-9939-1985-0770516-9/), and a sequel in vol. 94). The strongest result we found settles much more than the conjecture asks:

> For q an odd prime power with q > 169, there are always three consecutive primitive elements in the finite field F_q. Indeed, there are precisely eleven values of q ≤ 169 for which this is false.
> (S. D. Cohen, T. Oliveira e Silva, T. Trudgian, *On consecutive primitive elements in a finite field*, Bull. London Math. Soc. 47 (2015), 418–426; [arXiv:1410.6210](https://arxiv.org/abs/1410.6210))

## The math

For a prime p, F_p is the integers mod p, so "three consecutive primitive elements" means three consecutive integers n, n + 1, n + 2 that are all primitive roots mod p. In particular n and n + 1 are. So every prime p > 169 satisfies the conjecture.

The remaining primes are 11, 13, …, 167. For each of them a short search finds a consecutive pair; the first pair never starts later than n = 34. With p = 7 excluded, the conjecture is true for every prime p > 7. ∎

Two remarks:
- Szalay's bound (1975) needed p > 10¹⁹. Later work brought the bound down far enough for a direct check of the remaining cases to be trivial.
- The same 2015 paper conjectures how large q must be to guarantee four, five, … consecutive primitive roots. Those remain open.

## Verify it yourself

This script checks the 35 primes below 170 (and that 7 fails). As a sanity check it also confirms three consecutive primitive roots for every prime 170 < p < 20000, consistent with the theorem. It runs in a few seconds.

```python
from sympy import primerange, is_primitive_root

def first_run(p, length):
    # least n such that n, n+1, ..., n+length-1 are all primitive roots mod p
    for n in range(1, p - length + 1):
        if all(is_primitive_root(n + j, p) for j in range(length)):
            return n
    return None

assert first_run(7, 2) is None                                     # 7 is a real exception (roots 3, 5)
small = {p: first_run(p, 2) for p in primerange(11, 170)}
assert len(small) == 35 and all(n is not None for n in small.values()) and max(small.values()) == 34
assert all(first_run(p, 3) is not None for p in primerange(170, 20000))   # the 2015 theorem's range, sampled
print("consecutive primitive roots exist for every prime 7 < p < 170; theorem covers p > 169")
```
