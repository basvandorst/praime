# Conjecture 110: Why Underwood's cubic test has no known pseudoprimes (attempt)

> **In short:** Paul Underwood conjectured that no composite n passes his cubic test for *any* choice of a and b. He checked every n up to 200,000, trying all (a, b) pairs. We show that for squarefree n the test splits into independent conditions, one per prime factor p of n. Each asks for a cubic mod p on whose three roots "raise to the power n" acts as a cyclic shift, with a^(n−1) ≡ 1 (mod p). Counting those cubics explains the empty search: for a random n, a prime factor p has a suitable cubic in only about 0.16/p of the cases or fewer, and the largest prime factor of n is constrained much more. So pseudoprimes should be extremely rare, and none is in reach. That supports the conjecture but does not prove it.

**Conjecture:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_110.htm)

**Status:** open · local criterion proved for squarefree n · heuristic explanation (no counterexample, no proof)

## The question

For odd, composite, non-square n, integers a, b (b ≠ 0) and f = x³ − ax − b irreducible, Underwood asks whether all of these can hold at once:

```
gcd(a, n) = 1,   a^n ≡ a (mod n),   jacobi(4a³ − 27b², n) = 1,   x²(B² + B + 1) ≡ a (mod n, f)   with B = x^(n−1)
```

The conjecture is that no composite n satisfies them for any a, b. He checked n ≤ 200,000 with all a, b, and asks for a counterexample (Q1), a proof (Q2) or a larger verification (Q3).

## What the page had vs ours

| | On the page | Ours |
|---|---|---|
| Structure of the test | derivation for prime n | for squarefree n: a counterexample exists iff **every prime factor** has a suitable local cubic |
| Why no pseudoprimes | open | local solutions are rare (table below); the largest prime factor is the bottleneck |
| Q1, Q2, Q3 | | not settled |

## How we found it

For prime n the test says: x^n is another root of f, which is what Frobenius does. For composite n, the ring (ℤ/n)[x]/(f) splits by the Chinese remainder theorem into one piece per prime factor p. In each piece the congruence only involves the three roots of f mod p. All conditions are local except irreducibility over ℚ, and that one can always be arranged by changing a and b by multiples of n.

## The math

**Local criterion (squarefree n).** Write the roots of f mod p as r₁, r₂, r₃ (in F_p or in F_{p³}). The identity f(Y) = (Y − r)(Y² + rY + r² − a) shows that x²(B² + B + 1) − a = x² + x·xⁿ + x²ⁿ − a vanishes in the piece at r exactly when rⁿ is one of the *other two* roots. So the test holds mod p iff:

- f has three distinct roots mod p (otherwise the Jacobi symbol is 0 or −1, or the condition fails at the root in F_p);
- y ↦ yⁿ sends each root to a different root;
- a ≢ 0 and a^(n−1) ≡ 1 (mod p).

The Jacobi condition is then automatic: a cubic with 3 roots in F_p, or irreducible over F_p, has square discriminant. In the irreducible case the condition reads n ≡ p or p² modulo the order of a root in F_{p³}*. In the split case it needs n to permute the roots like a 3-cycle inside F_p*.

**How often a prime factor can pass.** For small p we counted, over all residues n mod p³ − 1 (which is all that matters at p), how many admit a suitable cubic:

```
p                      5       7       11      13      17      19      23
fraction of n          3.2%    2.3%    0.30%   1.1%    0.16%   0.35%   0.03%
p × fraction           0.16    0.16    0.03    0.14    0.03    0.07    0.01
```

Primes p ≡ 2 (mod 3) are much harder. There a 3-cycle in F_p* needs an element of order divisible by 3 under the n-th power map, which is rarely available. For the largest prime factor q of n the residue n mod q³ − 1 is not random: for n = pq it is pq itself. The condition then forces the order of a root of f to divide p − 1 or q − p, which a random cubic almost never satisfies.

**Non-squarefree n.** If p² | n, then xⁿ is a Teichmüller element mod p², so the root itself must be Teichmüller mod p². The coefficient a must also satisfy a^(p−1) ≡ 1 (mod p²). These are Wieferich-type conditions, rare but not impossible, so prime powers need separate care in any extended search.

**What would settle it.** A counterexample needs every prime factor of n to pass at once. The numbers above make that very unlikely for any n a computer can reach, so the empty search up to 200,000 is expected. A proof is out of reach: it would make this a deterministic primality test of the same strength as Baillie–PSW, for which no proof is known either.

## Validation

This script checks the local criterion against Underwood's congruence, computed literally in F_p[x]/(f), on random cases and on forced positive cases. It then recomputes the fractions in the table. It runs in a few seconds.

```python
import random
from sympy import factorint

def mul(u, v, p, a, b):                          # multiply in F_p[x]/(x^3 - a x - b)
    c = [0] * 5
    for i in range(3):
        for j in range(3):
            c[i + j] += u[i] * v[j]
    for k in (4, 3):
        t = c[k]; c[k] = 0; c[k - 2] += a * t; c[k - 3] += b * t
    return [c[0] % p, c[1] % p, c[2] % p]

def powx(p, a, b, e):
    r, base = [1, 0, 0], [0, 1, 0]
    while e:
        if e & 1: r = mul(r, base, p, a, b)
        base = mul(base, base, p, a, b); e >>= 1
    return r

def order_x(p, a, b):                            # order of x in F_{p^3}^* (f irreducible)
    M = o = p**3 - 1
    for q, e in factorint(M).items():
        for _ in range(e):
            if powx(p, a, b, o // q) == [1, 0, 0]: o //= q
            else: break
    return o

def underwood(p, a, b, n):                       # x^2 (B^2 + B + 1) == a in F_p[x]/(f), B = x^(n-1)
    B = powx(p, a, b, n - 1); B2 = mul(B, B, p, a, b)
    return mul([0, 0, 1], [(B2[0] + B[0] + 1) % p, (B2[1] + B[1]) % p, (B2[2] + B[2]) % p], p, a, b) == [a % p, 0, 0]

def local(p, a, b, n):                           # our criterion: n-th power moves every root to another root
    R = [r for r in range(p) if (r**3 - a*r - b) % p == 0]
    if len(R) == 3:
        img = [pow(r, n, p) for r in R]
        return all(img[i] in R and img[i] != R[i] for i in range(3))
    if len(R) == 0:
        o = order_x(p, a, b)
        return (n - p) % o == 0 or (n - p*p) % o == 0
    return None                                  # one root: excluded by the Jacobi condition

random.seed(110)
for p in [5, 7, 11, 13, 17, 19, 23, 29, 31]:
    for _ in range(3000):
        a, b = random.randrange(1, p), random.randrange(1, p)
        if (4 * a**3 - 27 * b * b) % p == 0: continue
        n = random.choice([random.randrange(2, 10**6), p * random.randrange(1, 10**5)])
        L = local(p, a, b, n)
        if L is not None: assert underwood(p, a, b, n) == L
for p in [7, 13, 19]:                            # forced positives in the irreducible case
    for a in range(1, p):
        for b in range(1, p):
            if (4 * a**3 - 27 * b * b) % p == 0 or any((r**3 - a*r - b) % p == 0 for r in range(p)): continue
            o = order_x(p, a, b)
            assert underwood(p, a, b, p + o) and underwood(p, a, b, p * p + 2 * o)

def fraction(p):                                 # share of n mod p^3-1 with a suitable local cubic
    M, good = p**3 - 1, bytearray(p**3 - 1)
    for a in range(1, p):
        oa = next(d for d in range(1, p) if pow(a, d, p) == 1)
        for b in range(1, p):
            if (4 * a**3 - 27 * b * b) % p == 0: continue
            R = [r for r in range(p) if (r**3 - a*r - b) % p == 0]
            if len(R) == 3:
                for n in range(1, M, oa):        # n = 1 (mod ord a)
                    img = [pow(r, n, p) for r in R]
                    if all(img[i] in R and img[i] != R[i] for i in range(3)): good[n] = 1
            elif not R:
                o = order_x(p, a, b)
                for n in range(1, M, oa):
                    if (n - p) % o == 0 or (n - p*p) % o == 0: good[n] = 1
    return sum(good) / M

table = {p: round(fraction(p), 4) for p in [5, 7, 11, 13]}
assert table == {5: 0.0323, 7: 0.0234, 11: 0.003, 13: 0.0109}, table
print("local criterion agrees with Underwood's congruence; local fractions:", table)
```
