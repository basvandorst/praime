# Conjecture 112: A Lucas–Lehmer type test for primes 2^k + 3

> **In short:** Davide Rotondo noticed that when N = 2^(n−1) + 3 is prime, the Lucas–Lehmer number s(n−2) leaves a fixed remainder mod N: s(n−2)/2 ≡ 7 for odd n and ≡ −2 for even n. We prove it. The proof is Fermat's little theorem in the numbers a + b√3, combined with the quadratic characters of 2 and 3 modulo N; the parity of n only decides whether √3 exists mod N. The test is a Lucas probable-prime test. No composite N passes it for n ≤ 10,000. It costs as much as an ordinary Fermat test, so it gives no speed-up for record hunting.

**Conjecture:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_112.htm)

**Status:** proved (Q1) · Q2: no pseudoprime for 3 ≤ n ≤ 10,000 (exhaustive) · Q3, Q4 answered

## The question

Start from the Lucas–Lehmer sequence s(0) = 4, s(k+1) = s(k)² − 2, and let N = 2^(n−1) + 3. Rotondo observed:

```
N prime, n odd:   s(n−2)/2 ≡  7  (mod N)
N prime, n even:  s(n−2)/2 ≡ −2  (mod N)
```

Q1 asks for a proof. Q2 asks for a pseudoprime (a composite N that gives the same remainder). Q3 and Q4 ask whether the test could beat the PRP record for 2^n + 3 (n = 2205444) or serve as a screening test.

## What the page had vs ours

| | On the page | Ours |
|---|---|---|
| Q1 | checked for the primes up to n = 785 | **proved** for every n ≥ 3 |
| Q2 | "probably pseudoprimes exist" | **none for 3 ≤ n ≤ 10,000**; heuristically none at all |
| Q3, Q4 | open | same cost as a Fermat test: no gain for record hunting; as a filter it is as strong as a Fermat test |

## How we found it

The Lucas–Lehmer numbers have a closed form. With ω = 2 + √3 (so ω · (2 − √3) = 1):

```
s(k) = ω^(2^k) + ω^(−2^k)
```

So s(n−2) is ω raised to the power 2^(n−2) = (N − 3)/2, plus its inverse. For a prime N, Fermat's little theorem in the ring of numbers a + b√3 says what ω^((N±1)/2) is. The only extra fact needed is a small identity:

```
2ω = 4 + 2√3 = (1 + √3)²
```

It turns ω^((N±1)/2) into a power of (1 + √3) divided by a power of 2. Euler's criterion handles the power of 2, and Frobenius handles the power of 1 + √3.

## The math

For n ≥ 4, 8 divides 2^(n−1), so N ≡ 3 (mod 8). Hence (2/N) = −1, and by quadratic reciprocity (3/N) = −(N/3). Since N ≡ 2^(n−1) (mod 3):

```
n odd:   N ≡ 1 (mod 3)  →  (3/N) = −1   (√3 does not exist mod N)
n even:  N ≡ 2 (mod 3)  →  (3/N) = +1   (√3 exists mod N)
```

**n odd.** In the field of N² elements built from √3, Frobenius sends √3 to −√3, so (1 + √3)^N = 1 − √3 and (1 + √3)^(N+1) = (1 + √3)(1 − √3) = −2. Then

```
ω^((N+1)/2) = (1 + √3)^(N+1) / 2^((N+1)/2) = −2 / 2^((N+1)/2) = −1 / 2^((N−1)/2) = −1/(2/N) = 1
s(n−2) = ω^((N+1)/2 − 2) + ω^(−(N+1)/2 + 2) = ω^(−2) + ω^2 = s(1) = 14
```

so s(n−2)/2 ≡ 7.

**n even.** Now √3 lies in the integers mod N, so (1 + √3)^(N−1) = 1 and

```
ω^((N−1)/2) = (1 + √3)^(N−1) / 2^((N−1)/2) = 1/(2/N) = −1
s(n−2) = ω^((N−1)/2 − 1) + ω^(−(N−1)/2 + 1) = −ω^(−1) − ω = −4
```

so s(n−2)/2 ≡ −2. ∎

The two edge cases also check out. For n = 3, N = 7 ≡ 7 (mod 8), and s(1)/2 = 7 ≡ 0 as stated. For n = 2 (outside the page's list), N = 5 ≡ 1 (mod 4) flips the sign of (3/N), so N = 5 follows the "odd" rule. The real rule is about (3/N), not about the parity of n.

**Why it is not a proof of primality.** The Lucas–Lehmer test proves primality because N + 1 = 2^p has only the prime factor 2. Here N + 1 = 4·(2^(n−3) + 1) and N − 1 = 2·(2^(n−2) + 1) have unknown factors. The test checks one Euler–Lucas condition with P = 4, Q = 1, so it is a probable-prime test.

**Pseudoprimes (Q2).** A composite N passes for odd n only if every prime q | N satisfies ω^(2^(n−2)) ≡ ω^(±2) (mod q), i.e. the order of ω mod q divides 2·(2^(n−3) ∓ 1). For a prime factor q of size about √N, that order is about q. The chance that it divides a given number is about 1/q ≈ 2^(−n/2). The sum over n converges, so we expect no pseudoprime at all, though we cannot rule one out. The search found none:

```
n = 3 … 10,000:   every N tested: prime N always give the predicted remainder,
                  composite N never do (0 pseudoprimes)
prime (PRP above 2^64) N found for n = 3, 4, 5, 7, 8, 13, 16, 17, 19, 29, 31, 56, 68, 85, 229, 391, 785, 1111,
                  1705, 2009, 2140, 2192, 2368, 2371, 4003, 4061, 4063, 4553, 5548, 8740
                  (= A057732 + 1)
```

**Speed (Q3, Q4).** One test is n − 2 modular squarings, the same work as a Fermat test, so it cannot find large primes 2^n + 3 faster than the usual sieve-plus-Fermat search. As a screening test it is as strong as a Fermat test (no false positives known). Numbers 2^n + 3 can only be proved prime with general methods such as ECPP, because N ± 1 does not factor.

## Verify it yourself

This script checks the proof's ingredients for many n (the characters of 2 and 3, and the values of ω^((N±1)/2) for prime N). It then runs the test on every N with 3 ≤ n ≤ 1500: prime N must give the predicted remainder, and composite N must not. It runs in a few seconds. We ran the same loop with gmpy2 up to n = 10,000 (about 10 minutes).

```python
from sympy import isprime, jacobi_symbol

LIMIT = 1500
for n in range(3, LIMIT + 1):
    N = 2**(n - 1) + 3
    s = 4
    for _ in range(n - 2):
        s = (s * s - 2) % N
    r = s * pow(2, -1, N) % N                     # s(n-2)/2 mod N
    target = 7 % N if n % 2 else N - 2
    if n >= 4:                                    # ingredients of the proof
        assert N % 8 == 3 and jacobi_symbol(2, N) == -1
        assert jacobi_symbol(3, N) == (-1 if n % 2 else 1)
    if isprime(N):
        assert r == target, n                     # the theorem
    else:
        assert r != target, n                     # no pseudoprime in this range

# omega^((N+1)/2) = 1 (n odd) and omega^((N-1)/2) = -1 (n even), computed in Z_N[sqrt 3]
def mul(a, b, N):                                 # (a0 + a1*r)(b0 + b1*r), r^2 = 3
    return ((a[0]*b[0] + 3*a[1]*b[1]) % N, (a[0]*b[1] + a[1]*b[0]) % N)
def power(a, e, N):
    out = (1, 0)
    while e:
        if e & 1: out = mul(out, a, N)
        a, e = mul(a, a, N), e >> 1
    return out
for n in [5, 7, 13, 17, 19, 29, 31, 85, 229, 391, 785, 4, 8, 16, 56, 68, 1111, 2140]:
    N = 2**(n - 1) + 3
    assert isprime(N)
    if n % 2: assert power((2, 1), (N + 1)//2, N) == (1, 0)
    else:     assert power((2, 1), (N - 1)//2, N) == (N - 1, 0)
print("Conjecture 112 verified for 3 <= n <=", LIMIT)
```
