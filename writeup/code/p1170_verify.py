# Puzzle 1170: perfect numbers and harmonic (Ore) numbers.
# Q1: an even perfect number has no harmonic proper divisor > 1 (proof in the write-up; checked here for the 8 even perfect numbers
#     below 6.5e19 and, by the formulas, for the first 20 Mersenne exponents).
# Q2: among the known harmonic numbers (OEIS A001599 b-file, 10000 terms to 6.4e19, complete below 1e14 by Goto-Shibata;
#     plus Borysiuk's conjectured list a001599_1.txt), the ones with no harmonic proper divisor > 1 are the even perfect numbers and 6200.
# Run from the repo root: .venv/bin/python3 writeup/code/p1170_verify.py
import os
from math import prod
from sympy import factorint, divisors, isprime
from sympy.functions.combinatorial.numbers import divisor_sigma

def harmonic(n, f=None):
    """n*tau(n)/sigma(n) is an integer (harmonic mean of the divisors of n)."""
    f = f or factorint(n)
    tau = prod(e + 1 for e in f.values()); sig = prod((p**(e + 1) - 1)//(p - 1) for p, e in f.items())
    return (n*tau) % sig == 0

here = os.path.dirname(__file__)
terms = set()
for fn in ("p1170_b001599.txt", "p1170_a001599_1.txt"):
    for line in open(os.path.join(here, "..", "data", fn)):
        parts = line.split()
        if parts: terms.add(int(parts[-1]))
terms = sorted(terms)
assert len(terms) == 10006 and terms[-1] == 64555935231675888000
assert all(harmonic(t) for t in terms), "a listed number is not harmonic"
assert terms[:10] == [1, 6, 28, 140, 270, 496, 672, 1638, 2970, 6200]

# Q2: seeds = harmonic numbers > 1 whose proper divisors > 1 are all non-harmonic.
seeds = [t for t in terms if t > 1 and not any(harmonic(d) for d in divisors(t) if 1 < d < t)]
perfect = [t for t in terms if divisor_sigma(t) == 2*t]
print("known harmonic numbers:", len(terms), "largest", terms[-1])
print("perfect numbers among them:", perfect)
print("harmonic numbers with no harmonic proper divisor > 1:", seeds)
assert perfect == [6, 28, 496, 8128, 33550336, 8589869056, 137438691328, 2305843008139952128]
assert seeds == sorted(perfect + [6200])
f = factorint(6200); print("6200 =", f, "harmonic mean", 6200*prod(e + 1 for e in f.values())//int(divisor_sigma(6200)))

# Q1 by the formulas: N = 2^(p-1) M, M = 2^p - 1 prime. Divisors 2^a (a <= p-1) and 2^a M (a <= p-2).
# H(2^a) = 2^a (a+1)/(2^(a+1) - 1): integer only for a = 0.  H(2^a M) = M (a+1) / (2^(p-a-1) (2^(a+1) - 1)).
from sympy import gcd
exps = [p for p in range(2, 128) if isprime(p) and isprime(2**p - 1)]
assert exps == [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127]
for p in exps:
    M = 2**p - 1
    for a in range(1, p):            # 2^a, a >= 1
        assert (2**a*(a + 1)) % (2**(a + 1) - 1) != 0
    for a in range(0, p - 1):        # 2^a M, proper divisor
        assert (M*(a + 1)) % (2**(p - a - 1)*(2**(a + 1) - 1)) != 0
        # the reason: gcd(2^(a+1) - 1, M) = 2^gcd(a+1, p) - 1 = 1 for 0 <= a < p - 1, and 2^(a+1) - 1 > a + 1 for a >= 1
        assert gcd(2**(a + 1) - 1, M) == 1
print("Q1 checked by the divisor formulas for the Mersenne exponents", exps)
print("all checks passed")
