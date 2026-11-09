# Puzzle 1185 / OEIS A247148: primes p that divide their look-and-say value LS(p).
# Checks the four known terms from scratch, searches all p < 1e7 in Python (independent of the
# C search in p1185_search.c, whose output for p < 1e11 is checked in data/p1185_1e11.txt),
# and computes the heuristic count.  Run from the repo root: .venv/bin/python3 writeup/code/p1185_verify.py
import re, sys
from math import log, exp, factorial
from sympy import isprime, factorint, primerange

def LS(n):
    s = str(n)
    return int("".join(str(len(m.group(0))) + m.group(2) for m in re.finditer(r"((\d)\2*)", s)))
assert LS(1211) == 111221 and LS(2) == 12 and LS(5) == 15 and LS(11111111111) == 111

known = [2, 5, 108947, 4940867]
for p in known:
    v = LS(p); assert isprime(p) and v % p == 0
    f = factorint(v); assert p in f
    print(p, v, f, "quotient", v // p)
# the quotient ends in 1: LS(p) ends in the last digit of p, so (m - 1) p = 0 (mod 10) for p coprime to 10
assert all((LS(p) // p) % 10 == 1 for p in known[2:])

# independent search below 1e7
sols = [p for p in primerange(2, 10**7) if LS(p) % p == 0]
print("solutions below 1e7:", sols); assert sols == known

# the C search output (p < 1e11)
lines = open("writeup/data/p1185_1e11.txt").read().split("\n")
csol = [(int(a), int(b)) for a, b in (l.split() for l in lines if l and l[0].isdigit())]
assert [p for p, _ in csol] == known and all(LS(p) == v for p, v in csol)
m = re.search(r"pi\((\d+)\) = (\d+), solutions = (\d+)", "\n".join(lines))
assert m.group(1) == "100000000000" and m.group(2) == "4118054813" and m.group(3) == "4"   # pi(1e11) (literature)
print("C search: pi(1e11) =", m.group(2), "solutions", m.group(3))

# heuristic: p | LS(p) with probability 1/p; expected count below X is sum 1/p = ln ln X + M
M = 0.2614972128
exact = sum(1.0 / p for p in primerange(2, 10**7))
print("sum 1/p for p < 1e7: %.4f, ln ln X + M: %.4f" % (exact, log(log(1e7)) + M))
def E(X): return log(log(X)) + M
for X in (1e7, 1e11, 4e12, 1e13, 1e15, 1e20, 1.6e25, 1e30):
    print("expected number of solutions with p <= %.1e: %.2f" % (X, E(X)))
lam = E(4e12)
print("Poisson probability of at least 4 below 4e12 (mean %.2f): %.2f" % (lam, 1 - sum(exp(-lam) * lam**k / factorial(k) for k in range(4))))
for X in (1e13, 1e15, 1e20, 1e25):
    print("probability that a(5) < %.0e given a(5) > 4e12: %.3f" % (X, 1 - exp(-(E(X) - E(4e12)))))
med = exp(2 * log(4e12))
print("median of a(5) given a(5) > 4e12: %.1e (ln ln X rises by ln 2, i.e. X = (4e12)^2)" % med)
print("all checks passed")
