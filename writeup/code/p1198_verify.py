# Puzzle 1198: ssd(n) = sum of sigma(d) over the divisors d of n; "Fred numbers" have ssd(n) = sigma(sigma(n)).
# Independent check with sympy (divisors listed explicitly, no multiplicative formula):
# (1) all n <= 2*10^5 by brute force: only 1, 2, 360;
# (2) every n given on the command line (default: the page's 15360620700) is a Fred number.
# Run from the repo root: .venv/bin/python3 writeup/code/p1198_verify.py [n ...]
import sys
from sympy import divisors, divisor_sigma, factorint
def ssd(n): return sum(int(divisor_sigma(d)) for d in divisors(n))
def fred(n): return ssd(n) == int(divisor_sigma(int(divisor_sigma(n))))
small = [n for n in range(1, 200001) if fred(n)]
print("Fred numbers up to 200000:", small); assert small == [1, 2, 360]
for n in [int(x) for x in sys.argv[1:]] or [15360620700]:
    assert fred(n), n
    print(f"{n} = {factorint(n)}: sigma = {divisor_sigma(n)}, ssd = sigma(sigma) = {ssd(n)}")
print("all checks passed")
