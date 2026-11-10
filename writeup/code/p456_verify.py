# Puzzle 456: independent check (sympy, plain Python lists, two-pointer).
# (1) For e = 2, 3, 5, 7, 13, 17, 19: all windows of k >= 2 consecutive primes with sum 2^e - 1, by a plain two-pointer
#     over sympy.primerange; compared with p456_search.windows.
# (2) The page's windows for 2^61 - 1 (175 primes) and 2^127 - 1 (3 primes): consecutive primes, exact sum.
# (3) 2^31 - 1: p456_smallk.py (a different method: next_prime around M/k) for every k up to 500 finds no window.
# Run from the repo root: (ulimit -v 1500000; .venv/bin/python3 writeup/code/p456_verify.py)
import sys, os, subprocess
from sympy import primerange, nextprime, isprime
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
from p456_search import windows
for e in (2, 3, 5, 7, 13, 17, 19):
    M = 2 ** e - 1; ps = list(primerange(2, M)); out = []; i = 0; s = 0
    for j, p in enumerate(ps):
        s += p
        while s > M: s -= ps[i]; i += 1
        if s == M and j > i: out.append((j - i + 1, ps[i], p))
    fast = windows(M)
    assert len(fast) == len(out) and [b for a, b in fast] == [w[2] for w in out], (e, fast, out)
    print(f"2^{e} - 1: {len(out)} windows (k, first, last): {out}")
for e, k, first, last in ((61, 175, 13176245766932173, 13176245766938521), (127, 3, 56713727820156410577229101238628035201, 56713727820156410577229101238628035283)):
    assert isprime(first); w = [first]
    while len(w) < k: w.append(nextprime(w[-1]))
    assert w[-1] == last and sum(w) == 2 ** e - 1
    print(f"2^{e} - 1 = sum of the {k} consecutive primes from {first} to {last}")
out = subprocess.run([sys.executable, os.path.join(here, "p456_smallk.py"), "31", "500"], capture_output=True, text=True, check=True).stdout
print(out.strip()); assert ": 0 windows" in out
print("all checks passed")
