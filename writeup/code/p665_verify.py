# Puzzle 665: independent checks (sympy / gmpy2, plain integers).
# (1) Direct: for all consecutive primes p < q < r with p < 2*10^7, (p + q)(q + r) is a perfect power only for p = 5.
# (2) The equation solver p665_power.py against a direct test of every A <= 200000, g <= 300 (run by the solver).
# (3) The largest solutions the solver reports are perfect powers; none of them is p + q of consecutive primes.
# Run from the repo root: (ulimit -v 1500000; .venv/bin/python3 writeup/code/p665_verify.py)   about 2 minutes
import subprocess, sys, os, gmpy2
from sympy import primerange, prevprime, nextprime
ps = list(primerange(2, 2 * 10 ** 7 + 1000))
hits = [(p, q, r) for p, q, r in zip(ps, ps[1:], ps[2:]) if p < 2 * 10 ** 7 and gmpy2.is_power((p + q) * (q + r))]
print("direct, p < 2e7:", hits); assert hits == [(5, 7, 11)]
here = os.path.dirname(os.path.abspath(__file__))
out = subprocess.run([sys.executable, os.path.join(here, "p665_power.py"), "300", "1e30", "200000"], capture_output=True, text=True, check=True).stdout
print(out.strip()); assert "the same 790 solutions" in out and "[(5, 7, 11, 3)]" in out
for A, g, n in ((990100714252, 1330, 3), (445971980288, 2900, 3), (349774825888, 1854, 3), (611085348, 30, 3)):
    root, exact = gmpy2.iroot(A * (A + g), n); assert exact
    p = prevprime(A // 2 + 1); q = nextprime(p); assert (p + q, q + nextprime(q)) != (A, A + g)
    print(f"A = {A}, g = {g}: A(A + g) = {root}^{n}; not from consecutive primes")
print("all checks passed")
