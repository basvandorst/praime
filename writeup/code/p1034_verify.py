# Puzzles 1034 / 1091: chains p, p+6k, p+6k+6(k+1), ... of primes = runs of consecutive n with 3n^2+3n+c prime.
# Checks: the record run (23 primes, c = 3789041) and the page's examples; the Hardy-Littlewood constant of the
# record polynomial against its actual prime density; every run listed in data/p1034_runs.txt and
# data/p1034_small.txt recomputed with sympy; the Python twin (p1034_search.py) rerun on |c| <= 300000 against
# the C output.   Run from the repo root: .venv/bin/python3 writeup/code/p1034_verify.py
import re, subprocess
from math import log
from sympy import isprime, primerange, legendre_symbol
def run_from(c, n0):
    k = 0
    while isprime(3 * (n0 + k)**2 + 3 * (n0 + k) + c): k += 1
    return k
def chain(p, step6):                 # the page's form: p, p + 6k, + 6(k+1), ...
    seq = [p]; s = step6
    while isprime(seq[-1] + s): seq.append(seq[-1] + s); s += 6
    return seq
assert len(chain(11, 6)) == 10 and len(chain(23, 6)) == 22 and len(chain(4745021, 3390)) == 23   # page: 11 (9 steps), 23 (21), Volpatti (22 steps)
assert chain(4745021, 3390)[-1] == 4820987 and not isprime(4824509)
print("page chains ok: 11 -> 10 primes, 23 -> 22 primes, Volpatti 4745021 -> 23 primes")
c = 4745021 - 3 * 564 * 565; assert c == 3789041 and run_from(c, 564) == 23 and not isprime(3*563*563 + 3*563 + c)
D = 9 - 12 * c; C = 3.0
for q in primerange(5, 100000):
    w = 1 + legendre_symbol(D % q, q) if D % q else 1
    C *= (1 - w / q) / (1 - 1 / q)
vals = [3*n*n + 3*n + c for n in range(1, 20001)]
ratio = sum(isprime(v) for v in vals) / sum(1 / log(v) for v in vals)
print(f"record polynomial: c = {c}, D = {D}, Hardy-Littlewood constant {C:.2f}, observed prime density / random = {ratio:.2f} (n <= 20000)")
for fn in ("writeup/data/p1034_runs.txt", "writeup/data/p1034_small.txt"):
    for m in re.finditer(r"run (\d+): c=(-?\d+) C=([\d.]+) n0=(\d+) p=(\d+)", open(fn).read()):
        r, cc, _, n0, p = int(m.group(1)), int(m.group(2)), m.group(3), int(m.group(4)), int(m.group(5))
        assert 3*n0*n0 + 3*n0 + cc == p and run_from(cc, n0) == r and (n0 == 0 or not isprime(3*(n0-1)**2 + 3*(n0-1) + cc)), (cc, n0)
summary = re.search(r"c in .*", open("writeup/data/p1034_runs.txt").read()).group(0)
print("all listed runs recomputed (exact lengths);", summary)
out = subprocess.run([".venv/bin/python3", "writeup/code/p1034_search.py", "300000", "2000", "6.5", "20"], capture_output=True, text=True).stdout
assert "run 22: c=-205609" in out and "longest run (22, -205609)" in out
print("Python twin on |c| <= 300000 agrees with the C search (longest run 22 at c = -205609)")
print("all checks passed")
