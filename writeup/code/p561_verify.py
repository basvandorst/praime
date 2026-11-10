# Puzzle 561: independent check (sympy, plain integers; every number is below 2^64, where sympy.isprime is a proof).
# (1) The page's runs: lengths recomputed.
# (2) The two new runs of 13: the primes are consecutive, the 13 sums are prime, and the run cannot be extended.
# (3) Completeness of the wheel search: for the runs of 7 starting at 1213 and 1231, every even k up to the period
#     2*5*13*17*29*37*41 is tested directly and compared with the output of p561_search.py.
# Run from the repo root: (ulimit -v 1500000; .venv/bin/python3 writeup/code/p561_verify.py)   about 4 minutes
import subprocess, sys, os
from sympy import isprime, nextprime, prevprime
import gmpy2

def run_length(p, k):
    """number of consecutive primes from p, with even numbers from k, that give prime sums"""
    n = 0
    while isprime(p * p + k * k): n += 1; p = nextprime(p); k += 2
    return n

page = [("puzzle text", 5, 16, 4), ("Alm", 10837, 1058, 5), ("Alm", 22717, 678, 5), ("Tomita", 29, 94, 5), ("Tomita", 59, 54, 6),
        ("Tomita", 1259, 306, 7), ("Tomita", 47, 7330, 8), ("Tomita", 35753, 97652, 9), ("Tomita", 136417, 91232, 10),
        ("Schalekamp", 53, 312, 7), ("Vantieghem", 3923, 300539662, 10), ("Heleen", 17581, 125130, 9),
        ("Vantieghem", 41549, 395498534, 11), ("Vantieghem", 648779, 2079679170, 11), ("Vantieghem", 827639, 987705646, 11),
        ("Vantieghem", 1048613, 2033153978, 12)]
for who, p, k, claimed in page:
    assert isprime(p)
    n = run_length(p, k)
    print(f"page ({who}): p = {p}, k = {k}: {n} primes (page: {claimed})")
    assert n >= claimed
assert 300539662 ** 2 == 90324088435074244

for p, k in ((236543033, 56455112), (242065997, 12859030)):
    assert isprime(p) and k % 2 == 0
    n = run_length(p, k)
    q = prevprime(p); before = isprime(q * q + (k - 2) ** 2)
    ps = [p]
    for _ in range(12): ps.append(nextprime(ps[-1]))
    sums = [x * x + (k + 2 * i) ** 2 for i, x in enumerate(ps)]
    assert all(s < 2 ** 64 and isprime(s) for s in sums) and n == 13 and not before
    print(f"new: p = {p}, k = {k}: {n} primes; the 14th sum and the sum before the run are composite")
    for i, x in enumerate(ps): print(f"   {x}^2 + {k + 2 * i}^2 = {sums[i]}")

M = 2 * 5 * 13 * 17 * 29 * 37 * 41
here = os.path.dirname(os.path.abspath(__file__))
out = subprocess.run([sys.executable, os.path.join(here, "p561_search.py"), "7", "1213", "1232"], capture_output=True, text=True, check=True).stdout
fast = sorted((int(l.split()[2]), int(l.split()[3])) for l in out.splitlines() if l.startswith("HIT") and int(l.split()[2]) in (1213, 1231))
slow = []
for p in (1213, 1231):
    ps = [p]
    for _ in range(6): ps.append(nextprime(ps[-1]))
    sq = [x * x for x in ps]
    for k in range(2, M + 1, 2):
        if gmpy2.is_prime(sq[0] + k * k) and all(gmpy2.is_prime(sq[i] + (k + 2 * i) ** 2) for i in range(1, 7)):
            assert all(isprime(sq[i] + (k + 2 * i) ** 2) for i in range(7)); slow.append((p, k))
assert fast == sorted(slow) and len(slow) > 0, (len(fast), len(slow))
print(f"wheel search = direct test for the runs of 7 from 1213 and 1231, even k <= {M}: {len(slow)} solutions each way")
print("all checks passed")
