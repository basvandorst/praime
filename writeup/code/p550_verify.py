# Puzzle 550: independent check (sympy, plain integers; no modulo-12 filter).
# (1) Direct: every run of >= 3 semiprimes S_i = p_i p_(i+1) + p_(i+2) p_(i+3) with p_i < 3*10^7, semiprime tested
#     by sympy.factorint; compared with the output of p550_search.py on the same range.
# (2) The runs of 5: the eight primes are consecutive, the five sums are 2 * prime (S/2 proved prime with PARI
#     isprime), the sums before and after are not semiprimes.
# Run from the repo root: (ulimit -v 1500000; .venv/bin/python3 writeup/code/p550_verify.py)   about 3 minutes
import subprocess, sys, os
from sympy import primerange, nextprime, prevprime, isprime, factorint
import cypari2
pari = cypari2.Pari()
X = 3 * 10 ** 7
ps = list(primerange(2, X + 1000))
def semi(n):
    if n % 2 == 0: return n > 4 and isprime(n // 2) or n == 4
    return sum(factorint(n).values()) == 2
flags = [semi(ps[i] * ps[i + 1] + ps[i + 2] * ps[i + 3]) for i in range(len(ps) - 3)]
runs = []; i = 0
while i < len(flags):
    if flags[i]:
        j = i
        while j < len(flags) and flags[j]: j += 1
        if j - i >= 3 and ps[i] < X: runs.append((j - i, ps[i]))
        i = j
    else: i += 1
here = os.path.dirname(os.path.abspath(__file__))
out = subprocess.run([sys.executable, os.path.join(here, "p550_search.py"), "3", "5", str(X)], capture_output=True, text=True, check=True).stdout
fast = sorted((int(l.split()[1]), int(l.split()[2])) for l in out.splitlines() if l.startswith("RUN"))
slow = sorted(r for r in runs if r[1] >= 5)
assert fast == slow, (len(fast), len(slow))
print(f"p < {X}: {len(slow)} runs of 3 or more both ways; small starts (p < 5): {[r for r in runs if r[1] < 5]}; runs of 4: {[r[1] for r in slow if r[0] == 4]}")
assert [r[1] for r in slow if r[0] >= 4] == [17017471]

for p0 in [int(x) for x in sys.argv[1:]] or [4748135713, 22304301167]:
    q = [prevprime(p0), p0]
    while len(q) < 10: q.append(nextprime(q[-1]))
    S = [q[i] * q[i + 1] + q[i + 2] * q[i + 3] for i in range(7)]
    ok = [s % 2 == 0 and bool(pari.isprime(s // 2)) for s in S]
    L = ok[1:].index(False)
    assert ok[0] is False and ok[1:1 + L] == [True] * L
    for s in (S[0], S[1 + L]): assert sum(factorint(s).values()) != 2
    print(f"run of {L} from {p0}:")
    for i in range(1, 1 + L): print(f"   {q[i]}*{q[i + 1]} + {q[i + 2]}*{q[i + 3]} = 2 * {S[i] // 2}")
print("all checks passed")
