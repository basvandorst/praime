# Puzzle 463: independent check (sympy; the index k from sympy.primepi / a plain list).
# (1) All runs of >= 5 consecutive k with T_k = k p_k + (k+1) p_(k+1) prime, for p_k < 2*10^7, from a plain list of
#     primes; compared with p463_search.py.
# (2) For each run given as "k length" on the command line (default: the page's 7 at k = 321150): p_k from
#     sympy.prime(k), the T values are prime (PARI isprime), and the values before and after the run are composite.
# Run from the repo root: (ulimit -v 1500000; .venv/bin/python3 writeup/code/p463_verify.py [k length ...])
import sys, os, subprocess
from sympy import primerange, isprime, prime, nextprime, prevprime
import cypari2
pari = cypari2.Pari()
X = 2 * 10 ** 7
ps = [0] + list(primerange(2, X + 2000))                             # ps[k] = k-th prime
ok = [False] + [isprime(k * ps[k] + (k + 1) * ps[k + 1]) for k in range(1, len(ps) - 1)]
runs = []; k = 1
while k < len(ok):
    if ok[k]:
        j = k
        while j < len(ok) and ok[j]: j += 1
        if j - k >= 5 and ps[k] < X and ps[k] >= 3: runs.append((j - k, k, ps[k]))
        k = j
    else: k += 1
here = os.path.dirname(os.path.abspath(__file__))
out = subprocess.run([sys.executable, os.path.join(here, "p463_search.py"), "5", "3", str(X)], capture_output=True, text=True, check=True).stdout
fast = sorted(tuple(int(x) for x in l.split()[1:]) for l in out.splitlines() if l.startswith("RUN"))
assert fast == sorted(runs), (fast[:5], runs[:5])
print(f"p_k < {X}: the same {len(runs)} runs of 5 or more both ways; the longest: {max(runs)}")
args = [int(x) for x in sys.argv[1:]] or [321150, 7]
for k0, L in zip(args[::2], args[1::2]):
    p = [prevprime(prime(k0)), prime(k0)]                            # p[j] = p_(k0 - 1 + j)
    while len(p) < L + 3: p.append(nextprime(p[-1]))
    T = [(k0 - 1 + j) * p[j] + (k0 + j) * p[j + 1] for j in range(L + 2)]
    flags = [bool(pari.isprime(t)) for t in T]
    assert flags == [False] + [True] * L + [False], flags
    print(f"k = {k0}, p_k = {p[1]}: {L} primes in a row; the values before and after are composite")
    for j in range(1, L + 1): print(f"   {k0 - 1 + j}*{p[j]} + {k0 + j}*{p[j + 1]} = {T[j]}")
print("all checks passed")
