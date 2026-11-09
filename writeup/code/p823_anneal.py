"""Puzzle 823, Python version of p823_anneal.cpp.  Strings over {1,3,7,9} of length L; K = number of distinct
primes among the substrings read left to right.  Simulated annealing over single-digit changes with a memo of
substring primality (gmpy2).  Usage: p823_anneal.py L iterations seed [start_string]"""
import sys, math, random
from gmpy2 import mpz, is_prime
L, iters, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
rng = random.Random(seed); D = "1379"
memo = {}
def prime(t):
    r = memo.get(t)
    if r is None: r = memo[t] = bool(is_prime(mpz(t)))
    return r
def K(s): return len({s[i:j] for i in range(len(s)) for j in range(i + 1, len(s) + 1) if prime(s[i:j])})
s = sys.argv[4] if len(sys.argv) > 4 else "".join(rng.choice(D) for _ in range(L))
cur = K(s); best, bests = cur, s
T0, T1 = 2.0, 0.2
for it in range(iters):
    T = T0 * (T1 / T0) ** (it / iters)
    t = list(s); mv = rng.random()
    if mv < 0.8:
        k = rng.randrange(L); t[k] = rng.choice(D.replace(t[k], ""))
    elif mv < 0.9:
        a, b = rng.randrange(L), rng.randrange(L); t[a], t[b] = t[b], t[a]
    else:
        for k in rng.sample(range(L), 2): t[k] = rng.choice(D)
    t = "".join(t); val = K(t)
    if val >= cur or math.exp((val - cur) / T) > rng.random():
        s, cur = t, val
        if cur > best: best, bests = cur, s; print(f"it {it}: K={best} {bests}", flush=True)
print(f"final best K={best} {bests}")
