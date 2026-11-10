# Puzzle 582 (follow-up question), local test.  Modulo q, can a product of n numbers c + k i (c^2 + k^2 coprime to q)
# be a +- k i?  If not, no m = a^2 + k^2 is a product of n primes c^2 + k^2 (all such primes exceed the prime
# factors of q).  usage: p582_local.py k [k ...]
import sys, math
def reach(k, q, nmax):
    gens = [(c % q, k % q) for c in range(q) if math.gcd(c * c + k * k, q) == 1]
    cur = set(gens); out = {}
    for n in range(1, nmax + 1):
        out[n] = any(((a, k % q) in cur or (a, -k % q) in cur) for a in range(q))
        cur = {((x * c - y * d) % q, (x * d + y * c) % q) for (x, y) in cur for (c, d) in gens}
    return out
for k in map(int, sys.argv[1:]):
    pmin = next(c * c + k * k for c in range(1, 10 ** 4) if all((c * c + k * k) % d for d in range(2, int((c * c + k * k) ** .5) + 1)))
    bad = {}
    for q in (16, 20, 80):
        if max(p for p in range(2, q + 1) if q % p == 0 and all(p % d for d in range(2, p))) >= pmin: continue
        r = reach(k, q, 11)
        blocked = [n for n in r if n >= 2 and not r[n]]
        if blocked: bad[q] = blocked
    print("k =", k, "smallest good prime", pmin, "numbers of factors (2..11) that are impossible, by modulus:", bad)
