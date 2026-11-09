"""Puzzle 917, constructive search.  P = A 10^d + B (2d digits, B may have leading zeros),
Q = P + A B, halved when even.  Solving for B:  Q + 10^d = (A+1)(B + 10^d)  in the odd branch,
2Q + 10^d = (A+1)(B + 10^d) in the even branch.  So the primes P that map to Q correspond to the
factorisations N = e f with f in [10^d, 2 10^d), N = Q + 10^d or 2Q + 10^d, P = e 10^d + f - 2 10^d
prime (and the parity of A B matching the branch, which is automatic: see the write-up).
Search: enumerate odd 5-free smooth N1 in [10^(2d-1) + 10^d, 10^(2d) + 10^d) with d(N1) >= T,
keep those with Q = N1 - 10^d prime, count prime P over both branches, report the best.
Usage: p917_construct.py d T [maxprime=61]        Needs gmpy2, sympy."""
import sys, time
from gmpy2 import mpz, is_prime
from sympy import primerange, divisors, factorint

d = int(sys.argv[1]); T = int(sys.argv[2]); PMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 61
M = 10**d; LO = 10**(2*d - 1) + M; HI = 10**(2*d) + M
primes = [p for p in primerange(3, PMAX + 1) if p != 5]
from math import log, floor

def preimages(N, branch):
    """primes P with (P + A B) = N - 10^d (branch 1) or (P + A B)/2 = (N - 10^d)/2 (branch 2)."""
    out = []
    for f in map(int, divisors(int(N))):
        if f < M or f >= 2*M: continue
        e = N // f; A = e - 1; B = f - M
        if A < M // 10 or A >= M: continue
        P = int(A * M + B)
        q = P + A * B
        if branch == 1 and (q % 2 == 1) and q == N - M and is_prime(P): out.append(P)
        if branch == 2 and (q % 2 == 0) and q // 2 == (N - M) // 2 and is_prime(P): out.append(P)
    return out

cands = {}   # N1 -> d(N1)
def dfs(i, prod, dn):
    if prod >= LO:
        if dn >= T: cands[prod] = dn
    if i == len(primes): return
    p = primes[i]
    # bound: remaining factor on d(N) at most 2^(log(HI/prod)/log p)
    if dn * 2 ** floor(log(HI / prod) / log(p)) < T: return
    dfs(i + 1, prod, dn)                      # skip p
    e = 1; q = prod * p
    while q < HI:
        dfs(i + 1, q, dn * (e + 1))
        e += 1; q *= p
t0 = time.time()
dfs(0, 1, 1)
print(f"d={d}: {len(cands)} candidates N1 with d(N1) >= {T}, {time.time()-t0:.1f}s", flush=True)
best = []; topk = 0; done = 0
for N1, dn in sorted(cands.items()):
    Q = N1 - M
    if not is_prime(Q): continue
    pre1 = preimages(N1, 1)
    pre2 = preimages(2 * Q + M, 2)
    k = len(pre1) + len(pre2)
    best.append((k, Q, N1, dn, sorted(pre1 + pre2))); done += 1
    if k > topk:                      # progress line, so that a killed run still leaves its best
        topk = k; print(f"  new best k={k} Q={Q} N1={N1} d(N1)={dn} after {done} prime Q, {time.time()-t0:.0f}s", flush=True)
best.sort(reverse=True)
print(f"{sum(1 for _ in best)} candidates with Q prime; top results:")
for k, Q, N1, dn, pre in best[:8]:
    strict = sum(1 for P in pre if P % M >= M // 10)
    print(f"k={k} (k={strict} with B >= 10^(d-1)) Q={Q} N1={N1}={dict(factorint(N1))} d(N1)={dn}")
    print("   P:", " ".join(map(str, pre)))
