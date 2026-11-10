# Puzzle 456: M = 2^e - 1 as a sum of k consecutive primes, for every k with 2 <= k <= KMAX (gmpy2.next_prime).
# For each k: take k consecutive primes around M/k, then slide the window up or down until the sum passes M.
# The window sum is strictly increasing in the position, so each k has at most one window and it is found.
# usage: p456_smallk.py e KMAX [KMIN]
import sys, gmpy2
from collections import deque
e, KMAX = int(sys.argv[1]), int(sys.argv[2]); KMIN = int(sys.argv[3]) if len(sys.argv) > 3 else 2
M = 2 ** e - 1
def prev_prime(n):
    n -= 1
    while not gmpy2.is_prime(n): n -= 1
    return n
hits = []
for k in range(KMIN, KMAX + 1):
    if k % 2 == 0 and k > 2: continue                     # an even number of odd primes has an even sum
    if k == 2:
        if gmpy2.is_prime(M - 2) and M - 2 == 3: hits.append((2, 2, 3))   # 2 + q with q the next prime: only 5 = 2 + 3
        continue
    c = M // k
    w = deque(); p = int(gmpy2.next_prime(c - 1)) if c > 2 else 2
    w.append(p)
    lo = p
    for _ in range((k - 1) // 2):
        if lo <= 2: break
        lo = prev_prime(lo); w.appendleft(lo)
    while len(w) < k: w.append(int(gmpy2.next_prime(w[-1])))
    s = sum(w)
    while s < M:
        nxt = int(gmpy2.next_prime(w[-1])); s += nxt - w.popleft(); w.append(nxt)
    while s > M and w[0] > 2:
        prv = prev_prime(w[0]); s += prv - w.pop(); w.appendleft(prv)
    if s == M: hits.append((k, w[0], w[-1])); print(f"2^{e} - 1 = sum of {k} consecutive primes from {w[0]} to {w[-1]}", flush=True)
print(f"2^{e} - 1: every k from {KMIN} to {KMAX} tested: {len(hits)} windows {[h[0] for h in hits]}")
