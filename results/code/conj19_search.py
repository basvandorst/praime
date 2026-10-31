# Conjecture 19 search: extreme triplets p = k 2^n + 1, q = p^2 + p - 1, r = (p^3 + p^2 - p + 1)/2 (PRP filter; proofs in conj19_verify.py)
# Run: .venv/bin/python results/code/conj19_search.py 560 1 30000000   (about 10 s, finds k = 25971283)
# extreme C3 triplets: p = k 2^n + 1, q = p^2 + p - 1, r = (p^3 + p^2 - p + 1)/2, all prime
import sys, time, numpy as np, gmpy2
from sympy import primerange
n = int(sys.argv[1]); K0, K1 = int(sys.argv[2]), int(sys.argv[3]); B = 100000
t0 = time.time()
two_n = pow(2, n)
ks = np.arange(K0, K1, dtype=np.int64)
ok = (ks % 2 == 1)                                    # k odd
for l in primerange(3, B):
    x = np.arange(l, dtype=np.int64)
    bad_x = np.nonzero((x == 0) | ((x * x + x - 1) % l == 0) | (((x * x % l) * x + x * x - x + 1) % l == 0))[0]
    inv = pow(two_n % l, -1, l)
    for bx in bad_x:
        kc = (int(bx) - 1) * inv % l               # p = k 2^n + 1 = bx (mod l)
        start = (kc - K0) % l
        ok[start::l] = False
surv = ks[ok]
print(f"n={n}: {len(surv)} survivors of {K1-K0} k (sieve to {B}) in {time.time()-t0:.1f}s", flush=True)
hits = []
for k in surv:
    p = gmpy2.mpz(int(k)) * two_n + 1
    if not gmpy2.is_prime(p, 2): continue
    q = p * p + p - 1
    if not gmpy2.is_prime(q, 2): continue
    r = (p**3 + p * p - p + 1) // 2
    if not gmpy2.is_prime(r, 2): continue
    hits.append(int(k)); print("HIT k =", int(k), "digits C3 =", len(str(p * q * r)), flush=True)
    break
print("tested", len(surv), "time", round(time.time() - t0, 1), "hits", hits)
