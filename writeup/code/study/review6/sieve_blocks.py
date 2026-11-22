# own segmented sieve: prime counts in [1e9+k*1e8, +1e8), k=0..9
import numpy as np, math, json, time
def small_primes(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]
L0, B, SEG = 10**9, 10**8, 10**7
P = small_primes(int(math.isqrt(2*10**9))+1)[1:]  # odd primes
counts = []
t0=time.time()
for k in range(10):
    n = 0
    for s in range(B//SEG):
        lo = L0 + k*B + s*SEG          # even
        # odd numbers lo+1, lo+3, ..., lo+SEG-1 : index i -> lo+1+2i
        m = np.ones(SEG//2, bool)
        for p in P:
            p = int(p)
            st = (-(lo+1)) % p        # first j with lo+1+j ≡ 0 mod p
            if st % 2: st += p        # need lo+1+j with j even (odd number)
            idx = st//2
            if lo+1+st == p: idx += p
            m[idx::p] = False
        n += int(m.sum())
    counts.append(n); print(k, n, round(time.time()-t0), flush=True)
json.dump(counts, open('block_counts.json','w'))
print('total', sum(counts))
