"""mersieve, main computation (numpy, vectorised over p).

For every prime p <= N and every k in 1..KMAX, decide whether q = 2kp+1
divides M_p = 2^p - 1 with q < M_p.  Writes
  ms_pairs.txt : every (p, k) with q | M_p, q < M_p, k <= KMAX (all factors, not only the smallest)
  ms_first.txt : p, smallest k (0 = none <= KMAX)
Arithmetic: q < 2*1000*500000+1 ~ 1e9 < 2^30, so r*r < 2^60 fits in uint64.
"""
import numpy as np, time

N = 500_000
KMAX = 1000

def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]

t0 = time.time()
P = primes_upto(N).astype(np.uint64)
print("primes <= N:", len(P))
nb = int(P.max()).bit_length()
pairs = []
first = np.zeros(len(P), dtype=np.int64)
for k in range(1, KMAX + 1):
    q = np.uint64(2 * k) * P + np.uint64(1)
    r = np.ones(len(P), dtype=np.uint64)
    for b in range(nb - 1, -1, -1):
        r = (r * r) % q
        bit = ((P >> np.uint64(b)) & np.uint64(1)).astype(bool)
        r = np.where(bit, (r * np.uint64(2)) % q, r)
    hit = (r == 1)
    # proper factor: q < 2^p - 1.  Only matters for small p.
    for i in np.nonzero(hit)[0]:
        p = int(P[i]); qq = 2 * k * p + 1
        if qq < (1 << p) - 1:
            pairs.append((p, k))
            if first[i] == 0:
                first[i] = k
print("pairs:", len(pairs), "time %.1fs" % (time.time() - t0))
with open("ms_pairs.txt", "w") as f:
    for p, k in sorted(pairs):
        f.write(f"{p} {k}\n")
with open("ms_first.txt", "w") as f:
    for p, k in zip(P.tolist(), first.tolist()):
        f.write(f"{p} {k}\n")
