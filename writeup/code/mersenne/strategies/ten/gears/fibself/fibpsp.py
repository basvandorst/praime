import sys, numpy as np, time
N = int(float(sys.argv[1])); CH = 4_000_000
s = np.ones(N, dtype=bool); s[:2] = False
for p in range(2, int(N ** 0.5) + 1):
    if s[p]: s[p * p::p] = False
BITS = int(N).bit_length()
def fib_at(m, n):
    a = np.zeros_like(n); b = np.ones_like(n) % n
    for i in range(BITS - 1, -1, -1):
        c = a * ((2 * b - a) % n) % n; d = (a * a + b * b) % n
        bit = ((m >> i) & 1).astype(bool)
        a, b = np.where(bit, d, c), np.where(bit, (c + d) % n, d)
    return a
def pow2(e, n):
    r = np.ones_like(n) % n; b = np.full_like(n, 2) % n
    for i in range(BITS):
        bit = ((e >> i) & 1).astype(bool); r = np.where(bit, r * b % n, r); b = b * b % n
    return r
fib, pair, byclass, first = [], 0, {1: 0, 2: 0, 3: 0, 4: 0}, None
t0 = time.time()
for lo in range(3, N, 2 * CH):
    n = np.arange(lo, min(N, lo + 2 * CH), 2, dtype=np.int64)
    n = n[(n % 5 != 0) & ~s[n]]
    r5 = n % 5; eps = np.where((r5 == 1) | (r5 == 4), 1, -1)
    hit = fib_at(n - eps, n) == 0
    h = n[hit]
    if h.size:
        fib.extend(h.tolist())
        f2 = pow2(h - 1, h) == 1
        pair += int(np.sum(f2 & ((h % 5 == 2) | (h % 5 == 3))))
        for k in (1, 2, 3, 4): byclass[k] += int(np.sum(h % 5 == k))
print(f"N={N:.0e}: odd composite Fibonacci pseudoprimes (F(n-(5/n)) = 0 mod n): {len(fib)}; first {fib[:12]}")
print(f"  by n mod 5: {byclass}; also base-2 Fermat and n = 2,3 mod 5: {pair}; {time.time()-t0:.0f}s")
open(f"fibpsp_{N:.0e}.txt", "w").write("\n".join(map(str, fib)) + "\n")
