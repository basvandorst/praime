"""mersieve, independent recheck (pure Python ints, sympy primes, loop over p then k).

Recomputes the smallest k <= 1000 with q = 2kp+1 | 2^p - 1 and q < 2^p - 1 for every
prime p <= 500000 using the built-in pow (no numpy, no gmpy2), and compares with
ms_first.txt.  Also recomputes the full pair list for k <= 1000 and compares with
ms_pairs.txt.
"""
import sympy, time

N, KMAX = 500_000, 1000
t0 = time.time()
first = {}
pairs = set()
for p in sympy.primerange(2, N + 1):
    M = (1 << p) - 1 if p < 64 else None
    fk = 0
    q = 1
    for k in range(1, KMAX + 1):
        q += 2 * p          # q = 2kp + 1
        if pow(2, p, q) == 1 and (M is None or q < M):
            pairs.add((p, k))
            if fk == 0:
                fk = k
    first[p] = fk
print("primes:", len(first), "pairs:", len(pairs), "time %.0fs" % (time.time() - t0))
with open("chk_first.txt", "w") as f:
    for p in sorted(first):
        f.write(f"{p} {first[p]}\n")
with open("chk_pairs.txt", "w") as f:
    for p, k in sorted(pairs):
        f.write(f"{p} {k}\n")

ref_first = {}
for line in open("ms_first.txt"):
    p, k = map(int, line.split()); ref_first[p] = k
ref_pairs = set()
for line in open("ms_pairs.txt"):
    p, k = map(int, line.split()); ref_pairs.add((p, k))
print("first-k table identical:", ref_first == first)
print("pair list identical:", ref_pairs == pairs)
