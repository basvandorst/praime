"""Primes k in [lo, hi) with no factor q = 2jk+1 of 2^k-1 for j <= J (J = 10^5).
Vectorised: for each k, 2^k mod q for all q = 2jk+1 with q = +-1 mod 8 (theorem:
every factor of 2^k-1 has this form), uint64 (q < 4.1e9 so q^2 < 2^64).
Output: survivors_<lo>_<hi>.txt  (k per line) and the count of primes < 10^4 left
for comparison with MERSENNE.md (532 of 1229 at K = 10^5, k >= 2).
"""
import sys, numpy as np, gmpy2

J = 10**5

def has_factor(k):
    j = np.arange(1, J + 1, dtype=np.uint64)
    q = np.uint64(2 * k) * j + np.uint64(1)
    m8 = q % np.uint64(8)
    q = q[(m8 == 1) | (m8 == 7)]
    if k < 64:
        q = q[q < np.uint64(2**k - 1)]   # proper factors only
    r = np.ones_like(q)
    for bit in bin(k)[2:]:
        r = (r * r) % q
        if bit == '1':
            r = (r * np.uint64(2)) % q
    hit = np.nonzero(r == 1)[0]
    return int(q[hit[0]]) if hit.size else 0

def run(lo, hi):
    out = []
    k = int(gmpy2.next_prime(lo - 1))
    while k < hi:
        if has_factor(k) == 0:
            out.append(k)
        k = int(gmpy2.next_prime(k))
    return out

if __name__ == '__main__':
    lo, hi = int(sys.argv[1]), int(sys.argv[2])
    s = run(lo, hi)
    with open(f'survivors_{lo}_{hi}.txt', 'w') as f:
        f.write('\n'.join(map(str, s)) + '\n')
    print(lo, hi, 'survivors', len(s))
