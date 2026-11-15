"""Independent check of the llres note with gmpy2: LL residues for odd primes p <= 3000, zeros = exponents, KS of the rest,
and the sign of s_(p-3) = +-2^((p+1)/2) at the exponents."""
import gmpy2, sympy, math
zeros, signs, ys = [], [], []
for p in sympy.primerange(3, 3001):
    M = (1 << p) - 1; s = gmpy2.mpz(4); prev = s
    for _ in range(p - 2): prev = s; s = (s * s - 2) % M
    if s == 0:
        h = pow(2, (p + 1) // 2, M); zeros.append(p); signs.append(1 if prev == h else -1 if prev == (M - h) % M else 0)
    else: ys.append(int(s) / M)
ys.sort(); n = len(ys); ks = max(max(abs(y - i / n), abs(y - (i + 1) / n)) for i, y in enumerate(ys))
print('zeros:', zeros); print('signs:', signs); print(f'nonzero {n}, KS {ks:.3f}, 5% critical {1.36/math.sqrt(n):.3f}')
