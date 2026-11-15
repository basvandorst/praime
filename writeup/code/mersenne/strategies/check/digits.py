"""Independent check of the digitrows note (plain Python, sympy): rows by digit sum in base B for n <= 500,000."""
import sympy, math
N = 500000
for B in (10, 7, 2):
    def ds(n):
        s = 0
        while n: s += n % B; n //= B
        return s
    rows = {}; prim = {}
    for n in range(1, N + 1):
        s = ds(n); rows[s] = rows.get(s, 0) + 1
        if sympy.isprime(n): prim[s] = prim.get(s, 0) + 1
    ruled = [s for s in rows if math.gcd(s, B - 1) > 1]
    leak = sorted(p for p in sympy.primerange(2, N + 1) if math.gcd(ds(p), B - 1) > 1)
    empty = sorted(s for s in rows if s not in prim)
    print(f'B={B}: rows {len(rows)} (1..{max(rows)}); ruled out {len(ruled)}; primes in ruled rows {leak}; rows without primes {empty}')
