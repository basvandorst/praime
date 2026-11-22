"""Python twin of qs.c / qs2.c on a reduced range (A <= 50, |C| <= 6000): longest run of consecutive y with |A y^2 + B y + C|
prime, 0 < B < A, A + B even, C odd, discriminant non-residue mod every odd p <= 23 not dividing A (p | A needs p | B).
Must find the Fung-Ruby 45 (36, 18, -1801) and 43 (47, 9, -5209). usage: qs_twin.py"""
from sympy import isprime, primerange
PR = list(primerange(3, 24)); Y = 80; out = []
for A in range(1, 51):
    for B in range(1, A):
        if (A + B) % 2 or any(A % p == 0 and B % p for p in PR): continue
        for C in range(-5999, 6000, 2):
            D = B * B - 4 * A * C
            if any(A % p and pow(D % p, (p - 1) // 2, p) != p - 1 for p in PR): continue
            if any(A % p == 0 and C % p == 0 for p in PR): continue
            cur = best = 0
            for y in range(-Y, Y + 1):
                cur = cur + 1 if isprime(abs(A * y * y + B * y + C)) else 0
                best = max(best, cur)
            if best >= 42: out.append((A, B, C, best))
print(out)
assert (36, 18, -1801, 45) in out and (47, 9, -5209, 43) in out and max(r[3] for r in out) == 45
print('twin agrees with qs.c on A <= 50, |C| <= 6000')
