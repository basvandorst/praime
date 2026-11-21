import gmpy2
from gmpy2 import mpz
def fibmod(m, n):
    a, b = mpz(0), mpz(1)
    for bit in bin(m)[2:]:
        c = a * ((2 * b - a) % n) % n; d = (a * a + b * b) % n
        a, b = (d, (c + d) % n) if bit == '1' else (c, d)
    return a
F = [mpz(0), mpz(1)]
for _ in range(2, 1502): F.append(F[-1] + F[-2])
rows = {"k%4!=0": [0, 0], "k%4==0": [0, 0]}; odd = []
for k in range(3, 1501):
    n = F[k]
    if k % 3 == 0 or k % 5 == 0 or gmpy2.is_prime(n): continue      # F_k odd iff 3 does not divide k; 5 | F_k iff 5 | k
    eps = 1 if n % 5 in (1, 4) else -1
    fpp = fibmod(n - eps, n) == 0
    key = "k%4==0" if k % 4 == 0 else "k%4!=0"
    rows[key][0] += 1; rows[key][1] += fpp
    if (k % 4 == 0) == fpp: odd.append(k)
    if k < 60: print(k, eps, fpp, "  k | n - eps:", (n - eps) % k == 0)
print(rows, "exceptions to 'pseudoprime iff 4 does not divide k':", odd)
from sympy import isprime as sp
bad = []; cnt = 0
for k in range(3, 1501):
    n = F[k]
    if k % 3 == 0 or k % 5 == 0 or gmpy2.is_prime(n): continue
    eps = 1 if n % 5 in (1, 4) else -1
    fpp = fibmod(n - eps, n) == 0
    pred = sp(k) or (k % 2 == 0 and sp(k // 2))
    cnt += 1
    if fpp != pred or fpp != ((n - eps) % k == 0): bad.append((k, fpp, pred))
print("checked", cnt, "odd composite F_k (3, 5 not dividing k), k <= 1500; pseudoprime iff k prime or twice a prime; exceptions:", bad)
