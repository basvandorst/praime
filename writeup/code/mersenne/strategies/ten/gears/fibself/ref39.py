"""Reference data for prime39.html, independent of the page: gmpy2 on the actual Fibonacci numbers."""
import json, os, gmpy2
from gmpy2 import mpz
from sympy import isprime, factorint
def fibmod(m, n):
    a, b = mpz(0), mpz(1) % n
    for bit in bin(m)[2:]:
        c = a * ((2 * b - a) % n) % n; d = (a * a + b * b) % n
        a, b = (d, (c + d) % n) if bit == '1' else (c, d)
    return a
def eps(n): return 1 if n % 5 in (1, 4) else -1
F = [mpz(0), mpz(1)]
for _ in range(2, 4200): F.append(F[-1] + F[-2])
grid = {}
for k in range(3, 1501):
    n = F[k]
    if k % 3 == 0: grid[k] = "even"; continue
    if k % 5 == 0: grid[k] = "five"; continue
    if gmpy2.is_prime(n): grid[k] = "prime"; continue
    grid[k] = "fools" if fibmod(n - eps(n), n) == 0 else "no"
fpp7 = [int(x) for x in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fibpsp_1e+08.txt")).read().split() if int(x) < 10**7]
ex = {}
for name, n in [("323", mpz(323)), ("4181", F[19]), ("17711", F[22]), ("F(31)", F[31]), ("F(4181)", F[4181]), ("F(49)", F[49]), ("1891", mpz(1891))]:
    ex[name] = {"digits": len(str(n)), "home": bool(fibmod(n - eps(n), n) == 0), "prime": bool(gmpy2.is_prime(n, 30))}
R = {"grid": grid, "fpp1e7": fpp7, "examples": ex,
     "fib_prime_k": [k for k, v in grid.items() if v == "prime"]}
json.dump(R, open("ref39.json", "w"))
from collections import Counter
print(Counter(grid.values())); print("Fibonacci-prime indices:", R["fib_prime_k"]); print(ex); print(len(fpp7), "FPP below 1e7")
