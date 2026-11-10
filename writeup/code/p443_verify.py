# Puzzle 443 / 393 Q1, independent check with Python integers only (no numpy, no floats): two-pointer scan over all
# windows of >= 2 consecutive primes with sum of k-th powers <= X, exact integer root test.
# Checks: (1) k = 3, X = 1e14: no cube.  (2) k = 2, X = 1e12: the 11 squares listed by p443_search.py, among them
# the page-128 example 586.  (3) Resta's and Fougeron's solutions of puzzle 128.  (4) p^2 = 1 (mod 24) for p > 3:
# a window of L primes > 3 has sum = L (mod 24), a square only for L = 0, 1, 4, 9, 12, 16 (mod 24).
# Run from the repo root: .venv/bin/python3 writeup/code/p443_verify.py
from sympy import primerange, integer_nthroot, isprime, nextprime
def scan(k, X):
    P = list(primerange(2, int((X / 2) ** (1 / k)) + 200)); out = []; n = 0
    for i in range(len(P) - 1):
        s = P[i] ** k; j = i
        while True:
            j += 1
            if j == len(P): break
            s += P[j] ** k
            if s > X: break
            n += 1
            r, ex = integer_nthroot(s, k)
            if ex: out.append((r, j - i + 1, P[i], P[j]))
        if j == i + 1 and s > X: break
    return n, sorted(out)
n, sol = scan(3, 10**14)
assert sol == [], sol
print(f"k = 3: {n} window sums <= 1e14, no cube")
n, sol = scan(2, 10**12)
roots = [r for r, L, a, b in sol]
assert roots == [586, 6088, 8174, 11585, 11707, 270106, 288818, 375661, 724909, 732910, 937423], roots
print(f"k = 2: {n} window sums <= 1e12, {len(sol)} squares, roots {roots}")
assert all(L % 24 in (0, 1, 4, 9, 12, 16) for r, L, a, b in sol if a > 3)
print("  window lengths mod 24:", [L % 24 for r, L, a, b in sol])
for A, first, last, L in ((2489647, 355363, 355951, 49), (27847739, 290209, 378193, 6913), (74930959, 654889, 797869, 10609)):
    s = 0; p = first; c = 0
    while p <= last: s += p * p; c += 1; p = nextprime(p)
    assert s == A * A and c == L and isprime(A) and L % 24 == 1
    print(f"  puzzle 128: {first}^2 + ... + {last}^2 ({L} primes, L = 1 mod 24) = {A}^2, {A} prime")
assert sorted({x * x % 24 for x in range(24)}) == [0, 1, 4, 9, 12, 16] and all(p * p % 24 == 1 for p in primerange(5, 10**4))
print("all checks passed")
