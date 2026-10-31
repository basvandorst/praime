# Build a Brier number with disjoint covering sets A (for k*2^n+1) and B (for k*2^n-1) by CRT, and prove it.
import json, sys
from sympy.ntheory.modular import crt
from sympy import n_order, isprime
d = json.load(open(sys.argv[1]))
plus = {int(p): r for p, r in d["plus"].items()}
minus = {int(p): tuple(v) for p, v in d["minus"].items()}
A, B = set(plus), set(minus)
assert not (A & B) and all(isprime(p) for p in A | B)
mods, res = [2], [1]                                   # k odd
for p, a in plus.items():
    m = n_order(2, p); mods.append(p); res.append((-pow(2, -a, p)) % p)
for q, (m, b) in minus.items():
    assert n_order(2, q) == m
    mods.append(q); res.append(pow(2, -b, q) % q)
k, M = crt(mods, res); k = int(k)
print("k has", len(str(k)), "digits; k mod 3 =", k % 3)
# proof: every n is covered on both sides
Lp = 24; Lm = d["L"]
for n in range(Lp):
    assert any(n % n_order(2, p) == a and (k * pow(2, n, p) + 1) % p == 0 for p, a in plus.items()), n
for n in range(Lm):
    assert any(n % m == b and (k * pow(2, n, q) - 1) % q == 0 for q, (m, b) in minus.items()), n
print(f"every n mod {Lp} has a prime of A dividing k*2^n+1, every n mod {Lm} a prime of B dividing k*2^n-1")
print("A =", sorted(A)); print("|B| =", len(B), " B =", sorted(B))
# k > max prime, so k*2^n +- 1 is a proper multiple, hence composite, for all n >= 0
assert k > max(A | B)
print("k =", k)
open("brier69_k.txt", "w").write(str(k))
