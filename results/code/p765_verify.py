# Puzzle 765: p^a = q^b + r^c s^d over {2,3,5,7}, exponents >= 1.
from itertools import permutations
# (1) every congruence fact used in the two impossibility proofs
assert {pow(2, a, 7) for a in range(6)} == {1, 2, 4}                     # powers of 2 mod 7 = quadratic residues
assert all((pow(5, b, 7) in {1, 2, 4}) == (b % 2 == 0) for b in range(12))  # 5^b is a QR mod 7 iff b even
assert all((pow(3, a, 7) in {1, 2, 4}) == (a % 2 == 0) for a in range(12))  # 3^a is a QR mod 7 iff a even
assert all((pow(3, a, 7) == 2) == (a % 6 == 2) for a in range(30))       # 3^a = 2 mod 7 iff a = 2 mod 6
assert all((pow(3, a, 5) == 2) == (a % 4 == 3) for a in range(30))       # 3^a = 2 mod 5 iff a = 3 mod 4
assert all(pow(3, a, 5) == pow(2, (3 * a) % 4, 5) for a in range(30))    # 3 = 2^3 mod 5
assert all((pow(5, c, 7) == 4) == (c % 6 == 2) for c in range(30))       # 5^c = 4 mod 7 iff c = 2 mod 6
# 2^m - 3^n = 1 (m >= 1, n >= 0): only (1,0), (2,1); 3^A - 2^B = 1: only (1,1), (2,3)  [elementary, checked here to 200]
assert [(m, n) for m in range(1, 200) for n in range(0, 130) if 2**m - 3**n == 1] == [(1, 0), (2, 1)]
assert [(A, B) for A in range(1, 130) for B in range(1, 200) if 3**A - 2**B == 1] == [(1, 1), (2, 3)]
print("all congruence steps hold")

# (2) bounded search of all 12 equations (each unordered {r, s} counted once), exponents 1..E
E = 150
sols = {}
for p, q, r, s in permutations([2, 3, 5, 7]):
    if r > s: continue
    key = (p, q, r, s); sols[key] = []
    P = {p**a: a for a in range(1, 4 * E)}
    for b in range(1, E):
        qb = q**b
        for c in range(1, E):
            rc = r**c
            if rc > 10**400: break
            for d in range(1, E):
                v = qb + rc * s**d
                if v in P: sols[key].append((P[v], b, c, d))
                if rc * s**d > 10**400: break
for k, v in sols.items():
    print(k, v if v else "geen oplossing")
assert sols[(2, 5, 3, 7)] == [] and sols[(3, 2, 5, 7)] == []
