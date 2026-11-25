"""check_arcs_indep.py: independent recomputation (plain Python integers, sympy) of three numbers of check_arcs.py.
  1. Lemma 2 for h = 5, 13, 25, 65 via sympy totient/mobius and a direct loop (no numpy).
  2. the zone sums (m = 4) for the axis and for 3 + 4i over Gaussian primes with N z <= 10^7, by a different
     enumeration (bounding box of the zone, sympy.isprime), compared with the values printed by check_arcs.py.
  3. the partial singular series S_1000 from Pythagorean triples found by brute force (x^2 + y^2 = h^2).
"""
import math, cmath
from math import gcd
from sympy import isprime, totient, mobius

X = 10 ** 7; M = 3.0; m = 4


def fpos(x):
    return math.exp(-1.0 / x) if x > 0 else 0.0


def step(x):
    a, b = fpos(x), fpos(1.0 - x)
    return a / (a + b)


def F(y):
    return step(2.0 - abs(y))


def V(v):
    return step((v - 0.05) / 0.25) * step((1.0 - v) / 0.25)


def e(x):
    return cmath.exp(2j * math.pi * x)


ok = True
# 1. Lemma 2
for h in (5, 13, 25, 65):
    reps = [(x, y) for x in range(-h, h + 1) for y in range(-h, h + 1) if x * x + y * y == h * h and gcd(x, y) == 1]
    units = [(a, b) for a in range(h) for b in range(h) if gcd(a * a + b * b, h) == 1]
    for (wx, wy) in reps:
        tot = sum(e(-((a * wx + b * wy) % h) / h) for (a, b) in units) / len(units)
        good = abs(tot - int(mobius(h)) / int(totient(h))) < 1e-12
        ok &= good
    print(f"1. h={h}: {len(reps)} primitive w, |units| = {len(units)}, mean e(-Re(z conj w)/h) = mu/phi = {int(mobius(h)) / int(totient(h)):+.6f}: {'ok' if good else 'WRONG'}")

# 2. zone sums by bounding-box enumeration
ref = {"axis": complex(66816, -65624) + 1645, "3+4i": complex(29464, -6410)}   # check_arcs.py output (axis: inert added back)
for name, (wx, wy) in (("axis", (1, 0)), ("3+4i", (3, 4))):
    h = math.isqrt(wx * wx + wy * wy)
    ux, uy = wx / h, wy / h
    Z = 0j; npr = 0
    rmax = math.isqrt(X)
    smax = 2 * M * X ** 0.25 + 2          # |s| <= 2 M sqrt t <= 2 M X^(1/4)
    for a in range(1, rmax + 1):          # both zones lie in Re z > 0
        blo = math.floor((a * uy - smax) / ux); bhi = math.ceil((a * uy + smax) / ux)
        for b in range(blo, bhi + 1):
            n = a * a + b * b
            if n < 0.05 * X or n > X:
                continue
            t = a * ux + b * uy
            if t <= 0:
                continue
            s = b * ux - a * uy
            if abs(s) > 2 * M * math.sqrt(t):
                continue
            if b != 0 and a != 0:
                prime = isprime(n)
            else:
                q = abs(a) + abs(b)
                prime = isprime(q) and q % 4 == 3
            if not prime:
                continue
            w = V(n / X) * F(s / (M * math.sqrt(t)))
            if w == 0.0:
                continue
            npr += 1
            r = math.sqrt(n)
            Z += math.log(n) * w * ((a + 1j * b) / r) ** m * e(-r)
    rel = abs(Z - ref[name]) / abs(ref[name])
    # check_arcs prints rounded integers: allow 1 unit in each component
    good = abs(Z.real - ref[name].real) <= 1.5 and abs(Z.imag - ref[name].imag) <= 1.5
    ok &= good
    print(f"2. zone {name}: {npr} primes, Z = {Z.real:.2f} {Z.imag:+.2f}i vs check_arcs {ref[name].real:.0f} {ref[name].imag:+.0f}i: {'ok' if good else 'WRONG'}")

# 3. S_1000 from Pythagorean triples
S = 1.0          # h = 1: the four axes, (1/4) * 4
H = 1000
tri = {}
for x in range(1, H + 1):
    for y in range(1, H + 1):
        hh = x * x + y * y
        h = math.isqrt(hh)
        if h * h == hh and h <= H and gcd(x, y) == 1:
            tri.setdefault(h, []).append((x, y))
for h, lst in tri.items():
    mu = int(mobius(h))
    if mu == 0:
        continue
    # all primitive w with |w| = h: the four sign patterns of each (x, y) with x, y > 0
    s4 = sum(complex(sx * x, sy * y) ** m for (x, y) in lst for sx in (1, -1) for sy in (1, -1)) / h ** m
    S += (mu / int(totient(h))) * s4.real / 4
good = abs(S - 1.587324) < 2e-6
ok &= good
print(f"3. S_1000 from {sum(len(v) for v in tri.values())} triples (x, y > 0): {S:.6f} vs check_arcs 1.587324: {'ok' if good else 'WRONG'}")
print("ALL INDEPENDENT CHECKS PASSED" if ok else "SOME INDEPENDENT CHECK FAILED")
assert ok
