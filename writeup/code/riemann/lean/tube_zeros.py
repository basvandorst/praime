import json, math, os, mpmath as mp
Z = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "weil_ref.json")))["zeros"]
gam = 0.57721566490153286061
def G(y):
    s = 0.0
    for n in range(2, int(y) + 1):
        m, p = n, None
        for q in range(2, int(n ** 0.5) + 1):
            if m % q == 0: p = q; break
        p = p or n
        while m % p == 0: m //= p
        if m == 1: s += math.log(p) * (y / n - 1)
    return s - y * math.log(y) + (1 + gam) * y
def triv(y): return sum(y ** (-2 * k) / (2 * k * (2 * k + 1)) for k in range(1, 40))
def zs(y):
    t = 0.0
    for g in Z:
        r = complex(0.5, g); t += 2 * (y ** r / (r * (1 - r))).real
    return t
T = Z[-1]; tail = math.log(T / (2 * math.pi)) / (math.pi * T)
print(f"{len(Z)} zeros up to {T:.1f}; truncation bound per sqrt(y) about {tail:.4f}")
for y in [1.5, 2, 2.5, 3, 5, 7.5, 10, 30, 100, 314.5, 1000, 5000]:
    g = G(y); lhs = g - math.log(2 * math.pi) + triv(y); rhs = zs(y)
    print(f"y={y:7}: (G - log2pi + trivial)/sqrt y = {lhs/math.sqrt(y):+.5f}  zero sum/sqrt y = {rhs/math.sqrt(y):+.5f}  diff {(lhs-rhs)/math.sqrt(y):+.5f}")
