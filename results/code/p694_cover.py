# Is there a covering set for "A = (4*10^n-1)/3 or B = (16*10^n-1)/3 is divisible by a small prime" for all n?
from sympy import primerange, n_order
from math import lcm
info = []
for p in primerange(7, 200000):
    o = n_order(10, p)
    if o > 2000: continue
    i4, i16 = pow(4, -1, p), pow(16, -1, p)
    res = set()
    x = 1
    for n in range(o):
        if x == i4 or x == i16: res.add(n)
        x = x * 10 % p
    if res: info.append((p, o, res))
# test candidate moduli L: which residues n mod L are not hit by any prime with o | L
best = None
for L in [2, 4, 6, 8, 12, 16, 18, 24, 30, 36, 48, 60, 72, 90, 120, 144, 180, 240, 360, 720, 1440, 2520, 5040, 10080, 55440]:
    hit = [False] * L
    for p, o, res in info:
        if L % o: continue
        for r in res:
            for n in range(r, L, o): hit[n] = True
    free = [n for n in range(L) if not hit[n]]
    print(f"L={L:6d}: uncovered residues {len(free):5d} / {L}  (fraction {len(free)/L:.4f})" + (f"  e.g. {free[:8]}" if free else "  COVERED"))
