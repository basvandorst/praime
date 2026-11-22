"""Rank real biquadratic fields Q(sqrt a, sqrt b) (a, b squarefree, coprime, a = b = 1 mod 4, a < b <= BMAX) by the
least prime that has a degree-1 prime above it (split completely, or ramified with the other character +1).
usage: biq_fields.py BMAX"""
import sys, json
from sympy import primerange, factorint
BMAX = int(sys.argv[1]); PR = list(primerange(3, 400))
sqf = [n for n in range(5, BMAX + 1, 4) if all(e == 1 for e in factorint(n).values())]
leg = {p: [0] + [1 if pow(x, (p - 1) // 2, p) == 1 else -1 for x in range(1, p)] for p in PR}
from math import gcd
out = []
for i, a in enumerate(sqf):
    for b in sqf[i + 1:]:
        if gcd(a, b) > 1: continue
        if a % 8 == 1 and b % 8 == 1: least = 2
        else:
            least = None
            for p in PR:
                la, lb = leg[p][a % p], leg[p][b % p]
                if (la == 1 and lb == 1) or (la == 0 and lb == 1) or (lb == 0 and la == 1) or (la == 0 and lb == 0):
                    least = p; break
            least = least or 999
        out.append((least, a, b))
out.sort(key=lambda t: (-t[0], t[1] * t[2]))
print(len(out), 'fields; best:', out[:15])
json.dump(out[:400], open('biq_fields.json', 'w'))
