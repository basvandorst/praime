"""Recheck of ms_extra.py: pure Python, sympy; reads chk_first.txt and ms_analyse_check.out."""
import sympy, math
first = {}
for line in open("chk_first.txt"):
    p, k = map(int, line.split()); first[p] = k
EX = [e for e in [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091] if e >= 100]
nb = sorted({sympy.prevprime(e) for e in EX} | {sympy.nextprime(e) for e in EX})
s = lambda p: first[p] == 0 or first[p] > 64
allp = sorted(first)
overall = sum(s(p) for p in allp) / len(allp)
def band_rate(p):
    lo, hi = (2, 10**4) if p < 10**4 else ((10**4, 10**5) if p < 10**5 else (10**5, 500001))
    m = [s(x) for x in allp if lo <= x < hi]; return sum(m) / len(m)
print(f"neighbours: {sum(s(p) for p in nb)}/{len(nb)} survivors; expected {len(nb)*overall:.1f} overall, {sum(band_rate(p) for p in nb):.1f} by band")
print("  non-survivors:", [(p, first[p]) for p in nb if not s(p)])
found = exp = var = 0.0; f2 = e2 = v2 = 0.0
for l in open("ms_analyse_check.out"):
    l = l.strip()
    if l.startswith("k ") and "admissible" in l:
        w = l.split()
        a, b = map(int, w[1].rstrip(":").split("-")); h = int(w[5]); e = float(w[7]); sd = float(w[9])
        found += h; exp += e; var += sd * sd
        if a >= 3 and b <= 64: f2 += h; e2 += e; v2 += sd * sd
print(f"all k: {int(found)} vs {exp:.1f}, z = {(found-exp)/math.sqrt(var):+.2f};  3..64: {int(f2)} vs {e2:.1f}, z = {(f2-e2)/math.sqrt(v2):+.2f}")
from fractions import Fraction
tail = 1 - sum(Fraction(1, math.factorial(j)) for j in range(6)) * Fraction(math.exp(-1))
print("Poisson >= 6 windows:", round(31848 * float(tail), 1))
# survivor fraction by p mod 8 at K = 64 (recheck of the mod-8 line in ms_analyse.out)
for r in (1, 3, 5, 7):
    m = [s(p) for p in allp if p % 8 == r]
    print(f"K=64 p={r} mod 8: n={len(m)} survivor frac {sum(m)/len(m):.4f}")
