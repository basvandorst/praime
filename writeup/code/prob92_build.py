# Problem 92 construction: core {3,5,7,13} (classes of Gurisatti's numbers) + disjoint greedy tails, glued by CRT.
# Run: .venv/bin/python writeup/code/prob92_build.py 1 720720 20000000   (s = 1 Sierpinski, -1 Riesel; a few minutes)
# Build k with many covering sets: fixed core {3,5,7,13} + disjoint greedy tails, glued by CRT.
import sys, numpy as np
from sympy import primerange, n_order
from sympy.ntheory.modular import crt
s = int(sys.argv[1]); L = int(sys.argv[2]); PMAX = int(sys.argv[3])
kref = 218164505156603 if s == 1 else 16107442173648029           # Gurisatti's numbers fix the core classes
def cls_of(k, p, s):
    e = n_order(2, p); t = (-s * pow(k, -1, p)) % p; x = 1
    for r in range(e):
        if x == t: return e, r
        x = x * 2 % p
core = {p: cls_of(kref, p, s) for p in (3, 5, 7, 13)}
cov = np.zeros(L, bool)
for p, (e, r) in core.items(): cov[r::e] = True
R = ~cov
print("core", core, "uncovered density", R.mean())
pool = [int(p) for p in primerange(11, PMAX) if p not in core and pow(2, L, int(p)) == 1]
pool = [(p, n_order(2, p)) for p in pool]
pool.sort(key=lambda t: t[1])
used = set(); tails = []
while True:
    unc = R.copy(); tail = []
    avail = [(p, e) for p, e in pool if p not in used]
    while unc.any():
        best = None
        for p, e in avail:
            if p in [t[0] for t in tail]: continue
            counts = unc.reshape(-1, e).sum(axis=0)        # L divisible by e
            r = int(np.argmax(counts)); c = int(counts[r])
            if c and (best is None or c * 1.0 > best[0]): best = (c, p, e, r)
        if best is None: break
        _, p, e, r = best; tail.append((p, e, r)); unc[r::e] = False
    if unc.any(): break
    # prune to minimal
    changed = True
    while changed:
        changed = False
        for t in list(tail):
            m = np.zeros(L, bool)
            for t2 in tail:
                if t2 is not t: m[t2[2]::t2[1]] = True
            if (m | ~R).all(): tail.remove(t); changed = True; break
    tails.append(tail); used |= {t[0] for t in tail}
    print("tail", len(tails), [(p, e, r) for p, e, r in tail])
print("number of disjoint tails:", len(tails))
# CRT: k = -s * 2^(-r) mod p for every prime used; k odd
mods = [2]; res = [1]
allp = dict(core)
for tail in tails:
    for p, e, r in tail: allp[p] = (e, r)
for p, (e, r) in allp.items():
    mods.append(p); res.append((-s * pow(pow(2, r, p), -1, p)) % p)
k, M = crt(mods, res); k = int(k)
while k <= max(allp): k += int(M)
print("k =", k, "digits", len(str(k)))
open(f"k_{'S' if s == 1 else 'R'}.txt", "w").write(f"{k}\n{allp}\n{tails}\n")
