"""Independent recount of the triples for p <= 10^6 with the clocks of check.py (sympy.n_order, matrix powers)."""
import json, sys
import sympy as sp
from check import clocks
js = json.load(open(sys.argv[1]))
lim = js["X"]
primes = set(sp.primerange(7, lim + 10))
ok = True
for key, v in js["res"].items():
    a, b = map(int, key.split(","))
    t = {k: 0 for k in v}
    for p in sorted(primes):
        if p > lim:
            break
        if p + a in primes and p + b in primes:
            fl = [clocks(q) for q in (p, p + a, p + b)]
            t["N"] += 1
            t["H%d" % sum(f[0] for f in fl)] += 1
            t["L%d" % sum(f[1] for f in fl)] += 1
    same = t == v
    ok &= same
    print("(p, p+%d, p+%d), p <= %d: %s %s" % (a, b, lim, t, "same" if same else "DIFFERENT %s" % v))
print("ALL CHECKS PASS" if ok else "CHECK FAILED")
