"""Triples: counted distribution of the number of even clocks against predict3.py (binomial z)."""
import json, math, sys
from predict3 import triple
js = json.load(open(sys.argv[1]))
zs = []
for key, v in js["res"].items():
    a, b = map(int, key.split(","))
    pr = triple(a, b, 10)
    n = v["N"]
    print("(p, p+%d, p+%d), p <= %d: %d triples" % (a, b, js["X"], n))
    for q in sorted(pr):
        p = float(pr[q])
        if p == 0:
            print("  %s predicted 0 (proved), counted %d %s" % (q, v[q], "exact" if v[q] == 0 else "FAIL"))
            continue
        z = (v[q] / n - p) / math.sqrt(p * (1 - p) / n)
        zs.append(z)
        print("  %s predicted %-6s = %.6f, counted %.6f, z %+.2f" % (q, pr[q], p, v[q] / n, z))
print("%d nonzero cells, max |z| = %.2f (cells of one pattern sum to 1, so they are correlated)" % (len(zs), max(map(abs, zs))))
