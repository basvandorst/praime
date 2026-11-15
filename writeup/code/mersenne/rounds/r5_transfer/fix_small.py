#!/usr/bin/env python3
"""Re-decide every record of data.jsonl with R < 2^64 by sympy.isprime (a proof below 2^64),
as gen.py now does. The first run marked n | b -/+ 1 composite ('ndiv'), wrong when R = n
(Wagstaff n = 3: R = 3). Records with log2 R >= 64 (the test universe) are untouched."""
import json, os, sympy
HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.join(HERE, "data.jsonl")
out = []; changed = []
for line in open(P):
    d = json.loads(line); kind, b = d["kind"], d["b"]
    for i, (n, st, x) in enumerate(d["recs"]):
        if st == "factor":
            continue
        R = (b ** n - 1) // (b - 1) if kind == "m" else (b ** n + 1) // (b + 1)
        if R < 2 ** 64:
            g = b - 1 if kind == "m" else b + 1
            new = "PRP" if sympy.isprime(R) else ("ndiv" if g % n == 0 else "fermat")
            if new != st:
                changed.append((kind, b, n, st, new)); d["recs"][i] = [n, new, x]
    out.append(json.dumps(d))
open(P, "w").write("\n".join(out) + "\n")
print("changed:", changed)
