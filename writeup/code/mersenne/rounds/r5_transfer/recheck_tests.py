#!/usr/bin/env python3
"""Independent recomputation of the six frozen tests (and the control) from data.jsonl:
own feature code (brute-force x^2+2y^2, own Collatz and continued-fraction loops, exact
rational check of frac(n/2pi) and cos n with mpmath), own nulls, and the Poisson-binomial
by polynomial products (np.convolve tree) instead of the truncated DP. Uses only the
discovery/holdout family lists from results.json."""
import json, math, os
import numpy as np
import mpmath

mpmath.mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
res = json.load(open(os.path.join(HERE, "results.json")))
DISC, HOLD = set(res["disc"]), set(res["hold"])


def is_prime(n):
    if n < 2: return False
    i = 2
    while i * i <= n:
        if n % i == 0: return False
        i += 1
    return True


def jacobi(a, n):
    a %= n; r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: r = -r
        a %= n
    return r if n == 1 else 0


def feats(n):
    f = {}
    f["H1"] = n % 7 == 1
    ang = None
    if n % 8 in (1, 3):
        for x in range(1, math.isqrt(n) + 1, 2):
            y2 = n - x * x
            if y2 > 0 and y2 % 2 == 0 and math.isqrt(y2 // 2) ** 2 == y2 // 2:
                ang = float(mpmath.atan2(mpmath.sqrt(2) * math.isqrt(y2 // 2), x) / (mpmath.pi / 2)); break
    f["H2"] = ang is not None and ang > 0.75
    fr = mpmath.frac(mpmath.mpf(n) / (2 * mpmath.pi))
    f["H3"] = fr <= mpmath.mpf("0.103")
    f["H4"] = mpmath.cos(n) > mpmath.mpf("0.803")
    v, t = n, 0
    while v >= n:
        v = 3 * v + 1 if v & 1 else v >> 1; t += 1
    f["H5"] = t > 21
    # CF of sqrt(n) via P, Q recurrences
    a0 = math.isqrt(n); P, Q = 0, 1; qs = []
    while True:
        P = a0 * Q - P if not qs else qs[-1] * Q - P
        Q = (n - P * P) // Q
        a = (a0 + P) // Q
        if a == 2 * a0: break
        qs.append(a)
    # note: first iteration uses a0 as the previous quotient
    f["H6"] = (max(qs) if qs else 0) / math.sqrt(n) <= 0.397
    return f


def pb_upper(w, O):
    polys = [np.array([1 - p, p]) for p in w if p > 0] or [np.array([1.0])]
    while len(polys) > 1:
        polys = [np.convolve(polys[i], polys[i + 1]) if i + 1 < len(polys) else polys[i]
                 for i in range(0, len(polys), 2)]
    pmf = polys[0]
    return float(pmf[int(O):].sum()), float(pmf[:int(O) + 1].sum())


out = {}
for line in open(os.path.join(HERE, "data.jsonl")):
    d = json.loads(line)
    kind, b, K = d["kind"], d["b"], d["K"]; g = b - 1 if kind == "m" else b + 1
    name = f"{kind}{b}"
    rows = []
    for n, st, x in d["recs"]:
        if st == "ndiv": continue
        lnR = math.log((b ** n - 1) // (b - 1) if kind == "m" else (b ** n + 1) // (b + 1))
        if lnR < 64 * math.log(2): continue
        rows.append((n, st == "PRP", st in ("PRP", "fermat"), lnR))
    O = sum(r[1] for r in rows)
    # null C: class survival rates mod 24
    A = {}
    for c in set(r[0] % 24 for r in rows):
        rs = [r for r in rows if r[0] % 24 == c]
        A[c] = sum(r[2] for r in rs) / sum(math.log(2 * r[0]) / math.log(2 * K * r[0] + 1) for r in rs)
    wC = [A[r[0] % 24] * math.log(2 * r[0]) / r[0] for r in rows]
    wS = [r[2] * math.log(2 * K * r[0] + 1) / r[3] for r in rows]
    wG = [math.log(r[0]) / r[0] for r in rows]
    W = {}
    for nm, w in (("G", wG), ("C", wC), ("S", wS)):
        s = sum(w); W[nm] = [min(0.95, O / s * x) if s else 0.0 for x in w]
    F = [feats(r[0]) for r in rows]
    for r, fe in zip(rows, F):
        q = 2 * r[0] + 1
        fe["P0"] = is_prime(q) and g % q != 0 and jacobi((b if kind == "m" else -b) % q, q) == 1
    half = "D" if name in DISC else "H"
    for key in ("H1", "H2", "H3", "H4", "H5", "H6", "P0"):
        for nm in ("G", "C", "S"):
            acc = out.setdefault((key, nm, half), [0, 0.0, []])
            for r, fe, w in zip(rows, F, W[nm]):
                if fe[key]:
                    acc[0] += r[1]; acc[1] += w; acc[2].append(w)
for key in ("H1", "H2", "H3", "H4", "H5", "H6", "P0"):
    s = []
    for nm in ("G", "C", "S"):
        for half in ("D", "H"):
            O, E, w = out[(key, nm, half)]
            pu, pl = pb_upper(w, O)
            s.append(f"{nm}{half} {O}/{E:.1f} p+{pu:.3f}" + (f" p-{pl:.1e}" if key == "P0" else ""))
    print(key, " | ".join(s))
