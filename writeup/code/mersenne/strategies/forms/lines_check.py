"""Question 5, method B: recheck every line reported in lines.json.

Independent route:
  - primes from sympy.primerange; the position of each prime from PARI
    qfbsolve (forms x^2 + P y^2 and x^2 + x y + y^2), moved into the drawn
    region R by the symmetries of the form (P=1: |x|,|y| sorted so b <= a;
    P=2,3: |x|,|y|; eisenstein: the 12 automorphisms, wedge a >= b >= 0);
  - primes per line counted in a Python dict over these positions;
  - lattice points on a line counted by a plain loop over the line;
  - Bateman-Horn constant with w(p) from the closed form
    (odd p not dividing the leading coefficient: 1 + (disc/p)), primes from
    sympy.primerange(2, 10001); sum of 1/ln n with math.log;
  - the top-5-by-count ranking over the four unit directions recomputed from
    the dict and compared with lines.json.
"""
import json
from math import gcd, isqrt, log

import cypari2
from sympy import primerange, legendre_symbol

pari = cypari2.Pari()
N = 500_000
SMALLP = list(primerange(2, 10_001))


def Qf(name, a, b):
    return a * a + a * b + b * b if name == "eis" else a * a + name * b * b


def in_R(name, a, b):
    if a < 0 or b < 0 or Qf(name, a, b) > N:
        return False
    if name == 1 or name == "eis":
        return b <= a
    return True


def eis_orbit(x, y):
    seen, todo = {(x, y)}, [(x, y)]
    while todo:
        u, v = todo.pop()
        for w in ((-u, -v), (v, u), (u + v, -v)):
            if w not in seen:
                seen.add(w)
                todo.append(w)
    return seen


def positions(name):
    Q = pari.Qfb(1, 1, 1) if name == "eis" else pari.Qfb(1, 0, name)
    pos = {}
    for p in primerange(2, N + 1):
        s = pari.qfbsolve(Q, p)
        if len(s) == 0:
            continue
        x, y = int(s[0]), int(s[1])
        if name == "eis":
            cand = [(a, b) for (a, b) in eis_orbit(x, y) if a >= b >= 0]
            assert len(cand) == 1
            a, b = cand[0]
        else:
            a, b = abs(x), abs(y)
            if name == 1 and b > a:
                a, b = b, a
        assert in_R(name, a, b) and Qf(name, a, b) == p
        pos[p] = (a, b)
    return pos


def omega(A, B, C, p):
    if p == 2:
        return sum(1 for t in (0, 1) if (A * t * t + B * t + C) % 2 == 0)
    if A % p == 0:
        return (1 if B % p else (p if C % p == 0 else 0))
    D = (B * B - 4 * A * C) % p
    return 1 if D == 0 else 1 + int(legendre_symbol(D, p))


def line_points(name, da, db, key):
    """All lattice points (a, b) of R with db*a - da*b = key (plain loop)."""
    pts = []
    M = isqrt(N)
    for a in range(0, M + 1):
        if da == 0:
            if db * a != key:
                continue
            bs = range(0, M + 1)
        elif db == 0:
            if key % da:
                continue
            bs = [-key // da]
        else:
            num = db * a - key
            if num % da:
                continue
            bs = [num // da]
        for b in bs:
            if in_R(name, a, b):
                pts.append((a, b))
    return pts


def main():
    ref = json.load(open("lines.json"))
    bad = 0
    out = {}
    for name in (1, 2, 3, "eis"):
        r = ref[str(name)]
        pos = positions(name)
        assert len(pos) == r["primes_in_R"], (len(pos), r["primes_in_R"])
        # ranking by count over the four unit directions
        counts = {}
        for p, (a, b) in pos.items():
            for (da, db) in ((1, 0), (0, 1), (1, 1), (-1, 1)):
                k = (da, db, db * a - da * b)
                counts[k] = counts.get(k, 0) + 1
        top = sorted(counts.items(), key=lambda kv: -kv[1])[:5]
        ref_top = [(d["dir"][0], d["dir"][1], d["key"], d["primes"]) for d in r["top_by_count"]]
        same_top = sorted(c for _, c in top) == sorted(x[3] for x in ref_top) and \
            all(counts.get((x[0], x[1], x[2])) == x[3] for x in ref_top)
        bad += not same_top
        res = {"primes_in_R": len(pos), "top5_counts_agree": same_top, "lines": []}
        allp = set(pos)
        seen = set()
        for tag in ("top_by_count", "top_by_contrast", "top_by_contrast_any_direction", "named"):
            items = r[tag].values() if tag == "named" else r[tag]
            for d in items:
                da, db = d["dir"]
                k = d["key"]
                if (da, db, k) in seen:
                    continue
                seen.add((da, db, k))
                pts = line_points(name, da, db, k)
                ns = [Qf(name, a, b) for (a, b) in pts]
                nprime = sum(1 for n in ns if n in allp)
                nprime2 = sum(1 for (a, b) in pos.values() if db * a - da * b == k)
                # polynomial along the line through its first point, step (da, db)
                a0, b0 = min(pts, key=lambda ab: (ab[1], ab[0]) if db else ab)
                A = Qf(name, da, db)
                C = Qf(name, a0, b0)
                B = Qf(name, a0 + da, b0 + db) - A - C
                Cbh = 1.0
                for p in SMALLP:
                    Cbh *= (1 - omega(A, B, C, p) / p) / (1 - 1 / p)
                S = sum(1 / log(n) for n in ns if n >= 2)
                ok = (len(pts) == d["points"] and nprime == d["primes"] == nprime2
                      and abs(Cbh - d["C_BH"]) < 1e-9 * max(1, Cbh)
                      and abs(S - d["sum_1_over_ln"]) < 1e-6)
                bad += not ok
                res["lines"].append({"line": d["line"], "points": len(pts), "primes": nprime,
                                     "C": Cbh, "pred": Cbh * S, "ok": ok})
                print(f"{name} {d['line']:>15s} pts={len(pts)} primes={nprime} "
                      f"C={Cbh:.4f} pred={Cbh*S:.1f} {'OK' if ok else 'MISMATCH'}", flush=True)
        out[str(name)] = res
    print("total mismatches:", bad)
    json.dump(out, open("lines_check.json", "w"), indent=1)


if __name__ == "__main__":
    main()
