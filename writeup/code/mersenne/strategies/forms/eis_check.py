"""Question 4, method B: every prime p <= N tested with PARI qfbsolve on the
form x^2 + x y + y^2, the solution moved into the wedge a >= b >= 0 by the
12 automorphisms of the form, and compared with eis.py (count, set, and the
drawn position of each prime).  Primes from sympy.primerange.
"""
import json

import cypari2
import numpy as np
from sympy import primerange

pari = cypari2.Pari()
N = 500_000


def orbit(x, y):
    """All images of (x, y) under the automorphism group of x^2+xy+y^2."""
    seen = {(x, y)}
    todo = [(x, y)]
    while todo:
        u, v = todo.pop()
        for w in ((-u, -v), (v, u), (u + v, -v)):
            if w not in seen:
                seen.add(w)
                todo.append(w)
    return seen


def main():
    ref = json.load(open("eis.json"))
    besta, bestb = np.load("eis_reps.npy")
    Q = pari.Qfb(1, 1, 1)
    primes = [int(p) for p in primerange(2, N + 1)]
    rep = set()
    pos_mismatch = 0
    orbit_sizes = {}
    for p in primes:
        s = pari.qfbsolve(Q, p)
        if len(s) == 0:
            continue
        x, y = int(s[0]), int(s[1])
        assert x * x + x * y + y * y == p
        orb = orbit(x, y)
        wedge = [(a, b) for (a, b) in orb if a >= b >= 0]
        orbit_sizes[len(orb)] = orbit_sizes.get(len(orb), 0) + 1
        assert len(wedge) == 1, (p, wedge)
        a, b = wedge[0]
        if (besta[p], bestb[p]) != (a, b):
            pos_mismatch += 1
        rep.add(p)
    res = {
        "primes_drawn": len(rep),
        "matches_eis.py_count": len(rep) == ref["primes_drawn"],
        "set_is_3_and_1mod3": rep == {p for p in primes if p == 3 or p % 3 == 1},
        "position_mismatches": pos_mismatch,
        "orbit_sizes": orbit_sizes,
    }
    print(res)
    json.dump(res, open("eis_check.json", "w"), indent=1)


if __name__ == "__main__":
    main()
