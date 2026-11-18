"""The page code of machine 23 (taken out of prime3.html by inherit_node.js) against inherit.py:
rounds, helpers, witnesses and the one division, for every wheel.
usage: inherit_js_check.py inherit_js.json"""
import json, sys
from sympy import primerange
from inherit import factor, helpers, one_cell, two_cells

J = json.load(open(sys.argv[1]))
tot = bad = 0
for N in (2, 3, 30, 250, 5000):
    rows = {r[0]: r for r in J[str(N)]['rows']}
    bad += sorted(rows) != list(primerange(2, N + 1))
    depth = {2: 0}
    lit = [2]
    for k in range(3, N + 1, 2):
        h = helpers(k, factor(k - 1, lit), factor(k + 1, lit), depth)
        if h is None:
            bad += k in rows
            continue
        g, F1, Q1, F2, Q2 = h
        w1 = one_cell(k, Q1) if Q1 else []
        w2 = two_cells(k, Q2) if Q2 else (0, [])
        ok = w1 not in (None, False) and w2 not in (None, False)
        x = 1
        while (x + 1) % F2:
            x += F1
        if ok and 1 < x < k and k % x == 0:
            ok = False
        tot += 1
        if not ok:
            bad += k in rows
            continue
        depth[k] = g + 1
        lit.append(k)
        want = [k, g + 1, F1, Q1, F2, Q2, [list(t) for t in w1], w2[0], [list(t) for t in w2[1]], x]
        if rows.get(k) != want:
            bad += 1
print('page code against Python for 2, 3, 30, 250 and 5,000 wheels: %d odd numbers, rounds, helpers, witnesses and division compared, %d differ' % (tot, bad))
