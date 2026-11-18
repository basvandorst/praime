"""The page code of machines 20, 21 and 22 (taken out of prime3.html by early_node.js) against early.py.
usage: early_js_check.py early_js.json"""
import json, sys
from math import isqrt
from sympy import primepi, primerange, isprime
from early import self_feed, lucy

J = json.load(open(sys.argv[1]))
bad = 0
n_checks = 0
for N in (2, 3, 30, 250, 1000, 5000):
    R = J[str(N)]
    # machine 21: gears, frontiers, new lamps, lamps, tick of every wheel
    rows = []
    done = False
    primes = list(primerange(2, N + 1))
    lit = set()
    dec = {}
    t = 0
    frontier = 2
    while frontier <= N:
        g = primes[t]
        t += 1
        q = g + 1
        while not isprime(q):
            q += 1
        new_frontier = min(q * q, N + 1)
        new = [m for m in range(frontier, new_frontier) if isprime(m)]
        for m in range(frontier, new_frontier):
            dec[m] = t
        lit.update(new)
        rows.append((g, q, q * q, len(new)))
        frontier = new_frontier
    sf = R['sf']
    n_checks += 6
    bad += sf['T'] != len(rows)
    bad += sf['g'][1:] != [r[0] for r in rows]
    bad += sf['fr'][1:] != [r[2] for r in rows]
    bad += sf['nw'][1:] != [r[3] for r in rows]
    bad += sf['lit'] != sorted(lit)
    bad += sf['dec'] != [dec[k] for k in range(2, N + 1)]
    # machine 20: tick of every wheel
    n_checks += 2
    bad += R['sq']['dec'] != [k.bit_length() - 1 for k in range(2, N + 1)]
    bad += R['sq']['lims'] != [1] + [2 ** (t + 1) - 1 for t in range(1, 14)]
    # machine 22: counters
    V, C, moves, pr = lucy(N * N)
    cn = R['cn']
    n_checks += 5
    bad += cn['V'] != V
    bad += cn['C'] != C
    bad += cn['moves'] != moves
    bad += cn['pr'] != list(primerange(2, N + 1))
    bad += cn['pi'] != int(primepi(N * N))
    for t, v, c, mv in cn['ex']:
        tt = min(t, N)
        want = N * N if t >= N else max(w for w in V if w < max(2, (t + 1) ** 2))
        n_checks += 2
        bad += v != want
        bad += c != int(primepi(v))
print('page code against Python for 2, 3, 30, 250, 1,000 and 5,000 wheels: %d comparisons, %d differ' % (n_checks, bad))
bc = 0
for x, r in J['count']:
    V, C, moves, pr = lucy(x) if x <= 10 ** 9 else (None, None, None, None)
    known = {10 ** 10: 455052511}
    pi = known[x] if x in known else int(primepi(x))
    if r['pi'] != pi or (moves is not None and (r['moves'] != moves or r['n'] != len(V) or r['gears'] != len(pr))):
        bc += 1
print('typed box of machine 22, x = %s: counts that differ from the true number of primes or from the moves of the Python twin: %d' % (', '.join('{:,}'.format(x) for x, r in J['count']), bc))
print('largest x: pi(%d) = %d in %d moves' % (J['count'][-1][0], J['count'][-1][1]['pi'], J['count'][-1][1]['moves']))
