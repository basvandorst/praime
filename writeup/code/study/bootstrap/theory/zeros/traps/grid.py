"""Rows g_j(u), j <= 24, on the scan grid u = i/1000 (i = 1..4000) and 4 + i/100 (i = 1..400) (numerics, arb).
Output grid.pkl: dict(u=[str], rows=[[str]*25]) with 100-digit midpoints (values are only used for signs
and ratios away from the double zeros).   Usage: python grid.py"""
import pickle, time
from fam import *
from flint import ctx
ctx.prec = 700
F = Fam(24, 0)
us = ["%d/1000" % i for i in range(1, 4001)] + ["%d/100" % (400 + i) for i in range(1, 401)]
rows = []
t0 = time.time()
for i, s in enumerate(us):
    p, q = s.split("/")
    u = arb(fmpq(int(p), int(q)))
    r = F.rows(u.exp())
    rows.append([r[j][0].mid().str(100, radius=False) for j in range(25)])
    if i % 500 == 0:
        print(i, s, "%.0fs" % (time.time() - t0), flush=True)
pickle.dump(dict(u=us, rows=rows), open("grid.pkl", "wb"))
print("done %.0fs" % (time.time() - t0))
