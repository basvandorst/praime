"""Independent reference for prime18.html: the five detective statistics of the true primes in fixed strips [lo, lo + L)
(sympy primes; same definitions as spot.py).   usage: ref18.py out.json"""
import sys, json
import numpy as np
from sympy import primerange
small = np.array(list(primerange(2, 32)))
def stats(S, lo, hi):
    S = np.array(S, dtype=np.int64)
    div = int(np.sum(np.any(S[:, None] % small[None, :] == 0, axis=1) & (S > 31)))
    g = np.diff(S).astype(float); gap2 = float((g * g).mean() / g.mean() ** 2)
    d = S % 10; dd = d[np.isin(d, (1, 3, 7, 9))]; rep = float(np.mean(dd[1:] == dd[:-1]))
    def bask(h):
        edges = np.arange(lo, hi + 1, h); c = np.histogram(S, bins=edges)[0].astype(float); mid = (edges[:-1] + edges[1:]) / 2; e = h / np.log(mid)
        return float(np.var(c - e) / e.mean())
    return dict(div=div, gap2=round(gap2, 9), rep=round(rep, 9), bask=round(bask(100), 9), bask1k=round(bask(1000), 9))
out = {}
for lo in (505852, 844887, 5000000):
    for L in (10000, 100000):
        out['%d_%d' % (lo, L)] = stats(list(primerange(lo, lo + L)), lo, lo + L)
json.dump(out, open(sys.argv[1], 'w')); print(out)
