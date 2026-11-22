"""Round 53 addendum 3: is the bias linear in g numerical? The order-2 law for g = 70, 300, 700 with the coarse lam
grid of r53.py and with a fine grid (step 0.25 from 30 to ln 4e18), and interpolation of ln share in lam."""
import os, json, math
import numpy as np
import common as C
import r47
from r53 import L2

data = json.load(open(os.path.join(r47.OUT, 'oes_gap_counts_4e18.json')))['N']
M = C.primes_mask(2, 10 ** 7 - 2); P = np.append(np.nonzero(M)[0] + 2, 10000019); g0 = np.diff(P)
u = np.linspace(math.log(1e7), math.log(4e18), 40001)
LF = np.concatenate([np.arange(16.0, 30.0, 1.0), np.arange(30.0, 42.84, 0.25), [math.log(4e18)]])
out = {}
for g in (70, 300, 700):
    st = r47.setup(g); base = int(np.sum(g0 == g)); d = data[str(g)]
    row = {}
    for name, grid in (('coarse', L2), ('fine', LF)):
        ls = [math.log(r47.share(st, lam, 2)) for lam in grid]
        law = float(np.trapezoid(np.exp(np.interp(u, grid, ls)) * np.exp(u) / u, u)) + base
        row[name] = (d - law) / d
    out[str(g)] = row
    print(g, {k: '%+.4f %%' % (100 * v) for k, v in row.items()}, 'shift %+.4f %%' % (100 * (row['fine'] - row['coarse'])), flush=True)
json.dump(out, open(os.path.join(r47.OUT, 'r53e.json'), 'w'), indent=1)
