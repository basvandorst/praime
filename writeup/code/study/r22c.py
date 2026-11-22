"""Round 22 check: the Bahadur orders and Aitken's acceleration against the exact Hardy-Littlewood share, at lam = 17.8.
Shows that the exponentiated expansion converges to a value above the exact one (an overlap bias), not just slowly."""
import os, json, math
import theory2 as T2
from r22 import U_sums
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
lam = 17.8
out = {}
for h in (30, 42, 52, 60, 70):
    U, _ = U_sums(h)
    ex = sum((-1) ** j * U[j] / lam ** (j + 1) for j in range(len(U)))
    pt = T2.Pattern([0, h], 4)
    L = [pt.log_avoid(lam, o) for o in (1, 2, 3, 4)]
    sh = lambda la: pt.SB / lam * math.exp(la)
    ait = lambda a, b, c: c - (c - b) ** 2 / ((c - b) - (b - a))
    out[h] = {'exact': ex, 'order3_rel': sh(L[2]) / ex - 1, 'order4_rel': sh(L[3]) / ex - 1,
              'aitken234_rel': sh(ait(L[1], L[2], L[3])) / ex - 1}
    print(h, {k: ('%.3e' % v) for k, v in out[h].items()})
json.dump(out, open(os.path.join(OUT, 'r22c.json'), 'w'), indent=1)
