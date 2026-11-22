"""First-order size of the interior terms dropped by the g >= 6 heuristic: same-offset interior pairs 210k apart
carry the same-class rigidity 2A/m (shifted pairs carry ~0, decomp.out). Extra deficit per cell
~ P * n_live * beta^2 * |2A/m|, beta = -rbar/(1 - rbar) (sensitivity of an empty-interior event to one live point)."""
import json, math, numpy as np
R = json.load(open('/home/bas/praime/writeup/data/study/r39.json'))
U = [a for a in range(210) if math.gcd(a, 210) == 1]
for w in ('train', 'A', 'B', 'C'):
    rbar = R[w]['m'] / R[w]['n']; twoA = abs(R[w]['vY'] - (1 - rbar))   # s1/m ~ 1 - rbar
    beta2 = (rbar / (1 - rbar)) ** 2
    num = den = extra = 0.0
    for g in range(2, 39, 2):
        rows = [a for a in U if math.gcd(a + g, 210) == 1]
        d = R[w][str(g)]; assert len(rows) == len(d['P'])
        for a, P, v in zip(rows, d['P'], d['per_row']):
            nl = sum(1 for dl in range(2, g, 2) if math.gcd(a + dl, 210) == 1)
            extra += P * nl * beta2 * twoA; num += v; den += 1 - P
    P40 = [p for p in R[w]['P_pooled40'] if p > 0]; num += sum(1 - p for p in P40); den += sum(1 - p for p in P40)
    print(w, 'rbar %.3f |2A/m| %.3f beta^2 %.3f implied %.3f -> with first-order interior terms %.3f' % (rbar, twoA, beta2, num / den, (num - extra) / den))
