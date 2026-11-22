"""Round 46 (preregister/round46.txt): the uncertainty about the next gap that remains when p mod q is known for all
primes q <= y. K = index of the next prime among the y-rough numbers after p; predicted geometric with rho = hit(y)
(r45b.json), so H = h2(rho)/rho bits. usage: r46.py"""
import os, json, math
import numpy as np
import theory_r3 as R3
from r45 import spf_window

h2 = lambda r: 0.0 if r >= 1 else -(r * math.log2(r) + (1 - r) * math.log2(1 - r))


def run(L, W, rows):
    spf = spf_window(L, W)
    pidx = np.nonzero(spf == 0)[0]
    starts = pidx[:-1] + 1
    out = []
    for row in rows:
        y = row['y']
        rough_comp = (spf > y).astype(np.int32)                    # composites with spf > y (primes have spf 0)
        K = 1 + np.add.reduceat(rough_comp, starts)               # the segment ends at the next prime (0)
        v, c = np.unique(K, return_counts=True); p = c / c.sum()
        H = float(-np.sum(p * np.log2(p)) + (len(c) - 1) / (2 * c.sum() * math.log(2)))
        rho = row['corrected']
        pk = dict(zip(v.tolist(), p.tolist()))
        out.append({'y': y, 'u': row['u'], 'rho': rho, 'H_measured': H, 'H_geometric': h2(rho) / rho,
                    'P1': pk.get(1, 0.0), 'P2': pk.get(2, 0.0), 'P2_over_P1': pk.get(2, 0.0) / pk.get(1, 1.0),
                    'P2_over_P1_geometric': 1 - rho, 'meanK': float(K.mean()), 'meanK_geometric': 1 / rho})
        print('y %7d u %5.2f rho %.4f H %.3f geo %.3f  P2/P1 %.3f geo %.3f  E K %.2f geo %.2f' % (
            y, row['u'], rho, H, h2(rho) / rho, out[-1]['P2_over_P1'], 1 - rho, out[-1]['meanK'], 1 / rho), flush=True)
    return out


if __name__ == '__main__':
    r45b = json.load(open(os.path.join(R3.OUT, 'r45b.json')))
    res = {}
    for name, (L, W) in (('1e9', (10 ** 9, 5 * 10 ** 7)), ('1e12', (10 ** 12, 10 ** 7))):
        print(name, flush=True)
        res[name] = run(L, W, r45b[name])
    json.dump(res, open(os.path.join(R3.OUT, 'r46.json'), 'w'), indent=1)
