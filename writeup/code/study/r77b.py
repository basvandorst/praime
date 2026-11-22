"""Round 77 analysis (after the run): per target prime q <= y, the class chains' mean reveal time T_q(y) against
W(y) = sum_{q'>y} 1/q'^2; least-squares line T_q = A_q + B_q W through the chains with y >= q (y = 7, 11, 13, 17),
its intercept A_q (the chain that holds every prime) against the data, and the line's predictions for y = 5 and 3."""
import os, json, math
import numpy as np
import common as C

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
p = C.small_primes(20_000_000).astype(float); X = p[-1]
W = lambda y: float(np.sum(1 / p[p > y] ** 2)) + 1 / (X * math.log(X))
runs = {7: 'r77_y7.json', 11: 'r77_y11.json', 13: 'r77_y13.json', 17: 'r76.json'}
R = {y: json.load(open(os.path.join(OUT, f))) for y, f in runs.items()}
d73 = json.load(open(os.path.join(OUT, 'r73.json')))['diag']
out = {'W': {str(y): W(y) for y in (3, 5, 7, 11, 13, 17)}, 'per_q': {}}
for q in (3, 5, 7, 11, 13):
    ys = [y for y in runs if y >= q]
    x = np.array([W(y) for y in ys]); t = np.array([R[y]['mean_T_q%d' % q] for y in ys])
    B, A = np.polyfit(x, t, 1) if len(ys) > 1 else (float('nan'), float('nan'))
    row = {'ys': ys, 'T': t.tolist(), 'A': float(A), 'B': float(B), 'data': d73['mean_T_q%d' % q],
           'max_residual': float(np.max(np.abs(t - (A + B * x)))) if len(ys) > 2 else None}
    for yp in (3, 5):
        if yp >= q:
            row['predicted_y%d' % yp] = float(A + B * W(yp))
    out['per_q'][str(q)] = row
    print(q, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()})
json.dump(out, open(os.path.join(OUT, 'r77b.json'), 'w'), indent=1)
