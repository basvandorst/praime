"""Round 11 check: exact Hardy-Littlewood inclusion-exclusion for the small discs (no truncation), against the
order-4 Bahadur values. P(no Gaussian prime in the disc | 0 prime) = sum_A (-1)^|A| S_G({0} u A)/S_G({0}) (kappa/ln N)^|A|."""
import math, json, os, time
from itertools import combinations
import numpy as np
import r11

out = {}
for R in (8, 10, 16):
    pts = [(x, y) for x in range(-10, 11) for y in range(-10, 11) if 0 < x * x + y * y < R and (x + y) % 2 == 0]
    P0 = max(4 * R, 50)
    zero = np.zeros((1, 1, 2), dtype=np.int64)
    I = np.array(pts, dtype=np.int64)
    s1 = r11.S_G(np.concatenate([np.repeat(zero, len(I), 0), I[:, None, :]], 1), P0)
    live = I[s1 > 0]; n = len(live)
    sums = []                                          # sums[k] = sum over k-subsets of S_G({0} u A)
    t0 = time.time()
    for k in range(0, n + 1):
        tot = 0.0
        it = combinations(range(n), k)
        while True:
            ch = np.array([x for _, x in zip(range(200000), it)], dtype=np.int64)
            if len(ch) == 0:
                break
            sets = np.concatenate([np.repeat(zero, len(ch), 0), live[ch]], 1) if k else np.repeat(zero, 1, 0)
            tot += float(r11.S_G(sets, P0).sum())
            if k == 0:
                break
        sums.append(tot)
    res = {}
    for box, (A, K) in r11.BOXES.items():
        lam_eff = math.log((A + K / 2) ** 2 * 2) / r11.KAPPA
        exact = sum((-1) ** k * sums[k] / lam_eff ** k for k in range(len(sums)))
        d = r11.Disc(R, 4)
        res[box] = {'exact': exact, 'order4': math.exp(d.log_avoid(lam_eff, 4)), 'order3': math.exp(d.log_avoid(lam_eff, 3))}
    out[R] = {'live': n, **res}
    print(R, n, round(time.time() - t0), res, flush=True)
json.dump(out, open(os.path.join(r11.OUT, 'r11_exact.json'), 'w'), indent=1)
