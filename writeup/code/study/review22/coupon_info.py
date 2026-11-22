"""Own coupon-collector information curve (lambda -> infinity, uniform residues), to compare with r72.json shape_coupon."""
import math, json, numpy as np
QS = (3, 5, 7, 11, 13)
def Hm(q, m):
    p1 = (q - m - 1) / ((q - m) * (q - 1)); h = -m * (1 / (q - 1)) * math.log(1 / (q - 1))
    if p1 > 0: h -= (q - m) * p1 * math.log(p1)
    return h
def curve(q, K):
    H0 = Hm(q, 1); dist = np.zeros(q); dist[1] = 1; out = []
    for k in range(K + 1):
        out.append(H0 - sum(dist[m] * Hm(q, m) for m in range(1, q)))
        new = np.zeros(q)
        for m in range(1, q):
            up = (q - 1 - m) / (q - 1); new[m] += dist[m] * (1 - up)
            if m + 1 < q: new[m + 1] += dist[m] * up
        dist = new
    return np.array(out)
tot = sum(curve(q, 60) for q in QS); full = sum((q - 2) / (q - 1) * math.log((q - 1) / (q - 2)) for q in QS)
sh = tot / full
r = json.load(open('/home/bas/praime/writeup/data/study/r72.json'))
print('full own', tot[-1], 'vs formula', full)
print('max |own - study coupon shape|', max(abs(a - b) for a, b in zip(sh, r['shape_coupon'])))
print([round(sh[k], 3) for k in (1, 2, 3, 5, 10, 20, 40, 60)])
print('data', [round(r['shape_data'][k], 3) for k in (1, 2, 3, 5, 10, 20, 40, 60)])
print('L_k', [round(r['L_k'][k], 3) for k in (0, 1, 2, 3, 5, 10, 20, 40, 60)], round(r['L_true'], 3))
