"""Hardy-Littlewood constants of the cluster shapes with the product over p < 10^7 plus the tail exp(-C(k,2) sum_{p>P} 1/p^2)
(approximated by C(k,2)/(P ln P)), and the parameter-free prediction C int_2^x dt / prod_i ln(t + h_i) against the counts
of clusters_1e9.txt at 10^8 and 10^9.   usage: clusters_check.py"""
import math
import numpy as np
from scipy.integrate import quad
P = 10 ** 7; s = np.ones(P + 1, bool); s[:2] = False
for i in range(2, int(P ** .5) + 1):
    if s[i]: s[i * i::i] = False
pr = np.nonzero(s)[0].astype(float)
PAT = {'twins': (0, 2), 'triplet': (0, 2, 6), 'quadruplet': (0, 2, 6, 8), 'quintuplet': (0, 2, 6, 8, 12), 'sextuplet': (0, 4, 6, 10, 12, 16)}
COUNT = {'twins': {8: 440312, 9: 3424506}, 'triplet': {8: 55600, 9: 379508}, 'quadruplet': {8: 4768, 9: 28388}, 'quintuplet': {8: 697, 9: 3633}, 'sextuplet': {8: 82, 9: 317}}
for n, pat in PAT.items():
    k = len(pat); small = pr[pr <= 20]; big = pr[pr > 20]; c = 1.0
    for p in small: c *= (1 - len(set(h % int(p) for h in pat)) / p) / (1 - 1 / p) ** k
    c *= float(np.exp(np.sum(np.log1p(-k / big) - k * np.log1p(-1 / big))))
    c *= math.exp(-(k * (k - 1) / 2) / (P * math.log(P)))
    out = []
    for e in (8, 9):
        x = 10.0 ** e
        i1 = quad(lambda t: 1 / math.log(t) ** k, 2, x, limit=500)[0]
        i2 = quad(lambda t: 1 / np.prod([math.log(t + h) for h in pat]), 2, x, limit=500)[0]
        out.append('10^%d: count %d, int from 2 %.4f, int of 1/prod ln(t+h) %.4f' % (e, COUNT[n][e], COUNT[n][e] / (c * i1), COUNT[n][e] / (c * i2)))
    print('%-11s C = %.6f | %s' % (n, c, ' | '.join(out)))
