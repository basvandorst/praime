"""Round 35 (preregister/round35.txt): prime counts in [1e9, 2e9) by blocks of 1e8 against Riemann's explicit formula
with the first 2419 zeros."""
import os, json, math
import numpy as np
import mpmath as mp
import common as C
from sympy import mobius

mp.mp.dps = 30
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', 'data', 'study')
Z = [mp.mpf(l.strip()) for l in open(os.path.join(HERE, '..', '..', 'data', 'riemann', 'zeros_hp2419.txt')) if l.strip()]
T = float(Z[-1])
L0, B = int(os.environ.get('R35_L0', 10 ** 9)), int(os.environ.get('R35_B', 10 ** 8))
edges = [L0 + k * B for k in range(11)]


def R(y):
    return sum(mp.mpf(int(mobius(n))) / n * mp.li(mp.mpf(y) ** (mp.mpf(1) / n)) for n in range(1, 41) if mobius(n) != 0)


def zsum(y):
    ly = mp.log(y)
    return sum(2 * mp.re(mp.ei(mp.mpc(0.5, g) * ly)) for g in Z)


Rv = {e: R(e) for e in edges}
Zv = {e: zsum(e) for e in edges}
liv = {e: mp.li(e) for e in edges}
counts = []
for k in range(10):
    lo = edges[k]
    n = 0
    for s in range(10):
        n += int(C.primes_mask(lo + s * (B // 10), B // 10).sum())
    counts.append(n)


def sig(x):
    return (2 * math.sqrt(x) / math.log(x)) * math.sqrt(math.log(T / (2 * math.pi)) / (2 * math.pi * T))


out = {'T': T, 'blocks': []}
for k in range(10):
    a, b = edges[k], edges[k + 1]
    pred = float((Rv[b] - Rv[a]) - (Zv[b] - Zv[a]))
    dli = float(liv[b] - liv[a]); s = math.sqrt(2) * sig(b)
    d = {'a': a, 'count': counts[k], 'count_minus_li': counts[k] - dli, 'prediction_minus_li': pred - dli,
         'count_minus_prediction': counts[k] - pred, 'sigma': s, 'z': (counts[k] - pred) / s}
    out['blocks'].append(d)
    print(k, counts[k], 'count-li %.0f' % d['count_minus_li'], 'pred-li %.0f' % d['prediction_minus_li'], 'diff %.0f' % d['count_minus_prediction'], 'z %.2f' % d['z'], flush=True)
a, b = edges[0], edges[-1]
pred = float((Rv[b] - Rv[a]) - (Zv[b] - Zv[a])); dli = float(liv[b] - liv[a]); N = sum(counts); s = math.sqrt(2) * sig(b)
out['whole'] = {'count': N, 'count_minus_li': N - dli, 'prediction_minus_li': pred - dli, 'count_minus_prediction': N - pred, 'sigma': s, 'z': (N - pred) / s}
print('whole', out['whole'])
json.dump(out, open(os.path.join(OUT, 'r35.json' if L0 == 10 ** 9 else 'r35_%d.json' % L0), 'w'), indent=1)
