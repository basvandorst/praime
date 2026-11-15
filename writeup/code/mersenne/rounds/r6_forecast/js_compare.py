#!/usr/bin/env python3
"""Check forecast.js against the exact prime sums of s1e9.npz.

Points: every integer 2..300, the 52 exponents, 142869328, and 400 log-spaced integers in [300, 2.9e8].
Exact quantile: smallest prime q > p with S(q) - S(p) >= -ln(1-l) (exact below 10^6; above 10^6 the
grid of S at multiples of 1000 is interpolated, so 'exact' is good to < 1000).
"""
import json, math, subprocess, sys
import numpy as np

LEVELS = [0.05, 0.25, 0.5, 0.75, 0.95]
EXPS = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253,
        4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091,
        756839, 859433, 1257787, 1398269, 2976221, 3021377, 6972593, 13466917, 20996011,
        24036583, 25964951, 30402457, 32582657, 37156667, 42643801, 43112609, 57885161,
        74207281, 77232917, 82589933, 136279841]
d = np.load('s1e9.npz')
sq, sW, gx, gW = d['small_q'], d['small_SW'], d['grid_x'], d['grid_SW']

def S(x):
    if x <= 10**6:
        return float(sW[np.searchsorted(sq, x, 'right') - 1])
    j = x // 1000
    return float(gW[j - 1] + (gW[j] - gW[j - 1]) * (x % 1000) / 1000)

def quant(p, l):
    T = S(p) - math.log1p(-l)
    if T <= sW[-1]:
        return int(sq[np.searchsorted(sW, T, 'left')])
    j = int(np.searchsorted(gW, T, 'left'))
    return float(gx[j - 1] + 1000 * (T - gW[j - 1]) / (gW[j] - gW[j - 1]))

pts = sorted(set(list(range(2, 301)) + EXPS + [142869328] +
                 [int(x) for x in np.round(np.logspace(math.log10(300), math.log10(2.9e8), 400))]))
js = r"""
const { nextExponentForecast } = require(process.argv[2]);
const pts = JSON.parse(require('fs').readFileSync(0, 'utf8'));
const t0 = Date.now();
const out = pts.map(p => { const r = nextExponentForecast(p); return [r.lo90, r.lo50, r.median, r.hi50, r.hi90, r.exact]; });
console.log(JSON.stringify({ out, ms: Date.now() - t0 }));
"""
open('node_check.js', 'w').write(js)
res = subprocess.run(['node', 'node_check.js', './forecast.js'], input=json.dumps(pts), capture_output=True,
                     text=True, check=True, timeout=600)
R = json.loads(res.stdout)
lines = [f"forecast.js on {len(pts)} start points, node time {R['ms']} ms"]
exact_mism = 0; exact_n = 0; worst = {}; worst_ex_big = 0.0
for p, row in zip(pts, R['out']):
    isexact = row[5]
    for l, v in zip(LEVELS, row[:5]):
        ex = quant(p, l)
        if isexact:
            exact_n += 1
            if isinstance(ex, int):
                if v != ex: exact_mism += 1; lines.append(f"MISMATCH exact p={p} l={l}: js {v} py {ex}")
            else:
                worst_ex_big = max(worst_ex_big, abs(v - ex))
        else:
            rel = v / ex - 1
            key = 'p<1e7' if p < 1e7 else ('p<1e8' if p < 1e8 else 'p>=1e8')
            if abs(rel) > abs(worst.get((key, l), (0, 0))[0]): worst[(key, l)] = (rel, p)
lines.append(f"exact branch (p < 1e6): {exact_n} quantiles, {exact_mism} mismatches with the Python prime sums "
             f"(where the Python value lies above 1e6 it is interpolated; max |js - py| there {worst_ex_big:.0f})")
for key in ['p<1e7', 'p<1e8', 'p>=1e8']:
    lines.append(f"closed branch {key}: max relative error by level " +
                 ', '.join(f"{l}: {worst[(key, l)][0]:+.1e}" for l in LEVELS if (key, l) in worst))
# the table rows for the report
for p in [2, 127, 4423, 216091, 82589933, 136279841, 142869328]:
    row = R['out'][pts.index(p)]
    lines.append(f"p = {p}: js {row[:5]}  py {[round(quant(p, l)) for l in LEVELS]}")
# closed form below 1e6 (if a page wants no sieve at all): its error there
EG = math.exp(0.5772156649015329); C12 = math.log2(math.sqrt(12))
def closed(p, L):
    y = p * 2 ** (L / EG)
    for _ in range(60):
        f = EG * (math.log2(y / p) + C12 * math.log(math.log(y) / math.log(p))) - L
        y -= f / (EG / math.log(2) * (1 + C12 * math.log(2) / math.log(y)) / y)
    return y
for lo, hi in [(1e3, 1e4), (1e4, 1e5), (1e5, 1e6)]:
    m = 0.0
    for p in [int(x) for x in np.logspace(math.log10(lo), math.log10(hi), 60)]:
        for l in LEVELS:
            m = max(m, abs(closed(p, -math.log1p(-l)) / quant(p, l) - 1))
    lines.append(f"(closed form alone for p in [{lo:.0e}, {hi:.0e}): max relative error {m:.1e})")
print('\n'.join(lines))
open('js_compare.txt', 'w').write('\n'.join(lines) + '\n')
