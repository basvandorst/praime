"""Bright and dark lines in the grid layouts of prime2.html, explained by Bateman-Horn; exponents on lines against a null.
For each layout: every lattice line in 24 directions (|dx|, dy <= 4, primitive), cut into stretches where the numbers along
the line follow one polynomial of degree <= 2 (constant second difference). For a stretch f(t), t = 0..m-1:
  observed primes; naive expectation sum 1/ln f(t) (random numbers of that size);
  Bateman-Horn prediction C(f) * naive, C(f) = prod_p (1 - delta_p)/(1 - 1/p) with delta_p the density of t with p | f(t)
  (p < 20000; delta_2 from t mod 8, odd p from the roots of 2f mod p).
Exponents: the largest number of exponents on one lattice line (any of the 24 directions) against 300 sets of 31
pseudo-exponents drawn from the primes <= N with weights w(p) (Wagstaff heuristic).
usage: lines.py layout [layout ...]  -> lines_out/<layout>.json and a printed summary"""
import sys, json, math, numpy as np
from math import gcd
D = 'data/'
EXP = [2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091]
N = 500000
comp = np.zeros(N + 1, dtype=bool); comp[:2] = True
for i in range(2, int(N ** 0.5) + 1):
    if not comp[i]: comp[i * i::i] = True
isp = ~comp
SP = [p for p in range(2, 20000) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
DIRS = [(dx, dy) for dy in range(0, 5) for dx in range(-4, 5) if (dy > 0 or dx > 0) and gcd(abs(dx), dy) == 1]
LMIN = 30
rng = np.random.default_rng(12345)

def load(lay):
    meta = json.load(open(D + lay + '.json'))
    raw = open(D + lay + '.pos.bin', 'rb').read()
    a = np.frombuffer(raw, dtype=np.uint16 if meta['pos16'] else np.int32).astype(np.int64).reshape(-1, 2)
    X, Y = a[1:, 0], a[1:, 1]
    n = np.arange(1, len(X) + 1)
    keep = (X < 60000) & (X >= 0)
    X, Y, n = X[keep], Y[keep], n[keep]
    if meta['dw'] == 2:
        j = Y // 2; r = X - j
        par = np.unique(r % 2); assert len(par) == 1, ('parity', par)
        X, Y = (r - par[0]) // 2, j
    return meta, X, Y, n

def bh_const(n0, d1, d2):
    """C(f) for f(t) = n0 + d1 t + d2 t(t-1)/2."""
    C = 1.0
    f = lambda t: n0 + d1 * t + d2 * t * (t - 1) // 2
    for p in SP:
        if p == 2:
            dens = sum(1 for t in range(8) if f(t) % 2 == 0) / 8
        else:
            A, B, Cc = d2 % p, (2 * d1 - d2) % p, (2 * n0) % p
            if A == 0:
                roots = (1 if B % p else (p if Cc % p == 0 else 0))
            else:
                disc = (B * B - 4 * A * Cc) % p
                roots = 1 if disc == 0 else (2 if pow(disc, (p - 1) // 2, p) == 1 else 0)
            dens = roots / p
        C *= (1 - dens) / (1 - 1 / p)
        if C == 0: return 0.0
    return C

def segments(X, Y, n):
    out = []
    for dx, dy in DIRS:
        L = dy * X - dx * Y; P = dx * X + dy * Y; s = dx * dx + dy * dy
        o = np.lexsort((P, L)); Ls, Ps, ns = L[o], P[o], n[o]
        step = (Ls[1:] == Ls[:-1]) & (Ps[1:] - Ps[:-1] == s)
        d1 = np.diff(ns)
        pair = step[1:] & step[:-1]
        d2 = d1[1:] - d1[:-1]
        if len(d2) == 0: continue
        newseg = np.ones(len(d2), dtype=bool)
        newseg[1:] = ~(pair[1:] & pair[:-1] & (d2[1:] == d2[:-1]))
        newseg |= ~pair
        sid = np.cumsum(newseg)
        valid = pair
        starts = np.flatnonzero(newseg)
        ends = np.append(starts[1:], len(d2)) - 1
        for a, b in zip(starts, ends):
            if not valid[a]: continue
            m = (b - a + 1) + 2
            if m < LMIN: continue
            vals = ns[a:a + m]
            out.append(((dx, dy), int(vals[0]), int(vals[1] - vals[0]), int(d2[a]), vals))
    return out

def analyse(lay):
    meta, X, Y, n = load(lay)
    segs = segments(X, Y, n)
    rows = []
    for (d, n0, d1, dd, vals) in segs:
        v = vals[vals >= 2]
        obs = int(isp[vals].sum()); naive = float((1 / np.log(np.maximum(v, 3))).sum())
        rows.append(dict(dir=d, n0=n0, d1=d1, d2=dd, m=len(vals), obs=obs, naive=naive, z=(obs - naive) / math.sqrt(max(naive, 1e-9))))
    rows.sort(key=lambda r: -r['z'])
    top = rows[:10]; bottom = sorted(rows, key=lambda r: r['z'])[:5]
    samp = [rows[i] for i in rng.choice(len(rows), size=min(150, len(rows)), replace=False)] if rows else []
    for r in top + bottom + samp:
        if 'C' not in r:
            r['C'] = bh_const(r['n0'], r['d1'], r['d2']); r['pred'] = r['C'] * r['naive']
            r['zbh'] = (r['obs'] - r['pred']) / math.sqrt(max(r['pred'], 1.0))
    # exponents on lattice lines
    E = np.array([k for k in EXP if k <= N]); pidx = np.flatnonzero(isp)
    w = np.array([1.0 if p <= 7 else min(1.0, math.exp(0.5772156649) * math.log2((6 if p % 4 == 1 else 2) * p) / p) for p in pidx])
    pos_of = np.full(N + 1, -1, dtype=np.int64); pos_of[n] = np.arange(len(n))
    def maxline(labels):
        best = 0; arg = None
        for dx, dy in DIRS:
            idx = pos_of[labels]; idx = idx[idx >= 0]
            L = dy * X[idx] - dx * Y[idx]
            u, c = np.unique(L, return_counts=True)
            if c.max() > best: best = int(c.max()); arg = ((dx, dy), int(u[c.argmax()]))
        return best, arg
    real, where = maxline(E)
    nulls = []
    for _ in range(300):
        lab = rng.choice(pidx, size=len(E), replace=False, p=w / w.sum())
        nulls.append(maxline(lab)[0])
    nulls = np.array(nulls); pval = float((1 + (nulls >= real).sum()) / (1 + len(nulls)))
    on = []
    if where:
        (dx, dy), Lv = where; idx = pos_of[E]; ok = idx >= 0
        on = [int(k) for k, i in zip(E[ok], idx[ok]) if dy * X[i] - dx * Y[i] == Lv]
    zs = [r['zbh'] for r in samp]
    res = dict(layout=lay, segments=len(rows), top=top, bottom=bottom,
               sample=dict(n=len(zs), mean=float(np.mean(zs)) if zs else None, sd=float(np.std(zs)) if zs else None,
                           frac_abs_gt3=float(np.mean(np.abs(zs) > 3)) if zs else None),
               exponents=dict(max_on_line=real, line=where, on=on, null_median=float(np.median(nulls)), null_95=float(np.quantile(nulls, .95)), p=pval))
    for r in res['top'] + res['bottom']:
        r['dir'] = list(r['dir'])
    json.dump(res, open(f'lines_out/{lay}.json', 'w'), default=float)
    return res

def fmt(r):
    a2 = r['d2'] / 2; b = r['d1'] - r['d2'] / 2
    poly = (f"{a2:g}t² + " if r['d2'] else "") + f"{b:g}t + {r['n0']}"
    return f"{poly:<28} dir {tuple(r['dir'])} len {r['m']:>4}: primes {r['obs']:>4}, random {r['naive']:7.1f}, Bateman-Horn {r['pred']:7.1f} (C {r['C']:.3f}), z vs BH {r['zbh']:+.1f}"

if __name__ == '__main__':
    import os; os.makedirs('lines_out', exist_ok=True)
    for lay in sys.argv[1:]:
        r = analyse(lay)
        print(f"== {lay}: {r['segments']} polynomial stretches of >= {LMIN} numbers")
        for x in r['top'][:5]: print('  bright ', fmt(x))
        for x in r['bottom'][:3]: print('  dark   ', fmt(x))
        s = r['sample']
        if s['n']: print(f"  random sample of {s['n']} stretches: z vs Bateman-Horn mean {s['mean']:+.2f}, sd {s['sd']:.2f}, |z| > 3: {100*s['frac_abs_gt3']:.1f}%")
        else: print(f"  no straight polynomial stretch of {LMIN} or more numbers in 24 directions")
        e = r['exponents']; print(f"  exponents: at most {e['max_on_line']} on one line {e['on']}; null median {e['null_median']:.0f}, 95% {e['null_95']:.0f}; p {e['p']:.3f}", flush=True)
