"""Round 100 (preregister/round100.txt): RMS of the swing of R(X, d) for the actual zeros against stationary GUE sets
(random origin), at d = 0.01, 0.005, 0.003. Uses r98.py's rms and GUE construction. usage: r100.py actual | gue s0 s1 | summary"""
import os, sys, json, math
import numpy as np
from scipy.linalg import eigh_tridiagonal
import r98

r98.DS = (0.01, 0.005, 0.003)
OUT = r98.OUT
PART = os.path.join(OUT, 'r100_parts.json')


def gue_stationary(seed, n=30_000):
    rng = np.random.default_rng(seed)
    dgl = rng.normal(0, math.sqrt(2), n)
    off = np.sqrt(rng.chisquare(2 * np.arange(n - 1, 0, -1)))
    ev = np.sort(eigh_tridiagonal(dgl, off, eigvals_only=True)) / math.sqrt(2 * n)
    R = float(np.max(np.abs(ev))) * (1 + 1e-9)
    x = np.clip(ev / R, -1, 1)
    u = (0.5 + (x * np.sqrt(1 - x * x) + np.arcsin(x)) / math.pi) * n
    o = u[10000] + rng.uniform(0, 1000)
    lo = float(np.interp(r98.Z[0], r98.t, r98.Nmean)); hi = float(np.interp(r98.Z[-1], r98.t, r98.Nmean))
    v = u - o; v = v[(v >= 0) & (v <= hi - lo)]
    return np.interp(lo + v, r98.Nmean, r98.t)


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode in ('actual', 'gue'):
        s0, s1 = (int(sys.argv[2]), int(sys.argv[3])) if mode == 'gue' else (0, 1)
        out = {}
        for s in range(s0, s1):
            g = r98.Z if mode == 'actual' else gue_stationary(s)
            out['%s_%d' % (mode, s)] = r98.rms_all(g)
        json.dump(out, open(PART + '.%s_%d' % (mode, s0), 'w')); print(mode, s0, s1, 'done', flush=True)
    else:
        res = {}
        for f in sorted(os.listdir(OUT)):
            if f.startswith('r100_parts.json.'):
                res.update(json.load(open(os.path.join(OUT, f))))
        summ = {}
        for d in map(str, r98.DS):
            g = np.array([v[d] for k, v in res.items() if k.startswith('gue')]); a = res['actual_0'][d]
            summ[d] = {'actual': a, 'n_gue': len(g), 'gue_mean': float(g.mean()), 'gue_sd': float(g.std(ddof=1)),
                       'z': float((a - g.mean()) / g.std(ddof=1)), 'rank_below': int(np.sum(g < a))}
            print(d, {k: (round(v, 5) if isinstance(v, float) else v) for k, v in summ[d].items()}, flush=True)
        json.dump({'sets': res, 'summary': summ}, open(os.path.join(OUT, 'r100.json'), 'w'), indent=1)
        for f in os.listdir(OUT):
            if f.startswith('r100_parts.json.'):
                os.remove(os.path.join(OUT, f))
