"""zeros 0 < gamma <= T2 + 1 of each L (PARI lfunzeros), completeness checks, gap re-search, deduplication. No phase sum here.
Checks: (1) drift of r_n = n - theta(gamma_n)/pi between the first and last 15 % of the zeros (a missed zero gives -1);
(2) the 30 largest unfolded gaps are searched again with divz = 64 (and gaps > 2.2 mean spacings with divz = 256);
(3) duplicates (|delta| < 1e-7) removed. usage: s3_zeros.py key [key ...]   (key 'Zeta' = Riemann zeta to 300)"""
import sys, os, time, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, cypari2
from ulib import *
pari = cypari2.Pari(); pari.allocatemem(700 * 10 ** 6)
def unfold(g, Q, mu): return theta(g, Q, mu) / math.pi
def drift(g, Q, mu):
    r = np.arange(1, len(g) + 1) - unfold(g, Q, mu); m = max(5, int(0.15 * len(g)))
    return float(r[-m:].mean() - r[:m].mean()), float(r.mean())
cat = {c[0]: c for c in CATALOGUE}
cat['Zeta'] = ('Zeta', ('zeta',), 300, 'zeta')
for key in sys.argv[1:]:
    t0 = time.time()
    _, spec, T2, fam = cat[key]
    L = pari.lfuncreate(1) if spec[0] == 'zeta' else build(pari, spec)
    Q, k, mu = params(pari, L)
    Tmax = T2 + 1.0
    L = pari.lfuninit(L, [1, int(Tmax) + 2])
    z = pari.lfunzeros(L, [0.01, Tmax])
    g = np.sort(np.array([float(x) for x in z]))
    t1 = time.time() - t0
    n0 = len(g); added = 0; log = []
    for it in range(3):
        u = unfold(g, Q, mu); gaps = np.diff(u)
        order = np.argsort(-gaps)[:30]
        new = []
        for i in order:
            dz = 256 if gaps[i] > 2.2 else 64
            zz = pari.lfunzeros(L, [g[i] + 1e-7, g[i + 1] - 1e-7], dz)
            new += [float(x) for x in zz]
        # also the first interval (0.01, g[0])
        zz = pari.lfunzeros(L, [0.01, g[0] - 1e-7], 64); new += [float(x) for x in zz]
        if not new: break
        g = np.sort(np.concatenate([g, new])); added += len(new); log.append(len(new))
    # dedupe
    keep = np.concatenate([[True], np.diff(g) > 1e-7]); ndup = int((~keep).sum()); g = g[keep]
    dr, rmean = drift(g, Q, mu)
    u = unfold(g, Q, mu); gaps = np.diff(u)
    info = dict(key=key, Q=Q, mu=mu, T2=T2, n_first=n0, added_by_gap_search=added, rounds=log, duplicates_removed=ndup,
                n=len(g), n_in_window=int(np.sum(g <= T2)), drift=dr, r_mean=rmean, max_gap=float(gaps.max()), min_gap=float(gaps.min()),
                gamma_min=float(g[0]), gamma_max=float(g[-1]), secs_lfunzeros=round(t1, 1), secs=round(time.time() - t0, 1))
    np.savez_compressed('zeros/zeros_%s.npz' % key, g=g, info=json.dumps(info))
    print('%-5s n=%d (first run %d, added %d, dup %d) drift=%+.3f r_mean=%+.3f max_gap=%.2f min_gap=%.4f g0=%.3f  lfunzeros %.0fs total %.0fs' % (
        key, len(g), n0, added, ndup, dr, rmean, gaps.max(), gaps.min(), g[0], t1, time.time() - t0), flush=True)
