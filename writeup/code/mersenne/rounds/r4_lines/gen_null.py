"""Draw Wagstaff-null sets conditioned on size, for three views:
full : primes < B, size 51        (all exponents)
page : primes < 500,000, size 31  (the page view)
h5   : primes < 100,000, size 28  (holdout selection set)
Output: null_<view>.npy, shape (N, K), each row sorted.
"""
import sys
import numpy as np
from common import load_primes, EXPS, B, PAGE, HERE
import os

N = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
p, pr = load_primes()
print("primes below B:", len(p), " expected exponents:", pr.sum())
for lim in (10**4, 10**5, PAGE, B):
    print(f"  E[count < {lim}] = {pr[p < lim].sum():.3f}  true = {sum(e < lim for e in EXPS)}")

rng = np.random.default_rng(20261115)
views = {"full": (B, 51), "page": (PAGE, 31), "h5": (10**5, 28)}
for name, (lim, K) in views.items():
    sel = p < lim
    pp, pq = p[sel], pr[sel]
    det = pp[pq >= 1.0]
    rnd = pq < 1.0
    rp = pp[rnd]
    mu = -np.log1p(-pq[rnd])
    cmu = np.cumsum(mu)
    tot = cmu[-1]
    need = K - len(det)
    out = np.empty((N, K), dtype=np.int64)
    got = 0
    tries = 0
    sizes = []
    while got < N:
        tries += 1
        M = rng.poisson(tot)
        idx = np.unique(np.searchsorted(cmu, rng.random(M) * tot, side="right"))
        idx = idx[idx < len(rp)]
        sizes.append(len(idx) + len(det))
        if len(idx) == need:
            out[got] = np.sort(np.concatenate([det, rp[idx]]))
            got += 1
    sizes = np.array(sizes)
    print(f"{name}: lim={lim} K={K} deterministic={det.tolist()} tries={tries} "
          f"unconditional size mean={sizes.mean():.2f} sd={sizes.std():.2f} "
          f"P(size=K)={N/tries:.4f}")
    np.save(os.path.join(HERE, f"null_{name}.npy"), out)
