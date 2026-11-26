"""Precompute g_j(u) = (-1)^j (Phi*Phi)^{(2j)}(u), j=0..J, on the uniform grid u = i*delta in [0, 6]
(delta = 1/400) and at u = log n (n = 2..403), plus the moments m_j = g_j(0) for j <= 2J.
Output: grid.pkl (mpmath values as strings).   Usage: python precompute.py [J] [dps] [outfile]   (2 processes)"""
import sys, pickle, time
from multiprocessing import Pool
import mpmath as mp
from core import Basis

J = int(sys.argv[1]) if len(sys.argv) > 1 else 16
dps = int(sys.argv[2]) if len(sys.argv) > 2 else 70
IMAX = 2400
NMAX = 403


def work(us):
    mp.mp.dps = dps
    bas = Basis(J, 0)
    out = []
    for u in us:
        v = bas.eval(mp.mpf(u) if not isinstance(u, tuple) else mp.log(u[1]))
        out.append([mp.nstr(v[j][0], dps - 5) for j in range(J + 1)])
    return out


if __name__ == "__main__":
    t0 = time.time()
    mp.mp.dps = dps
    pts = [mp.mpf(i) / 400 for i in range(IMAX + 1)] + [("log", n) for n in range(2, NMAX + 1)]
    chunks = [pts[k::8] for k in range(8)]
    with Pool(2) as pool:
        res = pool.map(work, chunks)
    vals = [None] * len(pts)
    for k in range(8):
        vals[k::8] = res[k]
    grid = vals[:IMAX + 1]
    logn = {n: vals[IMAX + 1 + n - 2] for n in range(2, NMAX + 1)}
    print("points done", time.time() - t0, flush=True)
    bm = Basis(2 * J, 0)
    v = bm.eval(0)
    mom = [mp.nstr(v[j][0], dps - 5) for j in range(2 * J + 1)]
    pickle.dump(dict(J=J, dps=dps, delta="1/400", imax=IMAX, grid=grid, logn=logn, mom=mom),
                open(sys.argv[3] if len(sys.argv) > 3 else "grid.pkl", "wb"))
    print("done", time.time() - t0)
