"""Referee's float evaluation of g (Lemma 10), J = 4 pi int g(r) r^-5 dr and the bound
c = (3/2) (2 J D^2)^(1/3), for mu = nu = ball (case A) or the K = 48 sphere mixture (case B).
Usage: ref_J.py ball|K48 hx Nb
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fractions import Fraction as Fr
from ref_core import Psi_ball, Mixture

case = sys.argv[1]; hx = float(sys.argv[2]); Nb = int(sys.argv[3])
HERE = os.path.dirname(os.path.abspath(__file__))
if case == 'ball':
    Psi = Psi_ball
    D = 0.6
    amin, amax = 2e-4, 1e4
    kappa = (5 * 5 ** 0.5 - 11) / 4
else:
    meas = json.load(open(os.path.join(HERE, '..', 'meas_K48.json')))
    print("json keys:", list(meas.keys()))
    MU = meas
    mr = [Fr(x) for x in meas['radii']]
    mw = [Fr(x) for x in meas['mu']]
    nwv = [Fr(x) for x in meas['nu']]
    print("sum mu =", sum(mw), " sum nu =", sum(nwv), " #radii", len(mr))
    nr = [r for r, w in zip(mr, nwv) if w != 0]; nw = [w for w in nwv if w != 0]
    mr2 = [r for r, w in zip(mr, mw) if w != 0]; mw2 = [w for w in mw if w != 0]
    print("nu atoms:", [(str(r), str(w)) for r, w in zip(nr, nw)], " mu atoms:", len(mr2), " min mu radius", min(mr2), " max", max(mr2))
    mix = Mixture([float(r) for r in mr2], [float(w) for w in mw2], [float(r) for r in nr], [float(w) for w in nw])
    print("pairs with c != 0:", len(mix.pairs))
    Psi = mix.Psi
    # exact D(mu,mu) = 1/2 sum mu_i mu_j / max(r_i, r_j)
    Dex = sum(wi * wj / max(ri, rj) for ri, wi in zip(mr2, mw2) for rj, wj in zip(mr2, mw2)) / 2
    print("D(mu,mu) exact =", Dex, "=", float(Dex))
    D = float(Dex)
    amin, amax = 1e-4, 1e4
    # kappa_mu = max_b b^3 (1 - b V_mu(b e1)), V_sphere_r(b) = 1/max(b,r)
    bb = np.linspace(0, 1, 2_000_001)
    Vm = sum(float(w) / np.maximum(bb, float(r)) for r, w in zip(mr2, mw2))
    Lam = bb ** 3 * (1 - bb * Vm)
    i = Lam.argmax()
    # refine
    from scipy.optimize import minimize_scalar
    f = lambda x: -(x ** 3 * (1 - x * sum(float(w) / max(x, float(r)) for r, w in zip(mr2, mw2))))
    res = minimize_scalar(f, bounds=(bb[max(i - 2, 0)], bb[min(i + 2, len(bb) - 1)]), method='bounded', options={'xatol': 1e-14})
    kappa = -res.fun
    print(f"kappa_mu = max Lambda = {kappa:.15f} at b = {res.x:.12f}")

def hfun(a, b):
    d = Psi(b, b)
    return Psi(a, b) - 0.5 * np.maximum(d, 0.0)

gr = (np.sqrt(5) - 1) / 2

def g_chunk(A):
    """g(a) for an array A of a values."""
    n = len(A)
    bmax = np.where(A <= 1, A, np.minimum(A, A / np.maximum(A - 1, 1e-300)))
    t = np.linspace(0, 1, Nb)
    B = bmax[:, None] * t[None, :]
    AA = np.broadcast_to(A[:, None], B.shape)
    H = hfun(AA, B)
    gmax = H.max(axis=1)
    # candidate local maxima on the grid (incl. endpoints)
    Hp = np.pad(H, ((0, 0), (1, 1)), constant_values=-np.inf)
    loc = (Hp[:, 1:-1] >= Hp[:, :-2]) & (Hp[:, 1:-1] >= Hp[:, 2:])
    Hl = np.where(loc, H, -np.inf)
    K = 4
    idx = np.argsort(-Hl, axis=1)[:, :K]
    rows = np.repeat(np.arange(n), K)
    ii = idx.ravel()
    ok = np.isfinite(Hl[rows, ii])
    rows = rows[ok]; ii = ii[ok]
    lo = B[rows, np.maximum(ii - 1, 0)]; hi = B[rows, np.minimum(ii + 1, Nb - 1)]
    a_r = A[rows]
    # golden section (maximise)
    x1 = hi - gr * (hi - lo); x2 = lo + gr * (hi - lo)
    f1 = hfun(a_r, x1); f2 = hfun(a_r, x2)
    for _ in range(48):
        m = f1 < f2
        lo = np.where(m, x1, lo); hi = np.where(m, hi, x2)
        x1n = np.where(m, x2, hi - gr * (hi - lo)); x2n = np.where(m, lo + gr * (hi - lo), x1)
        f1n = np.where(m, f2, np.nan); f2n = np.where(m, np.nan, f1)
        # evaluate new points
        need1 = ~m; need2 = m
        x1, x2 = x1n, x2n
        e1 = hfun(a_r[need1], x1[need1]); e2 = hfun(a_r[need2], x2[need2])
        f1n[need1] = e1; f2n[need2] = e2
        f1, f2 = f1n, f2n
    best = np.maximum(f1, f2)
    # also evaluate interval endpoints exactly (endpoint maxima)
    best = np.maximum(best, np.maximum(hfun(a_r, lo), hfun(a_r, hi)))
    gref = np.full(n, -np.inf)
    np.maximum.at(gref, rows, best)
    g = np.maximum(gmax, gref)
    g = np.maximum(g, 0.0)
    return g, gmax

t0 = time.time()
if case == 'ball':
    cuts = [amin, 1/1024, 7/8, 100.0, amax]
else:
    cuts = [amin, 2.0**-12, 600.0, amax]
CH = max(1, int(200_000 // (Nb * (len(Psi.__self__.pairs) if case != 'ball' else 1))))
segs = []
Aall = []; Gall = []
for lo, hi in zip(cuts[:-1], cuts[1:]):
    L = np.log(hi / lo)
    n = int(np.ceil(L / hx / 2)) * 2          # even number of steps -> Simpson
    X = np.linspace(np.log(lo), np.log(hi), n + 1)
    A = np.exp(X)
    G = np.empty_like(A)
    for s in range(0, len(A), CH):
        G[s:s + CH], _ = g_chunk(A[s:s + CH])
    F = G * np.exp(-4 * X)
    h = X[1] - X[0]
    T1 = h * (F.sum() - 0.5 * (F[0] + F[-1]))
    T2 = 2 * h * (F[::2].sum() - 0.5 * (F[0] + F[-1]))
    S = h / 3 * (F[0] + F[-1] + 4 * F[1:-1:2].sum() + 2 * F[2:-1:2].sum())
    segs.append((lo, hi, n, T1, T2, S))
    Aall.append(A); Gall.append(G)
    print(f"segment [{lo:.6g}, {hi:.6g}]: {n} steps  trap(h) {T1:.13f}  trap(2h) {T2:.13f}  Simpson {S:.13f}   t {time.time()-t0:.1f}s")
A = np.concatenate(Aall); G = np.concatenate(Gall)
c6 = Gall[0][0] / amin ** 6
small = c6 * amin ** 2 / 2
tail = kappa / amax
print(f"g(amin)/amin^6 = {c6:.10f}; g(amax)/amax^3 = {Gall[-1][-1] / amax ** 3:.15f} vs kappa {kappa:.15f}")
print(f"below amin: {small:.3e}; beyond amax: kappa/amax = {tail:.12e}")
for k, name in [(3, 'trapezoid h'), (4, 'trapezoid 2h'), (5, 'Simpson')]:
    I = sum(sg[k] for sg in segs)
    J = 4 * np.pi * (I + small + tail)
    c = 1.5 * (2 * J * D * D) ** (1 / 3)
    print(f"{name:13s}: int = {I:.13f}  J = {J:.10f}  c = {c:.10f}")
np.savez(os.path.join(HERE, f'ref_g_{case}.npz'), A=A, G=G)
if case != 'ball':
    m = A >= 600
    r = G[m] / A[m] ** 3
    print(f"max over a in [600, {amax:g}] of g(a)/a^3 = {r.max():.15f} at a = {A[m][r.argmax()]:.3f}; kappa_mu = {kappa:.15f}")
print(f"total time {time.time() - t0:.1f}s")
