"""Cell-by-cell attack on the certificate: for every a-interval I = [A0,A1] the certificate uses
   g(a) <= (G- + G+)/2 + (G+ - G-)/(2d) (a - am) + Haa+ d^2/2   on I.
We evaluate the referee's own g (refined maximisation of the referee's Psi) at the endpoints and at random interior
points of the intervals and look for a violation; and we compare a dense float sampling of Psi_aa (central second
differences of the referee's Psi) on every block with the certificate's Haa+.
Usage: ref_cells.py A|B [stride]"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fractions import Fraction as Fr
from ref_core import Psi_ball, Mixture

HERE = os.path.dirname(os.path.abspath(__file__))
part = sys.argv[1]; stride = int(sys.argv[2]) if len(sys.argv) > 2 else 1
rng = np.random.default_rng(2026)
if part == 'A':
    Psi = Psi_ball
    rows = json.load(open(os.path.join(HERE, 'dump_A.json')))
    A0 = np.array([r[0] for r in rows]); A1 = np.array([r[1] for r in rows])
    Gm = np.array([r[2] for r in rows]); Gp = np.array([r[3] for r in rows]); Haa = np.array([r[4] for r in rows])
    blk = np.arange(len(rows)) // 64
    Nb = 300; npair = 1
else:
    meas = json.load(open(os.path.join(HERE, '..', 'meas_K48.json')))
    mr = [Fr(x) for x in meas['radii']]; mw = [Fr(x) for x in meas['mu']]; nwv = [Fr(x) for x in meas['nu']]
    mix = Mixture([float(r) for r, w in zip(mr, mw) if w], [float(w) for w in mw if w],
                  [float(r) for r, w in zip(mr, nwv) if w], [float(w) for w in nwv if w])
    Psi = mix.Psi
    rows = json.load(open(os.path.join(HERE, 'dump_B.json')))
    A0 = np.array([r[0] for r in rows]); A1 = np.array([r[1] for r in rows])
    Gm = np.array([r[4] for r in rows]); Gp = np.array([r[5] for r in rows]); Haa = np.array([r[6] for r in rows])
    blk = np.array([r[7] for r in rows])
    Nb = 200; npair = len(mix.pairs)

def hfun(a, b):
    return Psi(a, b) - 0.5 * np.maximum(Psi(b, b), 0.0)

gr = (np.sqrt(5) - 1) / 2
def gvals(A):
    out = np.empty_like(A)
    CH = max(1, 150_000 // (Nb * npair))
    t = np.linspace(0, 1, Nb)
    for s in range(0, len(A), CH):
        a = A[s:s + CH]; n = len(a)
        bmax = np.where(a <= 1, a, np.minimum(a, a / np.maximum(a - 1, 1e-300)))
        B = bmax[:, None] * t[None, :]
        H = hfun(np.broadcast_to(a[:, None], B.shape), B)
        Hp = np.pad(H, ((0, 0), (1, 1)), constant_values=-np.inf)
        loc = (Hp[:, 1:-1] >= Hp[:, :-2]) & (Hp[:, 1:-1] >= Hp[:, 2:])
        Hl = np.where(loc, H, -np.inf)
        K = 3
        idx = np.argsort(-Hl, axis=1)[:, :K]
        rows_ = np.repeat(np.arange(n), K); ii = idx.ravel()
        ok = np.isfinite(Hl[rows_, ii]); rows_ = rows_[ok]; ii = ii[ok]
        lo = B[rows_, np.maximum(ii - 1, 0)]; hi = B[rows_, np.minimum(ii + 1, Nb - 1)]
        ar = a[rows_]
        x1 = hi - gr * (hi - lo); x2 = lo + gr * (hi - lo)
        f1 = hfun(ar, x1); f2 = hfun(ar, x2)
        for _ in range(45):
            m = f1 < f2
            lo = np.where(m, x1, lo); hi = np.where(m, hi, x2)
            nx1 = np.where(m, x2, hi - gr * (hi - lo)); nx2 = np.where(m, lo + gr * (hi - lo), x1)
            nf1 = np.where(m, f2, 0.0); nf2 = np.where(m, 0.0, f1)
            nf1[~m] = hfun(ar[~m], nx1[~m]); nf2[m] = hfun(ar[m], nx2[m])
            x1, x2, f1, f2 = nx1, nx2, nf1, nf2
        best = np.maximum(np.maximum(f1, f2), np.maximum(hfun(ar, lo), hfun(ar, hi)))
        gref = np.full(n, -np.inf); np.maximum.at(gref, rows_, best)
        out[s:s + CH] = np.maximum(np.maximum(H.max(axis=1), gref), 0.0)
    return out

t0 = time.time()
sel = np.arange(0, len(rows), stride)
am = (A0 + A1) / 2; d = A1 - am
def bound(i, a):
    return (Gm[i] + Gp[i]) / 2 + (Gp[i] - Gm[i]) / (2 * d[i]) * (a - am[i]) + Haa[i] * d[i] ** 2 / 2
# test points: left end, right end, two random interior
pts_i = np.concatenate([sel, sel, sel, sel])
u = np.concatenate([np.zeros(len(sel)), np.ones(len(sel)), rng.random(len(sel)), rng.random(len(sel))])
pts_a = A0[pts_i] + u * (A1[pts_i] - A0[pts_i])
gv = gvals(pts_a)
bd = bound(pts_i, pts_a)
marg = bd - gv
rel = marg / np.maximum(gv, 1e-300)
print(f"part {part}: {len(pts_a)} test points in {len(sel)} intervals, time {time.time() - t0:.0f}s")
print(f"  violations (bound < referee g): {np.sum(marg < 0)};  min relative margin {rel.min():.3e} at a = {pts_a[rel.argmin()]:.6g}")
print(f"  relative margin quantiles: 1% {np.quantile(rel, 0.01):.3e}, 50% {np.quantile(rel, 0.5):.3e}")
# margin relative to the Haa remainder term (how much of the remainder slack is actually used)
rem = Haa[pts_i] * d[pts_i] ** 2 / 2
print(f"  min (bound - g)/(Haa d^2/2) = {np.min(marg / np.maximum(rem, 1e-300)):.3e}")

# Haa+ vs dense sampling of Psi_aa on every block
t1 = time.time()
worst = -np.inf; wb = None; nviol = 0
nblk = blk.max() + 1
na_s, nb_s = (24, 300) if part == 'A' else (8, 120)
for k in range(nblk):
    m = np.where(blk == k)[0]
    if len(m) == 0: continue
    P, Q = A0[m[0]], A1[m[-1]]
    H = Haa[m[0]]
    a = rng.uniform(P, Q, na_s * nb_s)
    # Psi = 0 identically where the smeared charges are disjoint (Newton): ball: b >= a/(a-1); spheres: b >= a q/(a-r)
    # for every pair. Restrict sampling to the complement (beyond it the float second difference is pure rounding noise).
    if part == 'A':
        bsup = min(Q, P / (P - 1)) if P > 1 else Q
    else:
        bsup = min(Q, max(P * qq / (P - r) for r, qq, c in mix.pairs)) if P > 1 else Q
    b = rng.uniform(0, bsup, na_s * nb_s)
    # include points near the diagonal and near support boundaries by oversampling small b-ranges
    e = 1e-4 * a
    paa = (Psi(a + e, b) - 2 * Psi(a, b) + Psi(a - e, b)) / e ** 2
    r = paa.max() / H if H > 0 else (np.inf if paa.max() > 0 else 0)
    if paa.max() > H * (1 + 1e-6) + 1e-12:
        nviol += 1
    if r > worst:
        worst = r; wb = (P, Q, H, paa.max(), a[paa.argmax()], b[paa.argmax()])
print(f"  Haa+ test on {nblk} blocks: violations {nviol}; max sampled Psi_aa / Haa+ = {worst:.4f} (block {wb[0]:.4g}-{wb[1]:.4g}, Haa+ {wb[2]:.4g}, sampled {wb[3]:.4g} at a={wb[4]:.5g}, b={wb[5]:.5g}); {time.time() - t1:.0f}s")
