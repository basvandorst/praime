"""Exact (KKT) optimum of the semi-infinite LP at degree 2d = 4K, class A or E, by Newton in mpmath.

Structure (found by the LP, lp_opt.py): the optimum has K = d/2 touch points u_1..u_K where the constraint function
F (= g for class A, = S for class E) has a double zero, B > 0 strictly, b_d > 0.  For given u, the coefficients b
solve the square linear system   m.b = 1,  F(u_k) = 0,  F'(u_k) = 0   (2K+1 = d+1 equations).
Let y = A^{-T} c (A = that matrix, c = cost vector). Then c = y_0 m + sum_k y_k F(u_k) + sum_k z_k F'(u_k), and
d cost / d u_k = -z_k F''(u_k). Stationarity: z_k = 0, solved by Newton (central-difference Jacobian).
At the solution cost = lambda = y_0 and c = lambda m + sum mu_k F-row(u_k) with mu_k = y_k: if all mu_k >= 0 and the
primal is feasible (F >= 0 everywhere, B >= 0) this is the global optimum at this degree (LP duality); lambda is then
also a lower bound for every admissible g of degree <= 2d (up to the residual |z| ~ 1e-40).
Initial touch points: from results_<CLASS>_d<d>.json, or (option 'extend') from the degree-(d-2) optimum plus a scan
for one new point.   Usage: [GRID=grid24.pkl DPS=90] python polish.py CLASS d [extend|re]"""
import sys, json, time
import mpmath as mp
import numpy as np
import lp_opt as LO            # grid data, S rows
from core import Basis

import os
mp.mp.dps = int(os.environ.get("DPS", 60))
cls = sys.argv[1]
d = int(sys.argv[2])
extend = len(sys.argv) > 3 and sys.argv[3] == "extend"
K = d // 2
n = d + 1
bas = Basis(d, 2)
MOM = LO.MOM
cvec = [mp.mpf(0)] * n
for m in range(2, LO.NMAXL + 1):
    if LO.LAM[m]:
        w = 2 * LO.LAM[m] / mp.sqrt(m)
        for j in range(n):
            cvec[j] += w * LO.GL[m][j]

_cache = {}


def Frows(u):
    """rows F(u), F'(u) in the g_j basis (class A: g; class E: S)."""
    key = (cls, mp.nstr(u, 50))
    if key in _cache:
        return _cache[key]
    if cls == "A":
        v = bas.eval(u)
        r = ([v[j][0] for j in range(n)], [v[j][1] for j in range(n)])
    else:
        # terms with k u > log(e^u + 30) are below e^{-2 pi 30} (times polynomial factors) relative: dropped
        v = LO.S_point_row(u, bas, d, 1, vmax=mp.log(mp.exp(u) + 30))
        r = (v[0], v[1])
    _cache[key] = r
    return r


def system(us):
    A = mp.matrix(n, n)
    for j in range(n):
        A[0, j] = MOM[j]
    for k, u in enumerate(us):
        f, f1 = Frows(u)
        for j in range(n):
            A[1 + 2 * k, j] = f[j]
            A[2 + 2 * k, j] = f1[j]
    rhs = mp.matrix(n, 1); rhs[0] = 1
    b = mp.lu_solve(A, rhs)
    y = mp.lu_solve(A.T, mp.matrix(cvec))
    z = [y[2 + 2 * k] for k in range(len(us))]
    return b, y, z


def cost_of(b):
    return mp.fsum(b[j] * cvec[j] for j in range(n))


_ROWS = None


def feasibility(b):
    """min over the u-grid (i/400, i<=2400) of F/scale, plus Newton-refined local minima (class A)."""
    global _ROWS
    if _ROWS is None:
        if cls == "A":
            _ROWS = [LO.G[i][:n] for i in range(LO.IMAX + 1)]
        else:
            _ROWS = [None] + [LO.S_grid_row(i, d) for i in range(1, LO.IMAX + 1)]
    rows = _ROWS
    f = [mp.inf] + [mp.fsum(b[j] * rows[i][j] for j in range(n)) for i in range(1, LO.IMAX + 1)]
    sc = [mp.mpf(1)] + [max(abs(b[j] * rows[i][j]) for j in range(n)) for i in range(1, LO.IMAX + 1)]
    rel = [f[i] / sc[i] if sc[i] > 0 else mp.inf for i in range(LO.IMAX + 1)]
    i_min = min(range(1, LO.IMAX + 1), key=lambda i: rel[i])
    mins = [i for i in range(1, LO.IMAX) if f[i] <= f[i - 1] and f[i] <= f[i + 1]]
    return rel[i_min], float(i_min * LO.DELTA), [(float(i * LO.DELTA), mp.nstr(rel[i], 4)) for i in mins]


def newton(us, tol=None, maxit=30, verbose=True):
    tol = mp.mpf(10) ** -(mp.mp.dps - 20) if tol is None else tol
    us = [mp.mpf(u) for u in us]
    h = mp.mpf(10) ** -18
    for it in range(maxit):
        b, y, z = system(us)
        zn = max(abs(x) for x in z)
        if verbose:
            print("   newton it=%d cost=%s |z|=%s u=%s" % (it, mp.nstr(cost_of(b), 18), mp.nstr(zn, 3),
                                                       [mp.nstr(u, 12) for u in us]), flush=True)
        if zn < max(tol, mp.mpf(10) ** -20 * abs(cost_of(b))):
            break
        Jm = mp.matrix(K, K)
        for k in range(K):
            up = list(us); up[k] += h
            um = list(us); um[k] -= h
            zp = system(up)[2]; zm = system(um)[2]
            for i in range(K):
                Jm[i, k] = (zp[i] - zm[i]) / (2 * h)
        step = mp.lu_solve(Jm, -mp.matrix(z))
        mx = max(abs(step[k]) for k in range(K))
        lim = mp.mpf("0.02")
        if mx > lim:
            step = step * (lim / mx)
        us = [us[k] + step[k] for k in range(K)]
    return us, b, y, z


t0 = time.time()
if K == 0:
    raise SystemExit("d=0: nothing to polish")
if len(sys.argv) > 3 and sys.argv[3] == "re":      # re-polish from a stored optimum (e.g. at higher precision)
    P = json.load(open("polish_%s_d%d.json" % (cls, d)))
    us0 = [mp.mpf(u) for u in P["touch_points"]]
elif not extend:
    R = json.load(open("results_%s_d%d.json" % (cls, d)))
    us0 = sorted(t[0] for t in R["touch_points"])
    assert len(us0) == K, ("touch points", us0)
else:
    P = json.load(open("polish_%s_d%d.json" % (cls, d - 2)))
    old = [mp.mpf(u) for u in P["touch_points"]]
    best = None
    # the new touch point has always appeared beyond the old ones (checked for deg <= 32 with the full scan from 0.55)
    for unew in np.arange(float(max(old)) + 0.03, 3.6, 0.01):
        us = sorted(old + [mp.mpf(float(unew))])
        if min(abs(us[i + 1] - us[i]) for i in range(K - 1)) < 0.03:
            continue
        try:
            b, y, z = system(us)
        except ZeroDivisionError:
            continue
        cst = cost_of(b)
        fe = feasibility(b)[0] if cst > 0 else -1
        if fe > -mp.mpf("1e-3") and cst > 0 and (best is None or cst < best[0]):
            best = (cst, float(unew))
    print(" scan: best new point", best, flush=True)
    us0 = sorted(old + [mp.mpf(best[1])])
us, b, y, z = newton(us0, maxit=int(os.environ.get("MAXIT", 30)))
cost = cost_of(b)
lam = y[0]
mus = [y[1 + 2 * k] for k in range(K)]
fe = feasibility(b)
bb = list(b)
qs = mp.polyroots([bb[j] for j in range(n - 1, -1, -1)], maxsteps=800, extraprec=400)
real_pos = [r for r in qs if abs(mp.im(r)) < mp.mpf(10) ** -20 and mp.re(r) >= 0]
out = dict(cls=cls, d=d, deg=2 * d, K=K, cost=mp.nstr(cost, 25), lam=mp.nstr(lam, 25),
           mu=[mp.nstr(x, 10) for x in mus], all_mu_nonneg=bool(all(x >= 0 for x in mus)),
           stationarity_residual=mp.nstr(max(abs(x) for x in z), 3),
           touch_points=[mp.nstr(u, 30) for u in us], touch_points_exp=[mp.nstr(mp.exp(u), 12) for u in us],
           b=[mp.nstr(b[j], mp.mp.dps - 5) for j in range(n)],
           roots_t=[mp.nstr(mp.sqrt(r), 12) for r in qs], B_has_real_root=bool(len(real_pos) > 0),
           grid_min_rel=mp.nstr(fe[0], 5), grid_min_at=fe[1], grid_local_minima=fe[2],
           g_logn_over_sqrtn={m: mp.nstr(mp.fsum(b[j] * LO.GL[m][j] for j in range(n)) / mp.sqrt(m), 10)
                              for m in range(2, 13)},
           cost_share={m: mp.nstr(2 * LO.LAM[m] / mp.sqrt(m) * mp.fsum(b[j] * LO.GL[m][j] for j in range(n)) / cost, 6)
                       for m in range(2, 40) if LO.LAM[m]},
           time_s=time.time() - t0)
# arch-side check
import pickle
import os
XM = pickle.load(open(os.environ.get("XIMOM", "xi_mom.pkl"), "rb"))
aj = [mp.mpf(x) for x in XM["a"]]
mj = [mp.mpf(x) for x in XM["m"]]
arch = mp.fsum(b[j] * ((-mp.mpf(1) / 4) ** j / 2 + aj[j]) for j in range(n)) / mp.fsum(b[j] * mj[j] for j in range(n))
out["cost_arch"] = mp.nstr(arch, 25)
out["prime_minus_arch"] = mp.nstr(cost - arch, 4)
if cls == "E":
    out["S_logp"] = {p: mp.nstr(mp.fsum(b[j] * LO.S_logp_row(p, d)[j] for j in range(n)), 10) for p in [2, 3, 5, 7, 11, 13]}
json.dump(out, open("polish_%s_d%d.json" % (cls, d), "w"), indent=1)
print(json.dumps({k: out[k] for k in ["cls", "deg", "cost", "lam", "mu", "all_mu_nonneg", "stationarity_residual",
                                       "touch_points_exp", "grid_min_rel", "B_has_real_root", "prime_minus_arch"]}))
