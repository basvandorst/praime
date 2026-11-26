"""Semi-infinite LP: minimise cost(g)/g(0) over g = B(i d/du)(Phi*Phi), B(t) = sum_{j<=d} b_j t^{2j} >= 0,
class A: g(u) >= 0 on [0,inf);  class E: S(u) = sum_k e^{-ku/2} g(ku) >= 0 on (0,inf).
Variables: coefficients beta of B in the basis of orthonormal polynomials p_k(t^2) for (1/2pi) Xi(t)^2 dt
(Cholesky of the Hankel moment matrix m_{j+k} = g_{j+k}(0)); rows scaled to max-entry 1.
Constraints: u = i/400 (i=1..2400), u = log n (n<=403) [class E: S at these u and at log p];
B(t) >= 0 on t = i/50 (t<=80) and 80*1.02^i (i<400); tail sign: b_d >= 0 and (-1)^d b_d >= 0.
Cutting planes at the local minima of g (class A) / S (class E), located by Newton on g' (class A) or
on a fine local grid (class E), until the worst relative violation is > -5e-10.
Dual certificate: from the active rows, solve c = lam*m + sum mu_i row_i + sigma e_d in mpmath (60 digits);
if all mu_i >= 0 then lam is a LOWER bound for cost/g(0) over the whole class at this degree.
Usage: python lp_opt.py CLASS d1 d2 ...   (CLASS = A or E). Output: results_<CLASS>_d<d>.json"""
import sys, json, pickle, time
import numpy as np
import mpmath as mp
from scipy.optimize import linprog, nnls
from core import Basis, von_mangoldt

import os
mp.mp.dps = int(os.environ.get("DPS", 60))
DATA = pickle.load(open(os.environ.get("GRID", "grid.pkl"), "rb"))   # GRID=grid24.pkl for degrees up to 48
JMAX = DATA["J"]
DELTA = mp.mpf(1) / 400
IMAX = DATA["imax"]
G = [[mp.mpf(x) for x in row] for row in DATA["grid"]]          # G[i][j] = g_j(i*delta)
GL = {n: [mp.mpf(x) for x in row] for n, row in DATA["logn"].items()}
MOM = [mp.mpf(x) for x in DATA["mom"]]
NMAXL = max(GL)
LAM = {n: von_mangoldt(n) for n in range(2, NMAXL + 1)}
PRIMES = [p for p in range(2, NMAXL + 1) if all(p % q for q in range(2, int(p ** .5) + 1))]
EW = [mp.exp(-m * DELTA / 2) for m in range(IMAX + 1)]
UMAX = IMAX * DELTA


def S_grid_row(i, J):
    """class E sum at u = i*delta: sum_k e^{-k u/2} g_j(k u), terms with k u <= 6 (rest < e^{-2500} relative)."""
    v = [mp.mpf(0)] * (J + 1)
    k = 1
    while k * i <= IMAX:
        w = EW[k * i]
        row = G[k * i]
        for j in range(J + 1):
            v[j] += w * row[j]
        k += 1
    return v


def S_logp_row(p, J):
    v = [mp.mpf(0)] * (J + 1)
    q = p
    while q <= NMAXL:
        w = 1 / mp.sqrt(q)
        for j in range(J + 1):
            v[j] += w * GL[q][j]
        q *= p
    return v


def S_point_row(u, bas, J, r=0, vmax=None):
    """class E sum and its derivative at an arbitrary u (direct evaluation, k u <= vmax (default UMAX))."""
    v = [[mp.mpf(0)] * (J + 1) for _ in range(r + 1)]
    k = 1
    vmax = UMAX if vmax is None else vmax
    while k * u <= vmax:
        e = mp.exp(-k * u / 2)
        w = bas.eval(k * u)
        for j in range(J + 1):
            v[0][j] += e * w[j][0]
            if r >= 1:   # d/du [e^{-ku/2} g(ku)] = k e^{-ku/2} (g'(ku) - g(ku)/2)
                v[1][j] += k * e * (w[j][1] - w[j][0] / 2)
        k += 1
    return v


def to_float_row(vec):
    s = max(abs(x) for x in vec)
    if s == 0:
        return None
    return np.array([float(x / s) for x in vec])


def solve(cls, d, verbose=True):
    t0 = time.time()
    n = d + 1
    H = mp.matrix(n, n)
    for j in range(n):
        for k in range(n):
            H[j, k] = MOM[j + k]
    L = mp.cholesky(H)
    Li = mp.inverse(L)                      # p_k(s) = sum_j Li[k,j] s^j

    def tr(gvec):                           # transformed row: r_k = sum_j Li[k,j] g_j
        return [mp.fsum(Li[k, j] * gvec[j] for j in range(k + 1)) for k in range(n)]

    A, RM, LAB = [], [], []                 # float rows, mp rows (g_j basis), labels

    def add(row_mp, lab):
        fr = to_float_row(tr(row_mp))
        if fr is not None:
            A.append(fr); RM.append(row_mp); LAB.append(lab)

    if cls == "A":
        for i in range(1, IMAX + 1):
            add(G[i][:n], ("u", float(i * DELTA)))
        for m in range(2, NMAXL + 1):
            add(GL[m][:n], ("u", "log%d" % m))
    else:
        for i in range(1, IMAX + 1):
            add(S_grid_row(i, d), ("S", float(i * DELTA)))
        for p in PRIMES:
            add(S_logp_row(p, d), ("S", "log%d" % p))
    tg = [mp.mpf(i) / 50 for i in range(0, 50 * 80 + 1)] + [mp.mpf(80) * mp.mpf('1.02') ** i for i in range(1, 400)]
    for t in tg:
        s = t * t
        add([s ** j for j in range(n)], ("t", float(t)))
    cvec = [mp.mpf(0)] * n
    for m in range(2, NMAXL + 1):
        if LAM[m]:
            w = 2 * LAM[m] / mp.sqrt(m)
            for j in range(n):
                cvec[j] += w * GL[m][j]
    c_tr = tr(cvec)
    g0_tr = tr(G[0][:n])
    bounds = [(None, None)] * n
    if d > 0:
        bounds[d] = (0, None) if d % 2 == 0 else (0, 0)
    bas = Basis(d, 2)
    history = []
    for it in range(40):
        nsolved = len(A)
        res = linprog(np.array([float(x) for x in c_tr]), A_ub=-np.array(A), b_ub=np.zeros(nsolved),
                      A_eq=np.array([[float(x) for x in g0_tr]]), b_eq=[1.0], bounds=bounds,
                      method="highs", options=dict(primal_feasibility_tolerance=1e-10,
                                                    dual_feasibility_tolerance=1e-10))
        if res.status != 0:
            print("LP status", res.status, res.message, flush=True)
            if it == 0:
                return None
            res, nsolved = prev_res, prev_nsolved   # keep the last successful iterate
            history.append((it, "HiGHS failed: kept previous iterate", ""))
            break
        prev_res, prev_nsolved = res, nsolved
        beta = [mp.mpf(float(x)) for x in res.x]
        b = [mp.fsum(Li[k, j] * beta[k] for k in range(j, n)) for j in range(n)]
        g0 = mp.fsum(b[j] * MOM[j] for j in range(n))
        b = [x / g0 for x in b]
        cost = mp.fsum(b[j] * cvec[j] for j in range(n))
        if cls == "A":
            rows = [G[i][:n] for i in range(IMAX + 1)]
        else:
            rows = [None] + [S_grid_row(i, d) for i in range(1, IMAX + 1)]
        f = [mp.inf] + [mp.fsum(b[j] * rows[i][j] for j in range(n)) for i in range(1, IMAX + 1)]
        scale = [mp.mpf(1)] + [max(abs(b[j] * rows[i][j]) for j in range(n)) for i in range(1, IMAX + 1)]
        mins = [i for i in range(1, IMAX) if f[i] <= f[i - 1] and f[i] <= f[i + 1] and scale[i] > 0
                and f[i] / scale[i] < mp.mpf('1e-4')]
        new = 0
        worst = mp.mpf(0)
        touch = []
        for i in mins:
            u = i * DELTA
            if cls == "A":
                for _ in range(12):
                    v = bas.eval(u)
                    g1 = mp.fsum(b[j] * v[j][1] for j in range(n))
                    g2 = mp.fsum(b[j] * v[j][2] for j in range(n))
                    if g2 <= 0:
                        break
                    du = max(min(-g1 / g2, DELTA), -DELTA)
                    u += du
                    if abs(du) < mp.mpf(10) ** -22:
                        break
                v = bas.eval(u)
                row = [v[j][0] for j in range(n)]
            else:
                if u < mp.mpf('0.05'):
                    row = rows[i]
                else:   # secant on S' (few steps; each S evaluation costs ~6/u basis evaluations)
                    ua, ub = u - DELTA / 2, u + DELTA / 2
                    sa = S_point_row(ua, bas, d, 1); sb = S_point_row(ub, bas, d, 1)
                    fa = mp.fsum(b[j] * sa[1][j] for j in range(n)); fb = mp.fsum(b[j] * sb[1][j] for j in range(n))
                    for _ in range(6):
                        if fb == fa:
                            break
                        uc = ub - fb * (ub - ua) / (fb - fa)
                        if abs(uc - u) > DELTA:
                            break
                        sc = S_point_row(uc, bas, d, 1)
                        fc = mp.fsum(b[j] * sc[1][j] for j in range(n))
                        ua, fa, ub, fb = ub, fb, uc, fc
                        if abs(ub - ua) < mp.mpf(10) ** -18:
                            break
                    u = ub
                    row = S_point_row(u, bas, d, 0)[0]
            gv = mp.fsum(b[j] * row[j] for j in range(n))
            sc = max(abs(b[j] * row[j]) for j in range(n))
            rel = gv / sc
            worst = min(worst, rel)
            touch.append((float(u), float(rel), float(gv)))
            if rel < -mp.mpf(10) ** -13:
                add(row, ("u" if cls == "A" else "S", float(u))); new += 1
        history.append((it, mp.nstr(cost, 15), mp.nstr(worst, 3)))
        if verbose:
            print("  [%s d=%d it=%d] cost/g0=%s  worst rel=%s  new rows=%d  (%.1fs)" %
                  (cls, d, it, mp.nstr(cost, 15), mp.nstr(worst, 3), new, time.time() - t0), flush=True)
        if new == 0 or worst > -mp.mpf('5e-10'):
            break
        if len(history) >= 2 and abs(mp.mpf(history[-1][1]) - mp.mpf(history[-2][1])) < mp.mpf('1e-13') * abs(cost):
            break       # LP tolerance floor reached
    # ---------------- report ----------------
    out = dict(cls=cls, d=d, deg=2 * d, cost=mp.nstr(cost, 20), b=[mp.nstr(x, 30) for x in b],
               history=history, worst_rel_violation=mp.nstr(worst, 4), touch_points=touch)
    XM = pickle.load(open(os.environ.get("XIMOM", "xi_mom.pkl"), "rb"))
    aj = [mp.mpf(x) for x in XM["a"]]
    mj = [mp.mpf(x) for x in XM["m"]]
    arch = mp.fsum(b[j] * ((-mp.mpf(1) / 4) ** j / 2 + aj[j]) for j in range(n))
    g0xi = mp.fsum(b[j] * mj[j] for j in range(n))
    out["cost_arch"] = mp.nstr(arch / g0xi, 20)
    out["cost_prime_minus_arch"] = mp.nstr(cost - arch / g0xi, 5)
    bb = list(b)
    while len(bb) > 1 and abs(bb[-1]) < mp.mpf(10) ** -40:
        bb.pop()
    out["effective_deg"] = 2 * (len(bb) - 1)
    if len(bb) > 1:
        qs = mp.polyroots(list(reversed(bb)), maxsteps=800, extraprec=300)
        out["roots_t"] = [mp.nstr(mp.sqrt(r), 12) for r in qs]
    ts = np.linspace(0, 200, 40001)
    Bv = [(mp.fsum(b[j] * (mp.mpf(t) ** 2) ** j for j in range(n)), t) for t in ts]
    bm = min(Bv)
    out["Bmin_on_[0,200]"] = [mp.nstr(bm[0], 10), float(bm[1])]
    out["B(i/2)"] = mp.nstr(mp.fsum(b[j] * (-mp.mpf(1) / 4) ** j for j in range(n)), 20)
    out["g_logn_over_sqrtn"] = {m: mp.nstr(mp.fsum(b[j] * GL[m][j] for j in range(n)) / mp.sqrt(m), 10) for m in range(2, 13)}
    if cls == "E":
        out["S_logp"] = {p: mp.nstr(mp.fsum(b[j] * S_logp_row(p, d)[j] for j in range(n)), 10) for p in [2, 3, 5, 7, 11, 13]}
    out["cost_share"] = {m: mp.nstr(2 * LAM[m] / mp.sqrt(m) * mp.fsum(b[j] * GL[m][j] for j in range(n)) / cost, 6)
                         for m in range(2, 40) if LAM[m]}
    # ---------------- dual certificate (lower bound) ----------------
    marg = -res.ineqlin.marginals
    act = [k for k in range(nsolved) if marg[k] > 1e-14]
    out["active_rows"] = [(LAB[k][0], LAB[k][1], float(marg[k])) for k in act]
    # mp system: c = lam*m + sum_k mu_k RM[k] + sigma e_d ; support chosen by a small float LP (max lam, mu >= 0)
    cols = [MOM[:n]] + [RM[k] for k in act] + ([[mp.mpf(0)] * d + [mp.mpf(1)]] if d > 0 else [])
    kinds = ["lam"] + ["mu"] * len(act) + (["sig"] if d > 0 else [])
    M = mp.matrix(n, len(cols))
    for a_, col in enumerate(cols):
        for j in range(n):
            M[j, a_] = col[j]
    Mt = Li * M                              # orthonormal-basis coordinates (better conditioned)
    ct = Li * mp.matrix(cvec)
    cs = [max(abs(Mt[j, a_]) for j in range(n)) for a_ in range(len(cols))]
    Mf = np.array([[float(Mt[j, a_] / cs[a_]) for a_ in range(len(cols))] for j in range(n)])
    cf = np.array([float(ct[j]) for j in range(n)])
    bnds = []
    for kd in kinds:
        bnds.append((None, None) if kd == "lam" or (kd == "sig" and d % 2 == 1) else (0, None))
    obj = np.zeros(len(cols)); obj[0] = -1.0
    dl = linprog(obj, A_eq=Mf, b_eq=cf, bounds=bnds, method="highs")
    dual = dict(n=n, ncols=len(cols), lp_status=int(dl.status))
    if dl.status == 0:
        supp = [0] + [a_ for a_ in range(1, len(cols)) if abs(dl.x[a_]) > 1e-12]
        dual["support"] = [kinds[a_] for a_ in supp]
        Ms = mp.matrix(n, len(supp))
        for jj, a_ in enumerate(supp):
            for j in range(n):
                Ms[j, jj] = Mt[j, a_]
        if len(supp) == n:
            sol = mp.lu_solve(Ms, ct)
        else:   # least squares (over/under-determined); residual reported
            sol = mp.qr_solve(Ms, ct)[0] if len(supp) < n else mp.lu_solve(Ms.T * Ms, Ms.T * ct)
        resid = mp.norm(Ms * sol - ct)
        signs_ok = all((sol[jj] >= 0) or (kinds[a_] == "lam") or (kinds[a_] == "sig" and d % 2 == 1)
                       for jj, a_ in enumerate(supp))
        dual.update(lower_bound=mp.nstr(sol[0], 20), residual=mp.nstr(resid, 3), signs_ok=bool(signs_ok),
                    valid=bool(signs_ok and resid < mp.mpf(10) ** -40),
                    points=[LAB[act[a_ - 1]][1] for a_ in supp if kinds[a_] == "mu"])
    out["dual"] = dual
    json.dump(out, open("results_%s_d%d.json" % (cls, d), "w"), indent=1)
    print(json.dumps(dict((k, out[k]) for k in ["cls", "deg", "cost", "cost_prime_minus_arch", "dual", "touch_points",
                                                  "g_logn_over_sqrtn"]))[:3000], flush=True)
    return out


if __name__ == "__main__":
    cls = sys.argv[1]
    for d in map(int, sys.argv[2:]):
        solve(cls, d)
