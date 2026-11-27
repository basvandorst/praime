"""Second reviewer: exact moment formulation (no grid in the local parameter) and own barrier solver.

Unknown per prime p: c_{p,l} = E cos(l theta), l = 1..L_p, of ONE probability measure on [0, pi]
(moment set = { Toeplitz [c_{|i-j|}] PSD }, Caratheodory-Toeplitz; exact for every truncation).
Real characters: x = cos theta, E x^m = 2^-m sum_j C(m,j) c_{|m-2j|}.
Problem: maximise t subject to M_k(c) - k t I >= 0 for every form k (even and odd Legendre blocks).
Certificates: lower bound = min_k lambda_min(M_k(c))/k at the returned (Toeplitz-PSD) moments;
upper bound = sum_k <Z_k, M_k^0> + sum_p max_theta sum_l beta_{p,l} cos(l theta) for the dual matrices Z_k
(central path, rescaled so that sum_k k tr Z_k = 1); max over theta by a 20001-point grid + golden refinement.
"""
import sys
sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review2")
import numpy as np
from math import comb
from weil2 import Engine, prime_powers


def sym_weight(k):
    """Sym^k: coefficient of -L'/L at p^m = log p * sum_j cos((k-2j) m theta)."""
    def w(p, m):
        d = {}
        for j in range(k + 1):
            l = abs(k - 2 * j) * m
            d[l] = d.get(l, 0) + 1
        return d
    return w


def char_weight(p, m):
    """real character: log p * E x^m, x = cos theta."""
    d = {}
    for j in range(m + 1):
        l = abs(m - 2 * j)
        d[l] = d.get(l, 0) + comb(m, j) / 2 ** m
    return d


class Form:
    def __init__(self, kappas, k, weight, center=False, logq_extra=0.0, fixed=None):
        self.kappas, self.k, self.weight, self.center = kappas, k, weight, center
        self.logq_extra = logq_extra   # extra log-conductor (e.g. levels), added as + logq I
        self.fixed = fixed or {}       # p -> dict l -> value (known local data, e.g. Steinberg)


def build(eng, forms, decouple=False):
    pps = prime_powers(eng.X)
    # variables
    need = {}
    for F in forms:
        for p, m in pps:
            if p in F.fixed:
                continue
            for l in F.weight(p, m):
                if l > 0:
                    key = (p, m, l) if decouple else (p, l)
                    need.setdefault(key, True)
    if decouple:
        Lpm = {}
        for (p, m, l) in need:
            Lpm[(p, m)] = max(Lpm.get((p, m), 0), l)
        var = [(p, m, l) for (p, m) in sorted(Lpm) for l in range(1, Lpm[(p, m)] + 1)]
    else:
        Lp = {}
        for (p, l) in need:
            Lp[p] = max(Lp.get(p, 0), l)
        var = [(p, l) for p in sorted(Lp) for l in range(1, Lp[p] + 1)]
    vid = {v: i for i, v in enumerate(var)}
    nv = len(var) + 1  # last = t
    blocks = []
    meta = []
    for F in forms:
        M0 = eng.arch_sym(F.kappas) + F.logq_extra * np.eye(eng.N)
        if F.center:
            cv = eng.center_vec()
            M0 = M0 - np.outer(cv, cv)
        B = np.zeros((nv, eng.N, eng.N))
        for p, m in pps:
            S = eng.prime_sym(p ** m, 1.0) * np.log(p)   # = -2 log p p^{-m/2} S(m log p)
            for l, mult in F.weight(p, m).items():
                if p in F.fixed:
                    M0 = M0 + mult * F.fixed[p](m, l) * S
                elif l == 0:
                    M0 = M0 + mult * S
                else:
                    key = (p, m, l) if decouple else (p, l)
                    B[vid[key]] += mult * S
        B[-1] = -F.k * np.eye(eng.N)
        for ix in eng.blocks():
            blocks.append((M0[np.ix_(ix, ix)], B[:, ix][:, :, ix].copy()))
            meta.append(("form", F, ix, M0))
    # moment constraints
    if decouple:
        # each (p, m) its own measure: Toeplitz in its own l's (l multiples of m up to max)
        groups = {}
        for (p, m, l) in var:
            groups.setdefault((p, m), []).append(l)
        for (p, m), ls in groups.items():
            Lmax = max(ls)
            n = Lmax + 1
            C = np.eye(n)
            B = np.zeros((nv, n, n))
            # free cos-moments at l not used: add as extra variables? keep exact: unused l are free -> we
            # need them as variables; to stay simple, require ls to be all multiples present
            for l in range(1, n):
                key = (p, m, l)
                if key in vid:
                    B[vid[key]] = np.eye(n, k=l) + np.eye(n, k=-l)
                else:
                    raise ValueError("decouple needs full l range")
            blocks.append((C, B))
            meta.append(("toep", (p, m)))
    else:
        Lp = {}
        for (p, l) in var:
            Lp[p] = max(Lp.get(p, 0), l)
        for p in sorted(Lp):
            n = Lp[p] + 1
            C = np.eye(n)
            B = np.zeros((nv, n, n))
            for l in range(1, n):
                B[vid[(p, l)]] = np.eye(n, k=l) + np.eye(n, k=-l)
            blocks.append((C, B))
            meta.append(("toep", p))
    return var, blocks, meta


def Fx(blk, x):
    return blk[0] + np.tensordot(x, blk[1], axes=(0, 0))


def solve(var, blocks, meta, gap_end=1e-10, maxit=60, verbose=False):
    """log-barrier path following with damped Newton steps x += d/(1+lambda) (self-concordance: stays
    feasible and decreases, no function comparisons needed). Stops when mu*ntot < gap_end, or when
    centring fails (then the last centred point is returned)."""
    nv = len(var) + 1
    x = np.zeros(nv)
    tmin = min(np.linalg.eigvalsh(b[0])[0] / meta[i][1].k for i, b in enumerate(blocks) if meta[i][0] == "form")
    x[-1] = tmin - 1.0
    ntot = sum(b[0].shape[0] for b in blocks)
    mu = 1.0
    good = (x.copy(), mu)
    while True:
        ok = False
        for it in range(maxit):
            g = np.zeros(nv); g[-1] = -1.0 / mu
            H = np.zeros((nv, nv))
            for b in blocks:
                Fi = np.linalg.inv(Fx(b, x))
                Fi = 0.5 * (Fi + Fi.T)
                W = np.einsum("ab,ibc->iac", Fi, b[1], optimize=True)   # F^-1 A_i
                g -= np.einsum("iaa->i", W)
                H += np.einsum("iab,jba->ij", W, W, optimize=True)
            D = np.sqrt(np.diag(H))
            d = -np.linalg.solve(H / np.outer(D, D), g / D) / D
            lam = np.sqrt(max(-g @ d, 0.0))
            if lam < 1e-7:
                ok = True
                break
            step = 1.0 if lam < 0.2 else 1.0 / (1.0 + lam)
            xn = x + step * d
            try:
                for b in blocks:
                    np.linalg.cholesky(Fx(b, xn))
            except np.linalg.LinAlgError:
                break
            x = xn
        if not ok:
            break               # keep x (best primal so far, still feasible) for the lower bound
        good = (x.copy(), mu)
        if verbose:
            print("   mu %.1e t %.13f lam %.1e it %d" % (mu, x[-1], lam, it))
        if mu * ntot < gap_end:
            break
        mu /= 4.0
    xg, mug = good              # dual estimate from the last well-centred point (any PSD Z gives a valid bound)
    Z = [mug * np.linalg.inv(Fx(b, xg)) for b in blocks]
    return x, Z, mug * ntot


def theta_max(beta, n=20001):
    """max over theta in [0, pi] of sum_l beta[l] cos(l theta) (beta dict l -> value)."""
    th = np.linspace(0, np.pi, n)
    ls = np.array(sorted(beta)); bs = np.array([beta[l] for l in ls])
    f = lambda t: float(np.sum(bs * np.cos(ls * t)))
    vals = np.cos(np.outer(th, ls)) @ bs
    i = int(np.argmax(vals))
    lo, hi = th[max(i - 1, 0)], th[min(i + 1, n - 1)]
    gr = (np.sqrt(5) - 1) / 2
    for _ in range(80):
        c1 = hi - gr * (hi - lo); c2 = lo + gr * (hi - lo)
        if f(c1) > f(c2):
            hi = c2
        else:
            lo = c1
    tb = 0.5 * (lo + hi)
    return max(f(tb), vals[i]), tb


def certify(eng, var, blocks, meta, x, Z):
    nv = len(var) + 1
    # lower bound from the primal moments
    lam = {}
    for i, b in enumerate(blocks):
        if meta[i][0] == "form":
            k = meta[i][1].k
            M = Fx((b[0], b[1][:-1]), x[:-1])
            v = np.linalg.eigvalsh(M)[0] / k
            key = id(meta[i][1])
            lam[key] = min(lam.get(key, np.inf), v)
    LB = min(lam.values())
    toepmin = min(np.linalg.eigvalsh(Fx(b, x))[0] for i, b in enumerate(blocks) if meta[i][0] == "toep")
    # upper bound from the dual matrices of the form blocks
    Zf = [(Z[i], b, meta[i]) for i, b in enumerate(blocks) if meta[i][0] == "form"]
    scale = sum(m[1].k * np.trace(z) for z, b, m in Zf)
    UB = 0.0
    beta = np.zeros(nv - 1)
    for z, b, m in Zf:
        zz = z / scale
        UB += np.sum(zz * b[0])
        beta += np.einsum("ab,iab->i", zz, b[1][:-1])
    groups = {}
    for j, v in enumerate(var):
        key = v[:-1] if len(v) == 3 else v[0]
        groups.setdefault(key, {})[v[-1]] = beta[j]
    argmax = {}
    for key, bd in groups.items():
        mx, tb = theta_max(bd)
        UB += mx
        argmax[key] = tb
    return LB, UB, toepmin, argmax


def atoms(c):
    """atoms on [0, pi] of the symmetric circle measure with cos-moments c[0..L] (Prony on the Toeplitz kernel)."""
    L = len(c) - 1
    T = np.array([[c[abs(i - j)] for j in range(L + 1)] for i in range(L + 1)])
    ev = np.linalg.eigvalsh(T)
    r = int(np.sum(ev > 1e-6 * max(ev[-1], 1)))
    if r > L:
        return None, ev
    Tr = T[:r + 1, :r + 1]
    w, V = np.linalg.eigh(Tr)
    v = V[:, 0]
    z = np.roots(v[::-1])
    th = np.sort(np.abs(np.angle(z)))
    return th, ev


def moments_of(var, x):
    mom = {}
    for j, v in enumerate(var):
        mom.setdefault(v[0] if len(v) == 2 else v[:-1], {})[v[-1]] = x[j]
    out = {}
    for p, d in mom.items():
        Lm = max(d)
        out[p] = np.array([1.0] + [d.get(l, np.nan) for l in range(1, Lm + 1)])
    return out


def run(X, N, forms, decouple=False, nu=800, eng=None):
    eng = eng or Engine(X, N, nu=nu)
    var, blocks, meta = build(eng, forms, decouple=decouple)
    x, Z, gap = solve(var, blocks, meta)
    LB, UB, tmin, argmax = certify(eng, var, blocks, meta, x, Z)
    return dict(t=x[-1], LB=LB, UB=UB, gap=gap, toepmin=tmin, mom=moments_of(var, x), argmax=argmax,
                var=var, x=x, eng=eng)


def true_value(eng, forms, local):
    """min_k lambda_min(M_k)/k with known local data local[p] = array of cos-moments c_l (c_0 = 1)."""
    out = []
    for F in forms:
        M = eng.arch_sym(F.kappas) + F.logq_extra * np.eye(eng.N)
        if F.center:
            cv = eng.center_vec(); M = M - np.outer(cv, cv)
        for p, m in prime_powers(eng.X):
            S = eng.prime_sym(p ** m, 1.0) * np.log(p)
            for l, mult in F.weight(p, m).items():
                if p in F.fixed:
                    M = M + mult * F.fixed[p](m, l) * S
                else:
                    M = M + mult * local[p][l] * S
        out.append(min(np.linalg.eigvalsh(M[np.ix_(ix, ix)])[0] for ix in eng.blocks()) / F.k)
    return out
