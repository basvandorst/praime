"""Third reviewer's convex problem and solver (own code).

Unknown: for each prime p <= X one probability measure on theta in [0, pi], through its cosine moments
c_{p,l} = E cos(l theta), l = 1..L_p (Toeplitz [c_|i-j|] PSD: exact moment set, no angle grid).
Forms (coefficient of -L'/L at p^m, divided by log p):
  L(f):        2 E cos(m theta)                     Gamma_R(s + e +- i r), central zero for odd f (Q - |int phi|^2)
  Sym^2 f:     1 + 2 E cos(2 m theta)               Gamma_R(s) Gamma_R(s +- 2 i r)
  L(f x g):    b_m(g) * 2 E cos(m theta)            Gamma_R(s + d +- i(r + r_g)) Gamma_R(s + d +- i(r - r_g)), d = e_f + e_g mod 2
               b_m(g) = alpha_g^m + alpha_g^-m (given)
maximise t subject to M_F(c) - k_F t I >= 0 for every form F (k = 1, 2, krs), Toeplitz(c_p) >= 0.
Solver: log-barrier path following with damped Newton steps (own implementation).
Bounds: LB = min_F lambda_min(M_F(c))/k_F at the returned moments (a genuine measure exists: Toeplitz PD);
        UB = Lagrangian dual with Z_F = mu F^{-1}, normalised sum k tr Z = 1, max over each measure = max over theta.
"""
import sys
sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review3")
import numpy as np
from eng3 import Eng


def cheb_b(ap, M=8):
    """b_m = 2 cos(m theta) with 2 cos theta = ap: b_0 = 2, b_1 = ap, b_m = ap b_{m-1} - b_{m-2}."""
    b = [2.0, ap]
    for m in range(2, M + 1):
        b.append(ap * b[-1] - b[-2])
    return b


class Problem:
    def __init__(self, E, r, eps, center, known=(), sym2=True, krs=2.0, use_Lf=True):
        """known: list of dicts with keys r, eps, b (dict (p, m) -> b_m(g))."""
        self.E = E
        forms = []
        if use_Lf:
            forms.append(dict(name="L", kap=[eps + 1j * r, eps - 1j * r], k=1.0, center=center,
                              w=lambda p, m: (0.0, {m: 2.0})))
        if sym2:
            forms.append(dict(name="S2", kap=[0.0, 2j * r, -2j * r], k=2.0, center=False,
                              w=lambda p, m: (1.0, {2 * m: 2.0})))
        for g in known:
            d = (eps + g["eps"]) % 2
            kap = [d + 1j * (r + g["r"]), d - 1j * (r + g["r"]), d + 1j * (r - g["r"]), d - 1j * (r - g["r"])]
            b = g["b"]
            forms.append(dict(name="RS%.4f" % g["r"], kap=kap, k=krs, center=False,
                              w=(lambda bb: (lambda p, m: (0.0, {m: 2.0 * bb[(p, m)]})))(b)))
        self.forms = forms
        Lp = {}
        for F in forms:
            for p, m in E.pp:
                for l in F["w"](p, m)[1]:
                    Lp[p] = max(Lp.get(p, 0), l)
        self.Lp = Lp
        self.var = [(p, l) for p in sorted(Lp) for l in range(1, Lp[p] + 1)]
        vid = {v: i for i, v in enumerate(self.var)}
        nv = len(self.var) + 1
        self.nv = nv
        self.blocks = []   # (A0, A[nv], kind, form index, k)
        for fi, F in enumerate(forms):
            M0 = E.arch(F["kap"])
            if F["center"]:
                cv = E.center()
                M0 = M0 - np.outer(cv, cv)
            A = np.zeros((nv, E.N, E.N))
            for p, m in E.pp:
                P = E.prime(p, m)
                const, lin = F["w"](p, m)
                if const:
                    M0 = M0 + const * P
                for l, co in lin.items():
                    A[vid[(p, l)]] += co * P
            A[-1] = -F["k"] * np.eye(E.N)
            for ix in E.blocks:
                self.blocks.append((M0[np.ix_(ix, ix)].copy(), A[:, ix][:, :, ix].copy(), "form", fi, F["k"]))
        for p in sorted(Lp):
            n = Lp[p] + 1
            A = np.zeros((nv, n, n))
            for l in range(1, n):
                A[vid[(p, l)]] = np.eye(n, k=l) + np.eye(n, k=-l)
            self.blocks.append((np.eye(n), A, "toep", p, 0.0))

    def F(self, b, x):
        return b[0] + np.tensordot(x, b[1], axes=(0, 0))

    def solve(self, tol=1e-10, maxit=200):
        nv = self.nv
        x = np.zeros(nv)
        x[-1] = min(np.linalg.eigvalsh(b[0])[0] / b[4] for b in self.blocks if b[2] == "form") - 1.0
        nu = sum(b[0].shape[0] for b in self.blocks)
        mu = 1.0 / nu          # central point then lies about 1 below the optimum: short first centring
        last = None
        while True:
            centred = False
            for it in range(maxit):
                g = np.zeros(nv); g[-1] = -1.0 / mu
                H = np.zeros((nv, nv))
                for b in self.blocks:
                    Fi = np.linalg.inv(self.F(b, x))
                    W = np.einsum("ab,ibc->iac", Fi, b[1])
                    g -= np.einsum("iaa->i", W)
                    H += np.einsum("iab,jba->ij", W, W)
                Dg = np.sqrt(np.maximum(np.diag(H), 1e-300))
                dx = -np.linalg.solve(H / np.outer(Dg, Dg), g / Dg) / Dg
                lam = np.sqrt(max(-(g @ dx), 0.0))
                if lam < 1e-8:
                    centred = True
                    break
                st = 1.0 / (1.0 + lam) if lam > 0.25 else 1.0
                xn = x + st * dx
                ok = True
                for b in self.blocks:
                    try:
                        np.linalg.cholesky(self.F(b, xn))
                    except np.linalg.LinAlgError:
                        ok = False
                        break
                if not ok:
                    xn = x + dx / (1.0 + lam)
                    try:
                        for b in self.blocks:
                            np.linalg.cholesky(self.F(b, xn))
                    except np.linalg.LinAlgError:
                        break
                x = xn
            if not centred:
                break
            last = (x.copy(), mu)
            if mu * nu < tol:
                break
            mu /= 5.0
        self.x = x
        self.last = last
        return self.bounds()

    def bounds(self):
        x = self.x
        lb = np.inf
        for b in self.blocks:
            if b[2] == "form":
                M = b[0] + np.tensordot(x[:-1], b[1][:-1], axes=(0, 0))
                lb = min(lb, np.linalg.eigvalsh(M)[0] / b[4])
        toep = min(np.linalg.eigvalsh(self.F(b, x))[0] for b in self.blocks if b[2] == "toep")
        ub = np.nan
        if self.last is not None:
            xc, mu = self.last
            Zs = [(mu * np.linalg.inv(self.F(b, xc)), b) for b in self.blocks if b[2] == "form"]
            sc = sum(b[4] * np.trace(Z) for Z, b in Zs)
            ub = 0.0
            beta = np.zeros(self.nv - 1)
            for Z, b in Zs:
                Z = Z / sc
                ub += np.sum(Z * b[0])
                beta += np.einsum("ab,iab->i", Z, b[1][:-1])
            th = np.linspace(0, np.pi, 40001)
            for p in sorted(self.Lp):
                idx = [i for i, v in enumerate(self.var) if v[0] == p]
                ls = np.array([self.var[i][1] for i in idx])
                vals = np.cos(np.outer(th, ls)) @ beta[idx]
                ub += vals.max() + 1e-9 * np.abs(beta[idx]).sum()   # grid max + crude Lipschitz pad
        self.lb, self.ub, self.toep = lb, ub, toep
        return lb, ub

    def ap(self):
        """E 2 cos theta per prime, and the moments."""
        out = {}
        for i, (p, l) in enumerate(self.var):
            if l == 1:
                out[p] = 2 * self.x[i]
        return out

    def moments(self):
        mom = {}
        for i, (p, l) in enumerate(self.var):
            mom.setdefault(p, {})[l] = self.x[i]
        return mom


def true_min(E, r, eps, center, ap, known=(), sym2=True, krs=2.0):
    """smallest eigenvalue / k of every form with the given point-mass data ap[p] (2 cos theta_p)."""
    P = Problem(E, r, eps, center, known, sym2, krs)
    x = np.zeros(P.nv)
    for i, (p, l) in enumerate(P.var):
        x[i] = 0.5 * cheb_b(ap[p], l)[l]
    out = {}
    for b in P.blocks:
        if b[2] == "form":
            M = b[0] + np.tensordot(x[:-1], b[1][:-1], axes=(0, 0))
            nm = P.forms[b[3]]["name"]
            out[nm] = min(out.get(nm, np.inf), np.linalg.eigvalsh(M)[0] / b[4])
    return out
