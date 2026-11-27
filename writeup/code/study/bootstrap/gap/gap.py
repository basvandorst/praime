# The highest lowest zero: a convex problem with the primes given and everything else unknown ("degree-free").
# An L-function of degree d: Gamma factor prod_i Gamma_R(s + kappa_i) (the kappa_i closed under conjugation), conductor
# N >= 1, -L'/L = sum c(n) n^-s with c(p^m) = log p * sum_j alpha_j(p)^m, |alpha_j(p)| <= 1 (Ramanujan), real c(n).
# Divide the explicit formula by d.  Unknown: a probability measure nu on the shifts kappa, l = (log N)/d >= 0, and for
# every prime p <= X a probability measure sigma_p on the angle theta in [0, pi]: c(p^m)/(d log p) = <cos(m theta)>.
# (An alpha inside the disk is an average of boundary points: z -> z^m is harmonic.  So the circle is enough.)
# With H the truncated Weil form (test functions on [-a, a], X = e^(2a), basis exp(i pi k u/a)):
#   zeros real (GRH)                      =>  H >= 0,
#   and no zero with |gamma| < T          =>  G_T = (pi/a)^2 K H K - T^2 H >= 0 on {psi(+-a) = 0},  K = diag(k),
# because sum_gamma (gamma^2 - T^2) |psihat(gamma)|^2 = Q(psi') - T^2 Q(psi).
# t*(T) = max t with H >= t 2a and G_T >= t 2a T^2 over all data.  t*(T) < 0 at any support X proves (Rayleigh-Ritz, up
# to float64): every such L-function has a zero with |gamma| < T.  Kelley cutting planes + LP (HiGHS).
import os, sys, time, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'maass'))
from engine import Basis
from scipy.optimize import linprog
from sympy import primerange

class Gap:
    def __init__(self, X, shifts, N=None, nth=97, pole=0, margin=40.0):
        self.X = X; a = self.a = np.log(X)/2
        tmax = max(abs(complex(k).imag) for k in shifts)
        self.N = N or max(16, int(np.ceil(a*(tmax + margin)/np.pi)))
        B = self.B = Basis(X, self.N, maxfreq=tmax); n = self.n = B.n
        self.shifts = list(shifts)
        self.pp = [(p, m) for p in primerange(2, int(X) + 1) for m in range(1, 60) if m*np.log(p) <= 2*a - 1e-9]
        self.primes = sorted({p for p, m in self.pp})
        self.th = np.linspace(0, np.pi, nth)
        ks = np.arange(-self.N, self.N + 1); self.Kd = ks*np.pi/a
        e = np.where(ks % 2 == 0, 1.0, -1.0); q, _ = np.linalg.qr(np.column_stack([e, np.eye(n)[:, :n-1]])); self.Z = q[:, 1:]
        # elementary matrices: shifts, identity (conductor), prime powers
        self.Ak = np.array([B.arch([k], pole=pole) for k in self.shifts])
        self.S = np.array([B.S(m*np.log(p))[0]*np.log(p)*p**(-m/2) for p, m in self.pp])
        self.I = np.eye(n)
        # cosm[i_pp, j] = cos(m theta_j);  pidx[i_pp] = index of the prime
        self.cosm = np.array([np.cos(m*self.th) for p, m in self.pp]); self.pidx = np.array([self.primes.index(p) for p, m in self.pp])
    def G(self, M, T):
        return self.Z.T @ (self.Kd[:, None]*M*self.Kd[None, :] - T*T*M) @ self.Z
    def assemble(self, nu, l, sig):
        coef = np.einsum('ij,ij->i', self.cosm, sig[self.pidx])
        return np.einsum('i,ijk->jk', nu, self.Ak) + l*2*self.a*self.I - np.einsum('i,ijk->jk', coef, self.S)
    def row(self, v, T=None):
        """the linear functional (nu..., l, sigma...) -> v'Hv (T None) or w'G_T w (v = coordinates in the basis Z)."""
        if T is None: f = lambda M: np.einsum('j,ijk,k->i', v, M, v)
        else:
            w = self.Z @ v; kw = self.Kd*w
            f = lambda M: np.einsum('j,ijk,k->i', kw, M, kw) - T*T*np.einsum('j,ijk,k->i', w, M, w)
        s = f(self.S); blk = np.zeros((len(self.primes), self.cosm.shape[1]))
        np.add.at(blk, self.pidx, -s[:, None]*self.cosm)
        return np.concatenate([f(self.Ak), 2*self.a*f(self.I[None]), blk.ravel()])
    def solve(self, T, tol=1e-9, rounds=400, lmax=5.0, ncut=4, nu_mask=None):
        nk, nP, J = len(self.shifts), len(self.primes), self.cosm.shape[1]; nv = nk + 1 + nP*J + 1; two_a = 2*self.a
        Aeq = np.zeros((1 + nP, nv)); Aeq[0, :nk] = 1
        for i in range(nP): Aeq[1 + i, nk + 1 + i*J: nk + 1 + (i + 1)*J] = 1
        cost = np.zeros(nv); cost[-1] = -1
        bnd = [(0, 1 if nu_mask is None or nu_mask[i] else 0) for i in range(nk)] + [(0, lmax)] + [(0, 1)]*(nP*J) + [(-50, 50)]
        rows = []
        nu = np.array([b[1] for b in bnd[:nk]], float); nu /= nu.sum(); l = 0.0; sig = np.ones((nP, J))/J; t = np.inf
        for rd in range(rounds):
            H = self.assemble(nu, l, sig); ev, V = np.linalg.eigh(H); er, U = np.linalg.eigh(self.G(H, T))
            lam = min(ev[0]/two_a, er[0]/(two_a*T*T))
            if lam >= t - tol: break
            for i in range(ncut):
                if rd == 0 or ev[i]/two_a < t - tol: rows.append(np.concatenate([-self.row(V[:, i]), [two_a]]))
                if rd == 0 or er[i]/(two_a*T*T) < t - tol: rows.append(np.concatenate([-self.row(U[:, i], T), [two_a*T*T]]))
            res = linprog(cost, A_ub=np.array(rows), b_ub=np.zeros(len(rows)), A_eq=Aeq, b_eq=np.ones(1 + nP), bounds=bnd, method='highs')
            if res.status != 0: raise RuntimeError(res.message)
            x = res.x; nu = x[:nk]; l = x[nk]; sig = x[nk + 1:-1].reshape(nP, J); t = x[-1]
        self.nu, self.l, self.sig, self.rows = nu, l, sig, rows
        return lam, t, rd
    def describe(self):
        fmt = lambda k: ('%g' % k.real) if k.imag == 0 else ('%g%+gi' % (k.real, k.imag))
        s = 'shifts ' + ' '.join('%s:%.2f' % (fmt(complex(k)), w) for k, w in zip(self.shifts, self.nu) if w > 0.01)
        s += ' | mean shift %.3f | log N/d = %.4f | <cos theta_p>:' % (float(np.real(np.dot(self.nu, self.shifts))), self.l)
        return s + ' '.join(' %d:%+.4f' % (p, float(w @ np.cos(self.th))) for p, w in zip(self.primes, self.sig))

def tmax_fixed(g, H, lo=1.0, hi=60.0, eps=-1e-9):
    """largest T with G_T >= 0 for a fixed form H (bisection)."""
    for _ in range(50):
        mid = (lo + hi)/2
        if np.linalg.eigvalsh(g.G(H, mid))[0] >= eps: lo = mid
        else: hi = mid
    return lo

def shift_grid(cls):
    if cls == 'real':  return [float(v) for v in np.round(np.arange(0, 30.01, 0.25), 4)]
    if cls == 'complex': return [complex(d, t) for d in (0, 0.5, 1, 2, 4, 8) for t in np.arange(0, 30.01, 0.5)]
    raise ValueError(cls)

if __name__ == '__main__':
    # usage: gap.py <real|complex> <X> <T>   (one T per process)
    cls = sys.argv[1]; X = float(sys.argv[2]); T = float(sys.argv[3])
    t0 = time.time(); g = Gap(X, shift_grid(cls)); lam, t, rd = g.solve(T)
    print('%s X %g N %d T %.4f: t* = %+.3e (lam %+.3e, %d rounds, %.0fs)  %s' % (cls, X, g.N, T, t, lam, rd, time.time() - t0, g.describe()), flush=True)
