# Positivity bootstrap for Maass forms on SL(3,Z) (level 1).  Own float64 code.
# Data: spectral parameters (l1, l2), l3 = -l1 - l2; Gamma factor prod_j Gamma_R(s + i l_j) (FKL / FKLR convention).
# Unknown: for each prime p <= X a probability measure on the conjugacy classes of SU(3) (angles (phi1, phi2, phi3), sum 0).
# Forms (truncated Weil forms, complex test functions on [-a, a], basis exp(i pi k u / a), Hermitian matrices):
#   L(F):     Gamma_R(s + i l_j), weight at p^m: tau_m = sum_j e^(i m phi_j) (complex; the dual form has the conjugate)
#   L(Ad F):  Gamma_R(s)^2 prod_d Gamma_R(s + i d) Gamma_R(s - i d), d in {l1-l2, 2 l1+l2, l1+2 l2}; weight |tau_m|^2 - 1; no pole
# Q(f) = sum_j [ -(gamma + log pi) g(0) + int_0^inf (2 e^(-2u) g(0) - e^(-u/2) (e^(-i l_j u) g(u) + e^(i l_j u) g(-u))) / (1 - e^(-2u)) du ]
#        - sum_n n^(-1/2) (c(n) g(log n) + conj(c(n)) g(-log n)),     g(w) = int f(v + w) conj f(v) dv.
import sys, numpy as np
from scipy.optimize import linprog
from sympy import primerange
from engine import gl, EULER

def G_stack(a, N, w):
    """G(w)[k, j] = g_jk(w) = int phi_j(v + w) conj(phi_k(v)) dv for w >= 0; shape (len(w), n, n), complex."""
    w = np.atleast_1d(np.asarray(w, float)); ks = np.arange(-N, N + 1)
    E = np.exp(1j*np.pi*np.outer(w, ks)/a)                                   # E[w, k] = e^(i pi k w / a)
    d = ks[None, :] - ks[:, None]                                            # d[k, j] = j - k
    with np.errstate(divide='ignore', invalid='ignore'):
        fac = np.where(d != 0, np.where(d % 2 == 0, 1.0, -1.0)*a/(1j*np.pi*np.where(d == 0, 1, d)), 0.0)
    G = fac[None, :, :]*(E[:, :, None] - E[:, None, :])                      # (e^(i pi k w/a) - e^(i pi j w/a)), row k, column j
    idx = np.arange(len(ks)); G[:, idx, idx] = (2*a - w)[:, None]*E
    return G

class HBasis:
    def __init__(self, X, N, maxfreq):
        self.X = X; self.a = a = np.log(X)/2; self.N = N; self.n = 2*N + 1
        self.u, self.wt = gl(0.0, 2*a, int(40 + 2*a*(2*np.pi*N/a + maxfreq)/np.pi*1.5))
    def arch(self, lams):
        """Hermitian matrix of the archimedean part for prod_j Gamma_R(s + i lam_j), conductor 1, no pole."""
        a = self.a; n = self.n; ng = len(lams); I = np.eye(n); u = self.u
        A = np.zeros((n, n), complex)
        den = 1 - np.exp(-2*u)
        for i in range(0, len(u), 300):
            Gs = G_stack(a, self.N, u[i:i+300])
            ker = sum(np.exp(-0.5*u[i:i+300])*np.exp(-1j*l*u[i:i+300]) for l in lams)/den[i:i+300]*self.wt[i:i+300]
            A -= np.einsum('i,ijk->jk', ker, Gs)
        A = A + A.conj().T
        A += I*2*a*2*ng*np.sum(self.wt*np.exp(-2*u)/den)
        A += I*2*a*ng*(-(EULER + np.log(np.pi)) - np.log(1 - np.exp(-4*a)))
        return A

def torus_grid(ng):
    """Conjugacy classes of SU(3) on a grid: unordered triples of angles with sum 0 mod 2 pi."""
    seen = {}; h = 2*np.pi/ng
    for i in range(ng):
        for j in range(ng):
            k = (-i - j) % ng; key = tuple(sorted((i, j, k)))
            seen[key] = (i*h, j*h, k*h)
    return np.array(list(seen.values()))

class GL3:
    def __init__(self, l1, l2, X, ng=30, N=None, margin=12.0):
        self.lam = (l1, l2, -l1 - l2); a = np.log(X)/2
        ds = (l1 - l2, 2*l1 + l2, l1 + 2*l2); maxf = max(ds)
        self.N = N or max(16, int(np.ceil(a*(maxf + margin)/np.pi)))
        self.B = B = HBasis(X, self.N, maxf)
        self.names = ['F', 'Ad']; self.k = [1, 2]
        self.A = [B.arch(self.lam), B.arch([0, 0] + [s*d for d in ds for s in (1, -1)])]
        self.pp = [(p, m) for p in primerange(2, int(X) + 1) for m in range(1, 60) if m*np.log(p) <= 2*a - 1e-9]
        self.primes = sorted({p for p, m in self.pp}); self.pidx = np.array([self.primes.index(p) for p, m in self.pp])
        self.Gpm = np.array([G_stack(a, self.N, m*np.log(p))[0] for p, m in self.pp])
        self.grid = torus_grid(ng)
        self.fac = np.array([np.log(p)*p**(-m/2) for p, m in self.pp])
    def tau(self, ang, m):
        return np.exp(1j*m*np.asarray(ang)).sum(axis=-1)
    def weights(self, ang):
        """c_F[i_pp, g] (complex) and c_Ad[i_pp, g] (real) for angles ang[g] = (phi1, phi2, phi3)."""
        tm = np.array([self.tau(ang, m) for p, m in self.pp])
        return tm*self.fac[:, None], (np.abs(tm)**2 - 1)*self.fac[:, None]
    def form_mats(self, cF, cAd):
        """cF, cAd: per prime power coefficients (already averaged over the measures)."""
        PF = np.einsum('i,ijk->jk', cF, self.Gpm); PA = np.einsum('i,ijk->jk', cAd.astype(complex), self.Gpm)
        return [self.A[0] - (PF + PF.conj().T), self.A[1] - (PA + PA.conj().T)]
    def point(self, angles):
        """Smallest eigenvalues (per unit conductor exponent) at one conjugacy class per prime: angles[prime index] = (phi1, phi2, phi3)."""
        cF, cAd = self.weights(np.array([angles[i] for i in self.pidx])[:, None, :]) if False else (None, None)
        tm = np.array([self.tau(angles[i], m) for i, (p, m) in zip(self.pidx, self.pp)])
        Hs = self.form_mats(tm*self.fac, (np.abs(tm)**2 - 1)*self.fac)
        return [np.linalg.eigvalsh(H)[0]/(2*self.B.a*k) for H, k in zip(Hs, self.k)]
    def cut(self, iT, v, WF, WAd):
        s = np.einsum('j,ijk,k->i', v.conj(), self.Gpm, v)                   # s_pm = v^dagger G v (complex)
        per = 2*np.real((WF if iT == 0 else WAd)*s[:, None])                 # (n_pp, n_grid)
        row = np.zeros((len(self.primes), per.shape[1])); np.add.at(row, self.pidx, per)
        return float(np.real(v.conj() @ self.A[iT] @ v)), row.ravel()
    def solve(self, tol=1e-6, rounds=300, ncut=4):
        nP = len(self.primes); ng = len(self.grid); two_a = 2*self.B.a
        WF, WAd = self.weights(self.grid)
        rows, rhs, tk = [], [], []
        for iT in range(2):
            w, V = np.linalg.eigh(self.A[iT])
            for i in range(12):
                a_v, row = self.cut(iT, V[:, i], WF, WAd); rhs.append(a_v); rows.append(row); tk.append(two_a*self.k[iT])
        Aeq = np.zeros((nP, nP*ng + 1))
        for i in range(nP): Aeq[i, i*ng:(i + 1)*ng] = 1
        cost = np.zeros(nP*ng + 1); cost[-1] = -1
        for rd in range(rounds):
            res = linprog(cost, A_ub=np.hstack([np.array(rows), np.array(tk)[:, None]]), b_ub=np.array(rhs), A_eq=Aeq, b_eq=np.ones(nP),
                          bounds=[(0, 1)]*(nP*ng) + [(-50, 50)], method='highs')
            if res.status != 0: raise RuntimeError(res.message)
            w = res.x[:-1].reshape(nP, ng); t = res.x[-1]
            Hs = self.form_mats((WF[:, :] * w[self.pidx]).sum(axis=1), (WAd * w[self.pidx]).sum(axis=1)); added = 0; lam = np.inf; per_form = []
            for iT, H in enumerate(Hs):
                ev, V = np.linalg.eigh(H); per_form.append(ev[0]/(two_a*self.k[iT])); lam = min(lam, per_form[-1])
                for i in range(ncut):
                    if ev[i]/(two_a*self.k[iT]) < t - tol:
                        a_v, row = self.cut(iT, V[:, i], WF, WAd); rhs.append(a_v); rows.append(row); tk.append(two_a*self.k[iT]); added += 1
            if not added: break
        self.w = w; self.per_form = per_form
        self.a_p = [complex((w[i] @ self.tau(self.grid, 1))) for i in range(nP)]
        return lam, t, rd

if __name__ == '__main__':
    X = float(sys.argv[1])
    for pair in sys.argv[2:]:
        l1, l2 = map(float, pair.split(','))
        P = GL3(l1, l2, X); lam, t, rd = P.solve()
        print('(%.4f, %.4f) X %g N %d: t* = %.5f (%d rounds), per form %s, a_p %s' % (l1, l2, X, P.N, t, rd, ' '.join('%.4f' % v for v in P.per_form),
              ' '.join('%d:%+.3f%+.3fi' % (p, z.real, z.imag) for p, z in zip(P.primes, P.a_p))), flush=True)
