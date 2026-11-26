# Positivity bootstrap for Maass forms of level 1.  Input: spectral parameter r, parity eps, the primes.
# Unknown: for each prime p <= X a probability measure mu_p on the Satake angle theta in [0, pi] (grid of J angles).
# Constraints: the truncated Weil forms of L(f), L(Sym^2 f), L(Sym^3 f), L(Sym^4 f) (conductor q^k for Sym^k) are >= 0.
# Output: t*(r) = max t such that all forms are >= 0 at conductor q = e^(-t); q*(r) = e^(-t*) is the smallest
# admissible conductor.  t* < 0 at support X proves (Rayleigh-Ritz, up to float64) that no system with these
# Gamma factors, real Satake angles and conductor 1 satisfies the four positivities.
# Kelley cutting planes + LP (HiGHS).  Float64.
import sys, time, numpy as np
from scipy.optimize import linprog
from sympy import primerange
from engine import Basis

def types(r, eps, kmax=4):
    T = [('f',  [eps + 1j*r, eps - 1j*r], 1, lambda ph: 2*np.cos(ph))]
    if kmax >= 2: T.append(('S2', [0, 2j*r, -2j*r], 2, lambda ph: 1 + 2*np.cos(2*ph)))
    if kmax >= 3: T.append(('S3', [eps + 3j*r, eps - 3j*r, eps + 1j*r, eps - 1j*r], 3, lambda ph: 2*np.cos(ph) + 2*np.cos(3*ph)))
    if kmax >= 4: T.append(('S4', [0, 2j*r, -2j*r, 4j*r, -4j*r], 4, lambda ph: 1 + 2*np.cos(2*ph) + 2*np.cos(4*ph)))
    return T

class Problem:
    def __init__(self, r, eps, X, J=61, kmax=4, N=None, margin=12.0):
        self.r = r; self.eps = eps; self.X = X
        a = np.log(X)/2
        self.N = N or max(16, int(np.ceil(a*(kmax*r + margin)/np.pi)))
        self.B = B = Basis(X, self.N, maxfreq=kmax*r)
        self.T = types(r, eps, kmax)
        self.A = [B.arch(k) for (_, k, _, _) in self.T]
        self.th = np.linspace(0, np.pi, J)
        self.pp = [(p, m) for p in primerange(2, int(X) + 1) for m in range(1, 60) if m*np.log(p) <= 2*a - 1e-9]
        self.primes = sorted({p for p, m in self.pp})
        self.Spm = np.array([B.S(m*np.log(p))[0] for p, m in self.pp])          # (n_pp, n, n)
        # c[T][i_pp, j] = c_T(m theta_j) log p p^(-m/2)
        self.c = [np.array([cf(m*self.th)*np.log(p)*p**(-m/2) for p, m in self.pp]) for (_, _, _, cf) in self.T]
        self.pidx = np.array([self.primes.index(p) for p, m in self.pp])
    def forms(self, w):
        """w: (n_primes, J) measures.  Returns list of matrices H_T = A_T - P_T(w)."""
        out = []
        for A, c in zip(self.A, self.c):
            coef = np.einsum('ij,ij->i', c, w[self.pidx])
            out.append(A - np.einsum('i,ijk->jk', coef, self.Spm))
        return out
    def mineigs(self, w):
        two_a = 2*self.B.a
        return [np.linalg.eigvalsh(H)[0]/(two_a*k) for H, (_, _, k, _) in zip(self.forms(w), self.T)]
    def cut(self, iT, v):
        s = np.einsum('j,ijk,k->i', v, self.Spm, v)                               # (n_pp,)
        row = np.zeros((len(self.primes), len(self.th)))
        np.add.at(row, self.pidx, self.c[iT]*s[:, None])
        return float(v @ self.A[iT] @ v), row.ravel()
    def solve(self, tol=1e-7, rounds=400, ncut=4, verbose=False):
        nP, J = len(self.primes), len(self.th); two_a = 2*self.B.a
        rhs, rows, tk = [], [], []
        for iT, A in enumerate(self.A):
            w, V = np.linalg.eigh(A)
            for i in range(min(V.shape[1], 12)):
                a_v, row = self.cut(iT, V[:, i]); rhs.append(a_v); rows.append(row); tk.append(two_a*self.T[iT][2])
        Aeq = np.zeros((nP, nP*J + 1))
        for i in range(nP): Aeq[i, i*J:(i + 1)*J] = 1
        cost = np.zeros(nP*J + 1); cost[-1] = -1
        for rd in range(rounds):
            Aub = np.hstack([np.array(rows), np.array(tk)[:, None]])
            res = linprog(cost, A_ub=Aub, b_ub=np.array(rhs), A_eq=Aeq, b_eq=np.ones(nP),
                          bounds=[(0, 1)]*(nP*J) + [(-50, 50)], method='highs')
            if res.status != 0: raise RuntimeError(res.message)
            w = res.x[:-1].reshape(nP, J); t = res.x[-1]
            lam = np.inf; added = 0
            for iT, H in enumerate(self.forms(w)):
                ev, V = np.linalg.eigh(H); k = self.T[iT][2]
                lam = min(lam, ev[0]/(two_a*k))
                for i in range(ncut):
                    if ev[i]/(two_a*k) < t - tol:
                        a_v, row = self.cut(iT, V[:, i]); rhs.append(a_v); rows.append(row); tk.append(two_a*k); added += 1
            if verbose: print('   round', rd, 't_LP %.8f' % t, 'lam %.8f' % lam, 'cuts', len(rows), flush=True)
            if added == 0: break
        self.w = w
        return lam, t, rd
    def summary(self, w=None):
        w = self.w if w is None else w
        return ' '.join('%d:%+.3f' % (p, float(w[i] @ (2*np.cos(self.th)))) for i, p in enumerate(self.primes))

if __name__ == '__main__':
    eps = int(sys.argv[1]); X = float(sys.argv[2]); kmax = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    for r in [float(v) for v in sys.argv[3].split(',')]:
        t0 = time.time(); P = Problem(r, eps, X, kmax=kmax)
        lam, t, rd = P.solve()
        print('eps %d X %g kmax %d r %.4f N %d: t_low %.6f t_LP %.6f q* %.5f (%d rounds, %.0fs) per form %s  mean a_p %s' % (
            eps, X, kmax, r, P.N, lam, t, np.exp(-t), rd, time.time() - t0, ' '.join('%.4f' % v for v in P.mineigs(P.w)), P.summary()), flush=True)
