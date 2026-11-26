# Float64 bootstrap for a general degree-1 archimedean datum: conductor q (any real > 0), parity kappa,
# pole or not.  Unknowns: points x > 1 with a sign (real case).  No primes and no character are put in.
import sys, numpy as np
from scipy.optimize import minimize
from mpmath import mp, mpf, log
from weil import build_A
mp.dps = 30
class Form:
    def __init__(self, X, N, pole, kappa, q):
        self.X = X; self.N = N; self.a = float(log(mpf(X))/2)
        A = build_A(log(mpf(X))/2, N, pole=pole, kappa=kappa, logq=log(mpf(q)))
        n = 2*N + 1
        self.A = np.array([[float(A[i, j]) for j in range(n)] for i in range(n)])
        self.m = np.arange(-N, N+1); d = self.m[:, None] - self.m[None, :]
        self.pref = np.where(d != 0, (-1.0)**np.abs(d) * 2*self.a/(np.pi*np.where(d == 0, 1, d)), 0.0)
    def S(self, w):
        a = self.a; sn = np.sin(np.pi*self.m*w/a); cs = np.cos(np.pi*self.m*w/a)
        S = self.pref*(sn[None, :] - sn[:, None]); S[np.diag_indices_from(S)] = 2*(2*a - w)*cs
        return S
    def mineig(self, pos, sg):
        C = np.zeros_like(self.A)
        for x, s in zip(pos, sg):
            lx = np.log(x); k = 1
            while k*lx <= 2*self.a:
                C += s**k*lx*x**(-k/2)*self.S(k*lx); k += 1
        return np.linalg.eigvalsh(self.A - C)[0]/(2*self.a)
def polish(F, pos, sg):
    if not pos: return pos, F.mineig([], [])
    f = lambda v: 1e3 if (min(v) <= 1.001 or max(v) > F.X) else -F.mineig(list(v), sg)
    r = minimize(f, np.array(pos), method='Nelder-Mead', options=dict(xatol=1e-13, fatol=1e-16, maxiter=3000))
    r = minimize(f, r.x, method='Nelder-Mead', options=dict(xatol=1e-14, fatol=1e-17, maxiter=3000))
    return list(r.x), -r.fun
def run(pole, kappa, q, Xs, TOL=1e-10, verbose=True):
    pos, sg = [], []
    for X in Xs:
        N = max(12, int(8 + 4*np.log(X)/np.log(2)))
        F = Form(X, N, pole, kappa, q)
        pos, e = polish(F, pos, sg)
        for trial in range(3):
            if e > -TOL: break
            grid = np.arange(1.03, X, 0.01); best = None
            for s in (1, -1):
                vals = [F.mineig(pos + [y], sg + [s]) for y in grid]
                i = int(np.argmax(vals))
                if best is None or vals[i] > best[0]: best = (vals[i], grid[i], s)
            cand, e2 = polish(F, pos + [best[1]], sg + [best[2]])
            if e2 <= e + 1e-13: break
            pos, sg, e = cand, sg + [best[2]], e2
        for i in sorted(range(len(pos)), key=lambda i: -pos[i]):
            rp = pos[:i] + pos[i+1:]; rs = sg[:i] + sg[i+1:]
            if F.mineig(rp, rs) > -TOL: pos, sg = rp, rs; e = F.mineig(pos, sg)
        o = np.argsort(pos); pos = [pos[i] for i in o]; sg = [sg[i] for i in o]
        if verbose:
            print('  X=%-5g lam_min %9.2e  points: %s' % (X, e, '  '.join('%s%.10f' % ('+' if s > 0 else '-', x) for x, s in zip(pos, sg))), flush=True)
        if e < -1e-6:
            if verbose: print('  -> no positive configuration found at this support'); 
            return pos, sg, X, False
    return pos, sg, Xs[-1], True
if __name__ == '__main__':
    kappa = int(sys.argv[1]); q = float(sys.argv[2]); pole = len(sys.argv) > 4 and sys.argv[4] == 'pole'
    Xs = [float(v) for v in sys.argv[3].split(',')]
    print('datum: conductor q=%g parity %d pole %s' % (q, kappa, pole))
    run(pole, kappa, q, Xs)
