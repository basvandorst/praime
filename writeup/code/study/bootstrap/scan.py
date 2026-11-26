# Float64 scans: which multisets of generalized primes in (1, X] keep the form positive?
import sys, itertools, numpy as np
from mpmath import mp, mpf, log
from weil import build_A
mp.dps = 30
class Form:
    def __init__(self, X, N):
        self.X = X; self.N = N; self.a = float(log(mpf(X))/2)
        A = build_A(log(mpf(X))/2, N)
        n = 2*N + 1
        self.A = np.array([[float(A[i, j]) for j in range(n)] for i in range(n)])
        m = np.arange(-N, N+1)
        self.m = m
        d = m[:, None] - m[None, :]
        self.d = d
        with np.errstate(divide='ignore'):
            self.pref = np.where(d != 0, (-1.0)**np.abs(d) * 2*self.a/(np.pi*np.where(d == 0, 1, d)), 0.0)
    def S(self, w):
        a = self.a; sn = np.sin(np.pi*self.m*w/a); cs = np.cos(np.pi*self.m*w/a)
        S = self.pref*(sn[None, :] - sn[:, None])
        S[np.diag_indices_from(S)] = 2*(2*a - w)*cs
        return S
    def comb(self, pos):
        C = np.zeros_like(self.A)
        for x in pos:
            lx = np.log(x); k = 1
            while k*lx <= 2*self.a:
                C += lx*x**(-k/2)*self.S(k*lx); k += 1
        return C
    def mineig(self, pos):
        return np.linalg.eigvalsh(self.A - self.comb(pos))[0]/(2*self.a)
if __name__ == '__main__':
    X = float(sys.argv[1]); N = int(sys.argv[2])
    F = Form(X, N)
    print('X', X, 'N', N, 'empty', F.mineig([]), 'true primes', F.mineig([p for p in (2, 3, 5, 7, 11, 13) if p <= X]))
    g1 = np.arange(1.04, X + 1e-9, 0.01)
    e1 = np.array([F.mineig([x]) for x in g1])
    i = e1.argmax(); print('1 point : best', g1[i], e1[i], ' feasible grid points:', g1[e1 > -1e-9])
    best = []
    for i, x in enumerate(g1):
        for y in g1[i:]:
            best.append((F.mineig([x, y]), x, y))
    best.sort(reverse=True)
    print('2 points: top', [(round(e, 6), round(x, 2), round(y, 2)) for e, x, y in best[:8]])
    print('   feasible (> -1e-9):', [(round(x, 2), round(y, 2)) for e, x, y in best if e > -1e-9][:30])
    g3 = np.arange(1.06, X + 1e-9, 0.03)
    best = []
    for i, x in enumerate(g3):
        for j, y in enumerate(g3[i:]):
            for z in g3[i+j:]:
                best.append((F.mineig([x, y, z]), x, y, z))
    best.sort(reverse=True)
    print('3 points: top', [(round(e, 5), round(x, 2), round(y, 2), round(z, 2)) for e, x, y, z in best[:8]])
