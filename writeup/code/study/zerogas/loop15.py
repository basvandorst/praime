import sys, numpy as np
from fractions import Fraction as Fr
def H(n): return sum(Fr(1, i) for i in range(1, n + 1))
def c_exact(N, k):
    num = 1 - ((N + k)*(H(N + k) - H(N)) - (N - k)*(H(N) - H(N - k)))/N
    den = 1 - Fr(max(0, 2*k - N), N*k)
    return num/den
if __name__ == '__main__':
    from loop9 import cth
    rng = np.random.default_rng(15)
    for N, M in ((6, 40000), (10, 40000), (16, 40000)):
        ks = np.arange(1, N + 1); num = np.zeros((M, N)); den = np.zeros((M, N))
        for it in range(M):
            Z = (rng.standard_normal((N, N)) + 1j*rng.standard_normal((N, N)))/np.sqrt(2); Q, R = np.linalg.qr(Z)
            x = np.sort(np.angle(np.linalg.eigvals(Q*(np.diag(R)/np.abs(np.diag(R)))[None, :])) % (2*np.pi))*N/(2*np.pi)
            Sz = np.arange(N) - x
            E = np.exp(-2j*np.pi*np.outer(ks, x)/N); coef = E.sum(axis=1)/(2j*np.pi*ks); W = 2*np.real(coef[:, None]*np.conj(E))
            num[it] = W @ Sz; den[it] = (W*W).sum(axis=1)
        grp = np.array_split(np.arange(M), 20); worst = 0.0
        for i, k in enumerate(ks):
            c = num[:, i].mean()/den[:, i].mean(); e = np.array([num[g, i].mean()/den[g, i].mean() for g in grp]).std(ddof=1)/np.sqrt(20)
            ex = c_exact(N, int(k)); z = (c - float(ex))/e; worst = max(worst, abs(z))
            print('N %2d k %2d: simulation %+.5f +- %.5f | exact %+.5f = %s | limit curve %+.5f | z = %+.1f' % (N, k, c, e, float(ex), ex if N <= 6 else '...', float(cth(np.array([k/N]))[0]), z))
        print('N %d: largest |z| %.1f' % (N, worst), flush=True)
