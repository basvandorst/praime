# Reviewer's own law: c_E by the direct Euler product of F(s)zeta(s) at s=1, line factors rho(r) from the
# ratio of local Euler factors sum_k lam(ak) mu(bk) k^-s / sum_k lam(k) mu(k) k^-s at s=1, and the pair
# integral J12 = int H12(iu) K(u) du with H12 = prod_p sum_v rho1 rho2 p^-|v| p^{iuv} (exact local factors).
import os, sys, json, time
import numpy as np

SC = os.path.dirname(os.path.abspath(__file__))+'/'
D = np.load(SC + 'rv_curves.npz')
ids, N, w, primes, AP = D['ids'], D['N'], D['w'], D['primes'], D['ap'].astype(float)
nc = len(ids)
pr = primes.astype(float)
lam = AP / np.sqrt(pr)[None, :]                     # lambda_p (normalised), good and bad
bad = (N[:, None] % primes[None, :] == 0)           # p | N


def lam_pows(lp, isbad, V):
    """lambda(p^m), m = 0..V for arrays lp, isbad (good: Chebyshev recursion; bad: lp^m)."""
    L = np.zeros((V + 1,) + lp.shape)
    L[0] = 1.0
    L[1] = lp
    for m in range(2, V + 1):
        L[m] = np.where(isbad, L[m - 1] * lp, lp * L[m - 1] - L[m - 2])
    return L


def local_rho(lp, isbad, p, V):
    """rho(p^v) for v = -2..V (index v+2). Good: mu(p) = -lp, mu(p^2) = 1. Bad: mu(p) = -lp, mu(p^2)=0."""
    L = lam_pows(lp, isbad, V + 2)
    mu1 = -lp
    mu2 = np.where(isbad, 0.0, 1.0)
    F = 1.0 + lp * mu1 / p + L[2] * mu2 / p**2      # F_p(1/p) = sum_j lam(p^j) mu(p^j) p^-j
    R = np.zeros((V + 3,) + lp.shape)
    # v >= 0 (p^v || a): sum_j lam(p^{v+j}) mu(p^j) p^-j
    for v in range(0, V + 1):
        R[v + 2] = (L[v] + L[v + 1] * mu1 / p + L[v + 2] * mu2 / p**2) / F
    # v = -1 (p || b): sum_j lam(p^j) mu(p^{1+j}) p^-j = mu(p) + lam(p) mu(p^2)/p
    R[1] = (mu1 + lp * mu2 / p) / F
    # v = -2 (p^2 || b): mu(p^2)
    R[0] = mu2 / F
    return R, F


def cE(k, Pmax=None):
    m = slice(None) if Pmax is None else (primes <= Pmax)
    lp = lam[k, m]; ib = bad[k, m]; p = pr[m]
    L2 = lp * lp - 1.0
    F = np.where(ib, 1.0 - 1.0 / p**2, (1 - 1 / p) * (1 - L2 / p))
    return float(np.exp(np.sum(np.log(F / (1 - 1 / p)))))


def rho_r(k, a, b):
    val = 1.0
    for (n, sgn) in ((a, 1), (b, -1)):
        x = n
        for p in [2, 3, 5, 7, 11, 13]:
            v = 0
            while x % p == 0:
                x //= p; v += 1
            if v:
                i = int(np.searchsorted(primes, p))
                R, F = local_rho(np.array([lam[k, i]]), np.array([bad[k, i]]), float(p), 6)
                val *= R[sgn * v + 2][0]
        assert x == 1
    return val


LINES = [(1, 1), (2, 1), (1, 2), (3, 1), (1, 3), (4, 1), (1, 4), (5, 1), (1, 5), (7, 1), (1, 7), (3, 2), (2, 3)]

if __name__ == '__main__':
    t0 = time.time()
    c = np.array([cE(k) for k in range(nc)])
    c1e5 = np.array([cE(k, 100000) for k in range(nc)])
    PF = np.load(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/ecross/pred_fresh.npz')
    idx = {s: i for i, s in enumerate(ids)}
    perm = np.array([idx[s] for s in PF['ids']])
    cf = PF['c']
    print('c_E own (P=5e5) / frozen c: min %.5f max %.5f rms log %.5f' % ((c[perm] / cf).min(), (c[perm] / cf).max(),
                                                                          np.sqrt(np.mean(np.log(c[perm] / cf) ** 2))))
    print('c_E own (P=1e5) / own (P=5e5): rms log %.5f' % np.sqrt(np.mean(np.log(c1e5 / c) ** 2)))
    print('c range %.3f .. %.3f median %.3f' % (c.min(), c.max(), np.median(c)))
    print('frozen lines', PF['lines'].tolist())
    amp = np.zeros((nc, len(LINES)))
    for k in range(nc):
        for j, (a, b) in enumerate(LINES):
            amp[k, j] = -w[k] * c[k] * rho_r(k, a, b) / np.sqrt(a * b)      # per unit T/2pi
    # compare with frozen lineamp (find its normalisation)
    la = PF['lineamp']
    ratio = la / amp[perm][:, [LINES.index(tuple(x)) for x in PF['lines'].tolist()]]
    print('frozen lineamp / own amp per unit T/2pi: median %.6f, min %.6f, max %.6f' % (np.median(ratio), ratio.min(), ratio.max()))
    print('380/2pi = %.6f' % (380 / 2 / np.pi))
    np.savez(SC + 'rv_law.npz', ids=ids, c=c, amp=amp, lines=np.array(LINES))
    print('done', time.time() - t0)
