import numpy as np, sys, json
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
from readings import *
from mydata import unfold, sigma2, theta_count, LZ
from comparators import law, gls_chi2, P
import mypred

S, full, diag, X, primes, ks = load('B')
Ds = list(np.load(W + 'S_B.npz')['D'])
z = np.load(LZ + 'r134_B.npz')
n = len(Ds)

def gamma_of(S_, F_, D_, X_):
    Sp = S_[:, pidx]; Fp = F_[:, pidx]; Dp = D_[:, pidx]
    bf, _ = ols(X_, Sp); bF, _ = ols(X_, Fp); bD, _ = ols(X_, Dp)
    off = (bF - bD)[1:]; res = (bf - bD)[1:]
    return (res * off).sum() / (off ** 2).sum(), res, off

def boot_gamma(S_, F_, D_, X_, nb=500):
    rng = np.random.default_rng(5); g = []
    for b in range(nb):
        i = rng.integers(0, len(S_), len(S_))
        try: g.append(gamma_of(S_[i], F_[i], D_[i], X_[i])[0])
        except Exception: pass
    return np.std(g, ddof=1)

# 1. bias of np.var: (sample mean - l)^2 per set
bias = []
for D in Ds:
    g = np.sort(z['D%d' % D]); g = g[(g >= 100) & (g <= 600)]; u = unfold(g, abs(D))
    u0 = unfold(100.0, abs(D)); u1 = unfold(600.0, abs(D)); row = []
    for l in prim:
        x = np.arange(u0, u1 - l, 0.01); c = np.searchsorted(u, x + l) - np.searchsorted(u, x)
        row.append((c.mean() - l) ** 2)
    bias.append(row)
bias = np.array(bias)
print('1. np.var bias = -(mean - l)^2: mean over sets at l=0.5..8:', np.round(-bias.mean(0), 5))
bb, _ = ols(X, bias)
print('   its per-prime regression coefficients (max abs over cells): %.2e' % np.abs(bb[1:]).max())

# 2. finite height: split band into 100-300 and 300-600
for (a, b) in ((100, 300), (300, 600), (100, 350), (350, 600)):
    Ssub = []
    for D in Ds:
        g = np.sort(z['D%d' % D]); g = g[(g >= a) & (g <= b)]; u = unfold(g, abs(D))
        s, _ = sigma2(u, lgrid, u0=unfold(float(a), abs(D)), u1=unfold(float(b), abs(D))); Ssub.append(s)
    Ssub = np.array(Ssub)
    Fs = law([abs(D) for D in Ds], a, b, 'full'); Dg = law([abs(D) for D in Ds], a, b, 'diag')
    gm, res, off = gamma_of(Ssub, Fs, Dg, X)
    c2, _, _ = gls_chi2(Ssub, Fs, X)
    print('2. band %d-%d: gamma %.3f +- %.3f ; GLS chi2 full %.1f ; p=2 gamma %.3f' % (a, b, gm, boot_gamma(Ssub, Fs, Dg, X), c2,
          (res[0] * off[0]).sum() / (off[0] ** 2).sum()))

# 3. p = 2 split: 4||D vs 8|D ; and D>0 vs D<0
is2 = np.array([D % 4 == 0 for D in Ds]); is8 = np.array([D % 8 == 0 for D in Ds])
print('3. counts: 4||D %d, 8|D %d' % ((is2 & ~is8).sum(), is8.sum()))
X2 = np.column_stack([X[:, 0], (is2 & ~is8).astype(float), is8.astype(float), X[:, 2:]])
Sp = S[:, pidx]; Fp = full[:, pidx]; Dp = diag[:, pidx]
bf, se = ols(X2, Sp); bF, _ = ols(X2, Fp); bD, _ = ols(X2, Dp)
for j, nm in ((1, '4||D'), (2, '8|D')):
    off = bF[j] - bD[j]; res = bf[j] - bD[j]
    g = (res * off).sum() / (off ** 2).sum()
    # bootstrap
    rng = np.random.default_rng(3); gs = []
    for b in range(500):
        i = rng.integers(0, n, n)
        b1, _ = ols(X2[i], Sp[i]); b2, _ = ols(X2[i], Fp[i]); b3, _ = ols(X2[i], Dp[i])
        o_ = b2[j] - b3[j]; r_ = b1[j] - b3[j]; gs.append((r_ * o_).sum() / (o_ ** 2).sum())
    print('   %s: gamma_2 %.3f +- %.3f; resid(fit - full) l=6..8: %s (se %s)' % (nm, g, np.std(gs, ddof=1),
          ' '.join('%+.4f' % v for v in (bf[j] - bF[j])[11:]), ' '.join('%.4f' % v for v in se[j][11:])))
sgn = np.array([1.0 if D > 0 else 0.0 for D in Ds])
for nm, m in (('D>0', sgn == 1), ('D<0', sgn == 0)):
    gm, res, off = gamma_of(S[m], full[m], diag[m], X[m])
    print('   %s (n=%d): gamma %.3f +- %.3f' % (nm, m.sum(), gm, boot_gamma(S[m], full[m], diag[m], X[m])))

# 4. off-diagonal part only at shifted L (diag kept at true L)
ksB = [abs(D) for D in Ds]
for dL in (-0.3, 0.3):
    Fs = law(ksB, 100, 600, 'full', dL); Dg = law(ksB, 100, 600, 'diag', dL)
    M = diag + (Fs - Dg)
    c2, _, _ = gls_chi2(S, M, X)
    o = readings(S, M, diag, X, primes)
    print('4. off part at L%+.1f (diag at true L): GLS chi2 %.1f (BK 85.9); gamma vs this off %.3f; p2 pred %s' % (dL, c2, o['gamma'], ' '.join('%+.4f' % v for v in o['p2p'])))
