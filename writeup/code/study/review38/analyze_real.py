import numpy as np
from common import *
z = np.load('/home/bas/scratch/review-logic3/real.npz')
Ds, ls, L = z['Ds'], z['ls'], z['L']
n = len(Ds)
I = np.array([[1.0 if abs(D) % p == 0 else 0.0 for p in PR] for D in Ds])
X = np.column_stack([np.ones(n), I])
il1 = int(np.where(np.isclose(ls, 1))[0][0])

def pred_table(idx=None):
    idx = np.arange(n) if idx is None else idx
    P = np.zeros((len(ls), 6))
    for j, l in enumerate(ls):
        for i, p in enumerate(PR):
            sel = [k for k in idx if I[k, i]]
            P[j, i] = np.mean([dS(l, p, L[k]) for k in sel])
    return P

def fits(S, idx=None):
    idx = np.arange(n) if idx is None else idx
    F = np.zeros((len(ls), 6)); E = np.zeros((len(ls), 6))
    for j in range(len(ls)):
        b, se = ols(X[idx], S[idx, j]); F[j] = b[1:]; E[j] = se[1:]
    return F, E

P = pred_table()
print('pred l=1 at group mean L:', np.round(P[il1], 4))
print('pred at L=9.38:', [round(dS(1, p, 9.38), 4) for p in PR])
for key in ['dense', 'dense_lc', 'rand', 'trim', 'const']:
    S = z[key]; F, E = fits(S)
    sl1 = F[il1] @ P[il1] / (P[il1] @ P[il1])
    slall = (F * P).sum() / (P * P).sum()
    print(f'{key:9s} l=1 fit', np.round(F[il1], 4), '+-', np.round(E[il1], 4), f'slope1 {sl1:.3f} slope96 {slall:.3f}',
          f'pooled 7,11,13 {F[il1][3:].mean():.4f}')
S = z['dense']; F, E = fits(S)
for p_i, p in enumerate(PR):
    print(p, 'fit ', np.round(F[:, p_i], 4))
    print(p, 'pred', np.round(P[:, p_i], 4))
    print(p, 'se  ', np.round(E[:, p_i], 4))
j35 = int(np.where(np.isclose(ls, 3.5))[0][0]); j2 = int(np.where(np.isclose(ls, 2))[0][0])
print('p13 l=3.5 %.4f l=2 %.4f l=7.5 %.4f' % (F[j35, 5], F[j2, 5], F[-2, 5]))
print('p2 l=6.5/l=1', F[12, 0] / F[il1, 0])
# bootstrap over L-functions
rng = np.random.default_rng(7)
bs1, bs96, b13 = [], [], []
for b in range(1000):
    idx = rng.integers(0, n, n)
    Fb, _ = fits(S, idx); Pb = pred_table(idx)
    bs96.append((Fb * Pb).sum() / (Pb * Pb).sum())
    bs1.append(Fb[il1] @ Pb[il1] / (Pb[il1] @ Pb[il1]))
    b13.append(Fb[j35, 5])
print('bootstrap slope96 sd %.3f, 95%% [%.3f, %.3f]' % (np.std(bs96), *np.percentile(bs96, [2.5, 97.5])))
print('bootstrap slope(l=1,6 coeffs) sd %.3f, 95%% [%.3f, %.3f]' % (np.std(bs1), *np.percentile(bs1, [2.5, 97.5])))
print('bootstrap p13 l=3.5 sd %.4f' % np.std(b13))
# GLS-like slope: weight by OLS se
w = 1 / E**2
print('weighted slope96 %.3f' % ((w * F * P).sum() / (w * P * P).sum()))
# per-l slope
print('per-l slope', np.round([(F[j] @ P[j]) / (P[j] @ P[j]) for j in range(len(ls))], 2))
# residual correlation across l
R = np.array([S[:, j] - X @ ols(X, S[:, j])[0] for j in range(len(ls))])
C = np.corrcoef(R)
print('resid corr l=1 with l=0.5,1.5,2,4,8:', np.round(C[il1, [0, 2, 3, 7, 15]], 2))
ev = np.linalg.eigvalsh(C)[::-1]; print('effective number of independent l (sum ev)^2/sum ev^2: %.1f' % (ev.sum()**2 / (ev**2).sum()))
# confounds
sign = (Ds > 0).astype(float); logq = np.log(np.abs(Ds)); om = np.array([len(primes_of(D)) for D in Ds])
print('corr(2|D, D>0) %.3f; corr(2|D, log|D|) %.3f; corr(2|D, omega) %.3f; corr(2|D, L) %.3f' % (
    np.corrcoef(I[:, 0], sign)[0, 1], np.corrcoef(I[:, 0], logq)[0, 1], np.corrcoef(I[:, 0], om)[0, 1], np.corrcoef(I[:, 0], L)[0, 1]))
for extra, name in [(sign, 'sign'), (L, 'meanL'), (sign + 0, 'sign')]:
    pass
X2 = np.column_stack([X, sign, L])
for j in [il1, j2, j35, 12]:
    b, se = ols(X2, S[:, j]); b0, _ = ols(X, S[:, j])
    print(f'l={ls[j]}: base {np.round(b0[1:7],4)}  +sign+L {np.round(b[1:7],4)} sign {b[7]:.4f}+-{se[7]:.4f} L {b[8]:.4f}+-{se[8]:.4f}')
# with covariates: full slope96
F2 = np.zeros_like(F)
for j in range(len(ls)): F2[j] = ols(X2, S[:, j])[0][1:7]
print('slope96 with sign+L covariates %.3f' % ((F2 * P).sum() / (P * P).sum()))
# residual sign effect at l=1 for squarefree odd vs not
print('mean Sigma2(1) D>0 %.4f D<0 %.4f' % (S[sign == 1, il1].mean(), S[sign == 0, il1].mean()))
# without the constant shift: n zeros
print('nz range', z['nz'].min(), z['nz'].max())
