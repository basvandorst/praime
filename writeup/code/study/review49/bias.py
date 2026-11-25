# Estimator bias of the registered Sigma2 estimator on 831-point segments of the GUE bulk (sine process),
# Dumitriu-Edelman tridiagonal beta = 2 model, M = 2400, central 831 eigenvalues unfolded by the semicircle.
import sys, math, json, numpy as np
from scipy.linalg import eigvalsh_tridiagonal
LS = np.arange(0.5, 10.01, 0.5)
M = 2400; NPT = 831; nsamp = int(sys.argv[1]); seed = int(sys.argv[2])
rng = np.random.default_rng(seed)


def gue_s2(ls=LS, h=0.0005):
    out = []
    for l in ls:
        r = np.arange(0, l + h / 2, h); f = (l - r) * np.sinc(r) ** 2
        out.append(l - 2 * (h / 3) * (f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum()))
    return np.array(out)


def Fsc(lam):
    x = np.clip(lam / (2 * math.sqrt(M)), -1, 1)
    return M * (0.5 + (x * np.sqrt(1 - x * x) + np.arcsin(x)) / math.pi)


G = gue_s2(); acc = []
for s in range(nsamp):
    d = rng.normal(0, math.sqrt(2), M) / math.sqrt(2)
    e = np.sqrt(rng.chisquare(2 * np.arange(M - 1, 0, -1))) / math.sqrt(2)
    i0 = M // 2 - NPT // 2
    ev = eigvalsh_tridiagonal(d, e, select='i', select_range=(i0, i0 + NPT - 1))
    u = Fsc(np.sort(ev))
    v = []
    for l in LS:
        x = np.arange(u[0], u[-1] - l, 0.01); c = np.searchsorted(u, x + l) - np.searchsorted(u, x); v.append(c.var())
    acc.append(v)
acc = np.array(acc)
json.dump({'s2': acc.tolist(), 'gue': G.tolist()}, open('bias_%d.json' % seed, 'w'))
m = acc.mean(0) - G; se = acc.std(0, ddof=1) / math.sqrt(nsamp)
print('n', nsamp, 'bias (mean - GUE) at l=1,2,4,6,8,10:', ' '.join('%+.4f' % m[i] for i in (1, 3, 7, 11, 15, 19)))
print('se', ' '.join('%.4f' % se[i] for i in (1, 3, 7, 11, 15, 19)))
