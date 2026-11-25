"""Noise and bias of the windowed Sigma2 estimator (np.var of counts, window starts every 0.01, segment not periodic) on
CUE eigenvalues of size n (Haar unitary by QR with phase fix), against the exact CUE number variance.
usage: cue_noise.py n samples seed OUT.json"""
import sys, json, math
import numpy as np
n, K, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
LS = [0.5 * i for i in range(1, 21)]
rng = np.random.default_rng(seed)
def haar(n):
    Z = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / math.sqrt(2)
    Q, R = np.linalg.qr(Z); d = np.diag(R); return Q * (d / np.abs(d))
def sig(u):
    v = []
    for l in LS:
        x = np.arange(u[0], u[-1] - l, 0.01); c = np.searchsorted(u, x + l) - np.searchsorted(u, x); v.append(np.var(c))
    return np.array(v)
# exact CUE Sigma2 for an arc of l mean spacings: l + 2 int_0^l (l - r)(R(r) - 1) dr, R - 1 = -(sin(pi r)/(n sin(pi r/n)))^2
def exact(l):
    r = np.linspace(0, l, 20001); rr = np.where(r > 0, r, 1.0)
    f = np.where(r > 0, -(np.sin(np.pi * rr) / (n * np.sin(np.pi * rr / n))) ** 2, -1.0)
    return l + 2 * np.trapezoid((l - r) * f, r)
S = []
for k in range(K):
    th = np.sort(np.angle(np.linalg.eigvals(haar(n))) % (2 * np.pi)); u = th * n / (2 * np.pi)
    S.append(sig(u))
S = np.array(S); ex = np.array([exact(l) for l in LS])
out = {'n': n, 'samples': K, 'l': LS, 'mean': S.mean(0).tolist(), 'sd': S.std(0, ddof=1).tolist(), 'exact': ex.tolist(),
       'bias': (S.mean(0) - ex).tolist(), 'bias_se': (S.std(0, ddof=1) / math.sqrt(K)).tolist()}
json.dump(out, open(sys.argv[4], 'w'))
print('CUE n=%d, %d samples' % (n, K))
for i, l in enumerate(LS):
    if l in (0.5, 1, 2, 4, 6, 8, 10):
        print('l=%4.1f exact %.4f mean %.4f bias %+.4f (+- %.4f) per-sample sd %.4f' % (l, ex[i], S[:, i].mean(), S[:, i].mean() - ex[i], S[:, i].std(ddof=1) / math.sqrt(K), S[:, i].std(ddof=1)))
