"""Window over-dispersion before/after subtracting the 2,3,5-line leakage (windows (20,30], ..., (140,150] from r21_leak.npz),
with the curve-level sign-flip null."""
import math, numpy as np
Z = np.load('r21_leak.npz'); Ts = list(Z['Ts']); O = Z['OBS']; L = Z['LEAK']
idx = [Ts.index(t) for t in range(20, 151, 10)]
rng = np.random.default_rng(5)
for lab, X in (('observed', O), ('observed - leak', O - L), ('leak alone', L)):
    W = np.diff(X[:, idx], axis=1); n = len(W)
    for part in ('real', 'imag'):
        Y = W.real if part == 'real' else W.imag
        st = lambda Z_: np.sum((Z_.mean(0) / (Z_.std(0, ddof=1) / math.sqrt(n))) ** 2)
        obs = st(Y); null = np.array([st(Y * rng.choice([-1, 1], size=(n, 1))) for _ in range(3000)])
        print('%-16s %s: chi2 %.1f on %d windows, null p = %.4f' % (lab, part, obs, Y.shape[1], np.mean(null >= obs)))
