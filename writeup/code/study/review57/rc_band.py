"""Band split of my own fit (review 155): the twin (and cousin) zeta coefficient and the chi_5 coefficient per frequency band in u.
Method: weighted increments -> remove nuisance (WLS projection) -> real FFT -> in each band, complex least squares of the twin
transform on the transforms of (theta - x) and the 95 character regressors (all free); error from the residual power in the band.
usage: rc_band.py rc_full.npz rc_fit_out.npz
"""
import sys, math
import numpy as np
D = np.load(sys.argv[1]); F = np.load(sys.argv[2])
r, CH, NUI, w = F['r'], F['CH'], F['NUI'], F['w']
X = D['X'].astype(np.float64); K = np.arange(2, 4717); dx = X[K] - X[K - 1]; sq = np.sqrt(X[K]); u = np.log(X[K])
S = 2 * 0.6601618158468695739
sw = np.sqrt(w)
def proj_out(v):
    A = NUI * sw[:, None]; b = v * sw
    return b - A @ np.linalg.lstsq(A, b, rcond=None)[0]
du = 1 / 2048.0
n = len(u)
freqs = 2 * np.pi * np.fft.rfftfreq(n, d=du)          # angular frequency gamma
R = np.fft.rfft(proj_out(r)); C = np.fft.rfft(np.column_stack([proj_out(c) for c in CH.T]), axis=0)
k5 = 13   # placeholder, recomputed below
# locate chi_5 column: it is the real even character of conductor 5 -> find column with largest correlation to the Legendre-5 regressor
thc = D['thc'][K]
leg = np.array([0, 1, -1, -1, 1])
t5 = sum(leg[c] * thc[:, [k for k in range(420) if k % 5 == c]].sum(1) for c in (1, 2, 3, 4)) / sq
corr = [abs(np.corrcoef(t5, c)[0, 1]) for c in CH.T]
k5 = int(np.argmax(corr)); assert corr[k5] > 0.999
for label, h in (('twins', 2), ('cousins', 4)):
    y = (D['P%d' % h][K] - S * dx) / sq
    Y = np.fft.rfft(proj_out(y))
    print('\n%s: band       zeta coef            chi_5 coef' % label)
    for lo, hi in ((2, 100), (100, 300), (300, 1000), (1000, 2000), (2000, 4000), (4000, 6434)):
        m = (freqs >= lo) & (freqs < hi)
        Zb = np.column_stack([R[m], C[m]])
        # real-valued least squares on stacked real/imag parts
        A = np.vstack([Zb.real, Zb.imag]); b = np.concatenate([Y[m].real, Y[m].imag])
        if A.shape[0] <= A.shape[1] + 5:
            # too few bins for all characters: keep zeta and chi_5 only
            Zb = np.column_stack([R[m], C[m][:, k5]]); A = np.vstack([Zb.real, Zb.imag]); kk = 1
        else:
            kk = 1 + k5
        beta, *_ = np.linalg.lstsq(A, b, rcond=None)
        res = b - A @ beta; s2 = res @ res / (len(b) - A.shape[1]); cov = s2 * np.linalg.inv(A.T @ A)
        print('  [%4d,%4d) %4d bins  %7.3f +- %.3f      %7.3f +- %.3f' % (lo, hi, m.sum(), beta[0], math.sqrt(cov[0, 0]), beta[kk], math.sqrt(cov[kk, kk])))
