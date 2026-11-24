"""Validation stretch only: the lag covariance tau_j of the within-projected transition increments at long lags,
by FFT, summed over lag bins; errors from 10 sub-stretches (scatter of the per-sub-stretch bin sums / sqrt(10))."""
import numpy as np, math, json
P = np.load('data/val_primes.npy'); X0 = 150_000_000_000; SPAN = 100_000_000
NW = int(np.searchsorted(P, X0 + SPAN)); LAM = math.log(X0 + SPAN / 2)
BINS = [(1, 2), (2, 3), (3, 6), (6, 11), (11, 31), (31, 101), (101, 301), (301, 1001), (1001, 3001)]
JM = 3001
def proj_within(Pm, pi):
    C = len(pi); R = np.zeros((C * C, C * C)); V = np.zeros((C * C, C))
    for a in range(C):
        sp = np.sqrt(Pm[a]); R[a * C:(a + 1) * C, a * C:(a + 1) * C] = np.eye(C) - np.outer(sp, sp)
        for bb in range(C): V[a * C + bb, bb] = math.sqrt(pi[a] * Pm[a, bb])
    return R - R @ V @ np.linalg.pinv(V.T @ R @ V) @ V.T @ R
def autocov(x, jm):
    n = len(x); m = 1 << int(math.ceil(math.log2(2 * n)))
    f = np.fft.rfft(x, m); ac = np.fft.irfft(f * np.conj(f), m)[:jm]
    return ac / (n - np.arange(jm))
out = {'lam': LAM, 'bins': BINS}
for q in (3, 5, 7):
    r = [a for a in range(q) if math.gcd(a, q) == 1]; C = len(r); lut = np.full(q, -1); lut[r] = np.arange(C)
    c = lut[P % q]
    for k in (1, 2):
        a0 = c[:NW]; a1 = c[k:NW + k]
        X = np.zeros((C, C)); np.add.at(X, (a0, a1), 1); pi = X.sum(1) / X.sum(); Pm = X / X.sum(1, keepdims=True)
        PiW = proj_within(Pm, pi)
        w, U = np.linalg.eigh(PiW); U = U[:, w > 0.5]            # orthonormal basis of the within space
        Yv = np.zeros((C * C, C * C))
        for a in range(C):
            for bb in range(C):
                for b2 in range(C): Yv[a * C + bb, a * C + b2] -= Pm[a, b2] / math.sqrt(pi[a] * Pm[a, b2])
                Yv[a * C + bb, a * C + bb] += 1 / math.sqrt(pi[a] * Pm[a, bb])
        coord = Yv @ U                                            # within coordinates of each (a, b) increment
        idx = a0 * C + a1
        nsub = 10; L = NW // nsub
        sub = np.zeros((nsub, JM))
        for s in range(nsub):
            ii = idx[s * L:(s + 1) * L]
            for l in range(coord.shape[1]):
                x = coord[ii, l]; x = x - x.mean()
                sub[s] += autocov(x, JM)
        sub /= (C - 1) ** 2
        tot = sub.mean(0)
        res = []
        for lo, hi in BINS:
            v = sub[:, lo:hi].sum(1); res.append((float(v.mean()), float(v.std(ddof=1) / math.sqrt(nsub))))
        out['%d_%d' % (q, k)] = {'tau0': float(tot[0]), 'binsum': res}
        print('q=%d k=%d  tau0=%.4f  ' % (q, k, tot[0]) + '  '.join('[%d,%d):%+.5f(%.5f)' % (lo, hi, m, e) for (lo, hi), (m, e) in zip(BINS, res)), flush=True)
json.dump(out, open('data/val_tail.json', 'w'), indent=1)
