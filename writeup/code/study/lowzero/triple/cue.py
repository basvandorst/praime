"""GUE check of kappa3(l) by simulation: Haar-random U(N) (QR of a complex Ginibre matrix with the phase fix), N = 800,
eigenangles unfolded to unit density. (a) periodic windows on the whole circle, step 0.01: the mean over matrices
estimates CUE kappa3(l; N); compared with the exact CUE(N) value from the kernel S_N(x) = sin(pi x)/(N sin(pi x/N))
(T3 = 2 S S S, T2 = -S^2, same quadrature as tplib) and the sine-kernel limit. (b) the data estimator (windows from the
first to the last point of a segment of 750 points, no wrap): its bias and its scatter between matrices.
usage: cue.py NMAT SEED -> cue_out_SEED.txt, cue_SEED.npz"""
import sys, os, math, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import tplib
NMAT = int(sys.argv[1]); SEED = int(sys.argv[2]); N = 800
LS = [0.5 * i for i in range(1, 17)]
rng = np.random.default_rng(SEED)
def haar(n):
    Z = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / math.sqrt(2)
    Qm, R = np.linalg.qr(Z); d = np.diag(R); return Qm * (d / np.abs(d))
per, seg, segs2 = [], [], []
t0 = time.time()
for it in range(NMAT):
    th = np.sort(np.mod(np.angle(np.linalg.eigvals(haar(N))), 2 * np.pi)); x = th * N / (2 * np.pi)
    xx = np.concatenate([x, x + N]); a, b, c = [], [], []
    for l in LS:
        x0 = np.arange(0, N, 0.01); cnt = (np.searchsorted(xx, x0 + l) - np.searchsorted(xx, x0)).astype(float)
        d = cnt - cnt.mean(); a.append((d ** 3).mean())
        u = x[:750]; w0 = np.arange(u[0], u[-1] - l, 0.01); cs = (np.searchsorted(u, w0 + l) - np.searchsorted(u, w0)).astype(float)
        ds = cs - cs.mean(); b.append((ds ** 3).mean()); c.append((ds ** 2).mean())
    per.append(a); seg.append(b); segs2.append(c)
per = np.array(per); seg = np.array(seg); segs2 = np.array(segs2)
np.savez('cue_%d.npz' % SEED, per=per, seg=seg, segs2=segs2)
# exact CUE(N) and sine-kernel values
q = tplib.Quad()
SN = lambda x: np.where(np.abs(x) < 1e-12, 1.0, np.sin(np.pi * x) / (N * np.sin(np.pi * x / N) + (np.abs(x) < 1e-12)))
kf_cue = q.kf3(2 * SN(q.T) * SN(q.S - q.T) * SN(q.S)); kf_sin = q.kf3(tplib.T3_gue(0 * q.T, q.T, q.S))
xg, wg = np.polynomial.legendre.leggauss(40)
def s2(l, S):
    r = l / 2 * (xg + 1); return l + 2 * np.sum(l / 2 * wg * (l - r) * (-S(r) ** 2))
IX = [int(round(l / q.WS)) - 1 for l in LS]
k3_cue = np.array([kf_cue[i] + 3 * s2(l, SN) - 2 * l for i, l in zip(IX, LS)])
k3_sin = np.array([kf_sin[i] + 3 * s2(l, lambda r: np.sinc(r)) - 2 * l for i, l in zip(IX, LS)])
f = lambda v: ' '.join('%+.5f' % x for x in v)
print('%d matrices N=%d, %.0fs' % (NMAT, N, time.time() - t0))
print('l                        ' + ' '.join('%8.1f' % l for l in LS))
print('sim periodic mean        ' + f(per.mean(0)))
print('   its se                ' + ' '.join('%8.5f' % x for x in per.std(0, ddof=1) / math.sqrt(NMAT)))
print('exact CUE(N) kappa3      ' + f(k3_cue))
print('sine kernel kappa3       ' + f(k3_sin))
print('z (sim - exact CUE)      ' + ' '.join('%+8.1f' % x for x in (per.mean(0) - k3_cue) / (per.std(0, ddof=1) / math.sqrt(NMAT))))
print('segment estimator mean   ' + f(seg.mean(0)) + '   (bias vs exact: ' + f(seg.mean(0) - k3_cue) + ')')
print('segment estimator sd     ' + ' '.join('%8.5f' % x for x in seg.std(0, ddof=1)) + '   (750 points; per-L-function data sd for comparison in noise_out)')
print('segment Sigma2 sd        ' + ' '.join('%8.5f' % x for x in segs2.std(0, ddof=1)))
