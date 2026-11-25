"""Calibration of the completeness rule for one-sided sets stored only for |gamma| >= 100 (the complex characters):
centred count s_c = s - mean(s), s(t) = #{gamma in [lo, t]} - (theta(t) - theta(lo))/pi, and the step statistic
max_k |mean(s[:k]) - mean(s[k:])| (k = 20..180 of 200 heights). False flags on the USED complete sets of Fn (real, cut at
100) and the pilot complex sides; detection of one removed zero and of one removed close pair (random position)."""
import math, numpy as np
from scipy.special import loggamma

def theta(t, q, a): t = np.asarray(t, float); return t / 2 * math.log(q / math.pi) + np.imag(loggamma(0.25 + a / 2 + 0.5j * t))

def stats(g, q, a, lo, hi):
    t = np.linspace(lo, hi, 200); b = g[(g >= lo) & (g <= hi)]
    s = np.searchsorted(b, t, 'right') - (theta(t, q, a) - theta(lo, q, a)) / math.pi
    sc = s - s.mean(); cs = np.cumsum(s); tot = cs[-1]
    st = max(abs(cs[k - 1] / k - (tot - cs[k - 1]) / (200 - k)) for k in range(20, 181))
    u = b / (2 * math.pi) * np.log(q * b / (2 * math.pi * math.e))
    return np.abs(sc).max(), st, np.diff(u).max()

z = np.load(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'data', 'riemann', 'lzeros')) + '/r138_Fn.npz')
sets = [(np.sort(z[k]), abs(int(k[1:])), 0 if int(k[1:]) > 0 else 1, 100.0, 600.0) for k in z.files if k != 'D1153']
for N, f, lo, hi in ((7, 'q7_k1.txt', 100.0, 1090.0), (13, 'q13_k1.txt', 100.0, 790.0)):
    zz = np.loadtxt(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'data', 'riemann', 'lzeros')) + '/' + f)
    for side in (1, -1):
        g = np.sort(side * zz[side * zz > 0]); r = stats(g, N, 1, lo, hi)
        print('pilot q=%d side %+d: max|s_c| %.2f, step %.2f, gap %.2f' % (N, side, *r))
rng = np.random.default_rng(5)
A = np.array([stats(*s_) for s_ in sets])
print('Fn complete sets (n=%d): max|s_c| max %.2f (99%% %.2f), step max %.2f (99%% %.2f), gap max %.2f' % (
    len(A), A[:, 0].max(), np.quantile(A[:, 0], .99), A[:, 1].max(), np.quantile(A[:, 1], .99), A[:, 2].max()))
for nrem in (1, 2):
    R = []
    for g, q, a, lo, hi in sets:
        b = np.where((g > lo + 5) & (g < hi - 5))[0]; i = rng.choice(b[:-1]); g2 = np.delete(g, [i, i + 1][:nrem]); R.append(stats(g2, q, a, lo, hi))
    R = np.array(R)
    for th_c, th_s in ((1.5, 1.0), (1.5, 0.8), (1.2, 0.8), (1.5, 0.7)):
        flag = (R[:, 0] > th_c) | (R[:, 1] > th_s) | (R[:, 2] > 3.0); ff = ((A[:, 0] > th_c) | (A[:, 1] > th_s) | (A[:, 2] > 3.0)).mean()
        print('removed %d zero(s): thresholds s_c %.1f, step %.1f: detected %.3f; false flags on complete sets %.3f' % (nrem, th_c, th_s, flag.mean(), ff))
# the round-138 rule itself (count from 0) on the same injected sets, for comparison
for nrem in (1, 2):
    det = []
    for g, q, a, lo, hi in sets:
        b = np.where((g > lo + 5) & (g < hi - 5))[0]; i = rng.choice(b[:-1]); g2 = np.delete(g, [i, i + 1][:nrem])
        t = np.linspace(lo, hi, 200); s = np.searchsorted(g2, t) - theta(t, q, a) / math.pi
        det.append(np.abs(s).max() > 1.5)
    print('round-138 rule (count from 0), removed %d: detected %.3f' % (nrem, np.mean(det)))
