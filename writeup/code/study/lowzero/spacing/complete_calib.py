"""Calibrate the extra completeness checks on the USED set R120 (and on the zeros with one zero / one close pair removed
on purpose): pointwise |N - theta/pi| (registered 1.5), largest unfolded gap (3.0), the mean of N - theta/pi over each
quarter of the band, and the number of unfolded gaps below 0.25 against its expectation."""
import sys, os, math
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import numpy as np
from scipy.special import loggamma
import used_check
def theta(t, q, a): return np.imag(loggamma(0.25 + a / 2 + 0.5j * t)) + t / 2 * math.log(q / math.pi)
def checks(g, D, lo=100, hi=600):
    q = abs(D); a = 0 if D > 0 else 1; t = np.linspace(lo, hi, 200)
    s = np.searchsorted(g, t) - (theta(t, q, a) - theta(0.0, q, a)) / math.pi
    b = g[(g >= lo) & (g <= hi)]; u = b / (2 * math.pi) * np.log(q * b / (2 * math.pi * math.e)); d = np.diff(u)
    qm = s.reshape(4, 50).mean(1)
    return np.abs(s).max(), d.max(), qm, (d < 0.25).sum(), len(d), d.min()
Z = used_check.load_R120(); rng = np.random.default_rng(3)
A = []; Bq = []; Cq = []
for k in sorted(Z):
    D = int(k[1:]); g = np.sort(Z[k]); A.append(checks(g, D))
    # remove one zero at a random place in the band; and remove the closest pair in the band
    ib = np.nonzero((g >= 100) & (g <= 600))[0]; j = rng.choice(ib); Bq.append(checks(np.delete(g, j), D))
    d = np.diff(g[ib]); jj = ib[np.argmin(d)]; Cq.append(checks(np.delete(g, [jj, jj + 1]), D))
for name, R in (('as stored', A), ('one zero removed', Bq), ('closest pair removed', Cq)):
    mx = np.array([r[0] for r in R]); gp = np.array([r[1] for r in R]); qm = np.array([r[2] for r in R]); ns = np.array([r[3] for r in R]); n = np.array([r[4] for r in R])
    print('%-22s max|N-theta/pi|: median %.2f max %.2f (>1.5: %d) | max gap: median %.2f max %.2f (>3: %d) | quarter means: range %+.2f..%+.2f, flagged |m|>0.6: %d | gaps<0.25 per 800: %.1f'
          % (name, np.median(mx), mx.max(), (mx > 1.5).sum(), np.median(gp), gp.max(), (gp > 3).sum(), qm.min(), qm.max(), (np.abs(qm) > 0.6).any(1).sum(), (ns / n * 800).mean()))
qm = np.array([r[2] for r in A]); print('as stored: quarter-mean sd %.3f, mean %+.3f' % (qm.std(), qm.mean()))
