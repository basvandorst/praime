"""Exploratory (not registered): higher harmonics sum e^{2ik theta(gamma)}, k = 1..4, over the zeros of zeta and of 138 elliptic curves:
only k = 1 has a term linear in the height."""
import os, numpy as np, math, json
R = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))
from scipy.special import loggamma
z = np.loadtxt(os.path.join(R, 'data', 'riemann', 'zeros_40000_flint.txt'), usecols=-1); print('zeta zeros', len(z), z[-1])
th = lambda t: np.imag(loggamma(0.25 + 0.5j * t)) - t / 2 * math.log(math.pi)
for k in (1, 2, 3, 4):
    G = np.cumsum(np.exp(2j * k * th(z))); T = z
    row = ['T=%d: %8.1f%+8.1fi' % (T[i], G[i].real, G[i].imag) for i in (999, 4999, 19999, 39999)]
    print('zeta k=%d ' % k + '  '.join(row) + ('   law k=1: -T/2pi = %.1f' % (-T[-1] / (2 * math.pi)) if k == 1 else '   sqrt(T/2pi) = %.1f' % math.sqrt(T[-1] / (2 * math.pi))))
# elliptic curves of round 144: harmonics k = 1, 2, 3 pooled
Z = np.load(os.path.join(R, 'data', 'study', 'r144', 'fresh_zeros.npz')); P = np.load(os.path.join(R, 'code', 'study', 'lowzero', 'ecross', 'pred_fresh.npz'))
ids = [str(s) for s in P['ids']]; N = P['N']; w = P['w']; c = P['c']
thE = lambda t, n: t / 2 * math.log(n) - t * math.log(2 * math.pi) + np.imag(loggamma(1 + 1j * t))
for k in (1, 2, 3):
    S = []
    for i, s in enumerate(ids):
        g = np.sort(Z[s]); g = g[(g > 20) & (g <= 400)]; S.append(np.sum(np.exp(2j * k * thE(g, N[i]))) * (w[i] ** k))
    S = np.array(S); print('curves k=%d: mean of w^k * sum over 20<gamma<=400: %.2f%+.2fi +- %.2f (sd per curve %.1f); k=1 law mean %.2f' % (k, S.real.mean(), S.imag.mean(), S.real.std() / math.sqrt(len(S)), S.real.std(), (-c * 380 / (2 * math.pi)).mean()))
    if k > 1: print('   slope of Re on c: %.3f +- %.3f' % (np.polyfit(c, S.real, 1)[0], S.real.std() / (c.std() * math.sqrt(len(S)))))
