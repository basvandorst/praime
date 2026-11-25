# Completeness check per curve: D(t) = #{gamma <= t} - theta(t)/pi - 1 = S(t) (should average to ~0 on every stretch).
# Report per curve the window means of D over 40-wide windows in [0, 400]; a missing or extra pair shows as a +-2 shift.
import os, numpy as np
from scipy.special import loggamma
SC = os.path.dirname(os.path.abspath(__file__))+'/'
CU = np.load(SC + 'rv_curves.npz')
Z = np.load(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/data/study/r144/fresh_zeros.npz')
tg = np.arange(5, 400, 0.01)
worst = []
for k, s in enumerate(CU['ids']):
    N = int(CU['N'][k]); z = np.sort(Z[s]); r = int(CU['rank'][k])
    th = tg / 2 * np.log(N) - tg * np.log(2 * np.pi) + loggamma(1 + 1j * tg).imag
    # zeros stored from t > 0 only (one side); the central zero(s) of order r are not in the list:
    D = np.searchsorted(z, tg, side='right') - th / np.pi
    wm = np.array([D[(tg >= a) & (tg < a + 40)].mean() for a in range(0, 400, 40) if ((tg >= a) & (tg < a + 40)).any()])
    worst.append((np.ptp(wm), s, np.round(wm, 2)))
worst.sort(key=lambda x: -x[0])
allm = np.array([w[2] for w in worst])
print('window means of D(t) (40-wide windows, 0..400): overall mean %.3f; per-curve spread (max-min of window means): median %.2f, max %.2f' % (
    allm.mean(), np.median([w[0] for w in worst]), worst[0][0]))
for w in worst[:6]:
    print('  ', w[1], 'spread %.2f' % w[0], w[2])
idx = [i for i, w in enumerate(worst) if w[1].startswith('N1671')]
print('  N1671 rank in spread list:', idx, worst[idx[0]][2] if idx else '')
