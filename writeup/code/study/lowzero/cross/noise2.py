"""NOISE ONLY (used zero sets R120 + B + Dn): jackknife standard errors of the registered readings (readings.py, with the
swapped-term column from a crude predictor, primes <= 1e4) for old x old and for subsets that mimic the fresh designs.
Uses cache_C.npz from noise.py. Prints standard errors only; no reading value is printed or stored."""
import os, sys, json, glob, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from crosslib import *
from readings import readings
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', '..')); LOW = os.path.join(REPO, 'writeup/code/study/lowzero')
PJ = json.load(open(os.path.join(LOW, 'triple', 'predictions_kappa3.json')))
OLD = PJ['R120']['D'] + PJ['B']['D'] + PJ['Dn']['D']; N = len(OLD)
z = np.load(os.path.join(HERE, 'cache_C.npz')); C = z['C']; assert list(z['D']) == OLD
T1, T2, H = 100.0, 600.0, 0.02; W = 0.3 * np.arange(1, 17); PR = (2, 3, 5, 7)
def grid(w):
    M = int(math.floor((T2 - w - T1) / H + 1e-9)) + 1; return T1 + H * np.arange(M)
def AT(p, w):
    tg = grid(w); lp = math.log(p); A = 0.0; B = 0.0; m = 1
    while 1.0 / (m * m * p ** m) > 1e-14:
        base = 2 / math.pi ** 2 / (m * m * p ** m) * math.sin(m * w * lp / 2) ** 2 * (1 + float(np.mean(np.cos(2 * m * (tg + w / 2) * lp))))
        if m % 2: A += base
        else: B += base
        m += 1
    return A, B
TAB = np.array([[AT(p, w) for p in PR] for w in W]); A_T = TAB[:, :, 0]; B_T = TAB[:, :, 1]
P = primes_upto(200000); X = kron_table(OLD, P)
FULL = pair_cov_full(X, W, P); REST = pair_cov_full(X[:, 4:], W, P[4:])
# crude swapped-term predictor for the near pairs
SW = np.zeros((len(W), N, N)); idx = {D: i for i, D in enumerate(OLD)}; nsw = 0
for f in sorted(glob.glob(os.path.join(HERE, 'swap_OO_P10000_*.npz'))):
    s = np.load(f)
    for a, b, c in zip(s['D1'], s['D2'], s['C']):
        i, j = idx[int(a)], idx[int(b)]; SW[:, i, j] = c; SW[:, j, i] = c; nsw += 1
print('swap predictor (crude, primes <= 1e4) for %d near pairs' % nsw)
I, J = np.triu_indices(N, 1)
def run(name, sel):
    Ii, Jj = I[sel], J[sel]
    Cs = np.array([C[k][Ii, Jj] for k in range(len(W))]); F = np.array([FULL[k][Ii, Jj] for k in range(len(W))])
    R = np.array([REST[k][Ii, Jj] for k in range(len(W))]); S = np.array([SW[k][Ii, Jj] for k in range(len(W))])
    PSI = np.column_stack([X[Ii, a].astype(float) * X[Jj, a].astype(float) for a in range(4)])
    _, se, cov_s = readings(Cs, PSI, R, F, A_T, B_T, SWAP=S, I=Ii, J=Jj, jack=True)
    print('\n[%s] %d pairs, %d near pairs with a swap prediction (sum of squares of the swap predictor %.4f)' % (name, sel.sum(), int((S[0] != 0).sum()), float((S ** 2).sum())))
    print('  jackknife se: ' + ', '.join('%s %.4f' % (k, v) for k, v in se.items() if not isinstance(v, list)))
    print('  se of the intercepts per w: ' + ' '.join('%.5f' % v for v in se['icpt']))
    for a, p in enumerate(PR): print('  se of s_%d(w): ' % p + ' '.join('%.5f' % v for v in np.array(se['s'])[:, a]))
    return se, cov_s
out = {}
se, cov = run('old x old', np.ones(len(I), bool)); out['OO'] = {'se': se, 'cov_s': cov}
rng = np.random.default_rng(140)
for rep in range(2):
    new = np.zeros(N, bool); new[rng.choice(N, 186, replace=False)] = True
    se, cov = run('subset %d, 186 x 186' % rep, new[I] & new[J]); out['sub%d_NN' % rep] = {'se': se, 'cov_s': cov}
    se, cov = run('subset %d, 186 x 360' % rep, new[I] != new[J]); out['sub%d_NO' % rep] = {'se': se, 'cov_s': cov}
json.dump(out, open(os.path.join(HERE, 'noise2_se.json'), 'w'))
