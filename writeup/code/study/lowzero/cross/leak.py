"""Theory only (no zeros): the finite-range form of the diagonal law. On the registered grid the orbit cosines of
different prime powers are not exactly orthogonal (leakage) and the square of one has mean (1 + <cos 2x>)/2. Model:
dn_i(t) = sum_{n <= X0} c_n(w) chi_i(n) cos((t + w/2) ln n) + (orbits n > X0, T -> infinity diagonal), c_n = -(2/pi)
(Lambda(n)/ln n) n^-1/2 sin(w ln n/2). Prediction C_hyb = C_inf + [V K V' - sum_{n<=X0} c_n^2/2 psi(n)], K the exact grid
average of the cosine products. Prints the readings of the registered regression computed ON THE PREDICTIONS (C_inf
and C_hyb): the expected readings under H_law and the size of the finite-range effects. -> leak_out.txt"""
import os, sys, json, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from crosslib import *
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', '..')); LOW = os.path.join(REPO, 'writeup/code/study/lowzero')
PJ = json.load(open(os.path.join(LOW, 'triple', 'predictions_kappa3.json')))
OLD = PJ['R120']['D'] + PJ['B']['D'] + PJ['Dn']['D']; NEW = PJ['Fn']['D']
T1, T2, H = 100.0, 600.0, 0.02; W = 0.3 * np.arange(1, 17); PR = (2, 3, 5, 7)
X0 = int(sys.argv[1]) if len(sys.argv) > 1 else 300
def grid(w):
    M = int(math.floor((T2 - w - T1) / H + 1e-9)) + 1; return T1 + H * np.arange(M)
P = primes_upto(200000); PP = [p for p in P if p <= X0]
PW = []                                   # prime powers n = p^m <= X0: (n, p, m, index of p)
for j, p in enumerate(PP):
    m = 1
    while p ** m <= X0: PW.append((p ** m, p, m, j)); m += 1
lnn = np.array([math.log(n) for n, _, _, _ in PW]); mm = np.array([m for _, _, m, _ in PW]); jj = np.array([j for _, _, _, j in PW])
print('X0 = %d: %d prime powers' % (X0, len(PW)))
ALL = NEW + OLD; XA = kron_table(ALL, P); XA_small = XA[:, :len(PP)].astype(float)
CHIN = XA_small[:, jj] ** mm[None, :]                       # chi_i(n) = chi_i(p)^m
nN = len(NEW)
def AT(p, w):
    tg = grid(w); lp = math.log(p); A = 0.0; B = 0.0; m = 1
    while 1.0 / (m * m * p ** m) > 1e-14:
        base = 2 / math.pi ** 2 / (m * m * p ** m) * math.sin(m * w * lp / 2) ** 2 * (1 + float(np.mean(np.cos(2 * m * (tg + w / 2) * lp))))
        if m % 2: A += base
        else: B += base
        m += 1
    return A, B
ATAB = np.array([[AT(p, w) for p in PR] for w in W])
Cinf = pair_cov_full(XA, W, P)                               # (16, 732, 732), T -> infinity
DELTA = np.zeros_like(Cinf)
for k, w in enumerate(W):
    tg = grid(w); ph = np.outer(tg + w / 2, lnn); Cs = np.cos(ph)
    K = Cs.T @ Cs / len(tg)                                  # exact grid averages of cosine products
    c = -(2 / math.pi) / (mm * np.exp(lnn / 2)) * np.sin(w * lnn / 2)
    V = CHIN * c[None, :]
    DELTA[k] = V @ K @ V.T - (CHIN * (c ** 2 / 2)[None, :]) @ CHIN.T
REST = pair_cov_full(XA[:, 4:], W, P[4:])
def readings(Cv, I, J):
    psi = {p: XA[I, a].astype(float) * XA[J, a].astype(float) for a, p in enumerate(PR)}
    out = {'s': [], 'u': [], 'r': [], 'icp': []}; num = 0.0; den = 0.0; Fx = []; Fy = []
    for k in range(len(W)):
        Z = np.column_stack([np.ones(len(I))] + [psi[p] for p in PR] + [psi[p] ** 2 for p in PR] + [REST[k][I, J]])
        y = Cv[k][I, J]; c = np.linalg.lstsq(Z, y, rcond=None)[0]
        out['s'].append(c[1:5]); out['u'].append(c[5:9]); out['r'].append(c[9]); out['icp'].append(c[0])
        pred = sum(ATAB[k, a, 0] * psi[p] + ATAB[k, a, 1] * psi[p] ** 2 for a, p in enumerate(PR)) + REST[k][I, J]
        x = pred - pred.mean(); yy = y - y.mean(); num += (x * yy).sum(); den += (x * x).sum()
    s = np.array(out['s']); u = np.array(out['u'])
    g = {('gamma_%d' % p): float((s[:, a] * ATAB[:, a, 0]).sum() / (ATAB[:, a, 0] ** 2).sum()) for a, p in enumerate(PR)}
    g.update({('gammaB_%d' % p): float((u[:, a] * ATAB[:, a, 1]).sum() / (ATAB[:, a, 1] ** 2).sum()) for a, p in enumerate(PR)})
    g['gamma_rest'] = float(np.mean(out['r'])); g['gamma_full'] = float(num / den); g['max|intercept|'] = float(np.abs(out['icp']).max())
    return g
designs = {'NN (Fn x Fn)': np.triu_indices(nN, 1)}
I2, J2 = np.meshgrid(np.arange(nN), np.arange(nN, len(ALL)), indexing='ij'); designs['NO (Fn x old)'] = (I2.ravel(), J2.ravel())
I3, J3 = np.triu_indices(len(OLD), 1); designs['OO (old x old)'] = (I3 + nN, J3 + nN)
res = {}
for name, (I, J) in designs.items():
    a = readings(Cinf, I, J); b = readings(Cinf + DELTA, I, J)
    d = np.array([DELTA[k][I, J].std() for k in range(len(W))]); dm = np.array([DELTA[k][I, J].mean() for k in range(len(W))])
    print('\n%s: %d pairs' % (name, len(I)))
    print('  finite-range term Delta_ij: sd over pairs per w ' + ' '.join('%.4f' % v for v in d))
    print('                              mean over pairs per w ' + ' '.join('%+.4f' % v for v in dm))
    print('  readings on the T->inf law:     ' + ', '.join('%s %.4f' % (k, v) for k, v in a.items()))
    print('  readings on the finite-range law: ' + ', '.join('%s %.4f' % (k, v) for k, v in b.items()))
    res[name] = {'inf': a, 'hyb': b}
json.dump({'X0': X0, 'readings': res}, open(os.path.join(HERE, 'leak_%d.json' % X0), 'w'), indent=1)
