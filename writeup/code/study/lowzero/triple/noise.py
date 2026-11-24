"""NOISE ONLY. For the zero sets already used in rounds 120-136, compute per L-function the window-count cumulants
kappa3(l) (third central moment of the counts) and the third factorial cumulant k[3](l) = m3 - 3 m2 + 2 m1, plus the
spacing statistics, regress each on [1, I(2|D), ..., I(13|D)] and keep ONLY the residuals.  No mean, no coefficient and
no comparison with any prediction is printed or saved (blind rule of the brief).
Windows exactly as r134.py: u = (g/2pi) ln(q g/(2 pi e)), windows from the first to the last zero of the band, step 0.01.
usage: noise.py [clean]  -> noise_out.txt / noise_clean_out.txt (stdout), noise_resid.npz / noise_clean_resid.npz
'clean' drops the three sets with a missing pair found by gapcheck.py (C: D=-195; Dn: D=-939; A: N173)."""
import sys
CLEAN = len(sys.argv) > 1 and sys.argv[1] == 'clean'
DROP = {'C': [-195], 'Dn': [-939], 'A': [173]} if CLEAN else {}
import os, json, math
import numpy as np
REPO = __import__('os').path.abspath(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', '..', '..', '..', '..')); LZ = REPO + '/writeup/data/riemann/lzeros'; LOW = REPO + '/writeup/code/study/lowzero'
OUT = os.path.dirname(os.path.abspath(__file__))
LS = [0.5 * i for i in range(1, 17)]; PR = (2, 3, 5, 7, 11, 13)
P = json.load(open(LOW + '/bk/predictions_per_lfunction.json')); PA = json.load(open(LOW + '/amp/predictions_alt.json'))


def stats(g, q):
    """per-set: k3(l), kf3(l), m2(l) on LS; spacing stats [var s, m3 s, lag-1 corr]; number of zeros"""
    u = g / (2 * math.pi) * np.log(q * g / (2 * math.pi * math.e))
    k3, kf3, m2s = [], [], []
    for l in LS:
        x = np.arange(u[0], u[-1] - l, 0.01); c = (np.searchsorted(u, x + l) - np.searchsorted(u, x)).astype(float)
        m1 = c.mean(); d = c - m1; m2 = (d * d).mean(); m3 = (d ** 3).mean()
        k3.append(m3); kf3.append(m3 - 3 * m2 + 2 * m1); m2s.append(m2)
    s = np.diff(u); sd = s - s.mean()
    sp = [float((sd ** 2).mean()), float((sd ** 3).mean()), float((sd[1:] * sd[:-1]).mean() / (sd ** 2).mean())]
    return np.array(k3), np.array(kf3), np.array(m2s), np.array(sp), len(u)


def load(name):
    if name in ('R120', 'C'):
        Ds = P['C']['D']
        if name == 'R120':
            Z = {}
            for f in ('r117_part0.npz', 'r117_part1.npz', 'r120_part0.npz', 'r120_part1.npz'):
                z = np.load(os.path.join(LZ, f)); Z.update({k: z[k] for k in z.files})
            lo, hi = 100, 600
        else:
            Z = dict(np.load(LZ + '/r134_C.npz')); lo, hi = 600, 1000
    elif name == 'B':
        Ds = P['B']['D']; Z = dict(np.load(LZ + '/r134_B.npz')); lo, hi = 100, 600
    elif name == 'Dn':
        Ds = PA['Dn']['D']; Z = dict(np.load(LZ + '/r136_Dn.npz')); lo, hi = 100, 600
    rows = []
    for D in Ds:
        g = np.sort(Z['D%d' % D]); g = g[(g >= lo) & (g <= hi)]
        rows.append((D, abs(D), stats(g, abs(D)), float(np.mean(np.log(abs(D) * g / (2 * math.pi))))))
    return rows


def load_A():
    Z = {}
    for i in range(4):
        z = np.load(os.path.join(LZ, 'r121_part%d.npz' % i)); Z.update({k: z[k] for k in z.files})
    keys = [r['key'] for r in json.load(open(LOW + '/r121.json'))['rows']]
    lo, hi = P['A']['height']; rows = []
    for k in keys:
        q = int(k.split('_')[0][1:]); acc = None; ws = 0
        for side in (1, -1):
            g = np.sort(side * Z[k][side * Z[k] > 0]); g = g[(g >= lo) & (g <= hi)]
            st = stats(g, q); n = st[4]
            acc = [n * a for a in st[:4]] if acc is None else [a + n * b for a, b in zip(acc, st[:4])]; ws += n
        rows.append((q, q, tuple([a / ws for a in acc]) + (ws,), float(np.log(q * 300 / (2 * math.pi)))))
    return rows


def resid(Y, X):
    c = np.linalg.lstsq(X, Y, rcond=None)[0]; return Y - X @ c


res = {}
print('per-L-function scatter of the statistics after removing the regression on [1, I(2|D) .. I(13|D)] (residual sd, dof-corrected)')
for name in ('R120', 'B', 'Dn', 'C', 'A'):
    rows = load_A() if name == 'A' else load(name)
    rows = [r for r in rows if r[0] not in DROP.get(name, [])]
    Q = [r[1] for r in rows]; pr = (2, 3, 5, 7) if name == 'A' else PR
    X = np.array([[1] + [int(q % p == 0) for p in pr] for q in Q], float); n = len(Q); dof = n - X.shape[1]
    K3 = np.array([r[2][0] for r in rows]); KF = np.array([r[2][1] for r in rows]); M2 = np.array([r[2][2] for r in rows])
    SP = np.array([r[2][3] for r in rows]); NZ = np.array([r[2][4] for r in rows])
    rK3, rKF, rM2, rSP = resid(K3, X), resid(KF, X), resid(M2, X), resid(SP, X)
    res[name] = dict(K3=rK3, KF=rKF, M2=rM2, SP=rSP, NZ=NZ, Q=np.array(Q), L=np.array([r[3] for r in rows]), dof=dof)
    sd = lambda r: np.sqrt((r ** 2).sum(0) / dof)
    print('\n[%s] %d L-functions, mean zeros per set %.0f, mean L %.3f' % (name, n, NZ.mean(), np.mean([r[3] for r in rows])))
    print('  l        ' + ' '.join('%7.1f' % l for l in LS))
    print('  sd k3    ' + ' '.join('%7.4f' % v for v in sd(rK3)))
    print('  sd k[3]  ' + ' '.join('%7.4f' % v for v in sd(rKF)))
    print('  sd Sig2  ' + ' '.join('%7.4f' % v for v in sd(rM2)))
    C = np.corrcoef(rK3.T); print('  corr of k3 residuals across l (rows l=1,2,4,8 vs all):')
    for l in (1, 2, 4, 8):
        i = LS.index(l); print('   l=%-3g ' % l + ' '.join('%+5.2f' % v for v in C[i]))
    Ck = np.corrcoef(rKF.T); print('  corr k[3] l=1,4 rows: ' + ' '.join('%+5.2f' % v for v in Ck[1]) + ' | ' + ' '.join('%+5.2f' % v for v in Ck[7]))
    print('  corr(k3, Sigma2) residuals at same l: ' + ' '.join('%+5.2f' % np.corrcoef(rK3[:, j], rM2[:, j])[0, 1] for j in range(len(LS))))
    ev = np.linalg.eigvalsh(np.corrcoef(rK3.T)); neff = ev.sum() ** 2 / (ev ** 2).sum()
    print('  effective number of independent lengths (k3): %.1f' % neff)
    print('  spacing stats residual sd: var(s) %.4f, m3(s) %.4f, lag1 corr %.4f' % tuple(sd(rSP)))
    # dependence of noise on zeros count: scaled sd * sqrt(NZ)
    print('  sd(k3)*sqrt(zeros/800) at l=1,2,4,8: ' + ' '.join('%.4f' % (sd(rK3)[LS.index(l)] * math.sqrt(NZ.mean() / 800)) for l in (1, 2, 4, 8)))
np.savez(os.path.join(OUT, 'noise_clean_resid.npz' if CLEAN else 'noise_resid.npz'), **{'%s_%s' % (k, kk): v for k, d in res.items() for kk, v in d.items()})
