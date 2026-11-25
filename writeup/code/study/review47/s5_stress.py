"""stress tests on the OO design with own X: split-half check of the jackknife, node bootstrap, subsets (heights, psi(2),
conductor blocks, sigma), and complex design sides / exclusions."""
import sys, math, time
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import numpy as np
from rlib import *

D_ = __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/'
t0 = time.time()
P = np.load(FZ + 'pred_OO.npz'); Ds = np.array([int(x) for x in P['Ds']]); I, J = P['I'], P['J']
q = np.abs(Ds); a = (Ds < 0).astype(int)


def slope(x, p, m=None):
    if m is not None: x, p = x[m], p[m]
    A = np.column_stack([np.ones(len(p)), p]); return np.linalg.lstsq(A, x, rcond=None)[0][1]


def xgrid(funcs, sig, W1, W2, DT=0.02):
    tg = W1 + DT * np.arange(int(round((W2 - W1) / DT)) + 1)
    s = sig / math.sqrt(2); h = int(math.ceil(7 * s / DT))
    G = np.zeros((len(funcs), len(tg)), complex)
    for r, (qq, aa, g) in enumerate(funcs):
        b = g[(g > W1 - 7 * s) & (g < W2 + 7 * s)]
        ph = np.exp(2j * theta(b, qq, aa)); c0 = np.floor((b - W1) / DT).astype(int)
        ix = c0[:, None] + np.arange(-h, h + 2)[None, :]; ok = (ix >= 0) & (ix < len(tg))
        val = ph[:, None] * np.exp(-(W1 + DT * ix - b[:, None]) ** 2 / (2 * s * s))
        G[r] = np.bincount(ix[ok], weights=val[ok].real, minlength=len(tg)) + 1j * np.bincount(ix[ok], weights=val[ok].imag, minlength=len(tg))
    wt = np.full(len(tg), DT); wt[0] = wt[-1] = DT / 2
    return (G * wt) @ G.conj().T / (math.sqrt(math.pi) * s * (W2 - W1))


for sig in (0.5, 1.0):
    X = np.load(D_ + 'Xown_OO_%.1f.npy' % sig); x = X[I, J].real; p = P['pred_%.1f' % sig]; p1 = P['pred1_%.1f' % sig]
    one = np.ones(len(x))
    c, se, _, _ = ols_jack([one, p], x, I, J)
    print('sigma %.1f: slope %.4f, jackknife se %.4f' % (sig, c[1], se[1]))
    # split halves: independent function sets
    rng = np.random.default_rng(7); d = []; ses = []
    for rep in range(60):
        perm = rng.permutation(546); A_ = np.zeros(546, bool); A_[perm[:273]] = True
        mA = A_[I] & A_[J]; mB = ~A_[I] & ~A_[J]
        d.append(slope(x, p, mA) - slope(x, p, mB))
        if rep < 3:
            _, sA, _, _ = ols_jack([one, p], x, I, J, mA); _, sB, _, _ = ols_jack([one, p], x, I, J, mB); ses.append(math.hypot(sA[1], sB[1]))
    print('  split halves (60 random splits): sd(slope_A - slope_B) %.4f vs jackknife sqrt(seA^2 + seB^2) %.4f' % (np.std(d), np.mean(ses)))
    # node bootstrap
    bs = []
    for rep in range(200):
        wts = np.bincount(rng.integers(0, 546, 546), minlength=546).astype(float)
        w = wts[I] * wts[J]; m = w > 0
        A = np.column_stack([np.ones(m.sum()), p[m]]) * np.sqrt(w[m])[:, None]
        bs.append(np.linalg.lstsq(A, x[m] * np.sqrt(w[m]), rcond=None)[0][1])
    print('  node bootstrap (200): sd %.4f' % np.std(bs))
    # subsets
    chi2 = np.array([kron(int(D), 2) for D in Ds]); ps2 = chi2[I] * chi2[J]
    blk = np.where(q <= 400, 0, np.where(q <= 700, 1, 2)); bI, bJ = blk[I], blk[J]
    subs = [('psi(2)=+1', ps2 == 1), ('psi(2)=-1', ps2 == -1), ('psi(2)=0', ps2 == 0),
            ('both |D|<=400', (bI == 0) & (bJ == 0)), ('both 400-700', (bI == 1) & (bJ == 1)), ('both 700-1000', (bI == 2) & (bJ == 2)),
            ('mixed blocks', bI != bJ), ('law < 0', p < 0), ('law > 0.2', p > 0.2),
            ('gcd(q1,q2) > 1', np.array([math.gcd(int(q[i]), int(q[j])) > 1 for i, j in zip(I, J)]))]
    for nm, m in subs:
        c, se, _, _ = ols_jack([one, p], x, I, J, m)
        c2, se2, _, _ = ols_jack([one, p1, p - p1], x, I, J, m)
        print('  %-16s rows %6d: slope %.4f +- %.4f (z %+.1f); H-1 part %.4f +- %.4f (z %+.1f)' % (nm, m.sum(), c[1], se[1], (c[1] - 1) / se[1], c2[2], se2[2], (c2[2] - 1) / se2[2]))
    print('  (%.0fs)' % (time.time() - t0), flush=True)

# heights: two half windows
Z = load_group(R120F + ['r134_B.npz', 'r136_Dn.npz'])
funcs = [(abs(int(D)), 0 if D > 0 else 1, Z['D%d' % D]) for D in Ds]
for (W1, W2) in ((100.0, 347.0), (347.0, 594.0)):
    for sig in (0.5,):
        X = xgrid(funcs, sig, W1, W2); x = X[I, J].real; p = P['pred_%.1f' % sig]; p1 = P['pred1_%.1f' % sig]; one = np.ones(len(x))
        c, se, _, _ = ols_jack([one, p], x, I, J); c2, se2, _, _ = ols_jack([one, p1, p - p1], x, I, J)
        print('window [%.0f, %.0f] sigma %.1f: slope %.4f +- %.4f; [H=1, H-1] %.4f +- %.4f, %.4f +- %.4f' % (W1, W2, sig, c[1], se[1], c2[1], se2[1], c2[2], se2[2]), flush=True)
print('(%.0fs)' % (time.time() - t0))

# complex design: sides and exclusion
C = np.load(FZ + 'pred_C.npz'); labs = [str(s) for s in C['labels']]
unit = {}; uid = np.array([unit.setdefault(l.replace('_conj', ''), len(unit)) for l in labs])
cj = np.array([l.endswith('_conj') for l in labs]); isr = np.array([l.startswith('D') for l in labs])
for sig in (0.5,):
    X = np.load(D_ + 'Xown_C_%.1f.npy' % sig)
    I2 = np.r_[C['I_cc'], C['I_cr']]; J2 = np.r_[C['J_cc'], C['J_cr']]; p = np.r_[C['pred_cc_%.1f' % sig], C['pred_cr_%.1f' % sig]]
    x = X[I2, J2]; one = np.ones(len(x), complex)
    excl = np.array([l == 'N173_m2_conj' for l in labs]); keep = ~(excl[I2] | excl[J2])
    for nm, m in (('all incl. flagged N173 conj side', np.ones(len(x), bool)), ('registered (excl.)', keep),
                  ('chi x chi\' (both gamma>0)', keep & ~cj[I2] & ~cj[J2] & ~isr[J2] & (uid[I2] != uid[J2])),
                  ('conj x conj\' (both gamma<0)', keep & cj[I2] & cj[J2] & (uid[I2] != uid[J2])),
                  ('chi x conj chi\' (mixed sides)', keep & (cj[I2] != cj[J2]) & ~isr[J2] & (uid[I2] != uid[J2])),
                  ('complex(+) x real', keep & ~cj[I2] & isr[J2]), ('complex(-) x real', keep & cj[I2] & isr[J2])):
        c, se, sei, _ = ols_jack([one, p], x, uid[I2], uid[J2], m)
        print('C %-32s rows %5d: slope %.4f%+.4fi +- %.4f, %.4f' % (nm, m.sum(), c[1].real, c[1].imag, se[1], sei[1]))
print('done (%.0fs)' % (time.time() - t0))
