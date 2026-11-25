"""complex design: is the 0.8-1.3 % excess of C1 the finite-T value of the one-point sums? (pred1 -> pred1 x r_i conj r_j);
HI with the two flagged sets included."""
import sys, math
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import numpy as np
from rlib import *

D_ = __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/'
C = np.load(FZ + 'pred_C.npz'); labs = [str(s) for s in C['labels']]
Z = load_group(CPLXF + R120F)
W1, W2 = 106.0, 494.0
r = []
for k, l in enumerate(labs):
    if l.startswith('D'):
        g = Z[l]
    else:
        z = Z[l.replace('_conj', '')]; g = np.sort(-z[z < 0]) if l.endswith('_conj') else np.sort(z[z > 0])
    b = g[(g > W1) & (g <= W2)]
    r.append(np.sum(np.exp(2j * theta(b, int(C['q'][k]), int(C['a'][k])))) / (-C['eps'][k] * C['kappa'][k] * (W2 - W1) / (2 * math.pi)))
r = np.array(r)
cx = np.array([not l.startswith('D') for l in labs])
print('one-point ratio G/(-eps kappa T/2pi) on [106, 494]: complex sides mean %.4f%+.4fi sd %.3f; real mean %.4f' % (r[cx].real.mean(), r[cx].imag.mean(), np.abs(r[cx] - r[cx].mean()).std(), r[~cx].real.mean()))
unit = {}; uid = np.array([unit.setdefault(l.replace('_conj', ''), len(unit)) for l in labs])
excl = np.array([l == 'N173_m2_conj' for l in labs])
for sig in (0.5, 1.0):
    X = np.load(D_ + 'Xown_C_%.1f.npy' % sig)
    I = np.r_[C['I_cc'], C['I_cr']]; J = np.r_[C['J_cc'], C['J_cr']]
    p = np.r_[C['pred_cc_%.1f' % sig], C['pred_cr_%.1f' % sig]]; p1 = np.r_[C['pred1_cc_%.1f' % sig], C['pred1_cr_%.1f' % sig]]
    keep = ~(excl[I] | excl[J]); I, J, p, p1 = I[keep], J[keep], p[keep], p1[keep]
    x = X[I, J]; one = np.ones(len(x), complex)
    p1a = p1 * r[I] * np.conj(r[J])
    c, se, _, _ = ols_jack([one, p], x, uid[I], uid[J])
    c2, se2, _, _ = ols_jack([one, p1a + (p - p1)], x, uid[I], uid[J])
    c3, se3, _, _ = ols_jack([one, p1a, p - p1], x, uid[I], uid[J])
    print('sigma %.1f: C1 slope %.4f%+.4fi +- %.4f | with actual one-point sums in the H=1 part: %.4f%+.4fi +- %.4f | split [H=1 actual, H-1]: %.4f +- %.4f, %.4f +- %.4f' % (
        sig, c[1].real, c[1].imag, se[1], c2[1].real, c2[1].imag, se2[1], c3[1].real, se3[1], c3[2].real, se3[2]))
# HI incl. flagged
P = np.load(FZ + 'pred_OO.npz'); I, J = P['I'], P['J']; m = (I < 181) & (J < 181)
for sig in (0.5,):
    X = np.load(D_ + 'Xown_HI_%.1f.npy' % sig); x = X[I[m], J[m]].real; p = P['pred_%.1f' % sig][m]
    c, se, _, _ = ols_jack([np.ones(len(x)), p], x, I[m], J[m])
    print('HI incl. the 2 flagged sets: slope %.4f +- %.4f (registered 0.9931 +- 0.0069)' % (c[1], se[1]))
