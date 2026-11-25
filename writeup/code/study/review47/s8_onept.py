"""Is the 'H = 1 part' of X the product of the two one-point sums? Replace pred1 by pred1 * (G_i/law_i) conj(G_j/law_j)
with G the sharp one-point sum over the X window; also residual sd of Re X vs sd of Im X (all OO pairs)."""
import sys, math
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import numpy as np
from rlib import *

D_ = __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/'
P = np.load(FZ + 'pred_OO.npz'); Ds = [int(x) for x in P['Ds']]; I, J = P['I'], P['J']
Z = load_group(R120F + ['r134_B.npz', 'r136_Dn.npz'])
W1, W2 = 100.0, 594.0
g1 = []
for D in Ds:
    g = Z['D%d' % D]; b = g[(g > W1) & (g <= W2)]
    g1.append(np.sum(np.exp(2j * theta(b, abs(D), int(D < 0)))) / (-kappa(abs(D)) * (W2 - W1) / (2 * math.pi)))
g1 = np.array(g1)
print('one-point ratio G/law on [100, 594] (sharp): mean %.4f, sd %.4f, mean Im %.4f' % (g1.real.mean(), g1.real.std(), g1.imag.mean()))
for sig in (0.5, 1.0):
    x = np.load(D_ + 'Xown_OO_%.1f.npy' % sig)[I, J]; p = P['pred_%.1f' % sig]; p1 = P['pred1_%.1f' % sig]
    one = np.ones(len(p)); p1a = (p1 * g1[I] * np.conj(g1[J])).real
    c, se, _, _ = ols_jack([one, p1, p - p1], x.real, I, J)
    ca, sea, _, _ = ols_jack([one, p1a, p - p1], x.real, I, J)
    cb, seb, _, _ = ols_jack([one, p1a + (p - p1)], x.real, I, J)
    print('sigma %.1f: [H=1, H-1] with pred1: %.4f +- %.4f, %.4f +- %.4f | with pred1 x one-point ratios: %.4f +- %.4f, %.4f +- %.4f | slope on (pred - pred1 + pred1 x ratios): %.4f +- %.4f' % (
        sig, c[1], se[1], c[2], se[2], ca[1], sea[1], ca[2], sea[2], cb[1], seb[1]))
    r = x.real - (c[0] + c[1] * p1 + c[2] * (p - p1))
    print('   sd(Re X - fit) %.4f, sd(Im X) %.4f, sd(law) %.4f, R^2 %.3f' % (r.std(), x.imag.std(), p.std(), 1 - r.var() / x.real.var()))
