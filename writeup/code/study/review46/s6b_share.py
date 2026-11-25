"""share of the law's variance over the NN pairs carried by primes in ranges (T -> inf diagonal form)."""
import sys, math
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__file__))
from mylib import *
S2 = np.load(OUT + '/s2_pairs.npz'); I, J = S2['I'], S2['J']
X = np.load(OUT + '/kron.npy').astype(float); P = np.load(OUT + '/primes.npy'); lp = np.log(P.astype(float))
A = np.zeros((16, len(P))); B = np.zeros((16, len(P)))
for m in range(1, 40):
    c = 1.0 / (m * m * P.astype(float) ** m)
    t = c[None, :] * np.sin(m * W[:, None] * lp[None, :] / 2) ** 2
    if m % 2: A += t
    else: B += t
A *= 2 / math.pi ** 2; B *= 2 / math.pi ** 2
def part(mask, k):
    Xm = X[:, mask]; M = (Xm * A[k, mask]) @ Xm.T + (Xm ** 2 * B[k, mask]) @ (Xm ** 2).T
    return M[I, J]
for k in (0, 1, 4, 8, 15):
    tot = part(P > 0, k); rest = part(P > 7, k)
    out = []
    for lo, hi in ((2, 8), (8, 80), (80, 300), (300, 16000), (16000, 10 ** 7)):
        msk = (P >= lo) & (P < hi); v = part(msk, k)
        out.append('%d-%d: %.4f/%.4f' % (lo, hi, v.var() / tot.var(), v.var() / rest.var() if lo >= 8 else float('nan')))
    print('w=%.1f share of var over pairs (FULL / REST):  ' % W[k] + '  '.join(out))
