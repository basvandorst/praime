"""own X matrices for OO, HI, C (t-grid method, own code and own theta); saved to Xown_<des>_<sig>.npy.
usage: s1_X.py DES"""
import sys, math, time
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import numpy as np
from rlib import *

OUT = __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/'
des = sys.argv[1]
t0 = time.time()
P = np.load(FZ + 'pred_OO.npz'); Ds = [int(x) for x in P['Ds']]
C = np.load(FZ + 'pred_C.npz')
if des == 'OO':
    Z = load_group(R120F + ['r134_B.npz', 'r136_Dn.npz'])
    funcs = [(abs(D), 0 if D > 0 else 1, Z['D%d' % D]) for D in Ds]; W1, W2 = 100.0, 594.0
elif des == 'HI':
    Z = load_group(['r134_C.npz'])
    funcs = [(abs(D), 0 if D > 0 else 1, Z['D%d' % D]) for D in Ds[:181]]; W1, W2 = 600.0, 994.0
else:
    Z = load_group(CPLXF + R120F); funcs = []
    for lab in [str(x) for x in C['labels']]:
        if lab.startswith('D'):
            D = int(lab[1:]); funcs.append((abs(D), 0 if D > 0 else 1, Z[lab]))
        else:
            N, m = int(lab.split('_')[0][1:]), int(lab.split('_')[1][1:])
            tab = conrey_table(N, m); a = 0 if abs(tab[N - 1] - 1) < 1e-9 else 1
            z = Z['N%d_m%d' % (N, m)]
            g = np.sort(-z[z < 0]) if lab.endswith('_conj') else np.sort(z[z > 0])
            funcs.append((N, a, g))
    W1, W2 = 106.0, 494.0
DT = 0.02
tg = W1 + DT * np.arange(int(round((W2 - W1) / DT)) + 1)
for sig in (0.5, 1.0):
    s = sig / math.sqrt(2); h = int(math.ceil(7 * s / DT))
    G = np.zeros((len(funcs), len(tg)), complex)
    for r, (q, a, g) in enumerate(funcs):
        b = g[(g > W1 - 7 * s) & (g < W2 + 7 * s)]
        ph = np.exp(2j * theta(b, q, a))
        c0 = np.floor((b - W1) / DT).astype(int)
        idx = c0[:, None] + np.arange(-h, h + 2)[None, :]
        ok = (idx >= 0) & (idx < len(tg))
        val = ph[:, None] * np.exp(-(W1 + DT * idx - b[:, None]) ** 2 / (2 * s * s))
        G[r] = np.bincount(idx[ok], weights=val[ok].real, minlength=len(tg)) + 1j * np.bincount(idx[ok], weights=val[ok].imag, minlength=len(tg))
    # trapezoid in t
    wt = np.full(len(tg), DT); wt[0] = wt[-1] = DT / 2
    X = (G * wt) @ G.conj().T / (math.sqrt(math.pi) * s * (W2 - W1))
    np.save(OUT + 'Xown_%s_%.1f.npy' % (des, sig), X)
    del G
    print('%s sigma %.1f: %d functions, X done %.0fs; X mean offdiag %.4f' % (des, sig, len(funcs), time.time() - t0,
          (X.real.sum() - np.trace(X.real)) / (len(funcs) ** 2 - len(funcs))), flush=True)
