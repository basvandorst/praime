"""Round 62 diagnostic (after the failed reading): include q = 7 in the finite-lam hidden-residue model, gaps taken mod
210 (105 classes) and the information projected to residues mod 30, for k = 1, 2. usage: r62b.py"""
import os, json, math, itertools
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
r = json.load(open(os.path.join(OUT, 'r62.json'))); lam = r['lambda']


def Hb(p):
    p = p[p > 0]; return float(-np.sum(p * np.log2(p)))


g = np.arange(2, 6000, 2)
w = np.bincount((g // 2) % 105, weights=np.exp(-g / lam), minlength=105)       # g = 2c mod 210
out = {}
for k in (1, 2):
    n = k + 1
    J30 = np.zeros((15,) * n)
    for seq in itertools.product(range(105), repeat=n):
        pos = np.cumsum([0] + [2 * x for x in seq])
        f = 1.0
        for q in (3, 5, 7):
            f *= q - len(set((pos % q).tolist()))
            if f <= 0:
                break
        if f <= 0:
            continue
        J30[tuple(x % 15 for x in seq)] += np.prod(w[list(seq)]) * f
    J = J30 / J30.sum(); past = J.sum(axis=-1); nxt = J.reshape(-1, 15).sum(axis=0)
    I = Hb(past.ravel()) + Hb(nxt) - Hb(J.ravel())
    out[str(k)] = {'model_q3_5_7': I, 'data': r['k'][str(k)]['data'], 'model_q3_5': r['k'][str(k)]['model_finite_lambda']}
    print(k, {a: round(b, 5) for a, b in out[str(k)].items()}, flush=True)
json.dump(out, open(os.path.join(OUT, 'r62b.json'), 'w'), indent=1)
