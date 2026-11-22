"""Round 62 diagnostic 2: k = 1 with q = 3, 5, 7, 11 (gaps mod 2310, 1155 classes), projected to residues mod 30."""
import os, json, math
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
r = json.load(open(os.path.join(OUT, 'r62.json'))); lam = r['lambda']
g = np.arange(2, 8000, 2); M = 1155
w = np.bincount((g // 2) % M, weights=np.exp(-g / lam), minlength=M)
c1, c2 = np.meshgrid(np.arange(M), np.arange(M), indexing='ij'); c1 = c1.ravel(); c2 = c2.ravel()
x1 = 2 * c1; x2 = 2 * (c1 + c2)
f = np.ones(len(c1))
for q in (3, 5, 7, 11):
    R = np.stack([np.zeros_like(x1), x1 % q, x2 % q], 1); R.sort(1)
    nu = 1 + np.count_nonzero(np.diff(R, axis=1), axis=1)
    f *= q - nu
wt = w[c1] * w[c2] * f
J = np.zeros((15, 15)); np.add.at(J, (c1 % 15, c2 % 15), wt); J /= J.sum()
Hb = lambda p: float(-np.sum(p[p > 0] * np.log2(p[p > 0])))
I = Hb(J.sum(1)) + Hb(J.sum(0)) - Hb(J.ravel())
b = json.load(open(os.path.join(OUT, 'r62b.json')))
out = {'k1_model_q3_5_7_11': I, 'k1_model_q3_5_7': b['1']['model_q3_5_7'], 'k1_model_q3_5': r['k']['1']['model_finite_lambda'], 'k1_data': r['k']['1']['data']}
print({a: round(v, 5) for a, v in out.items()})
json.dump(out, open(os.path.join(OUT, 'r62c.json'), 'w'), indent=1)
