"""Round 113 (preregister/round113.txt): D(x) = pi(x;4,3) - pi(x;4,1) from sum_rho li(x^rho) (asymptotic Ei) plus the
squares (pi(sqrt x) - 1)/2, zeros of L(s, chi_-4) below 1500; regions with D < 0 and the share of primes to 1e9.
usage: r113.py"""
import os, sys, json, math
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'study'))
import common as C
H = os.path.dirname(os.path.abspath(__file__))
g = np.array([float(l) for l in open(os.path.join(H, '..', '..', '..', 'data', 'riemann', 'lzeros', 'D-4.txt')) if l.strip()]); g = g[g < float(os.environ.get('R113_T', 1500))]
rho = 0.5 + 1j * g
X = 1e9; du = 2 * math.pi / g.max() / 20
u = np.arange(math.log(100), math.log(X), du)
S = np.zeros_like(u)
for i0 in range(0, len(u), 2000):
    uu = u[i0:i0 + 2000]
    w = np.outer(uu, rho)                                   # rho ln x
    ser = 1 + 1 / w + 2 / w ** 2 + 6 / w ** 3 + 24 / w ** 4
    S[i0:i0 + 2000] = 2 * np.real(np.exp(w) / w * ser).sum(axis=1)     # both zeros of each conjugate pair
P = C.small_primes(int(math.sqrt(X)) + 10)
pis = np.searchsorted(P, np.floor(np.exp(u / 2)), side='right')
D = S + (pis - 1) / 2.0
neg = D < 0
idx = np.nonzero(np.diff(neg.astype(int)))[0]
if neg[0]: idx = np.concatenate([[-1], idx])
iv = [(float(np.exp(u[a + 1])), float(np.exp(u[b + 1]))) for a, b in zip(idx[::2], idx[1::2])]
cl = []
for a, b in iv:
    if cl and a < cl[-1][1] * 1.1: cl[-1][1] = b
    else: cl.append([a, b])
wt = np.exp(u) / u
share = float(wt[neg].sum() / wt.sum())
true = json.load(open(os.path.join(H, 'r112_clusters.json')))['clusters_true_to_1e9']
out = {'predicted_clusters': cl, 'true_clusters': true, 'predicted_share_to_1e9': share, 'true_share_to_1e9': 0.000150}
print('predicted clusters:', [(round(a), round(b)) for a, b in cl]); print('true clusters     :', true); print('share predicted %.5f, true 0.00015' % share)
json.dump(out, open(os.path.join(H, os.environ.get('R113_OUT', 'r113.json')), 'w'), indent=1)
