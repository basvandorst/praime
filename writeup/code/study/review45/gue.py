import numpy as np, json, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
from mylaw import kf3_P, nodes
ls = np.arange(0.5, 8.01, 0.5)
out = []
for l in ls:
    m = 120
    x, w = np.polynomial.legendre.leggauss(m)
    x = l*(x + 1)/2; w = l*w/2
    K = np.sinc(x[:, None] - x[None, :])
    A = np.sqrt(w)[:, None]*K*np.sqrt(w)[None, :]
    lam = np.linalg.eigvalsh(A)
    out.append((np.sum(lam*(1-lam)*(1-2*lam)), np.sum(lam*(1-lam)), 2*np.sum(lam**3)))
out = np.array(out)
d = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/triple/predictions_kappa3.json'))
k3q = kf3_P(ls, 10.0, [], 1.0, gue=True)
print(' l   k3_GUE(eig)   JSON         diff      | sig2 eig - JSON | kf3 eig - JSON | kf3 my quad - eig')
for i, l in enumerate(ls):
    print(f"{l:4.1f} {out[i,0]:+.6f} {d['gue_kappa3'][i]:+.6f} {out[i,0]-d['gue_kappa3'][i]:+.1e} | {out[i,1]-d['gue_sigma2'][i]:+.1e} | {out[i,2]-d['gue_kf3'][i]:+.1e} | {k3q[i]-out[i,2]:+.1e}")
