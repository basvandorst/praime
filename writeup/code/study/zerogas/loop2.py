import numpy as np
from sine import K
def counts(s, m=None, kmax=8):
    """E(j; s), j < kmax: number of points in (0, s) is a sum of independent Bernoulli(lambda_i) (determinantal)."""
    m = m or int(30 + 4*s); x, w = np.polynomial.legendre.leggauss(m); t = (x + 1)*s/2; w = w*s/2; sw = np.sqrt(w)
    lam = np.clip(np.linalg.eigvalsh(sw[:, None]*K(t, t)*sw[None, :]), 0, 1)
    poly = np.zeros(kmax); poly[0] = 1.0
    for l in lam:
        poly[1:] = poly[1:]*(1 - l) + poly[:-1]*l; poly[0] *= (1 - l)
    return poly
def V(kmax=8, smax=22.0, n=600):
    x, w = np.polynomial.legendre.leggauss(n); s = (x + 1)*smax/2; w = w*smax/2
    E = np.array([counts(v, kmax=kmax) for v in s])                       # (n, kmax)
    out = []
    for k in range(1, kmax + 1):
        F = sum((k - j)*E[:, j] for j in range(k)); out.append(2*(w @ F) - k*k)
    return np.array(out)
v = V(); c = np.zeros(len(v))                                             # V(k) = Var of the k-th neighbour distance
Vx = np.concatenate([[0.0], v]); c0 = v[0]
for j in range(1, len(v)): c[j] = 0.5*(Vx[j + 1] - 2*Vx[j] + Vx[j - 1])
r = c/c0; r[0] = 1
print('theory: Var of k-th neighbour distance', np.round(v, 5))
print('theory: lag correlations r_1..r_7', np.round(r[1:], 5), ' 1 + 2 sum r_j:', np.round(1 + 2*np.cumsum(r[1:]), 4))
np.save('theory_r.npy', r)
off = {3: 267653395647.0, 4: 144176897509546973000.0, 5: 1370919909931995300000.0}
import re
for f in (3, 4, 5):
    lines = open('/home/bas/praime/writeup/data/riemann/odlyzko_zeros%d.txt' % f).read().split('\n')
    g = np.array([float(x) for x in lines if re.fullmatch(r'\s*\d+\.\d+\s*', x)])
    s = np.diff(g)*np.log(off[f]/(2*np.pi))/(2*np.pi)
    d = s - s.mean(); n = len(s)
    rr = np.array([np.mean(d[:-k]*d[k:])/np.var(s) for k in range(1, 8)])
    z = (rr - r[1:8])*np.sqrt(n)
    print('zeros%d: %d spacings, mean %.5f var %.5f; r_1..r_7 %s; z = (data - theory) sqrt(n): %s' % (f, n, s.mean(), s.var(), np.round(rr, 4), np.round(z, 2)))
