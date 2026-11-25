import os, json, math, numpy as np
from scipy.special import loggamma
OUT = os.path.dirname(os.path.abspath(__file__))+'/'
def theta(t, Q, mu):
    t = np.asarray(t, float); d = len(mu)
    return 0.5 * t * math.log(Q / math.pi ** d) + sum(np.imag(loggamma((0.5 + m + 1j * t) / 2)) for m in mu)
a = json.load(open(OUT + 'r7f_n13_T400.json'))['res']; b = json.load(open(OUT + 'r7g_n13_400_600_fixed.json'))['res']
S = []
for k in ('N13', 'N13d'):
    S.append((k, np.sort(np.concatenate([a[k]['zeros'], b[k]['zeros']])), 600.0, 13.0, [0.5, 1.5], complex(*a[k]['eps']), 1.0215613, 1.1210163))
consts = {'F24a': (1.0, [11.5, 12.5], -0.0857724, 1.8846368, 'r10_f24a_T400.json'), 'F24b': (1.0, [11.5, 12.5], 0.5270993, 1.5755353, 'r10g_F24b_T400.json')}
pj = json.load(open(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/universal/predictions.json'))['Lfunctions']
for k in ('W23a', 'W23b', 'W29a', 'W29b', 'W31a', 'W31b'):
    cj = json.load(open(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/universal/constants_%s.json' % k))
    consts[k] = (pj[k]['Q'], pj[k]['mu'], pj[k]['c'], cj['A1'], 'r10g_%s_T400.json' % k)
for k, (Q, mu, c, A1, fn) in consts.items():
    r = json.load(open(OUT + fn)); S.append((k, np.array(r['zeros']), float(r['T']), Q, mu, 1 + 0j, c, A1))
rows = []
for k, g, T, Q, mu, eps, c, A1 in S:
    e = np.exp(2j * theta(g, Q, mu)); G = e.sum(); W = T / (2 * math.pi)
    vr = []
    for K in (4, 6, 8, 12, 20, 30):
        edges = np.linspace(0, T, K + 1); Gk = np.array([e[(g > edges[j]) & (g <= edges[j + 1])].sum() for j in range(K)])
        ph = -eps / abs(eps); vr.append((Gk * np.conj(ph)).real.var(ddof=1) * K)
    var = float(np.mean(vr))
    rows.append(dict(k=k, G=G, W=W, eps=eps, c=c, A1=A1, var=var))
    print('%-5s c=%+.4f  G=%8.2f%+8.2fi  law %8.2f%+8.2fi  slope %.3f +- %.3f' % (k, c, G.real, G.imag, (-eps * c * W).real, (-eps * c * W).imag, (G * np.conj(-eps * c * W)).real / abs(c * W) ** 2, math.sqrt(var) / abs(c * W)))
num = sum((r['G'] * np.conj(-r['eps'] * r['c'] * r['W'])).real / r['var'] for r in rows); den = sum(abs(r['c'] * r['W']) ** 2 / r['var'] for r in rows)
print('10 previously dropped: pooled inverse-variance slope %.3f +- %.3f' % (num / den, 1 / math.sqrt(den)))
# regression Re(G conj(-eps)/|eps|) = b1 * c W + b0 * W  (amplitude follows c: b1 = 1, b0 = 0 ; "c = 1": b1 = 0, b0 = 1)
Xm = np.array([[r['c'] * r['W'], r['W']] for r in rows]); y = np.array([(r['G'] * np.conj(-r['eps'] / abs(r['eps']))).real for r in rows]); w = np.array([1 / r['var'] for r in rows])
A = Xm.T @ (Xm * w[:, None]); bvec = np.linalg.solve(A, Xm.T @ (y * w)); cov = np.linalg.inv(A)
print('regression on [c W, W]: b_c = %.3f +- %.3f , b_1 = %.3f +- %.3f   (law: 1, 0 ; c=1 for all: 0, 1)' % (bvec[0], math.sqrt(cov[0, 0]), bvec[1], math.sqrt(cov[1, 1])))
Xm2 = np.array([[r['c'] * r['W'], r['W'] / r['A1']] for r in rows]); A2 = Xm2.T @ (Xm2 * w[:, None]); b2 = np.linalg.solve(A2, Xm2.T @ (y * w)); cov2 = np.linalg.inv(A2)
print('regression on [c W, W/L(1,Ad)]: b_c = %.3f +- %.3f , b_noB = %.3f +- %.3f' % (b2[0], math.sqrt(cov2[0, 0]), b2[1], math.sqrt(cov2[1, 1])))
