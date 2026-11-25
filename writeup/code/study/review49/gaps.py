import numpy as np, json, math
from scipy.special import loggamma
z=np.load(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/data/study/r144/fresh_zeros.npz')
fs={c['id']:c for c in json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/deg2/fresh_set.json'))['curves']}
def th(N,t): t=np.asarray(t,float); return t/2*math.log(N)-t*math.log(2*math.pi)+np.imag(loggamma(1+1j*t))
rows=[]; allg=[]
for k in z.keys():
    g=np.unique(np.sort(z[k])); b=g[(g>=20)&(g<=400)]; u=th(fs[k]['N'],b)/math.pi; d=np.diff(u)
    allg.append(d)
    for j in np.nonzero(d<0.06)[0]: rows.append((d[j], k, b[j], b[j+1]-b[j]))
rows.sort()
for r in rows: print('%.4f %-22s t=%.4f dt=%.5f'%r)
allg=np.concatenate(allg); print('n gaps', len(allg), 'mean', allg.mean())
# GUE expectation of gaps < s: approx (pi^2/9) s^3 for small s (Wigner surmise beta=2: (32/pi^2) s^2 exp(-4s^2/pi) -> cdf ~ 32 s^3/(3 pi^2))
for s in (0.03,0.06,0.1):
    print(s, (allg<s).sum(), 'Wigner-surmise expect', len(allg)*32*s**3/(3*math.pi**2))
