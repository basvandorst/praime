import numpy as np, json, math, cypari2
from scipy.special import loggamma
pari=cypari2.Pari(); pari.allocatemem(300*10**6)
z=np.load(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/data/study/r144/fresh_zeros.npz')
fs=json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/deg2/fresh_set.json'))['curves']
def th(N,t): t=np.asarray(t,float); return t/2*math.log(N)-t*math.log(2*math.pi)+np.imag(loggamma(1+1j*t))
for c in fs:
    if c['id'] not in ('N2357_0_m1_1_m3_1','N427_1_0_1_m8_7','N2699_0_0_1_2_2'): continue
    E=pari.ellinit(c['ainv']); gr=pari.ellglobalred(E); 
    print(c['id'], 'N file', c['N'], 'cond PARI', gr[0], 'rootno', pari.ellrootno(E), 'anrank', pari.ellanalyticrank(E))
    g=np.sort(z[c['id']])
    print(' n zeros', len(g), 'min', g[:3], 'dups', np.sum(np.diff(g)<1e-6))
    t=np.linspace(20,400,20)
    for r in (0,1):
        s=np.searchsorted(g,t,side='right')-th(c['N'],t)/math.pi+r/2
        print(' r=%d s:'%r, np.round(s,2))
    t=np.array([0.5,1,2,3,5,8,10,15,20])
    print(' small t s(r=1):', np.round(np.searchsorted(g,t,side='right')-th(c['N'],t)/math.pi+0.5,2))
