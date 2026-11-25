import numpy as np, json, math
from scipy.special import loggamma
Z=np.load(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/data/study/r144/fresh_zeros.npz')
def th(N,t): t=np.asarray(t,float); return t/2*math.log(N)-t*math.log(2*math.pi)+np.imag(loggamma(1+1j*t))
for cid,N in (('N2357_0_m1_1_m3_1',2357),('N427_1_0_1_m8_7',427)):
    for lab,g in (('stored',np.sort(Z[cid])),('dedup',np.unique(Z[cid]))):
        t=np.concatenate([[20],np.linspace(20,400,200),[400]])
        s=np.searchsorted(g,t,side='right')-th(N,t)/math.pi+0.5
        sm=np.convolve(s,np.ones(21)/21,mode='valid')
        b=g[(g>=20)&(g<=400)]; d=np.diff(th(N,b)/math.pi)
        print(cid,lab,'n band',len(b),'max|s| %.2f max|run.mean| %.2f mean s %+.2f max gap %.2f'%(np.abs(s).max(),np.abs(sm).max(),s.mean(),d.max()))
