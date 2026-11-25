# Sign/amplitude of p's wave in the real zeros: cov(S(t), dN_p(t)) against Var(dN_p), good-prime L-functions only.
import os, numpy as np, math
from scipy.special import loggamma
LZ=os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/data/riemann/lzeros/'
def theta(t,q,a):
    t=np.asarray(t,float); return np.imag(loggamma(0.25+a/2+0.5j*t))+t/2*math.log(q/math.pi)
def kron(D,p):
    if D%p==0: return 0
    if p==2: return 1 if D%8 in (1,7) else -1
    return 1 if pow(D%p,(p-1)//2,p)==1 else -1
def dN(t,p,chi): return np.angle(1-chi*p**-0.5*np.exp(1j*t*math.log(p)))/math.pi
for f,lo,hi in [('r134_B.npz',100,600),('r134_C.npz',600,1000)]:
    z=np.load(LZ+f)
    for p in (2,3,5,7,11,13):
        num=0;den=0;cnt=0; numw=0
        for k in z.files:
            D=int(k[1:]); chi=kron(D,p)
            if chi==0: continue
            g=z[k]; q=abs(D); a=0 if D>0 else 1
            tt=np.linspace(lo,hi,20001)
            S=np.searchsorted(g,tt)-theta(tt,q,a)/math.pi
            w=dN(tt,p,chi); w2=dN(tt,p,-chi)
            num+=np.mean(S*w); den+=np.mean(w*w); numw+=np.mean(S*w2); cnt+=1
        th=sum(p**-m/m**2 for m in range(1,60))/(2*math.pi**2)
        print(f,p,'n=%d  <S dN_p>/<dN_p^2> = %.3f  (wrong sign chi: %.3f)   Var dN_p sample %.5f theory %.5f'%(cnt,num/den,numw/den,den/cnt,th))
