import os, numpy as np, math
from sympy import primefactors
LZ=os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/data/riemann/lzeros/'
PR=(2,3,5,7,11,13)
def unfold(t,q): return t/(2*math.pi)*np.log(q*t/(2*math.pi*math.e))
def stt(t,q):
    s=np.diff(unfold(t,q)); r=np.minimum(s[1:],s[:-1])/np.maximum(s[1:],s[:-1]); return np.var(s),r.mean()
z=np.load(LZ+'r138_Bp.npz')
bands=[(100,350),(350,600),(600,800),(800,1000)]
Ds=[int(k[1:]) for k in z.files]; V=np.zeros((len(Ds),4)); R=np.zeros((len(Ds),4))
for i,D in enumerate(Ds):
    g=z['D%d'%D]
    for j,(lo,hi) in enumerate(bands):
        t=g[(g>=lo)&(g<=hi)]; V[i,j],R[i,j]=stt(t,abs(D))
X=np.array([[1.0]+[float(abs(d)%p==0) for p in PR]+[sum(1/p for p in primefactors(abs(d)) if p>13)] for d in Ds])
def res(y): return y-X@np.linalg.lstsq(X,y,rcond=None)[0]
for nm,M in (('Var',V),('r~',R)):
    Rr=np.array([res(M[:,j]) for j in range(4)])
    print(nm,'correlation matrix of indicator+largeprime residuals across bands',bands); print(np.round(np.corrcoef(Rr),2))
