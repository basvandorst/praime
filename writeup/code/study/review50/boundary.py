import os, numpy as np, math
from scipy.special import loggamma
from sympy import primefactors
LZ=os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/data/riemann/lzeros/'
PR=(2,3,5,7,11,13)
def unfold(t,q): return t/(2*math.pi)*np.log(q*t/(2*math.pi*math.e))
def theta(t,q,a): t=np.asarray(t,float); return np.imag(loggamma(0.25+a/2+0.5j*t))+t/2*math.log(q/math.pi)
z=np.load(LZ+'r138_Bp.npz'); Ds=[int(k[1:]) for k in z.files]
rows=[]
for D in Ds:
    g=z['D%d'%D]; q=abs(D); a=0 if D>0 else 1
    # boundary at a zero-free point: use t=350 exactly
    tb=350.0
    ta=g[(g>=100)&(g<=tb)]; tbb=g[(g>tb)&(g<=600)]
    sa=np.diff(unfold(ta,q)); sb=np.diff(unfold(tbb,q))
    # the boundary-straddling spacing
    sx=unfold(tbb[0],q)-unfold(ta[-1],q)
    Sb=np.searchsorted(g,tb)-theta(tb,q,a)/math.pi
    Slo=np.searchsorted(g,100.0)-theta(100.0,q,a)/math.pi; Shi=np.searchsorted(g,600.0)-theta(600.0,q,a)/math.pi
    rows.append((np.var(sa),np.var(sb),sa.mean(),sb.mean(),sx,Sb,Slo,Shi,len(sa),len(sb)))
A=np.array(rows)
X=np.array([[1.0]+[float(abs(d)%p==0) for p in PR]+[sum(1/p for p in primefactors(abs(d)) if p>13)] for d in Ds])
res=lambda y: y-X@np.linalg.lstsq(X,y,rcond=None)[0]
va,vb=res(A[:,0]),res(A[:,1])
print('corr(va,vb) %.2f'%np.corrcoef(va,vb)[0,1])
for k,nm in ((2,'mean s band a'),(3,'mean s band b'),(4,'straddling spacing'),(5,'S(350)'),(6,'S(100)'),(7,'S(600)')):
    print('%-20s corr with va %+.2f  vb %+.2f'%(nm,np.corrcoef(va,A[:,k])[0,1],np.corrcoef(vb,A[:,k])[0,1]))
# Var centred on 1 instead of band mean
