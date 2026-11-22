from rz import *
from tailchk import tail_num
import json
# off-diagonal O(L) = R(L) - diag on a fine L grid, via cumulative integral of x S(x)^2 in u = ln x
def R_curve(g,d,L0,L1,h=1e-4,span=14.0):
    c=coef(g,d)
    u=np.arange(L0-span,L1+h,h)
    S=S_grid(g,c,u)
    f=np.exp(2*(u-L1))*S**2 *np.exp(-u)  # x S^2 dx = e^{u} S^2 du, S already contains x^{1/2}; scale by e^{-2 L1}
    f=np.exp(u-2*L1)*S**2
    F=np.concatenate([[0],np.cumsum((f[1:]+f[:-1])*h/2)])
    idx=np.arange(int(span/h),len(u),200)
    L=u[idx]; R=2*F[idx]*np.exp(2*L1-2*L)/d
    return L,R
res={}
for N in (2419,10000):
    g=load(N); T=g[-1]
    for d in (0.1,0.01,0.001):
        if N==10000 and d==0.1: continue
        c=coef(g,d); diag=2*np.sum(np.abs(c)**2)/d
        L,R=R_curve(g,d,math.log(1e6),math.log(1e12), h=(1e-4 if N==2419 else 3e-5))
        off=R-diag
        res['%d_%g'%(N,d)]=dict(diag=diag,tail=tail_num(T,d),mean=float(off.mean()),std=float(off.std()),mn=float(off.min()),mx=float(off.max()),
            check={'%.0e'%X: float(np.interp(math.log(X),L,R)) for X in (1e8,1e9,1e10)})
        print(N,d,res['%d_%g'%(N,d)],flush=True)
        np.save('curve_%d_%g.npy'%(N,d),np.vstack([L,off]))
json.dump(res,open('avg.json','w'),indent=1)
