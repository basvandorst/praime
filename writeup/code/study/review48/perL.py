# per-L-function statistics from the stored values (own code)
import numpy as np, json, sympy, mpmath as mp, sys
sys.path.insert(0,'.')
z=np.load('../all.npz'); keys=sorted({k.rsplit('_',1)[0] for k in z.keys()})
k1=np.load('../k1_rows.npy'); k1d={(int(r[0])):r for r in k1}
out=[]
for k in keys:
    s,D=k.split('_'); D=int(D); q=abs(D)
    g=z[k+'_g']; ly=z[k+'_y']
    L=np.log(q*g/(2*np.pi)); y=ly-np.log(L)
    bad=sorted(sympy.factorint(q).keys())
    Aq=6/np.pi**2*np.prod([p/(p+1) for p in bad])
    X=np.exp(-2*ly)
    out.append(dict(set=s,D=D,q=q,bad=bad,n=len(g),mu=y.mean(),v=y.var(ddof=1),
        k3=((y-y.mean())**3).mean(), k4=((y-y.mean())**4).mean()-3*y.var()**2,
        lnL=np.log(L).mean(), invL=(1/L).mean(), S1=np.exp(2*ly).sum(), S1pred=k1d[D][2],
        Sm1=X.sum(), Sm1t=np.minimum(X,25.0).sum(), Xmax=X.max(), Pm1=Aq*500/(2*np.pi),
        Sm1_top5=np.sort(X)[-5:].sum(),
        lo=float(np.mean(g<350)), mu_lo=y[g<350].mean(), mu_hi=y[g>=350].mean(),
        v_lo=y[g<350].var(ddof=1), v_hi=y[g>=350].var(ddof=1),
        S1_lo=np.exp(2*ly[g<350]).sum(), S1_hi=np.exp(2*ly[g>=350]).sum(),
        Sm1_lo=X[g<350].sum(), Sm1_hi=X[g>=350].sum()))
json.dump(out,open('../perL.json','w'))
print(len(out))
