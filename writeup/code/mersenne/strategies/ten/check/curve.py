import numpy as np, math
from lib import sieve
import cypari2; pari=cypari2.Pari()
N=100000; s=sieve(N); P=np.nonzero(s)[0]; P=P[P>2]
CUR={11:([0,-1,1,0,0],0),37:([0,0,1,-1,0],1),389:([0,1,1,-2,0],2),5077:([0,0,1,-7,6],3),32:([0,0,0,-1,0],0),27:([0,0,1,0,0],0)}
def ap_all(c):
    a1,a2,a3,a4,a6=c; out=[]
    for p in P:
        p=int(p); x=np.arange(p,dtype=np.int64); f=(((x+a2)%p*x+a4)%p*x+a6)%p; D=(a3*a3+4*f)%p
        t=np.zeros(p,dtype=np.int8)-1; t[(x*x)%p]=1; t[0]=0; out.append(-int(t[D].sum()))
    return np.array(out)
lp=np.log(P); sq=np.sqrt(P)
for N0,(c,r) in CUR.items():
    a=ap_all(c); E=pari.ellinit(c); bad=sum(1 for p,x in zip(P[:2000],a[:2000]) if int(pari.ellap(E,int(p)))!=x)
    y=a/(2*sq); hist=np.histogram(y,bins=8,range=(-1,1))[0]/len(y)
    e=np.linspace(-1,1,9); F=lambda t:0.5+(t*np.sqrt(1-t*t)+np.arcsin(t))/math.pi; st=np.diff(F(e))
    nag=np.cumsum(a*lp/P); idx=[np.searchsorted(P,X)-1 for X in (1000,10000,100000)]
    ys=np.sort(y); ks=np.abs(F(ys)-np.arange(1,len(ys)+1)/len(ys)).max()
    print("conductor %d rank %d: PARI mismatches %d; a_p=0: %d of %d; mean y %.4f; KS to semicircle %.3f (5%% level %.3f)"%(N0,r,bad,int((a==0).sum()),len(a),y.mean(),ks,1.36/math.sqrt(len(y))))
    print("   hist",np.round(hist,3),"\n   semi",np.round(st,3))
    print("   sum a_p ln p / p at 1e3,1e4,1e5: %s ; divided by ln x: %s ; fitted slope between 1e3 and 1e5: %.2f"%([round(float(nag[i]),2) for i in idx],[round(float(nag[i]/math.log(P[i])),2) for i in idx],(nag[idx[2]]-nag[idx[0]])/(math.log(1e5)-math.log(1e3))))
    w=np.cumsum(a/sq); print("   sum a_p/sqrt p at 1e5: %.1f (coin sd %.0f)"%(w[-1],math.sqrt(len(a))))
