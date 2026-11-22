import numpy as np, json
from math import log, pi, e, cos
Z={}
for f in ['part0','part1']:
    d=np.load(f'/home/bas/praime/writeup/data/riemann/lzeros/r117_{f}.npz')
    for k in d.keys(): Z[int(k[1:])]=d[k]
print(len(Z))
def unfold(t,q): return t/(2*pi)*np.log(q*t/(2*pi*e))
rows=[]
for D,g in sorted(Z.items()):
    q=abs(D)
    g=np.sort(g); g=g[(g>=100)&(g<=600)]
    x=unfold(g,q)
    lo,hi=unfold(100.0,q),unfold(600.0,q)
    # continuous window positions
    s=np.arange(lo,hi-1,0.002)
    c=np.searchsorted(x,s+1)-np.searchsorted(x,s)
    v=c.var(); v1=((c-1.0)**2).mean()
    sp=np.diff(x); sv=sp.var()
    # mean L over zeros
    Lm=np.mean(np.log(q*g/(2*pi)))
    rows.append(dict(D=D,q=q,n=len(g),S2=v,S2m1=v1,spv=sv,Lm=Lm,c2=int(D%2==0),c3=int(D%3==0),c5=int(D%5==0),neg=int(D<0)))
json.dump(rows,open('/home/bas/scratch/review-logic2/rows.json','w'))
import numpy.linalg as la
R=rows
y=np.array([r['S2'] for r in R])
def fit(cols,yy=y,names=None):
    X=np.column_stack([np.ones(len(R))]+cols)
    b,res,rk,sv=la.lstsq(X,yy,rcond=None)
    r=yy-X@b; s2=r@r/(len(yy)-X.shape[1]); cov=s2*la.inv(X.T@X)
    return b,np.sqrt(np.diag(cov))
c2=np.array([r['c2'] for r in R]);c3=np.array([r['c3'] for r in R]);c5=np.array([r['c5'] for r in R])
lq=np.log([r['q'] for r in R]); Lm=np.array([r['Lm'] for r in R]); neg=np.array([r['neg'] for r in R])
print('n2,n3,n5',c2.sum(),c3.sum(),c5.sum(),'neg',neg.sum())
print('mean S2',y.mean())
for nm,cols in [('base',[c2,c3,c5]),('+lnq',[c2,c3,c5,lq]),('+Lm',[c2,c3,c5,Lm]),('+neg',[c2,c3,c5,neg]),('+lnq+neg',[c2,c3,c5,lq,neg])]:
    b,se=fit(cols); print(nm,np.round(b,4),np.round(se,4))
ym1=np.array([r['S2m1'] for r in R]); b,se=fit([c2,c3,c5],ym1); print('about 1',np.round(b,4),np.round(se,4))
ysp=np.array([r['spv'] for r in R]); b,se=fit([c2,c3,c5],ysp); print('spacing var',np.round(b,4),np.round(se,4))
print('corr c2 vs lnq',np.corrcoef(c2,lq)[0,1],'c2 vs Lm',np.corrcoef(c2,Lm)[0,1],'c3',np.corrcoef(c3,lq)[0,1],'c5',np.corrcoef(c5,lq)[0,1])
print('mean lnq c2=1/0',lq[c2==1].mean(),lq[c2==0].mean(), 'Lm', Lm[c2==1].mean(),Lm[c2==0].mean())
print('mean L overall',Lm.mean())
def berry(p,L,K=60):
    s=0
    for k in range(1,K):
        a=2*pi*k*log(p)/L; s+=log(p)**2*p**(-k)*(1-cos(a))/a**2
    return -4/L**2*s
for p in [2,3,5]:
    print(p,'at mean L',round(berry(p,Lm.mean()),5),'avg over L-fns',round(np.mean([berry(p,l) for l in Lm]),5))
# average over t within an L-function too, weighting by density
def berry_t(p,q):
    t=np.linspace(100,600,2001); w=np.log(q*t/(2*pi)); L=np.log(q*t/(2*pi))
    return np.sum(w*np.array([berry(p,l) for l in L]))/w.sum()
for p in [2,3,5]:
    print(p,'avg over t and fns',round(np.mean([berry_t(p,r['q']) for r in R]),5))
