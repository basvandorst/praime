import json, numpy as np
from sympy import primerange
z=json.load(open('zeros131.json'))
P=np.array(list(primerange(2,10**6+1)),dtype=np.int64)
lx=np.log(P.astype(float))
rows=[]
for qs,res in z.items():
    q=int(qs); g=res['g']; n=q-1
    dl={pow(g,j,q):j for j in range(n)}
    F={}
    for k in range(1,n):
        gp=np.array(res['chars'][str(k)]); gm=-np.array(res['chars'][str(n-k)])
        gam=np.concatenate([gp,gm]); rho=0.5+1j*gam
        acc=np.zeros(len(P),complex)
        for c in range(0,len(gam),64):
            gg=gam[c:c+64]; rr=rho[c:c+64]
            acc+= (np.exp(1j*np.outer(lx,gg))/rr).sum(1)
        F[k]=-acc
    # measured
    r=P%q
    cnt1=np.cumsum(r==1)
    for b in range(2,q):
        jb=dl[b]
        E=np.zeros(len(P),complex); E0=np.zeros(len(P),complex)
        for k in range(1,n):
            chib=np.exp(2j*np.pi*k*jb/n)
            eps=1.0 if 2*k==n else 0.0
            E+=(np.conj(chib)-1)*(F[k]-eps)
            E0+=(np.conj(chib)-1)*F[k]
        E/=n; E0/=n
        cb=np.cumsum(r==b)
        meas=np.mean(cb>cnt1)
        rows.append(dict(q=q,b=b,meas=float(meas),pred=float(np.mean(E.real>0)),pred0=float(np.mean(E0.real>0)),maximag=float(np.abs(E.imag).max()),
                        legendre=1 if jb%2==0 else -1))
json.dump(rows,open('race131_mine.json','w'))
m=np.array([r['meas'] for r in rows]); p=np.array([r['pred'] for r in rows]); p0=np.array([r['pred0'] for r in rows])
print(len(rows),'maximag',max(r['maximag'] for r in rows))
def rep(a,b,lab):
    print(lab,'n',len(a),'pearson %.4f'%np.corrcoef(a,b)[0,1],'MAE %.4f'%np.mean(abs(a-b)))
rep(m,p,'all eps'); rep(m,p0,'all noeps')
mid=(m>=0.1)&(m<=0.9); rep(m[mid],p[mid],'mid eps'); rep(m[mid],p0[mid],'mid noeps')
from scipy.stats import spearmanr
print('spearman all',spearmanr(m,p)[0],'mid',spearmanr(m[mid],p[mid])[0])
# baseline: legendre symbol only predictor
leg=np.array([r['legendre'] for r in rows])
print('corr meas vs (b nonresidue)', np.corrcoef(m,(leg==-1).astype(float))[0,1])
# within-class correlations
for s in [1,-1]:
    sel=leg==s; rep(m[sel],p[sel],'leg=%d'%s)
print('frac meas<0.1 or >0.9', np.mean(~mid))
