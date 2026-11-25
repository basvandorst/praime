import numpy as np, json, sys
P=[2,3,5,7,11,13]
d=json.load(open('../perL.json'))
excl = (sys.argv[1] if len(sys.argv)>1 else 'excl')
sel = (sys.argv[2] if len(sys.argv)>2 else 'all')
def keep(r):
    if excl=='excl' and r['q']==1153: return False
    if sel=='pos': return r['D']>0
    if sel=='neg': return r['D']<0
    if sel in ('R120','Dn','Fn'): return r['set']==sel
    if sel=='odd': return r['q']%2==1
    return True
d=[r for r in d if keep(r)]
n=len(d)
def design(rows, sets=True):
    cols=[np.ones(len(rows))]
    S=sorted({r['set'] for r in rows})
    if sets:
        for s in S[1:]: c=np.array([r['set']==s for r in rows],float); cols.append(c-c.mean())
    c=np.array([r['lnL'] for r in rows]); cols.append(c-c.mean())
    for p in P: c=np.array([p in r['bad'] for r in rows],float); cols.append(c-c.mean())
    return np.array(cols).T
X=design(d); ip=X.shape[1]-len(P)
def ols(X,y): return np.linalg.lstsq(X,y,rcond=None)[0]
mu=np.array([r['mu'] for r in d]); v=np.array([r['v'] for r in d])
invL=np.array([r['invL'] for r in d])
pm=np.array([sum(np.log(p)*np.log(p/(p-1)) for p in r['bad']) for r in d])*invL   # first-order mean shift, all bad primes
S1=np.array([r['S1'] for r in d]); S1p=np.array([r['S1pred'] for r in d])
Sm1=np.array([r['Sm1'] for r in d]); Sm1t=np.array([r['Sm1t'] for r in d]); Pm=np.array([r['Pm1'] for r in d])
Pm_none=np.full(n,6/np.pi**2*500/(2*np.pi))
Pm_k1=np.array([6/np.pi**2*np.prod([(1-1/p)/(1-1/p**2) for p in r['bad']]) for r in d])*500/(2*np.pi)  # factor 1-1/p instead of p/(p+1)
def readings(idx):
    Xi=X[idx]
    out={}
    b=ols(Xi,mu[idx]); bp=ols(Xi,pm[idx]); out['mu']=b[ip:]; out['mu_pred']=bp[ip:]
    out['gM']=np.sum(b[ip:ip+3]*bp[ip:ip+3])/np.sum(bp[ip:ip+3]**2)
    b=ols(Xi,v[idx]); out['v']=b[ip:]
    out['S1']=ols(Xi,np.log(S1[idx]/S1p[idx]))
    out['Sm1']=ols(Xi,np.log(Sm1[idx]/Pm[idx])); out['Sm1t']=ols(Xi,np.log(Sm1t[idx]/Pm[idx]))
    out['Sm1_none']=ols(Xi,np.log(Sm1[idx]/Pm_none[idx])); 
    out['Sm1_k1']=ols(Xi,np.log(Sm1[idx]/Pm_k1[idx]))
    # amplitude of p/(p+1) factor: log Sm1/Pnone = a + ... + g * log(prod p/(p+1)) 
    lf=np.log(Pm[idx]/Pm_none[idx]); X2=np.column_stack([Xi[:,:ip],lf-lf.mean()])
    out['gPC_direct']=ols(X2,np.log(Sm1[idx]/Pm_none[idx]))[-1]
    X2=np.column_stack([Xi[:,:ip],pm[idx]-pm[idx].mean()]); out['gM_direct']=ols(X2,mu[idx])[-1]
    out['S1level']=np.log(S1[idx].sum()/S1p[idx].sum()); out['Sm1level']=np.log(Sm1[idx].sum()/Pm[idx].sum())
    return out
full=readings(np.arange(n))
rng=np.random.default_rng(99); B=2000
boots=[readings(rng.integers(0,n,n)) for _ in range(B)]
def se(key): 
    a=np.array([np.atleast_1d(b[key]) for b in boots]); return a.std(axis=0,ddof=1)
np.set_printoptions(precision=4,suppress=True,linewidth=150)
print('n L-functions',n, excl, sel)
for key in ['mu','mu_pred','v','S1','Sm1','Sm1t','Sm1_none','Sm1_k1']:
    f=np.atleast_1d(full[key]); s=se(key)
    if key in ('S1','Sm1','Sm1t','Sm1_none','Sm1_k1'): f=f[ip:]; s=s[ip:]
    print('%-9s'%key,'fit',f,' se',s)
for key in ['gM','gM_direct','gPC_direct','S1level','Sm1level']:
    print('%-10s %.4f +- %.4f'%(key,full[key],se(key)[0]))
