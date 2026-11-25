# k = -1: heavy-tail analysis of sum 1/|L'|^2 per L-function and of the per-prime factor
import numpy as np, json
P=[2,3,5,7,11,13]
d=[r for r in json.load(open('../perL.json')) if r['q']!=1153]
z=np.load('../all.npz')
n=len(d)
def design(rows):
    cols=[np.ones(len(rows))]
    for s in ['Fn','R120']: c=np.array([r['set']==s for r in rows],float); cols.append(c-c.mean())
    c=np.array([r['lnL'] for r in rows]); cols.append(c-c.mean())
    for p in P: c=np.array([p in r['bad'] for r in rows],float); cols.append(c-c.mean())
    return np.array(cols).T
X=design(d); ip=4
Sm1=np.array([r['Sm1'] for r in d]); Sm1t=np.array([r['Sm1t'] for r in d]); Pm=np.array([r['Pm1'] for r in d])
Pnone=np.full(n,6/np.pi**2*500/(2*np.pi))
Pk1=np.array([6/np.pi**2*np.prod([(1-1/p) for p in r['bad']]) for r in d])*500/(2*np.pi)
print('pooled untruncated ratio sum Sm1/sum P_conj: %.4f ; truncated (25) %.4f'%(Sm1.sum()/Pm.sum(), Sm1t.sum()/Pm.sum()))
lr=np.log(Sm1/Pm)
o=np.argsort(lr)[::-1]
print('largest log(Sm1/P):',[(d[i]['D'],round(lr[i],3),round(d[i]['Xmax'],1)) for i in o[:6]])
print('smallest:',[(d[i]['D'],round(lr[i],3)) for i in o[-4:]])
print('per-L-function log ratio sd %.4f, skew %.2f'%(lr.std(), ((lr-lr.mean())**3).mean()/lr.std()**3))
# tail of the per-zero variable X = 1/|L'|^2 scaled by its mean (pooled)
allX=[];allXs=[]
for r in d:
    y=z['%s_%d_y'%(r['set'],r['D'])]; x=np.exp(-2*y); allX.append(x); allXs.append(x/x.mean())
allX=np.concatenate(allX)
xs=np.sort(allX)[::-1]
# Hill estimator of tail index of X
for k in [30,100,300,1000,3000]:
    h=1/np.mean(np.log(xs[:k]/xs[k])); print('Hill alpha(k=%d) = %.2f  (X_k=%.2f, |L\'|=%.3f)'%(k,h,xs[k],xs[k]**-0.5))
def ols(Xm,y): return np.linalg.lstsq(Xm,y,rcond=None)[0]
def coef(idx,stat='u'):
    S = Sm1 if stat=='u' else Sm1t
    b=ols(X[idx],np.log(S[idx]/Pm[idx]))[ip:]
    lf=np.log(Pm/Pnone); X2=np.column_stack([X[idx,:ip],lf[idx]-lf[idx].mean()])
    g=ols(X2,np.log(S[idx]/Pnone[idx]))[-1]
    return np.concatenate([b,[g]])
full=coef(np.arange(n))
print('full: beta_p (vs p/(p+1)) and amplitude g:',np.round(full,4))
rng=np.random.default_rng(5)
bs=np.array([coef(rng.integers(0,n,n)) for _ in range(2000)])
print('bootstrap se            :',np.round(bs.std(0,ddof=1),4))
# subsampling (Politis-Romano), without replacement, m = n/2, n/4, n/8; spread scaled to n by sqrt(m/n) [Gaussian rate]
# and by (m/n)^(1-1/alpha) with alpha = 1.5 [stable rate]
for m in [n//2,n//4,n//8]:
    ss=np.array([coef(rng.choice(n,m,replace=False)) for _ in range(2000)])
    sd=ss.std(0,ddof=1)
    # without-replacement finite-population correction for the Gaussian-rate scaling: sqrt(m/n * (1-m/n)^-1)...
    print('subsample m=%3d: se scaled gaussian'%m, np.round(sd*np.sqrt(m/n)/np.sqrt(1-m/n),4), ' stable(1.5)', np.round(sd*(m/n)**(1/3)/np.sqrt(1-m/n),4))
# remove the largest single terms: drop the top-k zeros overall
thr=np.sort(allX)[::-1]
for k in [2,8,20,50]:
    cut=thr[k-1]
    S=np.array([np.where(np.exp(-2*z['%s_%d_y'%(r['set'],r['D'])])>=cut,0,np.exp(-2*z['%s_%d_y'%(r['set'],r['D'])])).sum() for r in d])
    b=ols(X,np.log(S/Pm))[ip:]
    lf=np.log(Pm/Pnone); g=ols(np.column_stack([X[:,:ip],lf-lf.mean()]),np.log(S/Pnone))[-1]
    print('top %2d zeros removed (X>=%.1f): beta'%(k,cut),np.round(b,4),'g %.4f'%g, 'pooled %.4f'%(S.sum()/Pm.sum()))
# drop whole L-functions with largest log ratio
for k in [1,3,10]:
    idx=np.setdiff1d(np.arange(n),o[:k]); print('drop %d largest L-functions:'%k,np.round(coef(idx),4))
# alternatives per prime: z of data vs the three factor laws (using bootstrap se)
for j,p in enumerate(P[:4]):
    b=full[j]; s=bs.std(0,ddof=1)[j]
    alt_none=-np.log(p/(p+1)); alt_k1=np.log((1-1/p)*(p+1)/p)
    print('p=%d: beta %.4f+-%.4f; z vs p/(p+1) %.1f; vs factor 1 (pred %+.3f) %.1f; vs 1-1/p (pred %+.3f) %.1f'%(p,b,s,b/s,alt_none,(b-alt_none)/s,alt_k1,(b-alt_k1)/s))
# robust: Huber-type via median of log ratio per group is not regression; use LAD by IRLS
def lad(Xm,y,it=200):
    b=ols(Xm,y)
    for _ in range(it):
        w=1/np.maximum(np.abs(y-Xm@b),1e-4); b=np.linalg.lstsq(Xm*np.sqrt(w)[:,None],y*np.sqrt(w),rcond=None)[0]
    return b
print('LAD beta_p:',np.round(lad(X,np.log(Sm1/Pm))[ip:],4))
