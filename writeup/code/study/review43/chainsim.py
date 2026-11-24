# Null check: first-order 2-state chain (classes mod 3) with the primes' flip rate; do col and within come out at 1 with the 99 normalisation?
import numpy as np, math
rng=np.random.default_rng(7)
alpha=0.5430; nb=100; Nblk=37000; nwin=60
res={1:{'col':[],'win':[],'chi':[]},2:{'col':[],'win':[],'chi':[]}}
for w in range(nwin):
    n=nb*Nblk+10
    flips=rng.random(n)<alpha
    s=np.cumsum(flips)%2
    blk=np.arange(n-4)//Nblk; blk=blk[:nb*Nblk]
    for k in [1,2]:
        a=s[:nb*Nblk]; b=s[k:nb*Nblk+k]
        X=np.bincount(blk*4+a*2+b,minlength=nb*4).reshape(nb,2,2).astype(float)
        Na=X.sum(2); Ph=X.sum(0)/Na.sum(0)[:,None]
        E=Na[:,:,None]*Ph[None]; D=X-E
        chi=(D**2/E).sum((1,2)); c=D.sum(1)
        M=np.einsum('bi,ij->bj',Na,Ph)[:,:,None]*np.eye(2)[None]-np.einsum('bi,ij,ik->bjk',Na,Ph,Ph)
        col=np.einsum('bj,bjk,bk->b',c,np.linalg.pinv(M,hermitian=True),c)
        res[k]['chi'].append(chi.sum()/(99*2)); res[k]['col'].append(col.sum()/99); res[k]['win'].append((chi-col).sum()/99)
for k in [1,2]:
    for p in ['chi','col','win']:
        v=np.array(res[k][p]); print('chain k=%d %s %.4f +- %.4f'%(k,p,v.mean(),v.std(ddof=1)/math.sqrt(len(v))))
