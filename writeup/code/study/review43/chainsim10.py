import numpy as np, math, subprocess, sys
d=np.load('an_6e11.npz'); P1=d['Ph_10_1']
print('P1 rows sum', P1.sum(1).round(6), 'cols', P1.sum(0).round(4))
nb=100; Nblk=36000; nwin=40
out={k:{'chi':[],'col':[],'win':[]} for k in [1,2]}
for w in range(nwin):
    n=nb*Nblk+10
    raw=subprocess.run(['./chain10',str(n),str(w+1)]+['%.10f'%v for v in P1.ravel()],capture_output=True).stdout
    s=np.frombuffer(raw,np.uint8).astype(np.int64)
    blk=np.arange(nb*Nblk)//Nblk
    for k in [1,2]:
        a=s[:nb*Nblk]; b=s[k:nb*Nblk+k]
        X=np.bincount(blk*16+a*4+b,minlength=nb*16).reshape(nb,4,4).astype(float)
        Na=X.sum(2); Ph=X.sum(0)/Na.sum(0)[:,None]
        E=Na[:,:,None]*Ph[None]; D=X-E
        chi=(D**2/E).sum((1,2)); c=D.sum(1)
        M=np.einsum('bi,ij->bj',Na,Ph)[:,:,None]*np.eye(4)[None]-np.einsum('bi,ij,ik->bjk',Na,Ph,Ph)
        col=np.einsum('bj,bjk,bk->b',c,np.linalg.pinv(M,hermitian=True),c)
        out[k]['chi'].append(chi.sum()/(99*12)); out[k]['col'].append(col.sum()/(99*3)); out[k]['win'].append((chi-col).sum()/(99*9))
for k in [1,2]:
    print('chain mod 10, table k=%d: '%k+'  '.join('%s %.4f+-%.4f'%(p,np.mean(out[k][p]),np.std(out[k][p],ddof=1)/math.sqrt(nwin)) for p in ['chi','col','win']))
