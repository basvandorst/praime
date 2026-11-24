# Per-window per-block statistics for many (q,H) and (q,k,H). Own code.
import sys,os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, math
from stats import S1_blocks, red
x0=int(float(sys.argv[1])); L=int(float(sys.argv[2])); f=sys.argv[3]; out=sys.argv[4]
pr=np.load(f)
res={}
for q in [3,4,5,7,8,9,10,12,30,210]:
    for H in [10**4,10**5,10**6,10**7]:
        if L//H<2: continue
        v,_=S1_blocks(pr,x0,L,H,q); res['S1_%d_%d'%(q,H)]=v
def S2fast(pr,x0,L,H,q,k):
    nwin=np.searchsorted(pr,x0+L)
    a=pr[:nwin]%q; bb=pr[k:nwin+k]%q
    blk=(pr[:nwin]-x0)//H; nb=L//H
    R=red(q); phi=len(R); lut=-np.ones(q,np.int64); lut[R]=np.arange(phi)
    ia=lut[a]; ib=lut[bb]; m=(ia>=0)&(ib>=0)
    X=np.bincount(blk[m]*phi*phi+ia[m]*phi+ib[m],minlength=nb*phi*phi).reshape(nb,phi,phi).astype(float)
    Na=X.sum(2); Ph=X.sum(0)/Na.sum(0)[:,None]
    E=Na[:,:,None]*Ph[None]; D=X-E
    chi=(D**2/E).sum((1,2))
    c=D.sum(1)                       # nb x C
    # M per block: sum_i Na_i (diag P_i - P_i P_i')
    Mb=np.einsum('bi,ij->bj',Na,Ph)[:,:,None]*np.eye(phi)[None] - np.einsum('bi,ij,ik->bjk',Na,Ph,Ph)
    Mp=np.linalg.pinv(Mb, hermitian=True)
    col=np.einsum('bj,bjk,bk->b',c,Mp,c)
    return chi/(phi*(phi-1)), col/(phi-1), (chi-col)/((phi-1)**2), Ph
for q in [3,4,5,7,10]:
    for k in [1,2,3,4]:
        for H in [10**4,10**5,10**6]:
            if q>5 and H==10**4: continue
            chi,col,win,Ph=S2fast(pr,x0,L,H,q,k)
            res['S2chi_%d_%d_%d'%(q,k,H)]=chi; res['S2col_%d_%d_%d'%(q,k,H)]=col; res['S2win_%d_%d_%d'%(q,k,H)]=win
            if H==10**6: res['Ph_%d_%d'%(q,k)]=Ph
np.savez(out,**res)
print('ok',x0,len(res))
