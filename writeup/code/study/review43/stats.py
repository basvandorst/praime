# Own implementation of S1 and S2 (with column/within split). No repo code.
import numpy as np, math
def red(q): return [a for a in range(q) if math.gcd(a,q)==1]
def S1_blocks(pr, x0, L, H, q):
    """per-block X_beta/(phi-1) for equal-class chi2"""
    p=pr[pr<x0+L]; b=(p-x0)//H; nb=L//H
    R=red(q); phi=len(R); lut=-np.ones(q,np.int64); lut[R]=np.arange(phi)
    cl=lut[p%q]; m=cl>=0
    cnt=np.bincount(b[m]*phi+cl[m], minlength=nb*phi).reshape(nb,phi).astype(float)
    N=cnt.sum(1); E=N/phi
    X=((cnt-E[:,None])**2/E[:,None]).sum(1)
    return X/(phi-1), cnt
def S2_blocks(pr, x0, L, H, q, k):
    """per-block chi2, col, within (unnormalised) for transition tables; Phat per window"""
    nwin=np.searchsorted(pr, x0+L)
    a=pr[:nwin]%q; bb=pr[k:nwin+k]%q
    assert len(bb)==nwin
    blk=(pr[:nwin]-x0)//H; nb=L//H
    R=red(q); phi=len(R); lut=-np.ones(q,np.int64); lut[R]=np.arange(phi)
    ia=lut[a]; ib=lut[bb]; m=(ia>=0)&(ib>=0)
    X=np.bincount(blk[m]*phi*phi+ia[m]*phi+ib[m], minlength=nb*phi*phi).reshape(nb,phi,phi).astype(float)
    Na=X.sum(2)              # nb x R
    Ph=X.sum(0)/Na.sum(0)[:,None]   # R x C
    E=Na[:,:,None]*Ph[None]
    chi=((X-E)**2/E).sum((1,2))
    col=np.zeros(nb); win2=np.zeros(nb)
    for be in range(nb):
        d=X[be]-E[be]; c=d.sum(0)
        M=sum(Na[be,i]*(np.diag(Ph[i])-np.outer(Ph[i],Ph[i])) for i in range(phi))
        Mp=np.linalg.pinv(M)
        col[be]=c@Mp@c
        # independent route for within: residual after removing column projection
        w=0.0
        for i in range(phi):
            Si=Na[be,i]*(np.diag(Ph[i])-np.outer(Ph[i],Ph[i]))
            r=d[i]-Si@Mp@c
            w+=np.sum(r**2/(Na[be,i]*Ph[i]))
        win2[be]=w
    return chi, col, chi-col, win2, Ph, X
