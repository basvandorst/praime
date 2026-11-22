# rough cost of one law point at d=2400 with the same array shapes as r59.py's lnS_parts/share (own code, timing only)
import numpy as np, time, sympy, math
d=2400; h=np.arange(6,d,6); n=len(h); MID=list(sympy.primerange(17,3001))
t=time.time(); lnr=np.zeros((n,n))
for q in MID:
    r1=h%q; r2=(h+2)%q; o1=np.ones(n,bool); o2=o1
    delta=((r1[None,:]==r1[:,None])&o1[None,:]).astype(int)+((r1[None,:]==r2[:,None])&o1[None,:])+((r2[None,:]==r1[:,None])&o2[None,:])+((r2[None,:]==r2[:,None])&o2[None,:])
    nhk=4+delta
    lnr+=np.log1p(-nhk/q)
t1=time.time()-t
M=(np.random.rand(1485,n)<0.5).astype(float); R=np.random.rand(n,n)
t=time.time()
for _ in range(3): L=0.5*np.sum((M@R)*M,axis=1)
print('n',n,'lnS_parts-like',round(t1,2),'s; share-like',round(time.time()-t,2),'s')
