# Own segmented odd-only sieve: primes in [x0, x0+L+extra). Saves int64 array.
import numpy as np, sys, math, time
x0=int(float(sys.argv[1])); L=int(float(sys.argv[2])); extra=int(float(sys.argv[3])) if len(sys.argv)>3 else 10**5
out=sys.argv[4] if len(sys.argv)>4 else 'pr_%s.npy'%sys.argv[1]
hi=x0+L+extra
r=math.isqrt(hi)+1
s=np.ones(r+1,bool); s[:2]=False
for i in range(2,math.isqrt(r)+1):
    if s[i]: s[i*i::i]=False
base=np.nonzero(s)[0][1:].astype(np.int64)  # odd primes
SEG=10**7
res=[]
t=time.time()
lo=x0 if x0%2 else x0+1
while lo<hi:
    hi_s=min(lo+SEG, hi)
    n=(hi_s-lo+1)//2
    seg=np.ones(n,bool)
    m=-(-lo//base)*base
    m=np.maximum(m, base*base)
    m=np.where(m%2==0, m+base, m)
    st=(m-lo)//2
    for p,i0 in zip(base.tolist(), st.tolist()):
        if i0<n: seg[i0::p]=False
    idx=np.nonzero(seg)[0]
    pr=lo+2*idx.astype(np.int64)
    pr=pr[pr<hi_s]
    res.append(pr)
    lo+=2*n
pr=np.concatenate(res)
if x0<=2<hi: pr=np.concatenate([[2],pr])
pr=pr[pr>=x0]
np.save(out,pr)
print(x0,L,len(pr),'time %.1f'%(time.time()-t), 'in window', np.sum(pr<x0+L))
