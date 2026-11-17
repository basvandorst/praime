import numpy as np, math, sys
from lib import sieve
N=20_000_000; s=sieve(N); pr=np.nonzero(s[:4500])[0]; y0=871
sv=np.ones(N+1,dtype=bool); sv[:2]=False
for q in pr[pr<=y0]: sv[2*q::q]=False
spf=np.zeros(N+1,dtype=np.int32)
for q in pr[pr>y0][::-1]: spf[int(q)*int(q)::int(q)]=q
comp=sv&~s; comp[:y0*y0]&=False
nn=np.nonzero(comp)[0]; q=spf[nn].astype(np.int64); c=nn//q; r=np.floor(np.sqrt(nn)).astype(np.int64); r-=(r*r>nn); j=(q+c)//2-r
jf=np.full(N+1,10**9,dtype=np.int64); jf[nn]=j
print("survivors of the sieve to %d up to %d: %d; composite ones %d; with t - r <= 1, 2, 5, 10, 30, 100: %s"%(y0,N,int(sv.sum()),len(nn),[int((j<=J).sum()) for J in (1,2,5,10,30,100)]))
n=np.arange(10**6,N-40000,dtype=np.int64); rr=np.floor(np.sqrt(n)).astype(np.int64); rr-=(rr*rr>n); rr+=((rr+1)**2<=n); ks=list(range(-4,5))
for J in (0,1,2,5,10,30,100,10**8):
    L=s|(comp&(jf>J)) if J<10**8 else s.copy(); a=int(L[n].sum()); vals=[]
    for k in ks:
        m=n+2*rr+2*k; pp=int((L[n]&L[m]).sum()); vals.append(pp*len(n)/(a*int(L[m].sum())))
    mu=sum(vals)/len(vals); dev=[100*(v/mu-1) for v in vals]
    print("composites kept: t - r > %-9s  set size %8d  deviations %% (k=-4..4): %s   rms %.2f"%(J if J<10**8 else "none (primes)",a," ".join("%+.2f"%d for d in dev),math.sqrt(sum(d*d for d in dev)/len(dev))))
