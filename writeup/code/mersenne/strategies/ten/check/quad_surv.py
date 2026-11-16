import numpy as np, math, sys
from lib import sieve
N=20_000_000; s=sieve(N); pr=np.nonzero(s[:4500])[0]
n=np.arange(10**6,N-40000,dtype=np.int64); r=np.floor(np.sqrt(n)).astype(np.int64); r-=(r*r>n); r+=((r+1)**2<=n)
ks=list(range(-4,5)); rows=[]
for y in (13,100,871,2000,4472):
    sv=np.ones(N+1,dtype=bool); sv[:2]=False
    for q in pr[pr<=y]: sv[2*q::q]=False
    a=int(sv[n].sum()); vals=[]
    for k in ks:
        m=n+2*r+2*k; pp=int((sv[n]&sv[m]).sum()); vals.append(pp*len(n)/(a*int(sv[m].sum())))
    mu=sum(vals)/len(vals); print("sieve to %4d: pair ratio at h = 2r+2k, k=-4..4, as %% deviation from their mean (%.4f): %s"%(y,mu," ".join("%+.2f"%(100*(v/mu-1)) for v in vals)))
print("rows used: r from %d to %d"%(r.min(),r.max()))
# prediction from local factors, per row
for y in (13,100,871):
    qs=[int(q) for q in pr[(pr>2)&(pr<=y)]]; out=[]
    rr=np.arange(r.min(),r.max()+1); w=np.bincount(r-r.min(),minlength=len(rr)).astype(float)
    for k in ks:
        g=np.ones(len(rr))
        for q in qs: g*=np.where((rr+k)%q==0,(1-1/q),(1-2/q))/(1-1/q)**2
        out.append(float((w*g).sum()/w.sum()))
    mu=sum(out)/len(out); print("local-factor prediction, sieve to %4d: %s"%(y," ".join("%+.2f"%(100*(v/mu-1)) for v in out)))
