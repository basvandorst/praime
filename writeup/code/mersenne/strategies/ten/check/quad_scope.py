import numpy as np, math, sys
from lib import sieve
N=20_000_000; s=sieve(N); y0=871; sv=np.ones(N+1,dtype=bool); sv[:2]=False
for q in np.nonzero(s[:y0+1])[0]: sv[2*q::q]=False
n=np.arange(10**6,N-80000,dtype=np.int64)
def fl(x):
    r=np.floor(x).astype(np.int64); return r
r=fl(np.sqrt(n)); r-=(r*r>n); r+=((r+1)**2<=n)
rules={"2*floor(sqrt n)":2*r,"4*floor(sqrt n)":4*r,"6*floor(sqrt n)":6*r,"2*floor(sqrt(2n))":2*fl(np.sqrt(2*n)),"2*floor(sqrt(n/2))":2*fl(np.sqrt(n/2)),"2*floor(sqrt(3n))":2*fl(np.sqrt(3*n)),"2*floor(1.5 sqrt n)":2*fl(1.5*np.sqrt(n)),"2*round(sqrt n)":2*fl(np.sqrt(n)+0.5),"2*floor(n^(1/3))":2*fl(np.cbrt(n)),"2*floor(n^(2/3)/30)":2*fl(np.cbrt(n)**2/30),"2*floor(sqrt n * ln n / 10)":2*fl(np.sqrt(n)*np.log(n)/10),"fixed 6000":np.int64(6000)+0*r}
a_s=int(sv[n].sum()); a_p=int(s[n].sum()); noise=None
print("deviation of the pair ratio from its mean over k = -3..3, in %: survivors of the sieve to 871 | primes   (noise of the primes about 0.26)")
for name,h in rules.items():
    vs=[];vp=[]
    for k in range(-3,4):
        m=n+h+2*k; vs.append(int((sv[n]&sv[m]).sum())*len(n)/(a_s*int(sv[m].sum()))); vp.append(int((s[n]&s[m]).sum())*len(n)/(a_p*int(s[m].sum())))
    ms=sum(vs)/7; mp=sum(vp)/7; ds=[100*(v/ms-1) for v in vs]; dp=[100*(v/mp-1) for v in vp]
    print("  h = %-28s rms %.2f | %.2f    survivors: %s"%(name+" + 2k",math.sqrt(sum(d*d for d in ds)/7),math.sqrt(sum(d*d for d in dp)/7)," ".join("%+.2f"%d for d in ds)),flush=True)
