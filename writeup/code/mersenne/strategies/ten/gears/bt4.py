import numpy as np, math
T=100_000_000
s=np.ones(T+1,dtype=bool); s[:2]=False
for i in range(2,int(T**0.5)+1):
    if s[i]: s[i*i::i]=False
pr=np.nonzero(s)[0]; del s
g=np.diff(pr).astype(np.int64); sg=np.where(np.arange(len(g))%2==0,1,-1); alt=np.cumsum(g*sg)
print("alternating sum of the prime gaps (gap 1 added, gap 2 subtracted, ...), checkpoints:")
for X in (10**5,10**6,2*10**6,10**7,3*10**7,10**8):
    n=int(np.searchsorted(pr,X))-1; gg=g[:n].astype(float); sd=math.sqrt(n*gg.var())
    print("   up to %9d: value %7d, largest |value| so far %6d, spread expected for gaps in random order %6.0f, value/spread %+.2f"%(X,alt[n-1],np.abs(alt[:n]).max(),sd,alt[n-1]/sd))
ev=g[0::2]; od=g[1::2]; m=min(len(ev),len(od))
print("mean gap at odd places %.4f, at even places %.4f, difference %.4f +- %.4f"%(ev[:m].mean(),od[:m].mean(),ev[:m].mean()-od[:m].mean(),math.sqrt((ev[:m].var()+od[:m].var())/m)))
z=np.nonzero(alt==0)[0]; print("times the sum is exactly 0 up to 1e8:",len(z),"; last at prime",int(pr[z[-1]+1]) if len(z) else None)
