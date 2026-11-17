import numpy as np, math
K=22; N=1<<K
def sieve(n):
    s=np.ones(n,dtype=bool); s[:2]=False
    for i in range(2,int(n**0.5)+1):
        if s[i]: s[i*i::i]=False
    return s
pri=sieve(N)
irr=np.ones(N,dtype=bool); irr[:2]=False
def clm(f,g):
    r=np.zeros(len(g),dtype=np.int64); i=0
    while f>>i:
        if (f>>i)&1: r^=g<<i
        i+=1
    return r
lowirr=[]
for f in range(2,1<<(K//2+1)):
    if irr[f]:
        df=f.bit_length()-1; g=np.arange(2,1<<(K-df),dtype=np.int64); irr[clm(f,g)]=False
        if df<=4: lowirr.append(f)
n=np.arange(N,dtype=np.int64); pc=np.zeros(N,dtype=np.int64)
for j in range(K): pc+=(n>>j)&1
base=(n&1==1)&(pc%2==1)
print("numbers below 2^%d: primes %d, irreducible patterns %d (no carries), both %d"%(K,pri.sum(),irr.sum(),(pri&irr).sum()))
print("by length: among the odd numbers with an odd number of ones: primes, patterns, both, expected if independent, ratio, z")
for k in range(12,K+1):
    m=base&(n>=1<<(k-1))&(n<1<<k); a=(pri&m).sum(); b=(irr&m).sum(); c=(pri&irr&m).sum(); e=a*b/m.sum()
    print("  %2d bits: %6d %6d %6d  expected %8.1f  ratio %.3f  z %+.1f"%(k,a,b,c,e,c/e,(c-e)/math.sqrt(e)))
m=base&(n>=1<<15)
A=m.copy()
for p in (3,5,7,11,13): A&=(n%p!=0)
B=m.copy()
def pmodarr(f):
    r=n.copy(); df=f.bit_length()-1
    for b in range(K-1,df-1,-1):
        sel=((r>>b)&1)==1; r[sel]^=f<<(b-df)
    return r
for f in lowirr:
    B&=(pmodarr(f)!=0)
M=m.sum(); a=A.sum(); b=B.sum(); c=(A&B).sum()
print("layer of small factors (numbers of 16 to %d bits, odd, odd number of ones: %d): no prime factor up to 13: %d; no pattern factor of up to 5 bits (%s): %d; both: %d; expected if independent %.0f; ratio %.4f"%(K,M,a,",".join(map(str,lowirr)),b,c,a*b/M,c/(a*b/M)))
pa=(pri&m).sum(); ib=(irr&m).sum(); both=(pri&irr&m).sum()
print("primes %d (all inside the first set), patterns %d (all inside the second), both %d; expected from the overlap of the two sets alone: %.0f (ratio of found to that %.3f); expected with no link at all: %.0f"%(pa,ib,both,pa*ib*c/(a*b),both/(pa*ib*c/(a*b)),pa*ib/M))
for p,f in ((3,7),(3,3),(5,7),(7,11),(7,13),(5,31),(7,7)):
    x=(n%p==0)&(n>=1<<15); y=(pmodarr(f)==0)&(n>=1<<15); T=(n>=1<<15).sum()
    print("  %d divides n and the pattern %d divides n without carries: %d, expected if independent %.0f, ratio %.4f"%(p,f,(x&y).sum(),x.sum()*y.sum()/T,(x&y).sum()/(x.sum()*y.sum()/T)))
