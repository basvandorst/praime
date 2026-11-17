import numpy as np, math, sys
K=int(sys.argv[1]) if len(sys.argv)>1 else 21
N=1<<K
def sieve(n):
    s=np.ones(n+1,dtype=bool); s[:2]=False
    for i in range(2,int(n**0.5)+1):
        if s[i]: s[i*i::i]=False
    return s
s=sieve(N)[:N]
n=np.arange(N,dtype=np.int64); pc=np.zeros(N,dtype=np.int64)
for j in range(K): pc+=(n>>j)&1
top=(n>=N//2)&(n&1==1); P=top&s
print("numbers with %d bits (odd): %d, primes %d (share %.4f)"%(K,top.sum(),P.sum(),P.sum()/top.sum()))
print("LAYER 1: one bit layer at a time: share of primes with a 1 in the layer")
print("   "+" ".join("b%d %.4f"%(j,((n[P]>>j)&1).mean()) for j in range(1,K-1)))
print("LAYER 2: pairs of bit layers: largest information between two layers (bits)")
best=(0,0,0)
pp=n[P]
for i in range(1,K-1):
    for j in range(i+1,K-1):
        a=(pp>>i)&1; b=(pp>>j)&1; c=np.bincount(2*a+b,minlength=4)/len(pp); pa=c.reshape(2,2).sum(1); pb=c.reshape(2,2).sum(0)
        mi=sum(c[2*x+y]*math.log2(c[2*x+y]/(pa[x]*pb[y])) for x in (0,1) for y in (0,1) if c[2*x+y]>0)
        if mi>best[0]: best=(mi,i,j)
print("   largest: %.2e bits between layers %d and %d (one sample of %d primes gives about %.1e by chance)"%(best[0],best[1],best[2],len(pp),1/(2*len(pp)*math.log(2))))
print("LAYER 3: by the number of 1-bits: share of primes among the odd numbers with k ones")
for k in range(2,K):
    m=top&(pc==k)
    if m.sum()>=200: print("   k=%2d: %7d numbers, %.4f prime  (%s)"%(k,m.sum(),s[m].mean(),"odd" if k%2 else "even"))
ev=P&(pc%2==0); od=P&(pc%2==1)
print("   primes with an odd number of ones: %d, with an even number: %d, share odd %.4f"%(od.sum(),ev.sum(),od.sum()/P.sum()))
