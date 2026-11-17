import numpy as np, math, sys
N=100_000_000; y0=int(round(N**0.4027))
def sieve(n):
    s=np.ones(n+1,dtype=bool); s[:2]=False
    for i in range(2,int(n**0.5)+1):
        if s[i]: s[i*i::i]=False
    return s
s=sieve(N); small=np.nonzero(s[:y0+1])[0]; sv=np.ones(N+1,dtype=bool); sv[:2]=False
for q in small: sv[2*q::q]=False
print("N = %d, sieve to %d: primes %d, survivors %d"%(N,y0,int(s.sum()),int(sv.sum())),flush=True)
ks=list(range(-4,5)); CH=5_000_000
acc={k:np.zeros(6,dtype=np.int64) for k in ks}
for lo in range(10**7,N-100000,CH):
    n=np.arange(lo,min(lo+CH,N-100000),dtype=np.int64); r=np.floor(np.sqrt(n)).astype(np.int64); r-=(r*r>n); r+=((r+1)**2<=n); pn=s[n]; vn=sv[n]
    for k in ks:
        m=n+2*r+2*k; pm=s[m]; vm=sv[m]; acc[k]+=np.array([len(n),pn.sum(),pm.sum(),(pn&pm).sum(),vm.sum(),(vn&vm).sum()],dtype=np.int64)
    va=int(vn.sum())
    for k in ks: acc[k][1]+=0
A={k:acc[k] for k in ks}
tot_v=int(sv[10**7:N-100000].sum())
Rp=[A[k][3]*A[k][0]/(A[k][1]*A[k][2]) for k in ks]; Rs=[A[k][5]*A[k][0]/(tot_v*A[k][4]) for k in ks]
mp=sum(Rp)/9; ms=sum(Rs)/9
print("pairs n, n + 2 floor(sqrt n) + 2k, n from 1e7 to 1e8; prime pairs per k about %d (noise %.2f %%)"%(A[0][3],100/math.sqrt(A[0][3])))
print(" k:          "+" ".join("%+6d"%k for k in ks))
print(" primes    %: "+" ".join("%+6.2f"%(100*(v/mp-1)) for v in Rp)+"   rms %.2f"%math.sqrt(sum((100*(v/mp-1))**2 for v in Rp)/9))
print(" survivors %: "+" ".join("%+6.2f"%(100*(v/ms-1)) for v in Rs)+"   rms %.2f"%math.sqrt(sum((100*(v/ms-1))**2 for v in Rs)/9))
