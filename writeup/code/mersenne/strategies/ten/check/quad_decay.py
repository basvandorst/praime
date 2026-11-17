import numpy as np, math
NM=100_000_000
s=np.ones(NM+1,dtype=bool); s[:2]=False
for i in range(2,int(NM**0.5)+1):
    if s[i]: s[i*i::i]=False
pr=np.nonzero(s[:20000])[0]; ks=list(range(-4,5))
print("pairs n, n + 2 floor(sqrt n) + 2k for N/10 < n <= N; root mean square over k = -4..4 of the deviation of the pair ratio (%)")
for N in (10**6,3*10**6,10**7,3*10**7,10**8):
    y0=int(round(N**0.4027)); sv=np.ones(N+1,dtype=bool); sv[:2]=False
    for q in pr[pr<=y0]: sv[2*q::q]=False
    acc={k:np.zeros(6,dtype=np.int64) for k in ks}; tv=0; CH=5_000_000
    for lo in range(N//10,N-30000,CH):
        n=np.arange(lo,min(lo+CH,N-30000),dtype=np.int64); r=np.floor(np.sqrt(n)).astype(np.int64); r-=(r*r>n); r+=((r+1)**2<=n); pn=s[n]; vn=sv[n]; tv+=int(vn.sum())
        for k in ks:
            m=n+2*r+2*k; pm=s[m]; vm=sv[m]; acc[k]+=np.array([len(n),pn.sum(),pm.sum(),(pn&pm).sum(),vm.sum(),(vn&vm).sum()],dtype=np.int64)
    Rp=[acc[k][3]*acc[k][0]/(acc[k][1]*acc[k][2]) for k in ks]; Rs=[acc[k][5]*acc[k][0]/(tv*acc[k][4]) for k in ks]; mp=sum(Rp)/9; ms=sum(Rs)/9
    rp=math.sqrt(sum((100*(v/mp-1))**2 for v in Rp)/9); rs=math.sqrt(sum((100*(v/ms-1))**2 for v in Rs)/9)
    print("  N = %9d, sieve to %4d: survivors %.2f   primes %.2f (noise %.2f)   survivors at k = 0, 1: %+.2f %+.2f"%(N,y0,rs,rp,100/math.sqrt(acc[0][3]),100*(Rs[4]/ms-1),100*(Rs[5]/ms-1)),flush=True)
