import numpy as np, math, sys
from lib import sieve
N=int(float(sys.argv[1])); K=int(sys.argv[2]); JM=int(sys.argv[3]); QB=int(sys.argv[4])
s=sieve(N+K+2); allp=np.nonzero(s)[0]; P=allp[allp<=N]
lam=np.zeros(N+K+3,dtype=np.int8)
for p in allp[allp<=(N+K+2)//2]:
    q=int(p)
    while q<=N+K+2: lam[q::q]^=1; q*=int(p)
lam[allp[allp>(N+K+2)//2]]^=1
lam=(1-2*lam).astype(np.int8); lam[0]=0
ks=np.arange(1,K+1); S=np.array([int(lam[P+k].astype(np.int64).sum()) for k in ks],float); sd=math.sqrt(len(P))
small=[int(q) for q in allp[allp<=QB]]
LT={q:None for q in small}
for q in small:
    t=np.full(q,-1,dtype=np.int8); t[(np.arange(q,dtype=np.int64)**2)%q]=1; t[0]=0; LT[q]=t
mu=np.ones(JM+1,dtype=np.int64)
for q in allp[allp<=JM]:
    q=int(q); mu[q::q]*=-1; mu[q*q::q*q]=0
il=np.zeros(N+1); il[2:]=1/np.log(np.arange(2,N+1))
B=np.zeros(K); marks=[30,100,300,1000,3000,10000,30000,100000]; g=np.array([math.gcd(int(k),30) for k in ks])
for j in range(1,JM+1):
    if mu[j]:
        C=np.ones(K)
        for q in small:
            if q==2: w=np.ones(K) if j%2 else np.where(ks%2==0,2,0)
            elif j%q==0: w=np.where(ks%q==0,q,0)
            else: w=1+LT[q][(ks*pow(j,-1,q))%q]
            C*=(1-w/q)/(1-1/q)
        t=ks//j; C[(ks%j==0)&(np.rint(np.sqrt(t))**2==t)]=0
        m=np.arange(1,int(math.isqrt((N+K)//j))+2,dtype=np.int64); v=j*m[:,None]**2-ks[None,:]; ok=(v>=2)&(v<=N); li=np.where(ok,il[np.clip(v,0,N)],0).sum(axis=0)
        B+=mu[j]*C*li
    if j in marks and j<=JM:
        beta=(B*S).sum()/(B*B).sum(); print("J=%6d slope %.3f +- %.3f corr %.3f rms B %.0f sd resid %.0f (coin %.0f) | gcd30: S %.0f B %.0f | gcd6: S %.0f B %.0f | gcd1: S %.0f B %.0f | k=1: S %.0f B %.0f"%(j,beta,sd/math.sqrt((B*B).sum()),np.corrcoef(B,S)[0,1],math.sqrt((B*B).mean()),(S-B).std(),sd,S[g==30].mean(),B[g==30].mean(),S[g==6].mean(),B[g==6].mean(),S[g==1].mean(),B[g==1].mean(),S[0],B[0]),flush=True)
