import numpy as np, math, sys
from lib import sieve
N=5*10**7; K=120; PAD=K+2
s=sieve(N+PAD); P=np.nonzero(s[:N+1])[0]
lam=np.zeros(N+PAD+1,dtype=np.int8); allp=np.nonzero(s)[0]
for p in allp[allp<=(N+PAD)//2]:
    q=int(p)
    while q<=N+PAD: lam[q::q]^=1; q*=int(p)
lam[allp[allp>(N+PAD)//2]]^=1
lam=(1-2*lam).astype(np.int8); lam[0]=0
Q=P[P>K+1]
ks=[k for k in range(-K,K+1) if k!=0]
S={k:int(lam[Q+k].astype(np.int64).sum()) for k in ks}
np.save("lam_S.npy",np.array([[k,S[k]] for k in ks]))
sd=math.sqrt(len(Q)); print("primes used",len(Q),"sd of coin flips %.0f"%sd)
small=[q for q in range(2,20000) if s[q]]
leg={}
def legtab(q):
    t=np.full(q,-1,dtype=np.int8); t[(np.arange(q,dtype=np.int64)**2)%q]=1; t[0]=0; return t
LT={q:legtab(q) for q in small if q>2}
def sqfree_mu(J):
    mu=np.ones(J+1,dtype=np.int64)
    for q in range(2,J+1):
        if s[q]:
            mu[q::q]*=-1; mu[q*q::q*q]=0
    return mu
J=3000; mu=sqfree_mu(J); karr=np.array(ks,dtype=np.int64)
def Cvec(j):
    C=np.ones(len(ks))
    w2=np.where((j%2==0)&(karr%2==0),2,np.where((j%2==0)|(karr%2==0),0,1)) if j%2==0 else np.where(karr%2==0,0,1)
    if j%2==0: w2=np.where(karr%2==0,2,0)
    else: w2=np.ones(len(ks),dtype=np.int64)
    C*= (1-w2/2)/(1-1/2)
    for q in small[1:]:
        if j%q==0: w=np.where(karr%q==0,q,0)
        else:
            jinv=pow(j,-1,q); w=1+LT[q][(karr*jinv)%q]
        C*=(1-w/q)/(1-1/q)
    return C
def Li(j,k,X):
    m=np.arange(1,int(math.isqrt((X+k)//j))+2,dtype=np.float64); v=j*m*m-k; v=v[(v>K+1)&(v<=X)]
    return float((1/np.log(v)).sum())
def count(j,k,X):
    m=np.arange(1,int(math.isqrt((X+k)//j))+2,dtype=np.int64); v=j*m*m-k; v=v[(v>K+1)&(v<=X)]
    return int(s[v].sum())
B={30:np.zeros(len(ks)),300:np.zeros(len(ks)),3000:np.zeros(len(ks))}; A={30:np.zeros(len(ks)),300:np.zeros(len(ks)),3000:np.zeros(len(ks))}
chk=[]
for j in range(1,J+1):
    if mu[j]==0: continue
    C=Cvec(j)
    for i,k in enumerate(ks):
        if k>0:
            t=k//j
            if k%j==0 and math.isqrt(t)**2==t: C[i]=0
    li=np.array([Li(j,k,N) for k in ks]); pred=mu[j]*C*li
    cnt=np.array([count(j,k,N) for k in ks]) if j<=300 else None
    for JJ in B:
        if j<=JJ:
            B[JJ]+=pred
            if cnt is not None and JJ<=300: A[JJ]+=mu[j]*cnt
    if j<=6: chk.append((j,[(k,round(float(C[ks.index(k)]),3),round(float(C[ks.index(k)]*li[ks.index(k)])),int(cnt[ks.index(k)])) for k in (-1,1,-2,2,-67,3)]))
for c in chk: print(" j=%d (k, C, BH count, real count):"%c[0],c[1])
Sv=np.array([S[k] for k in ks],float)
for JJ in (30,300,3000):
    b=B[JJ]; beta=float((b*Sv).sum()/(b*b).sum()); z=float((b*Sv).sum()/(sd*math.sqrt((b*b).sum()))); r=float(np.corrcoef(b,Sv)[0,1])
    print("J=%d: rms prediction %.0f, slope of S on prediction %.2f +- %.2f, z = %.1f, correlation %.3f"%(JJ,math.sqrt((b*b).mean()),beta,sd/math.sqrt((b*b).sum()),z,r))
    if JJ<=300:
        a=A[JJ]; print("    exact structured part vs BH prediction: corr %.3f, rms diff %.0f; remainder S - A: mean %.0f, sd %.0f"%(np.corrcoef(a,b)[0,1],math.sqrt(((a-b)**2).mean()),(Sv-a).mean(),(Sv-a).std()))
print("S: mean %.0f, sd across k %.0f, share negative %.3f"%(Sv.mean(),Sv.std(),(Sv<0).mean()))
np.save("lam_B.npy",np.array([B[30],B[300],B[3000]]))
o=np.argsort(B[3000]); print("most negative predicted:",[(ks[i],round(B[3000][i]),S[ks[i]]) for i in o[:8]]); print("most positive predicted:",[(ks[i],round(B[3000][i]),S[ks[i]]) for i in o[-8:]])
