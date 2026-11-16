import numpy as np, math
from lib import sieve
N=5*10**7; s=sieve(N+2); P=np.nonzero(s[:N+1])[0]
lam=np.zeros(N+2,dtype=np.int8)
for p in P[P<=N//2+1]:
    q=int(p)
    while q<=N+1: lam[q::q]^=1; q*=int(p)
lam[P[P>N//2+1]]^=1
lam=(1-2*lam.astype(np.int8)); lam[0]=0
print("== (b) parity of the number of prime factors around the primes, N =",N)
L=int(lam[1:N+1].astype(np.int64).sum()); comp=~s[:N+1]; comp[:2]=False
up=int((lam[:N+1][comp]==1).sum()); dn=int((lam[:N+1][comp]==-1).sum()); print(" composites with an even number of prime factors %d, odd %d, difference %d; pi(N)+L(N)-1 = %d (L(N) = %d)"%(up,dn,up-dn,len(P)+L-1,L))
for X in (10**6,10**7,N):
    Q=P[(P>2)&(P<=X)]; a=int(lam[Q+1].astype(np.int64).sum()); b=int(lam[Q-1].astype(np.int64).sum()); c=int((lam[Q+1].astype(np.int64)*lam[Q-1]).sum()); sd=math.sqrt(len(Q))
    print(" x=%.0e: sum lambda(p+1) = %d, sum lambda(p-1) = %d, sum lambda(p-1)lambda(p+1) = %d   (one sd of coin flips %.0f)"%(X,a,b,c,sd))
mu=lambda k: 0 if any(k%(d*d)==0 for d in range(2,int(k**.5)+1)) else (-1)**sum(1 for d in range(2,k+1) if k%d==0 and all(d%e for e in range(2,int(d**.5)+1)))
def part(sign,K,X):
    tot=0; rows=[]
    for k in range(1,K+1):
        m=mu(k)
        if not m: continue
        mm=np.arange(1,int(math.isqrt((X+1)//k))+2,dtype=np.int64); v=k*mm*mm+sign; v=v[(v>2)&(v<=X)]; c=int(s[v].sum()); tot+=m*c; rows.append((k,m*c))
    return tot,rows
for sign,name in ((-1,"p+1 = k m^2  (primes k m^2 - 1)"),(1,"p-1 = k m^2  (primes k m^2 + 1)")):
    for K in (10,100,1000):
        t,rows=part(sign,K,N); 
        if K==10: print(" ",name,": signed counts for k=1..10:",rows)
        print("    part of the sum from square-free k <= %d: %d"%(K,t))
# log-averaged walks against the k<=1000 part
xs=np.unique(np.round(np.exp(np.linspace(math.log(1e5),math.log(N),200))).astype(np.int64))
Q=P[P>2]; cs_p=np.cumsum(lam[Q+1].astype(np.int64)); cs_m=np.cumsum(lam[Q-1].astype(np.int64)); idx=np.searchsorted(Q,xs,side='right')-1
for sign,cs,name in ((-1,cs_p,"lambda(p+1)"),(1,cs_m,"lambda(p-1)")):
    walk=cs[idx]; pred=np.array([part(sign,300,int(x))[0] for x in xs]); sc=np.sqrt(xs)/np.log(xs)
    print(" %s: mean of walk/(sqrt x/ln x) over ln x: %.2f; of the k<=300 part: %.2f; of the rest: %.2f (sd of the rest %.2f); share of log-time the walk is negative: %.2f"%(name,np.mean(walk/sc),np.mean(pred/sc),np.mean((walk-pred)/sc),np.std((walk-pred)/sc),np.mean(walk<0)))
print("== (e) Mersenne against Wagstaff exponents")
M=[2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377]
W=[3,5,7,11,13,17,19,23,31,43,61,79,101,127,167,191,199,313,347,701,1709,2617,3539,5807,10501,10691,11279,12391,14479,42737,83339,95369,117239,127031,138937,141079,267017,269987,374321,986191,4031399]
both=sorted(set(M)&set(W)); print(" both:",both)
eg=math.exp(0.5772156649); X=4*10**6
pm=lambda p: min(1,eg*math.log2((6 if p%4==1 else 2)*p)/p)
nM=sum(pm(int(p)) for p in P[(P>2)&(P<=X)]); cW=len([w for w in W if w<=X])/sum(eg*math.log2(int(p))/int(p) for p in P[(P>2)&(P<=X)])
exp_both=sum(pm(int(p))*min(1,cW*eg*math.log2(int(p))/int(p)) for p in P[(P>2)&(P<=X)])
print(" Mersenne exponents 3..4e6: %d (heuristic %.1f); Wagstaff: %d; expected in both if independent: %.1f; found %d, all <= 127"%(len([m for m in M if 2<m<=X]),nM,len(W),exp_both,len(both)))
form=lambda p: any(abs(p-2**k)==1 for k in range(1,30)) or any(abs(p-4**k)==3 for k in range(1,15))
print(" of the form 2^k+-1 or 4^k+-3 among primes <= 4e6:",[int(p) for p in P[P<=X] if form(int(p))])
