import numpy as np, math, gmpy2, sys
PMAX=60000; LB=30
def sieve(n):
    s=np.ones(n+1,dtype=bool); s[:2]=False
    for i in range(2,int(n**0.5)+1):
        if s[i]: s[i*i::i]=False
    return s
s=sieve(PMAX); P=[int(p) for p in np.nonzero(s)[0] if p>=1000]
SM=[3,5,7,11,13,17,19,23,29,31,37]
ff={}
for p in P:
    K=(1<<LB)//(2*p); found=0
    for k0 in range(1,K+1,200000):
        k=np.arange(k0,min(K,k0+199999)+1,dtype=np.int64); q=2*k*p+1; m=q%8; ok=(m==1)|(m==7)
        for r in SM: ok&=(q%r!=0)
        for qq in q[ok]:
            qq=int(qq)
            if pow(2,p,qq)==1: found=qq; break
        if found: break
    ff[p]=found
MEXP={1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497}
print("primes p from 1000 to %d: %d; smallest factor q = 2kp+1 of 2^p - 1 searched up to 2^%d"%(PMAX,len(P),LB))
print("exponents in range with a factor found (must be 0):",sum(1 for p in MEXP if ff[p]))
print("composite 'factors' (must be 0):",sum(1 for p in P if ff[p] and not gmpy2.is_prime(ff[p])))
print("share of p with no factor below B, measured | predicted ln(a p)/ln B with a = 2 (p = 3 mod 4), 6 (p = 1 mod 4)")
for cls,a in ((3,2),(1,6)):
    Q=[p for p in P if p%4==cls]
    for b in (18,20,22,24,26,28,30):
        B=1<<b; QQ=[p for p in Q if a*p<B/4]
        meas=sum(1 for p in QQ if ff[p]==0 or ff[p]>=B)/len(QQ); pred=sum(min(1.0,math.log(a*p)/math.log(B)) for p in QQ)/len(QQ)
        print("  p = %d mod 4, B = 2^%d: %5d primes, measured %.4f, predicted %.4f, ratio %.3f"%(cls,b,len(QQ),meas,pred,meas/pred))
for m in (12,):
    for c in (1,5,7,11):
        Q=[p for p in P if p%m==c]; print("  p = %2d mod 12: no factor below 2^%d for %.4f of %d primes"%(c,LB,sum(1 for p in Q if ff[p]==0)/len(Q),len(Q)))
near=sorted([p for p in P if ff[p]==0 and p not in MEXP])
print("primes with no factor below 2^%d that are not exponents: %d of %d; the exponents: %d"%(LB,len(near),len(P),len([p for p in MEXP if p>=1000 and p<=PMAX])))
import json; json.dump({str(p):int(f) for p,f in ff.items()},open('mers_law_ff.json','w'))
