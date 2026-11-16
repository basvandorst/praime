import math, sys
from sympy import primerange, isprime
MEXP=[2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,136279841]
PM=200000; K=3000
small=[3,5,7,11,13,17,19,23,29,31,37]
def firstk(p):
    for k in range(1,K+1):
        q=2*k*p+1
        if q%8 not in (1,7): continue
        if any(q%s==0 for s in small): continue
        if pow(2,p,q)==1 and isprime(q): return k
    return 0
res={}
for p in primerange(100,PM): res[p]=firstk(p)
def show(M):
    cls={}
    for p,k in res.items(): c=cls.setdefault(p%M,[0,0]); c[0]+=1; c[1]+=(k==0)
    tot=sum(c[0] for c in cls.values()); sv=sum(c[1] for c in cls.values()); E=[e for e in MEXP if e>100]
    rows=[]
    for r,(n,s) in sorted(cls.items()):
        if n<50: continue
        rows.append((r,n,s/n,sum(1 for e in E if e%M==r)))
    wsum=sum(x[2] for x in rows); 
    print("mod %d: survivors (no factor 2kp+1, k <= %d) %.1f%% of %d primes"%(M,K,100*sv/tot,tot))
    chi=0
    for r,n,sh,ne in rows:
        ex=len(E)*sh/wsum; chi+=(ne-ex)**2/ex
        print("   class %3d: survive %.1f%%   exponents %2d (expected by survival %.1f, by equal classes %.1f)"%(r,100*sh,ne,ex,len(E)/len(rows)))
    print("   chi2 against the survival weights: %.1f on %d classes"%(chi,len(rows)))
for M in (4,12,24,60): show(M)
