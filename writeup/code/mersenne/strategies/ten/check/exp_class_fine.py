import math, random
from sympy import primerange, nextprime, prevprime
MEXP=[2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,136279841]
PM=500000; K=2000
def surv(p):
    for k in range(1,K+1):
        q=2*k*p+1
        if q%8 not in (1,7) or q%3==0 or q%5==0 or q%7==0: continue
        if pow(2,p,q)==1: return 0
    return 1
res={p:surv(p) for p in primerange(1000,PM)}
print("primes 1000..%d: %d, share with no factor 2kp+1, k <= %d: %.3f"%(PM,len(res),K,sum(res.values())/len(res)))
E=[e for e in MEXP if e>500]; TR=E[:20]; TE=E[20:]; random.seed(8)
def rp(e):
    while True:
        p=nextprime(random.randint(int(e/1.25),int(e*1.25))); g=p-prevprime(p)
        if random.random()<2/g: return p
for M in (4,12,60,420):
    cl={}
    for p,s in res.items():
        if math.gcd(p,M)==1: c=cl.setdefault(p%M,[0,0]); c[0]+=1; c[1]+=s
    w={r:(v[1]+1)/(v[0]+2) for r,v in cl.items()}; wb=sum(w.values())/len(w)
    def llr(S): return sum(math.log(w[e%M]/wb) for e in S if e%M in w)
    v=(llr(TR),llr(TE),llr(E)); dec=sorted(llr([rp(e) for e in E]) for _ in range(300))
    print(" mod %3d (%2d classes, weights from %d to %d %%): gain on the first 20 exponents %+.2f, on the last 20 %+.2f, on all %+.2f; share of decoy sets with at least that gain %.3f"%(M,len(w),round(100*min(w.values())),round(100*max(w.values())),v[0],v[1],v[2],sum(d>=v[2] for d in dec)/300))
