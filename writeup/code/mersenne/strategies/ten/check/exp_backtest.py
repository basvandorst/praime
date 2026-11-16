import math, random, re
from sympy import nextprime, prevprime
MEXP=[2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,136279841]
E=[e for e in MEXP if e>500]; TR=E[:len(E)//2]; TE=E[len(E)//2:]
print("exponents above 500: %d; train on the first %d (up to %d), test on the last %d"%(len(E),len(TR),TR[-1],len(TE)))
print("== (a) mined rules: the classes r mod m (m = 3..500) that hold most training exponents")
rules=[]
for m in range(3,501):
    phi=sum(1 for r in range(m) if math.gcd(r,m)==1); cnt={}
    for e in TR: cnt[e%m]=cnt.get(e%m,0)+1
    for r,c in cnt.items():
        ex=len(TR)/phi; z=(c-ex)/math.sqrt(ex); rules.append((z,m,r,c,phi))
rules.sort(reverse=True)
for K in (10,50,200):
    hit=sum(sum(1 for e in TE if e%m==r) for z,m,r,c,phi in rules[:K]); ex=sum(len(TE)/phi for z,m,r,c,phi in rules[:K]); tin=sum(c for z,m,r,c,phi in rules[:K]); exin=sum(len(TR)/phi for z,m,r,c,phi in rules[:K])
    print(" best %3d rules: training hits %d where %.1f expected by chance; test hits %d where %.1f expected"%(K,tin,exin,hit,ex))
print(" the 5 strongest:",[(m,r,c) for z,m,r,c,phi in rules[:5]],"(modulus, class, training exponents)")
print("== (b) the factor rule: weight of a class = share of primes p in it whose 2^p - 1 has no factor 2kp+1, k <= 3000 (from exp_class.txt, p < 2e5)")
txt=open('../logs/exp_class.txt').read(); W={}
for M in (4,12,24,60):
    blk=txt.split("mod %d:"%M)[1].split("chi2")[0]; W[M]={int(a):float(b) for a,b in re.findall(r"class\s+(\d+): survive ([0-9.]+)%",blk)}
def llr(S,M): 
    w=W[M]; wb=sum(w.values())/len(w); return sum(math.log(w[e%M]/wb) for e in S if e%M in w)
random.seed(3)
def rp(e):
    while True:
        p=nextprime(random.randint(int(e/1.25),int(e*1.25))); g=p-prevprime(p)
        if random.random()<2/g: return p
for M in (4,12,24,60):
    for nm,S in (("train",TR),("test",TE),("all",E)):
        v=llr(S,M); dec=sorted(llr([rp(e) for e in S],M) for _ in range(300)); 
        print(" mod %2d %-5s: log-likelihood gain %+.2f; decoy primes: median %+.2f, 95%% point %+.2f; share of decoys >= real %.3f"%(M,nm,v,dec[150],dec[285],sum(d>=v for d in dec)/300))
