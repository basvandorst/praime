import math, random, sys
from sympy import isprime, prevprime, nextprime, primerange, randprime, primepi
MEXP=[2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,136279841]
random.seed(1); E=[e for e in MEXP if e>120]; print("exponents above 120:",len(E))
def rp(e):
    while True:
        p=nextprime(random.randint(int(e/1.25),int(e*1.25))); g=p-prevprime(p)
        if random.random()<2/g: return p
def decoy(): return [rp(e) for e in E]
MODS=list(range(3,121))
def scan(S):
    best=(0,0,0,0); chis={}
    for m in MODS:
        cls=[r for r in range(m) if math.gcd(r,m)==1]; ex=len(S)/len(cls); cnt={r:0 for r in cls}
        for e in S:
            if e%m in cnt: cnt[e%m]+=1
        chi=sum((c-ex)**2/ex for c in cnt.values()); chis[m]=chi
        for r,c in cnt.items():
            z=(c-ex)/math.sqrt(ex*(1-1/len(cls)))
            if z>best[0]: best=(z,m,r,c)
    return best,chis
rb,rch=scan(E); print("real: strongest class: z=%.2f, %d exponents = %d mod %d (expected %.1f)"%(rb[0],rb[3],rb[2],rb[1],len(E)/sum(1 for r in range(rb[1]) if math.gcd(r,rb[1])==1)))
T=300; zs=[]; chs={m:[] for m in MODS}
for t in range(T):
    b,ch=scan(decoy()); zs.append(b[0])
    for m in MODS: chs[m].append(ch[m])
print("decoys (%d sets of random primes of the same sizes): strongest z: median %.2f, 90%% %.2f, max %.2f; share of decoys with z >= real: %.3f"%(T,sorted(zs)[T//2],sorted(zs)[int(.9*T)],max(zs),sum(z>=rb[0] for z in zs)/T))
pm={m:sum(c>=rch[m] for c in chs[m])/T for m in MODS}; lo=sorted(pm.items(),key=lambda t:t[1])[:8]
print("moduli with the most uneven spread (share of decoys at least as uneven):",[(m,round(p,3)) for m,p in lo])
for m in (4,8,12,24,30,60):
    cnt={}
    for e in E: cnt[e%m]=cnt.get(e%m,0)+1
    print("  mod %d:"%m,dict(sorted(cnt.items())),"p=%.2f"%pm.get(m,float('nan')))
print("== neighbourhood of the exponents")
def nb(S):
    tw=0; g1=g2=0.0; dens=0.0; k=0
    for e in S:
        a=prevprime(e); b=nextprime(e); L=math.log(e); tw+=(a==e-2)+(b==e+2); g1+=(e-a)/L; g2+=(b-e)/L
        W=int(50*L); c=sum(1 for _ in primerange(e-W,e+W+1))-1; dens+=c/(2*W/L); k+=1
    return tw,g1/k,g2/k,dens/k
r=nb(E); print("real: twin sides %d, mean gap before %.3f after %.3f (in ln p), primes within 50 ln p on each side relative to 1/ln p: %.3f"%r)
D=[nb(decoy()) for _ in range(120)]
for i,name in enumerate(("twin sides","gap before","gap after","density")):
    v=sorted(d[i] for d in D); print("  decoys %s: median %.3f, 5%%..95%% %.3f..%.3f; share >= real %.3f"%(name,v[len(v)//2],v[int(.05*len(v))],v[int(.95*len(v))],sum(x>=r[i] for x in v)/len(v)))
