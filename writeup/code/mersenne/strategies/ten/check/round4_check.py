import math, random
import gmpy2
from sympy import primerange, n_order, isprime, primorial, prime
MEXP=[2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,136279841]
N=500000; P=list(primerange(2,N+1)); ip=lambda v:v>1 and gmpy2.is_prime(v,25)
print("== expclass: primes 3..500,000, 2^p - 1 has no divisor 2kp+1 with k <= 200")
cls={}
for p in P[1:]:
    Mp=(1<<p)-1 if p<29 else None; s=True
    for k in range(1,201):
        q=2*k*p+1
        if Mp is not None and q>=Mp: break
        if pow(2,p,q)==1: s=False; break
    for M in (12,24,60):
        if math.gcd(p,M)==1: c=cls.setdefault((M,p%M),[0,0]); c[0]+=1; c[1]+=s
for M in (12,24,60):
    rows=sorted(((r,100*v[1]/v[0]) for (m,r),v in cls.items() if m==M),key=lambda t:-t[1]); print(" mod %d:"%M," ".join("%d: %.1f%%"%t for t in (rows if M<60 else rows[:4]+rows[-4:])))
print("== per10 / per2: period of 1/p even")
for base,M in ((10,40),(2,24),(2,8)):
    up={};tot=ev=0
    for p in P[1:]:
        if p==5 and base==10: continue
        if math.gcd(p,M)>1: continue
        e=n_order(base,p)%2==0; c=up.setdefault(p%M,[0,0]); c[0]+=1; c[1]+=e; tot+=1; ev+=e
    print(" base %d mod %d: even %.1f%%; up only %s; down only %s; mixed %s"%(base,M,100*ev/tot,sorted(r for r,v in up.items() if v[1]==v[0]),sorted(r for r,v in up.items() if v[1]==0),{r:round(100*v[1]/v[0],1) for r,v in sorted(up.items()) if 0<v[1]<v[0]}))
print("== pow2: primes 2^n + d, 2^n - d, n <= 400, odd d <= 299 (gmpy2, 25 rounds)")
cnt={}; rows={}
for n in range(1,401):
    B=1<<n
    for d in range(1,300,2):
        for sg in (1,-1):
            if ip(B+sg*d): cnt[(sg,d)]=cnt.get((sg,d),0)+1; rows.setdefault((sg,d),[]).append(n)
print(" total %d; Mersenne row %s; Fermat row %s"%(sum(cnt.values()),rows[(-1,1)],rows[(1,1)]))
print(" empty rows:",[("+" if sg>0 else "-")+str(d) for sg in (1,-1) for d in range(1,300,2) if (sg,d) not in cnt],"; 2^n - 7:",rows.get((-1,7)),"; 2^n + 7:",rows.get((1,7)))
print(" fullest:",sorted(cnt.items(),key=lambda t:-t[1])[:4]," classes mod 60 open for 2^n - 7:",sum(1 for n in range(60,120) if all(((1<<n)-7)%q for q in (3,5,7,11,13))))
print("== riesel: k 2^n + 1, k 2^n - 1, n <= 400, odd k <= 299")
cnt={}
for n in range(1,401):
    B=1<<n
    for k in range(1,300,2):
        for sg in (1,-1):
            if ip(k*B+sg): cnt[(sg,k)]=cnt.get((sg,k),0)+1
print(" total %d; empty rows above %s below %s; fullest %s"%(sum(cnt.values()),[k for k in range(1,300,2) if (1,k) not in cnt],[k for k in range(1,300,2) if (-1,k) not in cnt],sorted(cnt.items(),key=lambda t:-t[1])[:4]))
print("== primorial: p# + d, p# - d, first 100 primes, odd d <= 499")
tot=0; B=1; plus=[]; minus=[]; fort=[]; comp=0
for n in range(1,101):
    p=prime(n); B*=p; f=0
    for d in range(1,500,2):
        for sg in (1,-1):
            if ip(B+sg*d):
                tot+=1
                if d==1: (plus if sg>0 else minus).append(p)
                elif sg>0 and not f: f=d
                if d>1 and n>=5 and not isprime(d): comp+=1
    fort.append(f)
print(" total %d; p#+1 prime for p = %s; p#-1 for p = %s; Fortune numbers %s; composite distances from column 5: %d"%(tot,plus,minus,fort[:14],comp))
print("== glue: column runs with step W")
S=set(P)
for W in (210,2310,708):
    at={}; best=(0,0)
    for p in P:
        l=1
        while p+l*W in S: l+=1
        for k in range(2,l+1): at[k]=at.get(k,0)+1
        if l>best[0] and (p-W) not in S: best=(l,p)
    print(" W=%d: runs of at least k:"%W,{k:at[k] for k in sorted(at) if k>=3},"longest %d from %d"%best)
print("== stairs: runs of rising gaps")
g=[P[i+1]-P[i] for i in range(len(P)-1)]
def cnt_r(g):
    c={}; r=0
    for i in range(1,len(g)):
        r=r+1 if g[i]>g[i-1] else 0
        for k in range(1,r+1): c[k]=c.get(k,0)+1
    return c
real=cnt_r(g); random.seed(5); sh=[]
for b in range(0,len(g),1000):
    blk=g[b:b+1000]; random.shuffle(blk); sh+=blk
print(" real:",{k:real[k] for k in sorted(real) if k>=2}," shuffled in blocks of 1000:",{k:v for k,v in sorted(cnt_r(sh).items()) if k>=2})
i=next(i for i in range(len(g)-6) if all(g[i+j+1]>g[i+j] for j in range(6))); print(" first 6 rising steps after",P[i],g[i:i+7])
print("== shield")
for X in (500000,1000000):
    Q=list(primerange(2,X+1)); gg=[Q[i+1]-Q[i] for i in range(len(Q)-1)]; al=[True]*len(Q)
    for i in sorted(range(len(gg)),key=lambda i:(gg[i],i)):
        if al[i] and al[i+1]: al[i]=al[i+1]=False
    print(" up to %d: %d of %d stay (%.2f%%); exponents shielded: %s"%(X,sum(al),len(Q),100*sum(al)/len(Q),[q for q,a in zip(Q,al) if a and q in MEXP]))
print("== lucky")
for X in (500000,1000000):
    L=list(range(1,X+1,2)); i=1
    while i<len(L) and L[i]<=len(L):
        st=L[i]; del L[st-1::st]; i+=1
    LS=set(L); Q=list(primerange(2,X+1)); lp=[q for q in Q if q in LS]; run=best=0; bs=0
    for q in Q:
        if q in LS: run=0
        else:
            if run==0: st0=q
            run+=1
            if run>best: best=run; bs=st0
    print(" up to %d: %d lucky numbers, %d lucky primes (all = 1 mod 6 after 3: %s); lucky = 1 mod 6: %d, primes = 1 mod 6: %d; longest run of unlucky primes %d from %d; lucky exponents %s"%(X,len(L),len(lp),all(q%6==1 for q in lp if q>3),sum(1 for v in L if v%6==1),sum(1 for q in Q if q%6==1),best,bs,[e for e in MEXP if e in LS]))
