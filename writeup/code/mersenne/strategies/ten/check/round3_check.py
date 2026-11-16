import math, sys
from sympy import factorint, primerange, isprime, mobius
import cypari2; pari=cypari2.Pari()
print("== A. lean by distance, N = 10,000, distances 1..250 (the page's numbers for n=10000)")
N=10000; K=250; primes=list(primerange(2,N+1))
lam={n:(-1)**sum(factorint(n).values()) for n in range(1,N+K+1)}
S=[sum(lam[p+k] for p in primes) for k in range(1,K+1)]
small=list(primerange(2,151)); J=3000
roots={q:[[sum(1 for m in range(q) if (a*m*m-b)%q==0) for b in range(q)] for a in range(q)] for q in small}
B=[0.0]*K
for j in range(1,J+1):
    mu=mobius(j)
    if mu==0: continue
    mm=math.isqrt((N+K)//j)+1
    for i in range(K):
        k=i+1
        if k%j==0 and math.isqrt(k//j)**2==k//j: continue
        C=1.0
        for q in small:
            w=roots[q][j%q][k%q]; C*=(1-w/q)/(1-1/q)
            if C==0: break
        if C==0: continue
        li=0.0
        for m in range(1,mm+1):
            v=j*m*m-k
            if 2<=v<=N: li+=1/math.log(v)
        B[i]+=mu*C*li
sb=sum(a*b for a,b in zip(S,B)); bb=sum(b*b for b in B); mS=sum(S)/K; mB=sum(B)/K
cs=sum((a-mS)**2 for a in S); cb=sum((b-mB)**2 for b in B); cc=sum((a-mS)*(b-mB) for a,b in zip(S,B)); rs=sum((a-b)**2 for a,b in zip(S,B))
print(" primes %d, coin spread %.0f, spread found %.0f; slope %.2f +- %.2f, correlation %.2f, spread left %.0f"%(len(primes),math.sqrt(len(primes)),math.sqrt(cs/K),sb/bb,math.sqrt(len(primes)/bb),cc/math.sqrt(cs*cb),math.sqrt(rs/K)))
for d in (30,6,3,1):
    idx=[i for i in range(K) if math.gcd(i+1,30)==d]; print("  gcd(n,30)=%d: mean %.0f expected %.0f"%(d,sum(S[i] for i in idx)/len(idx),sum(B[i] for i in idx)/len(idx)))
print("  n=1: %d expected %.0f;  p-1: %d"%(S[0],B[0],sum(lam[p-1] for p in primes if p>1 and p-1>=1)))
print("== B. square walks, all primes 100 < Q < 20000 (Euler's criterion, PARI class numbers)")
bad3=bad1=n3=n1=0
for Q in primerange(101,20000):
    st=[1 if pow(n,(Q-1)//2,Q)==1 else -1 for n in range(1,Q)]; w=0; W=[0]
    for x in st: w+=x; W.append(w)
    if Q%4==3:
        h=int(pari.qfbclassno(-Q)); two=1 if Q%8==7 else -1; n3+=1
        if W[(Q-1)//2]!=(2-two)*h or sum(W[1:Q])!=Q*h: bad3+=1
    else:
        h=int(pari.qfbclassno(-4*Q)); n1+=1
        if W[(Q-1)//2]!=0 or 2*W[Q//4]!=h: bad1+=1
print(" Q = 3 mod 4: %d primes, exceptions to (middle = (2 - (2/Q)) h, area = Q h): %d;  Q = 1 mod 4: %d primes, exceptions to (middle = 0, quarter = h(-4Q)/2): %d"%(n3,bad3,n1,bad1))
for Q in (4423,9941,216091,756839):
    t=[0]*Q
    for x in range(1,Q): t[x*x%Q]=1
    w=0; ar=0; mid=qt=0; mx=mn=0; pos=0
    for n in range(1,Q):
        w+=1 if t[n] else -1; ar+=w; mx=max(mx,w); mn=min(mn,w); pos+=w>=0
        if n==(Q-1)//2: mid=w
        if n==Q//4: qt=w
    print("  Q=%d: middle %d, quarter %d, area %d, h(-Q) %s, h(-4Q) %d, highest %d lowest %d, at or above the line %.1f%%"%(Q,mid,qt,ar,int(pari.qfbclassno(-Q)) if Q%4==3 else "-",int(pari.qfbclassno(-4*Q)),mx,mn,100*pos/(Q-1)))
print("== C. curves: counts by definition for p < 200, PARI for all odd p <= 100000")
CUR={11:[0,-1,1,0,0],37:[0,0,1,-1,0],389:[0,1,1,-2,0],5077:[0,0,1,-7,6],32:[0,0,0,-1,0],27:[0,0,1,0,0]}
for cd,c in CUR.items():
    E=pari.ellinit(c); a1,a2,a3,a4,a6=c; bad=0
    for p in primerange(3,200):
        cnt=sum(1 for x in range(p) for y in range(p) if (y*y+a3*y-(x**3+a2*x*x+a4*x+a6))%p==0)
        if p-cnt!=int(pari.ellap(E,p)) and int(pari.ellglobalred(E)[0])%p: bad+=1
    run=0.0; z=0; n=0; r3=0; hb=[0]*8
    for p in primerange(3,100001):
        a=int(pari.ellap(E,p)); run+=a*math.log(p)/p; z+=a==0; n+=1
        if p<=1000: r3=run
        y=a/(2*math.sqrt(p)); hb[max(0,min(7,int((y+1)*4)))]+=1
    print(" conductor %d (PARI: %d, rank %d): direct-count mismatches at good p < 200: %d; primes %d, a=0: %d; sum a ln p/p = %.2f = %.2f ln N, slope from 1000: %.2f; bands %s"%(cd,int(pari.ellglobalred(E)[0]),int(pari.ellanalyticrank(E)[0]),bad,n,z,run,run/math.log(100000),(run-r3)/(math.log(100000)-math.log(1000))," ".join("%.1f"%(100*v/n) for v in hb)))
print("== D. family trees up to 500,000 (smallest-prime-factor table, plain Python)")
N=500000; spf=list(range(N+2))
for i in range(2,int((N+1)**.5)+1):
    if spf[i]==i:
        for j in range(i*i,N+2,i):
            if spf[j]==j: spf[j]=i
def gpf(n):
    g=1
    while n>1: g=spf[n]; n//=spf[n]
    return g
pr=[q for q in range(3,N+1) if spf[q]==q]
for sh,name in ((1,"q + 1"),(-1,"q - 1")):
    root={}; dep={}
    for q in pr:
        pa=gpf(q+sh)
        if pa==2: root[q]=q; dep[q]=1
        else: root[q]=root[pa]; dep[q]=dep[pa]+1
    for X in (10**4,10**5,N):
        cnt={}; tot=0
        for q in pr:
            if q>X: break
            cnt[root[q]]=cnt.get(root[q],0)+1; tot+=1
        print("  %s, x=%d:"%(name,X)," ".join("%d: %.2f%%"%(r,100*c/tot) for r,c in sorted(cnt.items(),key=lambda t:-t[1])[:5]))
    dc={}
    for q in pr: dc[dep[q]]=dc.get(dep[q],0)+1
    md=max(dc); print("   primes per band:",[dc[d] for d in sorted(dc)],"; longest chain %d steps, first at %d"%(md,min(q for q in pr if dep[q]==md)))
