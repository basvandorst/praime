import math, numpy as np, random
def sieve(n):
    s=np.ones(n+1,dtype=bool); s[:2]=False
    for i in range(2,int(n**0.5)+1):
        if s[i]: s[i*i::i]=False
    return s
S=sieve(3_000_000); isp=lambda n:bool(S[n])
EXP=[2,3,5,7,13,17,19,31,61,89,107,127]
print("A. DIAL (turn = 120 ticks): best ray that holds only prime hands, for every number of hands N")
best_by_N={}
for N in range(2,251):
    best=(0,0,[])
    for t in range(1,120):
        rays={}
        for k in range(1,N+1): rays.setdefault((k*t)%120,[]).append(k)
        for v in rays.values():
            if len(v)>best[0] and all(isp(x) for x in v): best=(len(v),t,v)
    best_by_N[N]=best
mx=max(b[0] for b in best_by_N.values()); ns=[N for N,b in best_by_N.items() if b[0]==mx]
print("   most: %d primes alone on a ray, for N = %d to %d, at tick %d: %s"%(mx,ns[0],ns[-1],best_by_N[ns[0]][1],best_by_N[ns[0]][2]))
print("   five alone: N in",[N for N,b in best_by_N.items() if b[0]==5][:3],"...; by N (30,60,100,157,200,250):",[(N,best_by_N[N][0],best_by_N[N][2]) for N in (30,60,100,157,200,250)])
print("B. SLIPPING WHEELS: period of wheel k (ticks until position and phase repeat) against k^2/gcd(k, phi(k))")
def phi(k): return sum(1 for j in range(1,k+1) if math.gcd(j,k)==1)
ok=True; full=[]
for k in range(1,401):
    f=phi(k); pos=0; t=0
    while True:
        t+=1
        if math.gcd(t,k)==1: pos=(pos+1)%k
        if t%k==0 and pos==0: break
    if t!=k*k//math.gcd(k,f): ok=False
    if t==k*k: full.append(k)
cyc=[k for k in range(1,401) if math.gcd(k,phi(k))==1]
print("   formula holds for k <= 400:",ok,"; wheels that need the full k rounds:",full[:24],"...; equal to the numbers with gcd(k, phi(k)) = 1:",full==cyc,"; composite ones:",[k for k in full if k>1 and not isp(k)][:12])
print("C. FIBONACCI WHEELS: first tick wheel k is home, k <= 5000")
def alpha(k):
    a,b,t=0,1,0
    while True:
        a,b=b,(a+b)%k; t+=1
        if a==0: return t
al={k:alpha(k) for k in range(2,5001)}
r=sorted(((al[k]/k,k) for k in al),reverse=True)
print("   slowest (first tick / k):",[(k,al[k]) for _,k in r[:6]],"; wheels with first tick = 2k:",[k for k in al if al[k]==2*k],"; any above 2k:",any(al[k]>2*k for k in al))
print("   fastest:",[(k,al[k]) for _,k in r[-8:]])
print("D. ACCELERATING WHEELS: first tick wheel k is home, k <= 5000")
def acc(k):
    t=1
    while (t*(t+1)//2)%k: t+=1
    return t
ac={k:acc(k) for k in range(2,5001)}
r=sorted(((ac[k]/k,k) for k in ac),reverse=True)
print("   slowest:",[(k,ac[k]) for _,k in r[:8]],"; all powers of 2 have first tick 2k-1:",all(ac[2**a]==2**(a+1)-1 for a in range(1,13)),"; any other wheel above k:",[k for k in ac if ac[k]>k and (k&(k-1))][:5])
print("   odd primes at p-1:",all(ac[p]==p-1 for p in range(3,5001) if isp(p)),"; fastest:",[(k,ac[k]) for _,k in r[-6:]])
