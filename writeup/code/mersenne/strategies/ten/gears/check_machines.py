"""Independent check of machines.py: a plain tick-by-tick simulation that moves every wheel one tick at a time (no closed forms, no periods), as the page does."""
import math, sys
T=100_000
def isprime(n):
    if n<2: return False
    d=2
    while d*d<=n:
        if n%d==0: return False
        d+=1
    return True
PR=[isprime(i) for i in range(T+2)]
EX={2,3,5,7,13,17,19,31,61,89,107,127}
def run(rule,N):
    ks=list(range(2,N+1)); pos={k:(1%k if rule=='double' else 0) for k in ks}; prev={k:(1%k) for k in ks}  # fib state (a,b)
    fa={k:0 for k in ks}; fb={k:1%k for k in ks}
    pi=0; direction=1; first={}; quiet=qp=npr=0; recp=(0,0,[]); rece=(0,0,[]); moved=False
    for t in range(1,T+1):
        home=[]
        for k in ks:
            if rule=='div': pos[k]=(pos[k]+1)%k; h=pos[k]==0
            elif rule=='double': pos[k]=(2*pos[k])%k; h=pos[k]==1%k
            elif rule=='fib': fa[k],fb[k]=fb[k],(fa[k]+fb[k])%k; h=fa[k]==0
            elif rule=='pclock':
                if PR[t]: pos[k]=(pos[k]+1)%k
                h=pos[k]==0 and (pi>0 or PR[t])
            elif rule=='accel': pos[k]=(pos[k]+t)%k; h=pos[k]==0
            elif rule=='square': pos[k]=(pos[k]+2*t-1)%k; h=pos[k]==0
            elif rule=='flip': pos[k]=(pos[k]+direction)%k; h=pos[k]==0
            elif rule=='slip':
                if math.gcd(t,k)==1: pos[k]=(pos[k]+1)%k
                h=pos[k]==0
            if h:
                home.append(k)
                if PR[k] and k not in first: first[k]=t
        if PR[t]:
            pi+=1; direction=-direction if rule=='flip' else direction
        if t>=2:
            if PR[t]: npr+=1
            if not home:
                quiet+=1; qp+=PR[t]
            if home and all(PR[k] for k in home) and len(home)>recp[0]: recp=(len(home),t,home)
            if home and all(k in EX for k in home) and len(home)>rece[0]: rece=(len(home),t,home)
    return quiet,qp,npr,recp,rece,first
jobs=[(r,30) for r in ('div','double','fib','pclock','accel','square','flip','slip')]+[(r,250) for r in ('div','double','fib','accel')]
for rule,N in jobs:
    q,qp,npr,rp,re,first=run(rule,N)
    print("%-7s N=%3d: quiet %.2f%% of ticks, %.1f%% of them prime, holding %.1f%% of the prime ticks | only primes home: %d %s at tick %d | only exponents home: %d %s at tick %d | first home: %s"%(rule,N,100*q/(T-1),100*qp/max(1,q),100*qp/npr,rp[0],rp[2],rp[1],re[0],re[2],re[1],", ".join("%d (%d)"%(p,first.get(p,0)) for p in (2,3,5,7,11,13,17,19,23,29))),flush=True)
