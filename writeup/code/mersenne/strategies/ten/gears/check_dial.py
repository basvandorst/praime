"""Independent check of the dial (machine 7): plain loops."""
T=100_000
def isprime(n):
    if n<2: return False
    d=2
    while d*d<=n:
        if n%d==0: return False
        d+=1
    return True
PR=[isprime(i) for i in range(T+1)]
for N in (30,100,250):
    q=qp=npr=tot=0; first={}
    for t in range(1,T+1):
        home=[k for k in range(2,N+1) if (k*t)%120==0]
        for k in home:
            if PR[k] and k not in first: first[k]=t
        if t<2 or t%120==0: continue
        tot+=1; npr+=PR[t]
        if not home: q+=1; qp+=PR[t]
    print("dial N=%3d: quiet %.2f%% of ticks, %.1f%% of them prime, holding %.1f%% of the prime ticks; first home 2,3,5,7: %s"%(N,100*q/tot,100*qp/max(1,q),100*qp/npr,[first[p] for p in (2,3,5,7)]))
