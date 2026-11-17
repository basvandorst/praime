"""Independent check: quiet ticks of the accelerating wheels with 250 gears below 250^2, and the doubling sieve."""
def isprime(n):
    if n<2: return False
    d=2
    while d*d<=n:
        if n%d==0: return False
        d+=1
    return True
N=250; q=[]; s=0
for t in range(1,N*N):
    s+=t
    if t>=2 and all(s%k for k in range(2,N+1)): q.append(t)
form=[t for t in q if (t%2==0 and isprime(t//2) and isprime(t+1)) or (t%2==1 and isprime(t) and isprime((t+1)//2))]
print("accelerating, 250 gears: quiet ticks below 62,500: %d; of the two forms: %d; even ones (t = 2p, p and 2p+1 prime): %d; odd ones: %d"%(len(q),len(form),sum(1 for t in form if t%2==0),sum(1 for t in form if t%2)))
pr=[p for p in range(2,125) if isprime(p)]
el={p:[w for w in range(3,251,2) if pow(2,p,w)==1 and w!=2**p-1] for p in pr}
print("prime ticks up to 124: %d; ruled out by a wheel up to 250: %s; Mersenne exponents among them: %d"%(len(pr),{p:v for p,v in el.items() if v},sum(1 for p in pr if p in (2,3,5,7,13,17,19,31,61,89,107))))
