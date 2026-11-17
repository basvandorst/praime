from sympy import primerange, isprime, prime
ME={2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423}
T=5000; pr=list(primerange(2,T+1))
def out_plain(wheels):
    w=[q for q in wheels if q%2 and q>2]; r=[]
    for p in pr:
        if any(pow(2,p,q)==1 and q!=2**p-1 for q in w if (q-1)%(2*p)==0): r.append(p)
    return r
def out_fit(K):
    r=[]; tests=0; skipped=0
    for p in pr:
        hit=False
        for k in range(1,K+1):
            q=2*k*p+1
            if q%8 not in (1,7) or any(q%s==0 and q!=s for s in (3,5,7,11,13,17,19,23,29,31)): skipped+=1; continue
            tests+=1
            if pow(2,p,q)==1 and q!=2**p-1: hit=True
        if hit: r.append(p)
    return r,tests,skipped
print("prime ticks up to %d: %d, Mersenne exponents among them: %d"%(T,len(pr),len([p for p in pr if p in ME])))
for B in (250,5000):
    a=out_plain(range(1,B+1)); b=out_plain([prime(i) for i in range(1,B+1)]); c,t,s=out_fit(B)
    print("budget %4d wheels: numbers 1..%d rule out %d prime ticks | the first %d primes (up to %d) rule out %d | wheels that fit (2kp+1, k <= %d) rule out %d (%.1f%%), using %d of %d wheels (%.1f%% skipped as impossible); candidates left %d"%(B,B,len(a),B,prime(B),len(b),B,len(c),100*len(c)/len(pr),t,t+s,100*s/(t+s),len(pr)-len(c)))
    assert not (set(c)&ME)
for B in (30,250,5000):
    g=prime(B); print("budget %4d: pattern machine exact below %d^2 = %d with numbers 1..N, below %d^2 = %d with the first N primes (%.0f times further)"%(B,B,B*B,g,g*g,g*g/(B*B)))
