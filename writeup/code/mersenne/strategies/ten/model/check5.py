"""Independent recomputation (plain Python sieve in a bytearray, other code paths, other seeds) of the numbers quoted from ladder.py, goldbach.py and mers_law.py."""
import math, random, sys, gmpy2
N=20_000_000
def psieve(n):
    b=bytearray([1])*(n+1); b[0]=b[1]=0
    for i in range(2,int(n**0.5)+1):
        if b[i]: b[i*i::i]=bytearray(len(range(i*i,n+1,i)))
    return b
s=psieve(N); small=[i for i in range(2,20000) if s[i]]
def model(seed,z0=31):
    rnd=random.Random(seed); m=bytearray([1])*(N+1); m[0]=m[1]=0; P=1.0; zm=0
    for q in small:
        nq=math.exp(1/P) if P<1 else 2.0
        if nq>N: break
        P*=1-1/q; zm=q; lo=math.ceil(nq)
        if q<=z0: st=max(2*q,-(-lo//q)*q)
        else:
            a=rnd.randrange(q); st=lo+((a-lo)%q)
        m[st::q]=bytearray(len(range(st,N+1,q)))
    return m,zm
def li_inc(a,b,steps=64):
    h=(b-a)/steps; return h*(0.5/math.log(a)+0.5/math.log(b)+sum(1/math.log(a+i*h) for i in range(1,steps)))
def stats(L,name):
    A=N//2; idx=[i for i in range(2,N+1) if L[i]]
    cnt=len(idx); p2=sum(1 for i in range(A,N-2) if L[i] and L[i+2]); p6=sum(1 for i in range(A,N-6) if L[i] and L[i+6])
    q4=sum(1 for i in range(A,N-8) if L[i] and L[i+2] and L[i+6] and L[i+8])
    mg=0; g100=0; same=0; tot=0; d=0; pos=0; steps=0; prev=None
    for i in idx:
        if prev is not None:
            g=i-prev
            if g>mg: mg=g
            if g>=100: g100+=1
            if prev>3:
                tot+=1
                if i%3==prev%3: same+=1
        if i>3 and i%2==1:
            d+=1 if i%4==3 else -1
        if i>3:
            steps+=1
            if d>0: pos+=1
        prev=i
    H=8192; c=[0]*(N-A+1); t=0
    for i in range(A,N):
        t+=L[i]; c[i-A+1]=t
    res=[]; mean=0
    for j in range(0,N-A-H+1,H//2):
        k=c[j+H]-c[j]; res.append(k-li_inc(A+j,A+j+H)); mean+=k
    mean/=len(res); mu=sum(res)/len(res); var=sum((x-mu)**2 for x in res)/(len(res)-1)/mean
    v=[]; lo=10**6; hi=N-2*math.isqrt(N)-20; a=sum(1 for n in range(lo,hi) if L[n])
    for k in range(-3,4):
        both=0; b=0
        for n in range(lo,hi):
            m=n+2*math.isqrt(n)+2*k
            if L[m]:
                b+=1
                if L[n]: both+=1
        v.append(both*(hi-lo)/(a*b))
    mu2=sum(v)/7; rms=math.sqrt(sum((100*(x/mu2-1))**2 for x in v)/7)
    print("%-14s count %d, pairs h=2: %d, h=6: %d, quadruples: %d, largest gap %d, gaps >= 100: %d, same remainder mod 3: %.4f, 3 mod 4 ahead: %.4f, end %d, variance ratio H=8192 (li removed): %.3f, square rms %.2f%%"%(name,cnt,p2,p6,q4,mg,g100,same/tot,pos/steps,d,var,rms),flush=True)
stats(s,"primes")
for seed in (101,102):
    m,zm=model(seed); stats(m,"class seed %d"%seed)
print("random classes up to",zm)
