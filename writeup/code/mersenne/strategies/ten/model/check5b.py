"""Independent checks of goldbach.py (direct counts for 1,500 random even n, no FFT) and mers_law.py (every 8th prime, candidates filtered by a primality test instead of remainders)."""
import numpy as np, math, random, gmpy2
N=2_000_000
b=bytearray([1])*(N+1); b[0]=b[1]=0
for i in range(2,1415):
    if b[i]: b[i*i::i]=bytearray(len(range(i*i,N+1,i)))
s=np.frombuffer(bytes(b),dtype=np.uint8).astype(bool); P=np.nonzero(s)[0]; small=[int(p) for p in P[P<5000]]
def model(seed,z0=31):
    rnd=random.Random(seed); m=bytearray([1])*(N+1); m[0]=m[1]=0; Pm=1.0; CL={}
    for q in small:
        nq=math.exp(1/Pm) if Pm<1 else 2.0
        if nq>N: break
        Pm*=1-1/q; lo=math.ceil(nq)
        if q<=z0: st=max(2*q,-(-lo//q)*q)
        else:
            a=rnd.randrange(q); st=lo+((a-lo)%q); CL[q]=a
        m[st::q]=bytearray(len(range(st,N+1,q)))
    return np.frombuffer(bytes(m),dtype=np.uint8).astype(bool),CL
rnd=random.Random(9); ns=sorted(rnd.sample(range(N//2,N+1,2),1500))
lg=np.zeros(N+1); lg[2:]=1/np.log(np.arange(2,N+1))
def run(L,sing,name):
    idx=np.nonzero(L)[0]; z=[]; t=[]
    R=[];E=[]
    for n in ns:
        r=int(L[n-idx[idx<=n-2]].sum()); I=float((lg[2:n-1]*lg[n-2:1:-1]).sum()); R.append(r); E.append(sing(n)*I)
    R=np.array(R,float); E=np.array(E); c=(R/E).mean(); zz=(R-c*E)/np.sqrt(c*E)
    print("%-22s 1,500 random even n: mean r/E %.4f, variance of (r - E)/sqrt(E) %.3f, as a share of a random count %.3f"%(name,c,zz.var(),zz.var()/2),flush=True)
def sing_true(n):
    v=2.0; m=n
    for q in small[1:]:
        if q*q>m and m%q: 
            pass
        v*=1-1/(q-1)**2
        if n%q==0: v*=(q-1)/(q-2)
    return v
run(s,sing_true,"primes")
for seed in (201,202):
    L,CL=model(seed)
    def sing_m(n,CL=CL):
        v=2.0
        for q in small[1:]:
            if q<=31:
                v*=1-1/(q-1)**2
                if n%q==0: v*=(q-1)/(q-2)
            elif q in CL:
                v*=1-1/(q-1)**2
                if (n-2*CL[q])%q==0: v*=(q-1)/(q-2)
        return v
    run(L,sing_m,"class model seed %d"%seed)
PMAX=60000; ps=[int(p) for p in P[(P>=1000)&(P<=PMAX)]][::8]; LB=30; ff={}
for p in ps:
    K=(1<<LB)//(2*p); f=0
    for k in range(1,K+1):
        q=2*k*p+1
        if gmpy2.is_prime(q) and gmpy2.powmod(2,p,q)==1: f=q; break
    ff[p]=f
print("Mersenne: %d primes (every 8th from 1000 to %d)"%(len(ps),PMAX))
for cls,a in ((3,2),(1,6)):
    Q=[p for p in ps if p%4==cls]
    for bb in (20,24,28,30):
        B=1<<bb; QQ=[p for p in Q if a*p<B/4]; meas=sum(1 for p in QQ if ff[p]==0 or ff[p]>=B)/len(QQ); pred=sum(min(1.0,math.log(a*p)/math.log(B)) for p in QQ)/len(QQ)
        print("  p = %d mod 4, B = 2^%d: %4d primes, measured %.4f, predicted %.4f, ratio %.3f"%(cls,bb,len(QQ),meas,pred,meas/pred))
import json; json.dump({str(p):int(f) for p,f in ff.items()},open('check5b_ff.json','w'))
