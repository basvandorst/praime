import numpy as np, math, sys
N=2_000_000
def sieve(n):
    s=np.ones(n+1,dtype=bool); s[:2]=False
    for i in range(2,int(n**0.5)+1):
        if s[i]: s[i*i::i]=False
    return s
s=sieve(N); PR=np.nonzero(s)[0]; Z1=int(N**0.4027)
def m1(rng):
    sv=np.ones(N+1,dtype=bool); sv[:2]=False
    for q in PR[PR<=Z1]: sv[2*q::q]=False
    B=4096; pc=np.add.reduceat(s.astype(np.int64),np.arange(0,N+1,B)); sc=np.add.reduceat(sv.astype(np.int64),np.arange(0,N+1,B))
    return sv&(rng.random(N+1)<np.repeat(pc/np.maximum(sc,1),B)[:N+1])
def m3(rng,z0):
    out=np.ones(N+1,dtype=bool); out[:2]=False; P=1.0; CL={}
    for q in PR:
        q=int(q); nq=math.exp(1/P) if P<1 else 2.0
        if nq>N: break
        P*=1-1/q; lo=int(math.ceil(nq))
        if q<=z0: st=max(2*q,((lo+q-1)//q)*q)
        else: a=int(rng.integers(q)); st=lo+((a-lo)%q); CL[q]=a
        out[st::q]=False
    return out,CL
M=1<<22
def conv(L):
    f=np.fft.rfft(L.astype(float),M); return np.rint(np.fft.irfft(f*f,M)[:N+1]).astype(np.int64)
w=np.zeros(N+1); w[2:]=1/np.log(np.arange(2,N+1)); f=np.fft.rfft(w,M); I=np.fft.irfft(f*f,M)[:N+1]
ns=np.arange(N//2,N+1,2)
def S(zmax):
    out=np.ones(N+1); 
    for q in PR[(PR>2)&(PR<=zmax)]:
        q=int(q); out*=1-1/(q-1)**2; out[::q]*=(q-1)/(q-2)/(1-1/(q-1)**2)*(1-1/(q-1)**2)
    return 2*out
def S2(zmax):
    out=np.full(N+1,2.0)
    for q in PR[(PR>2)&(PR<=zmax)]:
        q=int(q); out*=(1-1/(q-1)**2); out[::q]*=(q-1)/(q-2)
    return out
Sfull=S2(int(N**0.5)+1); S31=S2(31)
def stat(L,Sx):
    r=conv(L)[ns]; E=Sx[ns]*I[ns]; t=r/E; z=(r-E*t.mean())/np.sqrt(E*t.mean()); return t.mean(),float(z.var()),int(r.min())
print("even n from %d to %d; r(n) = ordered pairs of members with sum n; E(n) = singular series x sum 1/(ln m ln(n-m))"%(N//2,N))
a=stat(s,Sfull); print("primes:          mean r/E %.4f, variance of (r - E)/sqrt(E) %.3f"%a[:2])
rng=np.random.default_rng(5)
def Sc(z0,CL):
    out=np.full(N+1,2.0)
    for q in PR[(PR>2)&(PR<=z0)]:
        q=int(q); out*=(1-1/(q-1)**2); out[::q]*=(q-1)/(q-2)
    for q,a in CL.items():
        out*=(1-1/(q-1)**2); out[(2*a)%q::q]*=(q-1)/(q-2)
    return out
v=[stat(m1(rng),Sfull) for _ in range(4)]
print("coin:            mean r/E %.4f, variance %.3f (values %s)"%(np.mean([x[0] for x in v]),np.mean([x[1] for x in v])," ".join("%.3f"%x[1] for x in v)))
for z0 in (31,Z1):
    v=[]
    for _ in range(4):
        L,CL=m3(rng,z0); v.append(stat(L,Sc(z0,CL)))
    print("class z0=%-4d   mean r/E %.4f, variance %.3f (values %s)"%(z0,np.mean([x[0] for x in v]),np.mean([x[1] for x in v])," ".join("%.3f"%x[1] for x in v)))
print("a Poisson count of unordered pairs gives variance 2 in these units")
