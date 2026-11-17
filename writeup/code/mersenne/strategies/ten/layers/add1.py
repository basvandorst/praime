import numpy as np, math, sys
N=1_000_000
def sieve(n):
    s=np.ones(n+1,dtype=bool); s[:2]=False
    for i in range(2,int(n**0.5)+1):
        if s[i]: s[i*i::i]=False
    return s
s=sieve(N); odd=s.copy(); odd[2]=False
M=1<<22
def conv(a,b): return np.rint(np.fft.irfft(np.fft.rfft(a.astype(float),M)*np.fft.rfft(b.astype(float),M),M)).astype(np.int64)
def lpf(n):
    for p in (3,5,7,11,13,17,19,23,29,31,37,41,43,47):
        if n%p==0: return p
    return 0
def score(c,name,low=True,lo=N//2,hi=N):
    n=np.arange(lo|1,hi,2); v=c[n].astype(float)
    B=2048; out=np.zeros(len(n),bool); nb=len(n)//B
    hit=0; tot=0
    for b in range(nb):
        sl=slice(b*B,(b+1)*B); k=int(s[n[sl]].sum()); o=np.argsort(v[sl]); o=o[:k] if low else o[::-1][:k]
        hit+=int(s[n[sl]][o].sum()); tot+=k
    base=s[n].mean()
    print("%-34s odd n in [%d, %d): the %s counts, as many as there are primes: %.1f%% are prime (share of primes among odd numbers %.1f%%)"%(name,lo,hi,"lowest" if low else "highest",100*hit/tot,100*base))
    for p in (3,5,7,11,13,31,47):
        sel=(n%p==0)
        for q in (3,5,7,11,13,17,19,23,29,31,37,41,43):
            if q<p: sel&=(n%q!=0)
        rr=[]
        for bb in range(nb):
            sl=slice(bb*B,(bb+1)*B); pm=s[n[sl]]; sm=sel[sl]
            if pm.sum()>=5 and sm.sum()>=5: rr.append(np.median(v[sl][sm])/np.median(v[sl][pm]))
        print("    smallest factor %2d: median count / median for primes, inside blocks of %d odd numbers, mean over %d blocks = %.3f"%(p,B,len(rr),np.mean(rr)),end="")
        print("   (rule (p-1)/(p-2) = %.3f)"%((p-1)/(p-2)) if low else "")
    return
two=np.zeros(N+1,bool); q=np.nonzero(odd[:N//2+1])[0]; two[2*q]=True
c1=conv(odd,two)[:N+1]; score(c1,"n = p + 2q")
A=np.fft.rfft(odd.astype(float),M); t=np.rint(np.fft.irfft(A*A*A,M)[:N+1]); score(t,"n = p + q + r",low=False)
X=N; f=np.fft.rfft(s.astype(float),M); ac=np.rint(np.fft.irfft(f*np.conj(f),M)).astype(np.int64)
g=np.zeros(N+1); nn=np.arange(1,N//20); g[nn]=ac[2*nn]
score(g,"pairs of primes 2n apart",lo=3,hi=N//20-1)
