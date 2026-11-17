import numpy as np, math, sys
from lib import sieve
N=20_000_000; s=sieve(N); y0=int(round(N**0.4027)); sv=np.ones(N+1,dtype=bool); sv[:2]=False
for q in np.nonzero(s[:y0+1])[0]: sv[2*q::q]=False
BL=12; nb=(N>>BL)+1; idx=np.arange(N+1)>>BL; bp=np.bincount(idx[y0+1:],weights=s[y0+1:],minlength=nb); bs=np.bincount(idx[y0+1:],weights=sv[y0+1:],minlength=nb); pb=np.where(bs>0,bp/np.maximum(bs,1),0)
rng=np.random.default_rng(21); ND=12; D=[]
for i in range(ND):
    d=(sv&(rng.random(N+1)<pb[idx])); d[:y0+1]=s[:y0+1]; D.append(d)
n=np.arange(10**6+1,N-200000,2,dtype=np.int64); fl=lambda x: np.floor(x).astype(np.int64)
def rev10(a):
    r=np.zeros_like(a); t=a.copy()
    while t.any(): r=np.where(t>0,r*10+t%10,r); t//=10
    return r
def rev2(a):
    r=np.zeros_like(a); t=a.copy()
    while t.any(): r=np.where(t>0,r*2+(t&1),r); t>>=1
    return r
dsum=lambda a: sum((a//10**i)%10 for i in range(8))
R={"n + 2":n+2,"n + 6":n+6,"n + 30":n+30,"n + 2310":n+2310,"2n + 1":2*n+1,"2n - 1":2*n-1,"4n + 1":4*n+1,"6n + 1":6*n+1,"n + 2 floor(sqrt n)":n+2*fl(np.sqrt(n)),"n + 2 floor(sqrt n) + 2":n+2*fl(np.sqrt(n))+2,"n + 2 floor(n^(1/3))":n+2*fl(np.cbrt(n)),"n + 2 floor(n^(1/4))":n+2*fl(n**0.25),"n + 2 floor(n^(2/3))":n+2*fl(np.cbrt(n)**2),"n + 2 floor(ln^2 n)":n+2*fl(np.log(n)**2),"digits reversed (base 10)":rev10(n),"bits reversed":rev2(n),"n xor 2":n^2,"n xor 4":n^4,"n xor 1024":n^1024,"n xor 2^20":n^(1<<20),"n + 2 (n mod 1000)":n+2*(n%1000),"n + 2 (n div 1000 mod 1000)":n+2*((n//1000)%1000),"2 floor(n phi / 2) + 1":2*fl(n*0.8090169943749475)+1,"N - n":N-n,"n + 2 floor(sqrt(n) phi)":n+2*fl(np.sqrt(n)*1.618033988749895)}
print("relation n -> m, n odd from 1e6 to 2e7: prime pairs, ratio for the primes, for the decoys (sieve to %d, then a coin), z"%y0)
for name,m in R.items():
    ok=(m>1)&(m<=N)&(m!=n); a=n[ok]; b=m[ok]; np_=len(a)
    def Rr(L):
        u=L[a]; v=L[b]; pp=int((u&v).sum()); return pp,(pp*np_/(int(u.sum())*max(1,int(v.sum()))) if pp else 0)
    pp,Rp=Rr(s); rd=[Rr(d) for d in D]; mean=np.mean([t[0] for t in rd]); Rd=np.mean([t[1] for t in rd]); sdR=np.std([t[1] for t in rd],ddof=1)
    if pp<30 or mean<30: print("  %-30s too few pairs (%d)"%(name,pp)); continue
    print("  %-30s %8d  %.4f  %.4f  z %+5.2f  (by decoy spread %+5.2f)"%(name,pp,Rp,Rd,(Rp-Rd)/(Rd*math.sqrt(1/pp+1/(ND*mean))),(Rp-Rd)/(sdR*math.sqrt(1+1/ND))),flush=True)
