import math, numpy as np, random, sys
from sympy import isprime, nextprime, prevprime
def sieve(n):
    s=np.ones(n+1,dtype=bool); s[:2]=False
    for i in range(2,int(n**0.5)+1):
        if s[i]: s[i*i::i]=False
    return s
print("E. HOW FAR A MACHINE CAN BE TRUSTED: share of quiet ticks that are real finds, at t near N^u")
SM=sieve(1000); small=[p for p in range(2,1000) if SM[p]]
def rough_window(N,lo,hi):
    r=np.ones(hi-lo,dtype=bool)
    for p in small:
        if p>N: break
        st=(-lo)%p; r[st::p]=False
    return r
def primes_window(lo,hi):
    r=np.ones(hi-lo,dtype=bool)
    for p in range(2,int(hi**0.5)+1):
        if p<1000 and not SM[p]: continue
        st=max(p*p,((lo+p-1)//p)*p)-lo
        r[st::p]=False
    return r
BIG=sieve(6000); bigp=[p for p in range(2,6000) if BIG[p]]
def pw(lo,hi):
    r=np.ones(hi-lo,dtype=bool)
    for p in bigp:
        if p*p>hi: break
        st=max(p*p,((lo+p-1)//p)*p)-lo; r[st::p]=False
    return r
w=lambda u:(1+math.log(u-1))/u
for N in (30,60,100,250):
    row=[]
    for u in (2.25,2.5,2.75,3.0):
        x=int(N**u); lo=int(x*0.9); hi=int(x*1.1)+4
        if hi>3.3e7: row.append("u=%.2f: too far"%u); continue
        r=rough_window(N,lo,hi); p=pw(lo,hi)
        one=(r[:-2]&p[:-2]).sum()/r[:-2].sum(); tw=r[:-2]&r[2:]; twin=(tw&p[:-2]&p[2:]).sum()/max(1,tw.sum())
        row.append("u=%.2f: primes %.3f (Buchstab %.3f), twin pairs %.3f"%(u,one,1/(u*w(u)),twin))
    print("   N=%3d: "%N+"; ".join(row),flush=True)
print("G. BACK AND FORTH: the distance D(t) = prime gaps added and subtracted in turn")
T=2_000_000; s=sieve(T); pr=np.nonzero(s)[0]; g=np.diff(pr); sg=np.where(np.arange(len(g))%2==0,1,-1)
A=np.cumsum(g*sg*(-1)); 
alt=np.cumsum(g*sg)
print("   primes to %d: %d gaps; alternating sum of gaps: final %d, largest |value| %d; sqrt(sum of gap^2) = %.0f"%(T,len(g),alt[-1],np.abs(alt).max(),math.sqrt((g.astype(float)**2).sum())))
rnd=np.random.default_rng(3); mx=[];fin=[]
for _ in range(300):
    gg=rnd.permutation(g); a=np.cumsum(gg*sg); mx.append(np.abs(a).max()); fin.append(a[-1])
print("   same gaps in random order, 300 times: largest |value| mean %.0f (sd %.0f), final sd %.0f; real largest is at percentile %.1f"%(np.mean(mx),np.std(mx),np.std(fin),100*np.mean(np.array(mx)<np.abs(alt).max())))
for X in (10**4,10**5,10**6,2*10**6):
    n=np.searchsorted(pr,X)-1; print("   up to %d: largest |alternating sum| %d = %.2f x sqrt(number of primes)"%(X,np.abs(alt[:n]).max(),np.abs(alt[:n]).max()/math.sqrt(n)))
