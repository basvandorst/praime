import numpy as np, math
N=20_000_000
def sieve(n):
    s=np.ones(n+1,dtype=bool); s[:2]=False
    for i in range(2,int(n**0.5)+1):
        if s[i]: s[i*i::i]=False
    return s
s=sieve(N); PR=np.nonzero(s)[0]
def m3(rng,z0):
    out=np.ones(N+1,dtype=bool); out[:2]=False; P=1.0
    for q in PR:
        q=int(q); nq=math.exp(1/P) if P<1 else 2.0
        if nq>N: break
        P*=1-1/q; lo=int(math.ceil(nq))
        if q<=z0: st=max(2*q,((lo+q-1)//q)*q)
        else: a=int(rng.integers(q)); st=lo+((a-lo)%q)
        out[st::q]=False
    return out
def share(L):
    od=np.nonzero(L)[0]; od=od[od>3]; d=np.cumsum(np.where(od%4==3,1,-1)); return float((d>0).mean()),int(d[-1])
print("primes: remainder 3 ahead at %.4f of the primes up to %d, final lead %d"%(*share(s),)[:0] if False else "primes: share %.4f, final lead %d"%share(s))
rng=np.random.default_rng(77); v=[share(m3(rng,31)) for _ in range(40)]
sh=np.array([x[0] for x in v]); en=np.array([x[1] for x in v])
print("class model, 40 draws: mean share %.3f, sd %.3f (arcsine law: 0.5 and 0.354); draws above 0.99: %d, below 0.01: %d, between 0.25 and 0.75: %d; final lead mean %.0f, sd %.0f (sqrt of the count: %.0f)"%(sh.mean(),sh.std(ddof=1),(sh>0.99).sum(),(sh<0.01).sum(),((sh>0.25)&(sh<0.75)).sum(),en.mean(),en.std(ddof=1),math.sqrt(s.sum())))
