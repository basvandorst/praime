import numpy as np, math, random
N=20_000_000
def sieve(n):
    s=np.ones(n+1,dtype=bool); s[:2]=False
    for i in range(2,int(n**0.5)+1):
        if s[i]: s[i*i::i]=False
    return s
s=sieve(N); PR=[int(p) for p in np.nonzero(s[:20000])[0]]
def model(pick,z0=31):
    out=np.ones(N+1,dtype=bool); out[:2]=False; P=1.0
    for q in PR:
        nq=math.exp(1/P) if P<1 else 2.0
        if nq>N: break
        P*=1-1/q; lo=int(math.ceil(nq))
        if q<=z0: st=max(2*q,((lo+q-1)//q)*q)
        else: a=pick(q); st=lo+((a-lo)%q)
        out[st::q]=False
    return out
def share(L):
    od=np.nonzero(L)[0]; od=od[od>3]; d=np.cumsum(np.where(od%4==3,1,-1)); return round(float((d>0).mean()),3),int(d[-1])
for seed in (101,102,103,104,105,106):
    r=random.Random(seed); print("python random seed",seed,share(model(r.randrange)))
rng=np.random.default_rng(77); v=sorted(share(model(lambda q:int(rng.integers(q))))[0] for _ in range(24))
print("numpy seed 77, 24 draws sorted:",v)
