import numpy as np, math, sys, time
N=int(sys.argv[1]) if len(sys.argv)>1 else 10_000_000
NR=int(sys.argv[2]) if len(sys.argv)>2 else 4
SEED=int(sys.argv[3]) if len(sys.argv)>3 else 1
def sieve(n):
    s=np.ones(n+1,dtype=bool); s[:2]=False
    for i in range(2,int(n**0.5)+1):
        if s[i]: s[i*i::i]=False
    return s
s=sieve(N); PR=np.nonzero(s)[0]
Z1=int(N**0.4027)
def m1(rng):
    sv=np.ones(N+1,dtype=bool); sv[:2]=False
    for q in PR[PR<=Z1]: sv[2*q::q]=False
    B=4096; nb=(N+1+B-1)//B
    pc=np.add.reduceat(s.astype(np.int64),np.arange(0,N+1,B)); sc=np.add.reduceat(sv.astype(np.int64),np.arange(0,N+1,B))
    rate=np.repeat(pc/np.maximum(sc,1),B)[:N+1]
    return sv&(rng.random(N+1)<rate)
def m3(rng,z0):
    out=np.ones(N+1,dtype=bool); out[:2]=False; P=1.0; zmax=0
    for q in PR:
        q=int(q); nq=math.exp(1/P) if P<1 else 2.0
        if nq>N: break
        P*=1-1/q; zmax=q; lo=int(math.ceil(nq))
        if q<=z0: st=max(2*q,((lo+q-1)//q)*q)
        else:
            a=int(rng.integers(q)); st=lo+((a-lo)%q)
        out[st::q]=False
    return out,zmax
A0=N//2
def tests(L,name):
    r={}
    r['count']=int(L.sum()); r['countA']=int(L[A0:].sum())
    for h in (2,4,6,30,210): r['pair%d'%h]=int((L[A0:N-h]&L[A0+h:N]).sum())
    r['t026']=int((L[A0:N-6]&L[A0+2:N-4]&L[A0+6:N]).sum()); r['t046']=int((L[A0:N-6]&L[A0+4:N-2]&L[A0+6:N]).sum())
    r['q0268']=int((L[A0:N-8]&L[A0+2:N-6]&L[A0+6:N-2]&L[A0+8:N]).sum())
    r['ap6']=int((L[A0:N-18]&L[A0+6:N-12]&L[A0+12:N-6]&L[A0+18:N]).sum())
    n=np.arange(10**6,N-2*math.isqrt(N)-20,dtype=np.int64); rr=np.floor(np.sqrt(n)).astype(np.int64); rr-=(rr*rr>n); rr+=((rr+1)**2<=n)
    a=int(L[n].sum()); v=[]
    for k in range(-3,4):
        m=n+2*rr+2*k; v.append(int((L[n]&L[m]).sum())*len(n)/(a*int(L[m].sum())))
    mu=sum(v)/7; r['sq_rms']=math.sqrt(sum((100*(x/mu-1))**2 for x in v)/7); r['sq_k0']=100*(v[3]/mu-1)
    c=np.concatenate([[0],np.cumsum(L[A0:N].astype(np.int64))]); M=N-A0
    for e in range(4,19):
        H=1<<e; st=max(1,H//2); i=np.arange(0,M-H+1,st); cnt=(c[i+H]-c[i]).astype(float); x=(i-i.mean())/max(1.0,i.std())
        co=np.polyfit(x,cnt,2); res=cnt-np.polyval(co,x); r['var%d'%e]=float((res**2).sum()/(len(cnt)-3)/cnt.mean())
    idx=np.nonzero(L)[0]; g=np.diff(idx); r['maxgap']=int(g.max()); r['gap100']=int((g>=100).sum())
    od=idx[idx>3]
    d=np.cumsum(np.where(od%4==3,1,-1)*(od%2==1)); r['race4_end']=int(d[-1]); r['race4_pos']=float((d>0).mean()); r['race4_max']=float(np.abs(d).max()/math.sqrt(len(d)))
    m3_=od%3; r['same3']=float((m3_[1:]==m3_[:-1]).mean())
    g2=np.diff(od).astype(float); r['gapcorr']=float(np.corrcoef(g2[1:],g2[:-1])[0,1])
    return r
t0=time.time(); RP=tests(s,'primes'); keys=list(RP.keys())
rows={'primes':[RP]}
rng=np.random.default_rng(SEED)
for name,gen in (('coin',lambda: m1(rng)),('class z0=1',lambda: m3(rng,1)[0]),('class z0=31',lambda: m3(rng,31)[0]),('class z0=%d'%Z1,lambda: m3(rng,Z1)[0])):
    rows[name]=[tests(gen(),name) for _ in range(NR)]
    print('#',name,'done',round(time.time()-t0),file=sys.stderr,flush=True)
print("N = %d, statistics on [N/2, N] unless stated, %d realisations per model, seed %d, coin sieve to %d, random classes up to %d"%(N,NR,SEED,Z1,m3(np.random.default_rng(0),1)[1]))
x=0.75*N
print("%-12s %14s"%("statistic","primes")+"".join(" %22s"%k for k in list(rows)[1:])+"   (mean +- sd over realisations)")
for k in keys:
    line="%-12s %14.6g"%(k,RP[k])
    for name in list(rows)[1:]:
        v=np.array([r[k] for r in rows[name]],float); line+=" %12.6g +- %-7.3g"%(v.mean(),v.std(ddof=1) if len(v)>1 else 0)
    if k.startswith('var'):
        H=1<<int(k[3:]); line+="   MS %.4f"%((math.log(x/H)+1-0.5772156649-math.log(2*math.pi))/math.log(x))
    print(line)
