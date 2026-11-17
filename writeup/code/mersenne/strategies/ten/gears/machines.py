import numpy as np, math, json, random
T=100_000; NS=(30,100,250); NMAX=250
t=np.arange(0,T+1,dtype=np.int64)
def sieve(n):
    s=np.ones(n+1,dtype=bool); s[:2]=False
    for i in range(2,int(n**0.5)+1):
        if s[i]: s[i*i::i]=False
    return s
isp=sieve(T); PI=np.cumsum(isp).astype(np.int64)
D=np.zeros(T+1,dtype=np.int64); sign=np.where(PI[:-1]%2==1,-1,1); D[1:]=np.cumsum(sign)
PRIMES=[k for k in range(2,NMAX+1) if isp[k]]; EXPS=[2,3,5,7,13,17,19,31,61,89,107,127]
def periodic(seq_start,step,k):
    seen={}; seq=[]; x=seq_start; i=0
    while x not in seen:
        seen[x]=i; seq.append(x); x=step(x); i+=1
    mu=seen[x]; lam=i-mu; idx=np.where(t<mu,t,mu+(t-mu)%lam); return np.array([s if isinstance(s,int) else s[0] for s in seq],dtype=np.int64)[idx]
def home(rule,k):
    if rule=='div': return t%k==0
    if rule=='double': return periodic(1%k,lambda x:(2*x)%k,k)==(1%k)
    if rule=='fib': return periodic((0,1%k),lambda p:(p[1],(p[0]+p[1])%k),k)==0
    if rule=='pclock': return PI%k==0
    if rule=='accel': return (t*(t+1)//2)%k==0
    if rule=='square': return (t*t)%k==0
    if rule=='flip': return D%k==0
    if rule=='slip':
        c=np.cumsum([0]+[1 if math.gcd(j,k)==1 else 0 for j in range(1,k+1)]); return ((t//k)*int(c[k])+c[t%k])%k==0
    if rule=='fast': return (k*t)%120==0
RULES=[('div','Wheels, one tooth per tick (machines 1 to 6, 8 to 10, 18)'),('fast','One dial, wheel k turns k times as fast (7)'),('double','Doubling wheels (11)'),('fib','Fibonacci wheels (12)'),('pclock','Prime clock (13)'),('accel','Accelerating wheels (14)'),('square','Square jumps (15)'),('flip','Back and forth (16)'),('slip','Slipping wheels (17)')]
out={}
for rule,name in RULES:
    H=np.zeros((NMAX+1,T+1),dtype=bool)
    for k in range(1,NMAX+1): H[k]=home(rule,k)
    H[:,0]=False
    if rule=='pclock': H[:,PI==0]=False
    first={p:(int(np.argmax(H[p])) if H[p].any() else 0) for p in PRIMES}
    res={'name':name,'first':first,'N':{}}
    for N in NS:
        ks=np.arange(2,N+1); pk=[k for k in ks if isp[k]]; ck=[k for k in ks if not isp[k]]; ek=[k for k in EXPS if k<=N]; nek=[k for k in ks if k not in ek]
        ph=H[pk].sum(0); chh=H[ck].sum(0); eh=H[ek].sum(0); neh=H[nek].sum(0); allh=ph+chh
        rng=np.arange(T+1)>=2
        if rule=='fast': rng&=(t%120!=0)
        quiet=(allh==0)&rng; prime_t=isp&rng
        pure=np.where((chh==0)&rng,ph,0); pe=np.where((neh==0)&rng,eh,0)
        r={'mean_home':float(allh[rng].mean()),'quiet_share':float(quiet.sum()/rng.sum()),'precision':float((quiet&prime_t).sum()/max(1,quiet.sum())),'recall':float((quiet&prime_t).sum()/prime_t.sum()),'base':float(prime_t.sum()/rng.sum()),
           'rec_p':[int(pure.max()),int(np.argmax(pure))],'rec_e':[int(pe.max()),int(np.argmax(pe))],'pure_share':float((pure>0).sum()/rng.sum()),'never':sum(1 for p in pk if first[p]==0)}
        r['rec_p_list']=[int(k) for k in pk if H[k][r['rec_p'][1]]] if r['rec_p'][0] else []
        r['rec_e_list']=[int(k) for k in ek if H[k][r['rec_e'][1]]] if r['rec_e'][0] else []
        if rule=='accel' and N==250:
            q=np.nonzero(quiet)[0]; q=q[q<N*N]; sg=[int(x) for x in q if (x%2==0 and isp[x//2] and isp[x+1]) or (x%2==1 and isp[x] and isp[(x+1)//2])]
            r['accel_pairs']=[len(q),len(sg),[int(x) for x in q[:12]]]
        res['N'][N]=r
    rat_e=[first[p]/p for p in EXPS if first[p]]; rat_o=[first[p]/p for p in PRIMES if p not in EXPS and first[p]]
    obs=np.median(rat_e)-np.median(rat_o); allr=rat_e+rat_o; rnd=random.Random(1); cnt=0
    for _ in range(20000):
        rnd.shuffle(allr)
        if abs(np.median(allr[:len(rat_e)])-np.median(allr[len(rat_e):]))>=abs(obs)-1e-12: cnt+=1
    res['exp_med']=float(np.median(rat_e)); res['oth_med']=float(np.median(rat_o)); res['perm_p']=cnt/20000
    out[rule]=res
    print("\n== %s =="%name)
    print("   first tick home for the prime gears p <= 47: "+", ".join("%d (%d)"%(p,first[p]) for p in PRIMES if p<=47))
    print("   first tick / p: median for the 12 Mersenne-exponent gears %.3f, for the other 41 prime gears %.3f, permutation p = %.3f; prime gears never home in %d ticks: %d"%(res['exp_med'],res['oth_med'],res['perm_p'],T,res['N'][250]['never']))
    for N in NS:
        r=res['N'][N]; print("   N=%3d: gears home per tick %.2f | quiet ticks (no gear but 1 home) %.2f%% of ticks; of them %.1f%% are prime (all ticks: %.1f%%); they hold %.1f%% of the prime ticks | only primes home: record %d %s at tick %d (such ticks: %.2f%%) | only exponents home: record %d %s at tick %d"%(N,r['mean_home'],100*r['quiet_share'],100*r['precision'],100*r['base'],100*r['recall'],r['rec_p'][0],r['rec_p_list'],r['rec_p'][1],100*r['pure_share'],r['rec_e'][0],r['rec_e_list'],r['rec_e'][1]))
    if 'accel_pairs' in res['N'][250]: a=res['N'][250]['accel_pairs']; print("   accelerating, N=250, quiet ticks below 250^2: %d; of the form t = 2p with p and 2p+1 prime, or t prime with (t+1)/2 prime: %d; first: %s"%(a[0],a[1],a[2]))
json.dump(out,open('machines.json','w'))
print("\nLucas-Lehmer wheels: p <= 127 with tooth 0 at tick p-2:",[p for p in range(3,128) if isp[p] and (lambda M:(__import__('functools').reduce(lambda s,_:(s*s-2)%M,range(p-2),4)==0))(2**p-1)])
print("doubling machine as a Mersenne sieve: prime ticks p whose 2^p-1 has a factor up to N")
for N in NS:
    el=[p for p in range(2,2000) if isp[p] and any(pow(2,p,q)==1 for q in range(3,N+1,2) if q!=2**p-1)]
    print("   N=%3d: eliminated %s; a wheel can only help for p <= %d"%(N,el,(N-1)//2))
