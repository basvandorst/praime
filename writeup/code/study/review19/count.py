# own recount of twin starts p in [L0,L1) (p, p+2 prime) and gaps between consecutive starts
import sys, math, json, numpy as np, time
L0=int(float(sys.argv[1])); L1=int(float(sys.argv[2])); SEG=int(float(sys.argv[3])); out=sys.argv[4]
lim=int(math.isqrt(L1+2))+1
s=np.ones(lim+1,bool); s[:2]=False
for i in range(2,int(lim**0.5)+1):
    if s[i]: s[i*i::i]=False
ps=np.nonzero(s)[0][1:].astype(np.int64)  # odd primes
cnt={}; ntw=0; prev=None; first=None; t0=time.time()
lo=L0 if L0%2 else L0+1   # odd numbers lo, lo+2, ...
# each segment covers odd n in [a, a+SEG) plus 2 extra for p+2 check (overlap)
a=lo
while a<L1:
    b=min(a+SEG,L1)          # starts p in [a,b)
    n=(b+2-a+1)//2+1          # odd numbers a..b+2
    m=np.ones(n,bool)
    for p in ps:
        p=int(p)
        if p*p>a+2*n: break
        st=max(p*p, ((a+p-1)//p)*p)
        if st%2==0: st+=p
        m[(st-a)//2::p]=False
    idx=np.nonzero(m[:-1]&m[1:])[0]
    tw=a+2*idx
    tw=tw[tw<b]
    if len(tw):
        if first is None: first=int(tw[0])
        ntw+=len(tw)
        if prev is not None: tw=np.concatenate([[prev],tw])
        prev=int(tw[-1])
        v,c=np.unique(np.diff(tw),return_counts=True)
        for x,y in zip(v.tolist(),c.tolist()): cnt[x]=cnt.get(x,0)+y
    a=b if b%2 else b+1
json.dump({'twins':ntw,'gaps':sum(cnt.values()),'first':first,'last':prev,'data':{str(k):v for k,v in sorted(cnt.items())}},open(out,'w'))
print('twins',ntw,'gaps',sum(cnt.values()),'first',first,'last',prev,round(time.time()-t0),'s')
