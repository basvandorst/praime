# own segmented sieve: twin pairs (p, p+2) with p in [5e9, 6e9), counted per class mod 2310 and mod 210
import numpy as np, math, json
L0, W = 5*10**9, 10**9
lim = int(math.isqrt(L0+W+2))+1
s = np.ones(lim+1, bool); s[:2]=False
for i in range(2, int(lim**0.5)+1):
    if s[i]: s[i*i::i]=False
base = np.nonzero(s)[0].astype(np.int64)
seg = 10**7
cnt2310 = np.zeros(2310, np.int64); total=0
for st in range(L0, L0+W, seg):
    n = seg + 2  # include p+2 beyond segment end
    m = np.ones(n, bool)
    for p in base:
        p=int(p)
        f = max(p*p, ((st+p-1)//p)*p)
        if f >= st+n: continue
        m[f-st::p]=False
    idx = np.nonzero(m[:seg] & m[2:seg+2])[0].astype(np.int64) + st
    cnt2310 += np.bincount(idx % 2310, minlength=2310); total += len(idx)
adm2 = [a for a in range(2310) if math.gcd(a*(a+2),2310)==1]
cnt210 = np.zeros(210,np.int64)
for a in range(2310): cnt210[a%210]+=cnt2310[a]
adm1 = [a for a in range(210) if math.gcd(a*(a+2),210)==1]
X2=cnt2310[adm2].astype(float); X1=cnt210[adm1].astype(float)
print('total', total, 'adm sums', X2.sum(), X1.sum(), len(adm2), len(adm1))
print('v2310 meas', X2.var(ddof=1)/X2.mean(), 'v210 meas', X1.var(ddof=1)/X1.mean())
# subclasses per class mod 210
sub = {a:[b for b in adm2 if b%210==a] for a in adm1}
print('subclass counts', set(len(v) for v in sub.values()))
m2=X2.mean(); ss2=((X2-m2)**2).sum(); ss1=((X1-X1.mean())**2).sum()
ssw = sum(((cnt2310[v]-cnt2310[v].mean())**2).sum() for v in sub.values())
print('SS2310', ss2, 'SSwithin', ssw, 'SS210/9', ss1/9, 'check', ssw+ss1/9)
json.dump({'cnt2310':cnt2310.tolist()}, open('/home/bas/scratch/study-review10/count5.json','w'))
