# independent recount of twin starts in [1e9,2e9) and consecutive distances
import numpy as np, sys
L0, L1 = 10**9, 2*10**9
lim = int(L1**0.5)+2
sm = np.ones(lim+1, bool); sm[:2]=False
for i in range(2,int(lim**0.5)+1):
    if sm[i]: sm[i*i::i]=False
ps = np.nonzero(sm)[0][1:]  # odd primes
SEG = 10**8
starts=[]; prev_last=None
for s in range(L0, L1, SEG):
    e = min(s+SEG, L1)
    # odd numbers s+1, s+3, ... (s even)
    n = (e - s)//2
    m = np.ones(n, bool)
    for p in ps:
        p=int(p)
        st = max(p*p, ((s+1+p-1)//p)*p)
        if st % 2 == 0: st += p
        if st >= e: continue
        m[(st-(s+1))//2::p] = False
    pr = (np.nonzero(m)[0]*2 + s+1).astype(np.int64)
    if prev_last is not None: pr = np.concatenate([[prev_last], pr])
    tw = pr[:-1][np.diff(pr)==2]
    starts.append(tw); prev_last = int(pr[-1])
st = np.concatenate(starts)
st = st[(st>=L0)]
print('twins', len(st), 'max', st[-1])
dd = np.diff(st)
print('gaps', len(dd), 'non-mult6', int(np.sum(dd%6)))
v,c = np.unique(dd, return_counts=True)
np.save('/home/bas/scratch/study-review15/gaps_hist.npy', np.vstack([v,c]))
D = dict(zip(v.tolist(), c.tolist()))
for lo,hi in [(6,30),(36,60),(96,120),(576,600),(606,660),(966,1020),(1146,1200)]:
    print(lo,hi,sum(x for d,x in D.items() if lo<=d<=hi))
print('tail>1200', sum(x for d,x in D.items() if d>1200), 'd6',D[6],'d30',D[30],'d300',D[300])
