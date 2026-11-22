import numpy as np, json, math, sys, time
sys.path.insert(0,'/home/bas/scratch/study-review20'); from sv import primes_in
L0=10**13; L1=L0+10**10; SEG=10**8
cnt={}; lastp=None; lastw=None; ntw=0; t=time.time()
for s in range(L0,L1,SEG):
    P=primes_in(s,s+SEG)
    if lastp is not None: P=np.concatenate([[lastp],P])
    lastp=int(P[-1])
    tw=P[:-1][np.diff(P)==2]; ntw+=len(tw)
    if lastw is not None: tw=np.concatenate([[lastw],tw])
    lastw=int(tw[-1]); v,c=np.unique(np.diff(tw),return_counts=True)
    for a,b in zip(v.tolist(),c.tolist()): cnt[a]=cnt.get(a,0)+b
print('twins',ntw,'gaps',sum(cnt.values()),round(time.time()-t),'s')
# note: a twin (p,p+2) with p<L1<=p+2 would be missed/included? check boundary
d=json.load(open('/home/bas/praime/writeup/data/study/r67.json'))
rep={int(k):v for k,v in d['data'].items()}
print('histogram equal to r67.json:', rep==cnt, 'twins equal', d['twins']==ntw)
json.dump({str(k):v for k,v in sorted(cnt.items())},open('/home/bas/scratch/study-review20/tw13_counts.json','w'))
