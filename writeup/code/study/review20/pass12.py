# full recount of consecutive gap pairs and triples in [1e12, 1e12+1e10), own sieve; checkpoints for sub-windows
import numpy as np, math, time, sys
sys.path.insert(0,'/home/bas/scratch/study-review20')
from sv import primes_in
L0=10**12; L1=L0+10**10; SEG=10**8; S=2048
CK={3*10**8,10**9,3*10**9,10**10}
acc2={}; acc3={}; prev=None; nprimes=0; t=time.time(); maxg=0
def add(acc,v,c):
    for a,b in zip(v.tolist(),c.tolist()): acc[a]=acc.get(a,0)+b
for s in range(L0,L1,SEG):
    P=primes_in(s,s+SEG); nprimes+=len(P)
    if prev is not None: P=np.concatenate([prev,P])
    g=np.diff(P); maxg=max(maxg,int(g.max())); assert g.max()<S
    if prev is None:
        k2=g[:-1]*S+g[1:]; k3=(g[:-2]*S+g[1:-1])*S+g[2:]
    else:  # prev holds 3 primes -> 2 old gaps; pairs: skip first old pair
        k2=(g[:-1]*S+g[1:])[1:]; k3=(g[:-2]*S+g[1:-1])*S+g[2:]
    add(acc2,*np.unique(k2,return_counts=True)); add(acc3,*np.unique(k3,return_counts=True))
    prev=P[-3:]
    w=s+SEG-L0
    if w in CK:
        k=np.array(list(acc3.keys()),np.int64); c=np.array(list(acc3.values()),np.int64)
        np.savez('/home/bas/scratch/study-review20/tri_%d.npz'%w,k=k,c=c)
        print('ck',w,'primes',nprimes,'triples',c.sum(),'pairs',sum(acc2.values()),'cells3',len(k),round(time.time()-t),flush=True)
k=np.array(list(acc2.keys()),np.int64); c=np.array(list(acc2.values()),np.int64)
np.savez('/home/bas/scratch/study-review20/pairs_full.npz',k=k,c=c)
print('done primes',nprimes,'maxgap',maxg)
