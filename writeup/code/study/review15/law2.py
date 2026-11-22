import math, sys, itertools, numpy as np
import law
from sympy import primerange
QS = list(primerange(2, 140)); GEN2 = np.array([q for q in law.PR if q > 139])
GT2 = {}
def gt2(k):
    if k not in GT2: GT2[k] = float(np.sum(np.log1p(-k/GEN2) - k*np.log1p(-1/GEN2))) - k*(k-1)/2/(law.Qe*math.log(law.Qe))
    return GT2[k]
def lnSf(pts):
    k=len(pts); s=0.0
    for q in QS:
        nu=len({x%q for x in pts})
        if nu>=q: return -math.inf
        s+=math.log(1-nu/q)-k*math.log(1-1/q)
    return s+gt2(k)
def exact(d):
    lam=law.lam; B=[0,2,d,d+2]; H=list(range(6,d,6)); l02=lnSf([0,2]); tot=0.0
    # prune: DFS over subsets, stop extending if S==0 (adding points cannot make it admissible)
    def rec(i, pts, m):
        nonlocal tot
        l=lnSf(pts)
        if l==-math.inf: return
        tot += (-1)**m*math.exp(l-l02-(2+2*m)*math.log(lam))
        for j in range(i,len(H)): rec(j+1, pts+[H[j],H[j]+2], m+1)
    rec(0,B,0); return tot
import json
ref=json.load(open('/home/bas/praime/writeup/data/study/r59.json'))['law']
for d in map(int, sys.argv[1].split(',')):
    o2,_=law.order2(d); ex=exact(d)
    print('d=%d exact=%.9f order2=%.9f repo=%.9f o2/exact-1=%.2e'%(d,ex,o2,ref[str(d)][1],o2/ex-1),flush=True)
