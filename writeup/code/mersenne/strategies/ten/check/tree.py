import numpy as np, math
from lib import sieve
MEXP=[2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,136279841]
N=10**7; s=sieve(N+1); gp=np.zeros(N+2,dtype=np.int32)
for p in np.nonzero(s)[0]: gp[p::p]=p
P=np.nonzero(s[:N+1])[0]
for name,sh in (("largest prime factor of q - 1 (Fermat founders)",-1),("largest prime factor of q + 1 (Mersenne founders)",1)):
    par=gp[P+sh]; root=np.zeros(N+1,dtype=np.int32); dep=np.zeros(N+1,dtype=np.int16)
    for q,pa in zip(P.tolist()[1:],par.tolist()[1:]):
        if pa==2: root[q]=q; dep[q]=1
        else: root[q]=root[pa]; dep[q]=dep[pa]+1
    print("==",name,"up to",N)
    for X in (10**4,10**5,5*10**5,10**6,10**7):
        Q=P[(P>2)&(P<=X)]; rr,cc=np.unique(root[Q],return_counts=True); oo=np.argsort(-cc)[:6]; print("  x=%8d:"%X," ".join("%d: %.2f%%"%(rr[i],100*cc[i]/len(Q)) for i in oo))
    Q=P[P>2]; d=dep[Q]; print("  depth: mean %.3f max %d (first at %d); shares by depth:"%(d.mean(),d.max(),Q[d.argmax()]),np.round(np.bincount(d)/len(Q),4).tolist())
    ex=[e for e in MEXP if 2<e<=N]; rr=[int(root[e]) for e in ex]; sh_={int(r):float((root[Q]==r).mean()) for r in set(rr)}
    print("  the %d odd exponents by founder:"%len(ex),{r:(rr.count(r),round(len(ex)*sh_[r],1)) for r in sorted(set(rr))},"(found, expected)")
