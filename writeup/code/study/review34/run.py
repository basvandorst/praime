import numpy as np, math, sys, json
from rq import *
DS=(0.01,0.005,0.003)
EDGES=[0,50,100,200,300,500,1000,2000,5000,1e9]
def stats(g):
    W=rows(g,DS); out={}
    for d in DS:
        out[str(d)]=2/d*math.sqrt(4*W[d].sum())
        out['bins_'+str(d)]=[float(W[d][(g>=a)&(g<b)].sum()*2) for a,b in zip(EDGES[:-1],EDGES[1:])]
    for H in (5000,2000):
        gg=g[g<=H]; r=rms(gg,DS)
        for d in DS: out['H%d_%s'%(H,d)]=r[d]
    out['n']=len(g)
    return out
s0,s1=int(sys.argv[1]),int(sys.argv[2]); res={}
if s0<0:
    res['actual']=stats(Z)
else:
    lo1=float(Nbar(Z[0])); hi=float(Nbar(Z[-1]))
    for s in range(s0,s1):
        u=gue_unf(s); rng=np.random.default_rng(10_000+s); o=rng.uniform(12000,14000)
        for name,lo in (('v1',lo1),('v2',-0.125)):
            v=u-o; v=v[(v>=0)&(v<=hi-lo)]; g=Ninv(lo+v)
            res['%s_%d'%(name,s)]=stats(g)
        print(s,flush=True)
json.dump(res,open('out_%d.json'%s0,'w'))
