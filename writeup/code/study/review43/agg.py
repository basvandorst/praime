import sys,os,json,math, numpy as np
# usage: agg.py L x0a:file x0b:file ...   prints table
L=int(float(sys.argv[1])); wins=[(int(float(s.split(':')[0])), s.split(':')[1]) for s in sys.argv[2:]]
hl=json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'hlpred.json')))
D=[(x0,np.load(f)) for x0,f in wins]
keys=sorted(D[0][1].files)
def agg(key, s2):
    per=[]; allv=[]; law=[]
    for x0,d in D:
        v=d[key]; nb=len(v)
        fac = nb/(nb-1) if s2 else 1.0
        per.append(v.mean()*fac); allv.append(v*fac)
    allv=np.concatenate(allv)
    return per, np.mean(per), allv.std(ddof=1)/math.sqrt(len(allv)), len(allv)
print('S1: measured chi2/dof vs HL law (exact finite-H sums)')
for q in [3,4,5,7,8,9,10,12,30,210]:
    for H in [10**4,10**5,10**6,10**7]:
        key='S1_%d_%d'%(q,H)
        if key not in keys: continue
        per,pool,err,n=agg(key,False)
        a=hl['%d_%d'%(q,H)]['a_exact']
        law=np.mean([1+a/math.log(x0+L/2) for x0,_ in D])
        print('q=%3d H=%.0e  meas %.4f +- %.4f  law %.4f  z %+.2f  nblk %d  per-window %s'%(q,H,pool,err,law,(pool-law)/err,n,' '.join('%.3f'%p for p in per)))
