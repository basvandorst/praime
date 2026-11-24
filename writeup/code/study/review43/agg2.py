import sys,os,json,math, numpy as np
L=int(float(sys.argv[1])); wins=[(int(float(s.split(':')[0])), s.split(':')[1]) for s in sys.argv[2:]]
D=[(x0,np.load(f)) for x0,f in wins]
def agg(key):
    per=[]; allv=[]
    for x0,d in D:
        v=d[key]; nb=len(v); fac=nb/(nb-1)
        per.append(v.mean()*fac); allv.append(v*fac)
    allv=np.concatenate(allv)
    return per,np.mean(per),allv.std(ddof=1)/math.sqrt(len(allv)),len(allv)
for q in [3,4,5,7,10]:
    for k in [1,2,3,4]:
        for H in [10**4,10**5,10**6]:
            kk='S2win_%d_%d_%d'%(q,k,H)
            if kk not in D[0][1].files: continue
            s=[]
            for part in ['chi','col','win']:
                per,pool,err,n=agg('S2%s_%d_%d_%d'%(part,q,k,H))
                s.append('%s %.3f+-%.3f'%(part,pool,err))
            perw=agg(kk)[0]
            print('(%2d,%d) H=%.0e  %s   within per window %s'%(q,k,H,'  '.join(s),' '.join('%.3f'%p for p in perw)))
