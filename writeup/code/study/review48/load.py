# own loader: per-L-function arrays from the stored vals files -> one compact npz
import numpy as np, sys
R=__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/data/study/r142/'
out={}
rows=[]
for s in ['R120','Dn','Fn']:
    z=np.load(R+'vals_'+s+'.npz')
    Ds=sorted({int(k.split('__')[0][1:]) for k in z.keys()})
    for D in Ds:
        g=z['D%d__gamma'%D]; la=z['D%d__logabs'%D]; aL=z['D%d__absL'%D]
        m=(g>=100)&(g<=600)
        rows.append((s,D,len(g),int(m.sum()),g.min(),g.max(),aL.max()))
        out['%s_%d_g'%(s,D)]=g[m]; out['%s_%d_y'%(s,D)]=la[m]
np.savez_compressed(sys.argv[1],**out)
import collections
c=collections.Counter(r[0] for r in rows); print(c)
print('total zeros', sum(r[3] for r in rows), 'outside range', sum(r[2]-r[3] for r in rows))
print('1153 present:', [r for r in rows if abs(r[1])==1153])
print('max absL', max(r[6] for r in rows))
print('gmin/gmax', min(r[4] for r in rows), max(r[5] for r in rows))
