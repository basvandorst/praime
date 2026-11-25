import numpy as np, sympy, json, sys
z=np.load(sys.argv[1]); keys=sorted({k.rsplit('_',1)[0] for k in z.keys()})
rng=np.random.default_rng(7)
info=[]
for k in keys:
    s,D=k.split('_'); D=int(D)
    f=sorted(sympy.factorint(abs(D)).keys())
    info.append((s,D,f))
# patterns: by set and by smallest-prime pattern
want=[]
def add(cond,n):
    c=[x for x in info if cond(x) and x not in want]
    for i in rng.choice(len(c),size=min(n,len(c)),replace=False): want.append(c[i])
for s in ['R120','Dn','Fn']:
    add(lambda x:x[0]==s and 2 in x[2] and abs(x[1])%8==0,2)
    add(lambda x:x[0]==s and 2 in x[2] and abs(x[1])%8==4,1)
    add(lambda x:x[0]==s and 3 in x[2] and 2 not in x[2],1)
    add(lambda x:x[0]==s and 5 in x[2] and 2 not in x[2] and 3 not in x[2],1)
    add(lambda x:x[0]==s and 7 in x[2] and min(x[2])==7,1)
    add(lambda x:x[0]==s and len(x[2])==1 and x[2][0]>13,2)
    add(lambda x:x[0]==s and len(x[2])>=3,1)
# 1153 explicitly
want.append([x for x in info if x[1]==1153][0])
jobs=[]
for (s,D,f) in want:
    g=z['%s_%d_g'%(s,D)]; y=z['%s_%d_y'%(s,D)]
    idx=sorted({0,len(g)-1,int(rng.integers(1,len(g)-1))})
    for i in idx: jobs.append([s,D,float(g[i]),float(y[i])])
# the smallest |L'| zeros overall (k=-1 tail)
allv=[]
for (s,D,f) in info:
    y=z['%s_%d_y'%(s,D)]; g=z['%s_%d_g'%(s,D)]
    for i in np.argsort(y)[:3]: allv.append((float(y[i]),s,D,float(g[i])))
allv.sort()
for (yy,s,D,gg) in allv[:10]: jobs.append([s,D,gg,yy])
print(len(want),'L-functions;',len(jobs),'zeros'); print([ (w[1],w[2]) for w in want])
print('smallest:',[(a[2],round(np.exp(a[0]),4)) for a in allv[:10]])
json.dump(jobs,open(sys.argv[2],'w'))
