import numpy as np, json, math
B=json.load(open('bk_small.json'))
PR=(2,3,5,7,11,13)
for s0 in ('0.2','0.3'):
    X=[];Y=[];E=[];W=[]
    for name in ['B','Dn','Fn','C','Bp']:
        for r in B[name]:
            X.append([1.0]+[float(abs(r['D'])%p==0) for p in PR]+[1.0 if name in ('C','Bp') else 0.0])
            Y.append(r['obs'][s0]/r['n']); E.append(r['pred'][s0]); W.append(r['n'])
    X,Y,E=map(np.array,(X,Y,E))
    c=np.linalg.lstsq(X,Y/E,rcond=None)[0]   # ratio obs/BK on indicators: BK correct per prime -> coefficients 0
    cs=np.linalg.lstsq(X,Y,rcond=None)[0]; ce=np.linalg.lstsq(X,E,rcond=None)[0]
    rng=np.random.default_rng(5); n=len(Y); bs=[]
    for _ in range(2000):
        b=rng.integers(0,n,n); bs.append(np.linalg.lstsq(X[b],Y[b],rcond=None)[0])
    se=np.array(bs).std(0)
    print('s<%s per-prime change of the small-gap rate (absolute, per spacing): data / BK pred / se / data in %% of the base rate'%s0)
    for j,p in enumerate(PR):
        print('   p=%2d  data %+.5f  BK %+.5f  se %.5f   (data %+.1f%%, BK %+.1f%%)'%(p,cs[j+1],ce[j+1],se[j+1],100*cs[j+1]/cs[0],100*ce[j+1]/ce[0]))
    w=1/se[1:4]**2; g=(w*cs[1:4]*ce[1:4]).sum()/(w*ce[1:4]**2).sum()
    gb=[(w*x[1:4]*ce[1:4]).sum()/(w*ce[1:4]**2).sum() for x in bs]
    print('   amplitude on BK per-prime pattern p=2,3,5: %.2f +- %.2f'%(g,np.std(gb)))
