import json, numpy as np, math
from sympy import jacobi_symbol, primerange
from sympy.ntheory import legendre_symbol
D='/home/bas/praime/writeup/code/study/lowzero/'
rows=json.load(open(D+'r132_part0.json'))+json.load(open(D+'r132_part1.json'))
print(len(rows), len(set(r['d'] for r in rows)))
PR=list(primerange(2,2000))
def chi(d,p): return 0 if d%p==0 else int(legendre_symbol(d%p,p)) if p>2 else (1 if d%8 in (1,7) else -1)
for r in rows:
    d=r['d']
    assert d%8==5 and all(chi(d,p)==-1 for p in (3,5,7,11))
    ni=[p for p in PR if p>=13 and chi(d,p)!=-1]
    assert ni[0]==r['p13']
    r['p2']=ni[1]; r['p3']=ni[2]
    r['eul100']=sum(chi(d,p)*math.log(p)/(p-1) for p in PR if 13<=p<=100)
    r['eul1000']=sum(chi(d,p)*math.log(p)/(p-1) for p in PR if 13<=p<=1000)
    r['ninert50']=sum(1 for p in PR if 13<=p<=50 and chi(d,p)==-1)
y=np.array([r['run'] for r in rows],float)
z=lambda a:(lambda b:(b-b.mean())/b.std())(np.asarray([float(t) for t in a]))
V={'lp13':z(np.log([r['p13'] for r in rows])),'S':z([r['S'] for r in rows]),'lp2':z(np.log([r['p2'] for r in rows])),
   'lp3':z(np.log([r['p3'] for r in rows])),'eul100':z([r['eul100'] for r in rows]),'eul1000':z([r['eul1000'] for r in rows]),
   'lnd':z(np.log([r['d'] for r in rows])),'ni50':z([r['ninert50'] for r in rows])}
def fit(names):
    X=np.column_stack([np.ones(len(y))]+[V[n] for n in names]); c=np.linalg.lstsq(X,y,rcond=None)[0]; res=y-X@c
    s2=res@res/(len(y)-X.shape[1]); se=np.sqrt(np.diag(s2*np.linalg.inv(X.T@X))); r2=1-res@res/np.sum((y-y.mean())**2)
    print('%-28s R2 %.3f  '%('+'.join(names),r2)+'  '.join('%s %+.2f(t %.1f)'%(n,c[i+1],c[i+1]/se[i+1]) for i,n in enumerate(names)))
for names in [["lp13","ni50","lnd"],["lp13","ni50","lnd","S"],["lp13","eul100","lnd"],["lp13","eul100","lnd","S"],["lp13","lp2","lp3","lnd","S"]]:
    fit(names)
print('corr S eul1000',np.corrcoef(V['S'],V['eul1000'])[0,1],'corr S lp13',np.corrcoef(V['S'],V['lp13'])[0,1])
print('run distribution',np.unique(y,return_counts=True))
