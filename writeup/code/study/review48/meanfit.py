import numpy as np, json, sys
P=[2,3,5,7,11,13]
d=json.load(open('../perL.json')); ms={m['D']:m for m in json.load(open(sys.argv[1]))}
def run(sel,label):
    rows=[r for r in d if sel(r)]
    n=len(rows)
    cols=[np.ones(n)]
    S=sorted({r['set'] for r in rows})
    for s in S[1:]: c=np.array([r['set']==s for r in rows],float); cols.append(c-c.mean())
    c=np.array([r['lnL'] for r in rows]); cols.append(c-c.mean())
    for p in P: c=np.array([p in r['bad'] for r in rows],float); cols.append(c-c.mean())
    X=np.array(cols).T; ip=X.shape[1]-len(P)
    mu=np.array([r['mu'] for r in rows]); f1=np.array([ms[r['D']]['first'] for r in rows]); fo=np.array([ms[r['D']]['off'][-1] for r in rows])
    ols=lambda y: np.linalg.lstsq(X,y,rcond=None)[0]
    b=ols(mu); b1=ols(f1); b2=ols(f1+fo)
    rng=np.random.default_rng(11); bs=[]
    for _ in range(1000):
        i=rng.integers(0,n,n); Xi=X[i]
        bb=np.linalg.lstsq(Xi,mu[i],rcond=None)[0]; p1=np.linalg.lstsq(Xi,f1[i],rcond=None)[0]; p2=np.linalg.lstsq(Xi,(f1+fo)[i],rcond=None)[0]
        g1=np.sum(bb[ip:ip+3]*p1[ip:ip+3])/np.sum(p1[ip:ip+3]**2); g2=np.sum(bb[ip:ip+3]*p2[ip:ip+3])/np.sum(p2[ip:ip+3]**2)
        # direct amplitude: mu on nuisance + predicted shift
        X2=np.column_stack([Xi[:,:ip],(f1+fo)[i]-(f1+fo)[i].mean()]); gd=np.linalg.lstsq(X2,mu[i],rcond=None)[0][-1]
        bs.append(np.concatenate([bb[ip:],[g1,g2,gd]]))
    bs=np.array(bs); se=bs.std(0,ddof=1)
    g1=np.sum(b[ip:ip+3]*b1[ip:ip+3])/np.sum(b1[ip:ip+3]**2); g2=np.sum(b[ip:ip+3]*b2[ip:ip+3])/np.sum(b2[ip:ip+3]**2)
    X2=np.column_stack([X[:,:ip],(f1+fo)-(f1+fo).mean()]); gd=np.linalg.lstsq(X2,mu,rcond=None)[0][-1]
    np.set_printoptions(precision=4,suppress=True,linewidth=160)
    print('==',label,'n',n)
    print(' data fit      ',b[ip:],' se',se[:6])
    print(' first order   ',b1[ip:],'  z',(b[ip:]-b1[ip:])/se[:6])
    print(' ratios full   ',b2[ip:],'  z',(b[ip:]-b2[ip:])/se[:6])
    print(' amplitude (p=2,3,5) vs first order %.4f +- %.4f ; vs full ratios %.4f +- %.4f ; direct (all D, full) %.4f +- %.4f'%(g1,se[6],g2,se[7],gd,se[8]))
run(lambda r:r['q']!=1153,'F3 (547)')
run(lambda r:True,'F3 + 1153')
for s in ['R120','Dn','Fn']: run(lambda r,s=s:r['set']==s and r['q']!=1153,s)
run(lambda r:r['D']>0 and r['q']!=1153,'D>0'); run(lambda r:r['D']<0,'D<0')
