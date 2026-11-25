# stress tests of readings (1)-(4) on splits; predictions recomputed per split (own code)
import numpy as np, json, sys, sympy
sys.path.insert(0,'.')
from gonek_k1 import poly_coeffs
P=[2,3,5,7,11,13]
z=np.load('../all.npz'); keys=sorted({k.rsplit('_',1)[0] for k in z.keys()})
G=np.load('../zgrid_2e5.npz'); r=G['r']; wr=G['w']; Z=G['Z']*np.exp(-(r/30.)**2); lnr=np.log(r)
pp=json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/lprime/predictions_lprime.json'))['per_prime']
def li2(x): return sum(x**k/k**2 for k in range(1,200))
def rho_m1(bad):
    s=np.zeros(len(r),complex)
    for p in bad:
        lp=np.log(p); s+=np.log1p((4/p)*np.sin(r*lp/2)**2/((1-2/p)+np.exp(-1j*r*lp)/p))
    return np.expm1(s)
def Soff(bad,Ls):
    h=Z*rho_m1(bad)
    return np.array([2*np.sum(wr*lnr*np.real(np.exp(-1j*r*L)*h)) for L in Ls])/(np.pi*Ls)
cache={}
def polyS(q,T):
    if q not in cache: cache[q]=[float(x) for x in poly_coeffs(q)]
    Pq=cache[q]; x=np.log(q*T/(2*np.pi)); return T/(2*np.pi)*sum(Pq[m]*x**m for m in range(5))
rows=[]
for k in keys:
    s,D=k.split('_'); D=int(D); q=abs(D); bad=sorted(sympy.factorint(q))
    g=z[k+'_g']; ly=z[k+'_y']
    for band,(a,b) in {'all':(100,600),'lo':(100,350),'hi':(350,600)}.items():
        m=(g>=a)&(g<b) if band!='all' else np.ones(len(g),bool)
        L=np.log(q*g[m]/(2*np.pi)); y=ly[m]-np.log(L)
        # predicted mean shifts at the zeros (exact ratios) : evaluate on 8 quantiles of L in the band
        Lq=np.quantile(L,(np.arange(8)+0.5)/8)
        first=sum(np.log(p)*np.log(p/(p-1)) for p in bad)*np.mean(1/L)
        off=np.mean(Soff(bad,Lq))
        varlead=-sum(0.5*li2(1/p) for p in bad)
        varldm=varlead-sum(pp[str(p)]['kappa_ldm_per_lnp_over_L'][1]*np.log(p)*np.mean(1/L) if str(p) in pp else 0 for p in bad)
        X=np.exp(-2*ly[m])
        rows.append(dict(set=s,D=D,q=q,bad=bad,band=band,mu=y.mean(),v=y.var(ddof=1),lnL=np.log(L).mean(),
            first=first,full=first+off,varlead=varlead,varldm=varldm,
            S1=np.exp(2*ly[m]).sum(),S1p=polyS(q,b)-polyS(q,a),Sm1=X.sum(),
            Pm=6/np.pi**2*np.prod([p/(p+1) for p in bad])*(b-a)/(2*np.pi), Pn=6/np.pi**2*(b-a)/(2*np.pi)))
json.dump(rows,open('../stress_rows.json','w'))
def analyse(sel,label):
    R=[x for x in rows if sel(x)]; n=len(R)
    cols=[np.ones(n)]
    S=sorted({x['set'] for x in R})
    for s in S[1:]: c=np.array([x['set']==s for x in R],float); cols.append(c-c.mean())
    c=np.array([x['lnL'] for x in R]); cols.append(c-c.mean())
    for p in P: c=np.array([p in x['bad'] for x in R],float); cols.append(c-c.mean())
    X=np.array(cols).T; ip=X.shape[1]-len(P)
    A=lambda key: np.array([x[key] for x in R])
    def est(i):
        Xi=X[i]; o=lambda y: np.linalg.lstsq(Xi,y[i],rcond=None)[0][ip:]
        bm=o(A('mu')); bf=o(A('full')); b1=o(A('first'))
        gM=np.sum(bm[:3]*bf[:3])/np.sum(bf[:3]**2); g1=np.sum(bm[:3]*b1[:3])/np.sum(b1[:3]**2)
        bv=o(A('v')); bl=o(A('varlead')); bd=o(A('varldm'))
        gV=np.sum((bv-bl)[:3]*(bd-bl)[:3])/np.sum(((bd-bl)[:3])**2)
        bs=o(np.log(A('S1')/A('S1p')))
        lv=np.log(A('S1')[i].sum()/A('S1p')[i].sum())
        bn=o(np.log(A('Sm1')/A('Pm')))
        lf=np.log(A('Pm')/A('Pn')); X2=np.column_stack([Xi[:,:ip],lf[i]-lf[i].mean()])
        gPC=np.linalg.lstsq(X2,np.log(A('Sm1')/A('Pn'))[i],rcond=None)[0][-1]
        return np.concatenate([bm[:4],[gM,g1],bv[:3],[gV],bs[:3],[lv],bn[:3],[gPC]])
    full=est(np.arange(n)); rng=np.random.default_rng(3)
    se=np.array([est(rng.integers(0,n,n)) for _ in range(400)]).std(0,ddof=1)
    names=['mu2','mu3','mu5','mu7','gM(full)','gM(1st)','v2','v3','v5','g1/L(var)','S1_2','S1_3','S1_5','S1lev','Sm1_2','Sm1_3','Sm1_5','gPC']
    print('== %-22s n=%d'%(label,n)); print('   '+'  '.join('%s %+.4f(%.4f)'%(a,b,c) for a,b,c in zip(names,full,se)))
    return full,se
res={}
res['all']=analyse(lambda x:x['band']=='all' and x['q']!=1153,'all 547')
analyse(lambda x:x['band']=='all','all + 1153')
res['lo']=analyse(lambda x:x['band']=='lo' and x['q']!=1153,'heights 100-350')
res['hi']=analyse(lambda x:x['band']=='hi' and x['q']!=1153,'heights 350-600')
for s in ['R120','Dn','Fn']: res[s]=analyse(lambda x,s=s:x['band']=='all' and x['set']==s and x['q']!=1153,s)
res['pos']=analyse(lambda x:x['band']=='all' and x['D']>0 and x['q']!=1153,'D>0')
res['neg']=analyse(lambda x:x['band']=='all' and x['D']<0,'D<0')
f0,s0=res['all']
print('largest |split - all|/se_split:')
for k,(f,s) in res.items():
    if k=='all': continue
    zz=(f-f0)/s; j=np.argmax(np.abs(zz)); print('  %-5s max at item %d: %.1f'%(k,j,zz[j]))
