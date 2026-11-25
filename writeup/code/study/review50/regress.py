import os, numpy as np, json, math, sys
st=json.load(open('base_stats.json'))
P=json.load(open(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/spacing/predictions_spacing.json'))
PR=(2,3,5,7,11,13)
EXCL={'Fn':[1153],'C':[197,-395]}
def design(D): return np.array([[1.0]+[float(abs(d)%p==0) for p in PR] for d in D])
def ols(X,y): return np.linalg.lstsq(X,y,rcond=None)[0]
res={}
mode=sys.argv[1] if len(sys.argv)>1 else 'registered'
for name in ['B','Dn','Fn','C','Bp','R120']:
    S=P[name]; Dp=S['D']; idx={d:i for i,d in enumerate(Dp)}
    recs={x['D']:x for x in st[name]}
    keepD=[d for d in Dp if not (mode=='registered' and d in EXCL.get(name,[]))]
    ii=np.array([idx[d] for d in keepD])
    V=np.array([recs[d]['var'] for d in keepD]); R=np.array([recs[d]['r'] for d in keepD])
    Lm=np.array([recs[d]['Lmean'] for d in keepD])
    X=design(keepD); n=len(keepD)
    rng=np.random.default_rng(7); BI=[rng.integers(0,n,n) for _ in range(2000)]
    out={}
    for stn,y in (('var',V),('r',R)):
        ign=np.array(S['pred']['ignored'][stn])[ii]; law=np.array(S['pred']['law'][stn])[ii]; law2=np.array(S['pred']['law2'][stn])[ii]
        c=ols(X,y-ign); cl=ols(X,law-ign); cl2=ols(X,law2-ign)
        # my own regression: stat on [1, L-10.5, I_p] (no use of the level formula)
        X2=np.column_stack([X[:,0],Lm-10.5,X[:,1:]]); c2=ols(X2,y)
        cb=np.array([ols(X[b],(y-ign)[b]) for b in BI]); se=cb.std(0)
        c2b=np.array([ols(X2[b],y[b]) for b in BI]); se2=c2b.std(0)
        w=1/se[1:4]**2; l=cl[1:4]
        gam=(w*c[1:4]*l).sum()/(w*l*l).sum(); gb=np.array([(w*x[1:4]*l).sum()/(w*l*l).sum() for x in cb])
        gam2=(w*c2[2:5]*l).sum()/(w*l*l).sum(); gb2=np.array([(w*x[2:5]*l).sum()/(w*l*l).sum() for x in c2b])
        glaw2=(w*cl2[1:4]*l).sum()/(w*l*l).sum()
        # gamma over p=2..13 (all six) for p=7,11,13 check
        w6=1/se[1:]**2; l6=cl[1:]
        g713=(w6[3:]*c[4:]*l6[3:]).sum()/(w6[3:]*l6[3:]**2).sum(); g713b=np.array([(w6[3:]*x[4:]*l6[3:]).sum()/(w6[3:]*l6[3:]**2).sum() for x in cb])
        lev=np.mean(y-law); levse=np.std([np.mean((y-law)[b]) for b in BI])
        out[stn]=dict(c=c[1:].tolist(),se=se[1:].tolist(),law=cl[1:].tolist(),gamma=gam,gse=gb.std(),gamma_ownreg=gam2,gse_own=gb2.std(),
                      c_own=c2[2:].tolist(),Lslope=c2[1],gamma_law2=glaw2,g713=g713,g713se=g713b.std(),lev=lev,levse=levse,meanL=Lm.mean(),
                      gb=gb.tolist(), g713b=g713b.tolist())
        print('%-4s %-3s n=%d  c(2..13)='%(name,stn,n)+' '.join('%+.5f'%x for x in c[1:])+'  se '+' '.join('%.5f'%x for x in se[1:]))
        print('          law  '+' '.join('%+.5f'%x for x in cl[1:])+'   gamma %.3f +- %.3f (own L-regression %.3f +- %.3f; law2 %.3f)  gamma(7,11,13) %.3f +- %.3f  level-law %+.5f +- %.5f  L-slope %+.5f'%(gam,gb.std(),gam2,gb2.std(),glaw2,g713,g713b.std(),lev,levse,c2[1]))
    res[name]=out
json.dump(res,open('regress_%s.json'%mode,'w'))
def pool(names,stn,key='gamma',sek='gse'):
    g=np.array([res[n][stn][key] for n in names]); s=np.array([res[n][stn][sek] for n in names]); w=1/s**2
    m=(w*g).sum()/w.sum(); return m,1/math.sqrt(w.sum()),((g-m)**2*w).sum()
for names in (['B','Dn','Fn'],['C','Bp'],['B','Dn','Fn','C','Bp']):
    for stn in ('var','r'):
        m,s,chi=pool(names,stn); m2,s2,_=pool(names,stn,'gamma_ownreg','gse_own'); m3,s3,chi3=pool(names,stn,'g713','g713se')
        print('POOL',names,stn,'gamma %.3f +- %.3f (het chi2 %.1f)  own-reg %.3f +- %.3f   p=7,11,13: %.3f +- %.3f (het %.1f)'%(m,s,chi,m2,s2,m3,s3,chi3))
