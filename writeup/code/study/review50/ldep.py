import os, numpy as np, json, math
st=json.load(open('base_stats.json'))
P=json.load(open(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/spacing/predictions_spacing.json'))
PR=(2,3,5,7,11,13); EXCL={'Fn':[1153],'C':[197,-395]}
def design(D): return np.array([[1.0]+[float(abs(d)%p==0) for p in PR] for d in D])
def ols(X,y): return np.linalg.lstsq(X,y,rcond=None)[0]
def get(name,Ds,stn):
    S=P[name]; idx={d:i for i,d in enumerate(S['D'])}; ii=[idx[d] for d in Ds]; rec={x['D']:x for x in st[name]}
    y=np.array([rec[d][stn] for d in Ds]); ign=np.array(S['pred']['ignored'][stn])[ii]; law=np.array(S['pred']['law'][stn])[ii]
    law2=np.array(S['pred']['law2'][stn])[ii]
    return y,ign,law,law2
for a,b in (('B','Bp'),('R120','C')):
    common=sorted(set(P[a]['D'])&set(P[b]['D'])-set(EXCL.get(a,[]))-set(EXCL.get(b,[])))
    print(a,b,'common discriminants',len(common),'of',len(P[a]['D']),len(P[b]['D']))
    X=design(common); n=len(common)
    for stn in ('var','r'):
        ya,ia,la,_=get(a,common,stn); yb,ib,lb,_=get(b,common,stn)
        ca=ols(X,ya-ia); cla=ols(X,la-ia); cb=ols(X,yb-ib); clb=ols(X,lb-ib)
        rng=np.random.default_rng(9); BI=[rng.integers(0,n,n) for _ in range(3000)]
        da=np.array([ols(X[i],(ya-la)[i]) for i in BI]); db=np.array([ols(X[i],(yb-lb)[i]) for i in BI])
        wa=1/da.std(0)[1:4]**2; wb=1/db.std(0)[1:4]**2
        ga=(wa*ca[1:4]*cla[1:4]).sum()/(wa*cla[1:4]**2).sum(); gb=(wb*cb[1:4]*clb[1:4]).sum()/(wb*clb[1:4]**2).sum()
        gab=np.array([(wa*(x[1:4]+cla[1:4])*cla[1:4]).sum()/(wa*cla[1:4]**2).sum() for x in da])
        gbb=np.array([(wb*(x[1:4]+clb[1:4])*clb[1:4]).sum()/(wb*clb[1:4]**2).sum() for x in db])
        d=gb-ga; dsd=(gbb-gab).std(); unp=math.hypot(gab.std(),gbb.std())
        # residual correlation per L-function
        ra=(ya-la)-X@ols(X,ya-la); rb=(yb-lb)-X@ols(X,yb-lb)
        print('  %-3s gamma(%s)=%.3f gamma(%s)=%.3f diff %+.3f  paired se %.3f  unpaired se %.3f  corr(gamma_a,gamma_b) %.2f  corr of per-L residuals %.2f'%(stn,a,ga,b,gb,d,dsd,unp,np.corrcoef(gab,gbb)[0,1],np.corrcoef(ra,rb)[0,1]))
        if stn=='var': print('     "no L-dependence" alternative +0.16: (0.16 - diff)/se = %.1f (paired), %.1f (unpaired)'%((0.16-d)/dsd,(0.16-d)/unp))
        # per-prime ratio data/law at the two heights
        print('     per-prime c_b/c_a data: %s   law: %s'%(' '.join('%.3f'%(cb[j]/ca[j]) for j in (1,2,3)),' '.join('%.3f'%(clb[j]/cla[j]) for j in (1,2,3))))
