import os, numpy as np, json, math
from sympy import primefactors
R={}; R.update(json.load(open('remap_B.json'))); R.update(json.load(open('remap_Dn_Fn_C_Bp_R120.json')))
st=json.load(open('base_stats.json'))
P=json.load(open(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/spacing/predictions_spacing.json'))
PR=(2,3,5,7,11,13); EXCL={'Fn':[1153],'C':[197,-395]}
def design(D,extra=False):
    rows=[]
    for d in D:
        r=[1.0]+[float(abs(d)%p==0) for p in PR]
        if extra: r.append(sum(1.0/p for p in primefactors(abs(d)) if p>13))
        rows.append(r)
    return np.array(rows)
def ols(X,y): return np.linalg.lstsq(X,y,rcond=None)[0]
G={}
for name in ['B','Dn','Fn','C','Bp','R120']:
    S=P[name]; Ds=[d for d in S['D'] if d not in EXCL.get(name,[])]; idx={d:i for i,d in enumerate(S['D'])}; ii=[idx[d] for d in Ds]
    n=len(Ds); rec=[R[name][str(d)] for d in Ds]
    rng=np.random.default_rng(11); BI=[rng.integers(0,n,n) for _ in range(2000)]
    G[name]={}
    for extra in (False,True):
        X=design(Ds,extra)
        for k,stn in ((0,'var'),(1,'r')):
            y=np.array([r['v0' if k==0 else 'r0'] for r in rec])
            ign=np.array(S['pred']['ignored'][stn])[ii]; law=np.array(S['pred']['law'][stn])[ii]
            yy=y-ign; c=ols(X,yy); cl=ols(X,law-ign)
            # surrogate data set: every p|D replaced, for p in 2,3,5 separately -> predicted coefficient a_p = c_p - c_p(after)
            A=np.zeros(3); Tt=np.zeros(3); addv={}; remv={}
            for j,p in enumerate((2,3,5)):
                Am=np.array([abs(d)%p==0 for d in Ds])
                addv[p]=np.array([rec[i]['add%d'%p][k] if Am[i] else 0 for i in range(n)])
                remv[p]=np.array([rec[i]['remt%d'%p][k] if not Am[i] else 0 for i in range(n)])
                A[j]=c[j+1]-ols(X,yy+addv[p])[j+1]; Tt[j]=ols(X,yy+remv[p])[j+1]*0-(ols(X,yy+remv[p])[j+1]-c[j+1])
            # bootstrap: data coefficients and the real-base predictions jointly
            cb=[];Ab=[];Tb=[]; lb=[]
            for b in BI:
                Xb=X[b]; cc=ols(Xb,yy[b]); cb.append(cc[1:4]); lb.append(ols(Xb,(law-ign)[b])[1:4])
                Ab.append([cc[j+1]-ols(Xb,(yy+addv[p])[b])[j+1] for j,p in enumerate((2,3,5))])
                Tb.append([-(ols(Xb,(yy+remv[p])[b])[j+1]-cc[j+1]) for j,p in enumerate((2,3,5))])
            cb,Ab,Tb,lb=map(np.array,(cb,Ab,Tb,lb))
            se=cb.std(0); w=1/se**2
            def gam(cv,pat): return (w*cv*pat).sum()/(w*pat*pat).sum()
            gA=gam(c[1:4],A); gAb=np.array([gam(cb[i],Ab[i]) for i in range(len(BI))])
            gT=gam(c[1:4],Tt); gTb=np.array([gam(cb[i],Tb[i]) for i in range(len(BI))])
            gL=gam(c[1:4],cl[1:4]); gLb=np.array([gam(cb[i],lb[i]) for i in range(len(BI))])   # joint bootstrap (law matched recomputed)
            aL=gam(A,cl[1:4]); tL=gam(Tt,cl[1:4])
            G[name][(stn,extra)]=dict(gA=gA,gAse=gAb.std(),gT=gT,gTse=gTb.std(),gL=gL,gLse=gLb.std(),aL=aL,tL=tL,c=c[1:4].tolist(),A=A.tolist(),T=Tt.tolist(),law=cl[1:4].tolist(),
                                      gAb=gAb,gTb=gTb,gLb=gLb)
            print('%-4s %-3s extra=%d  data/ADD %.3f +- %.3f   data/REMT %.3f +- %.3f   data/law(joint boot) %.3f +- %.3f   ADD/law %.3f  REMT/law %.3f   largeprime coef data %s law %s'%(
                name,stn,extra,gA,gAb.std(),gT,gTb.std(),gL,gLb.std(),aL,tL, ('%+.4f'%c[-1]) if extra else '-', ('%+.4f'%cl[-1]) if extra else '-'))
for extra in (False,True):
    for stn in ('var','r'):
        for grp in (['B','Dn','Fn'],['C','Bp'],['B','Dn','Fn','C','Bp']):
            out=[]
            for key in ('gA','gT','gL'):
                g=np.array([G[n][(stn,extra)][key] for n in grp]); s=np.array([G[n][(stn,extra)][key+'se'] for n in grp]); w=1/s**2
                out.append('%s %.3f +- %.3f'%(key,(w*g).sum()/w.sum(),1/math.sqrt(w.sum())))
            aL=np.mean([G[n][(stn,extra)]['aL'] for n in grp]); tL=np.mean([G[n][(stn,extra)]['tL'] for n in grp])
            print('POOL extra=%d %-3s %-20s'%(extra,stn,'+'.join(grp)),' | '.join(out),'  mean ADD/law %.3f REMT/law %.3f'%(aL,tL))
