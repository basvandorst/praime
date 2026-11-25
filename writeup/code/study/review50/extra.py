import os, numpy as np, json, math
st=json.load(open('base_stats.json'))
P=json.load(open(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/spacing/predictions_spacing.json'))
PR=(2,3,5,7,11,13); EXCL={'Fn':[1153],'C':[197,-395]}
def ols(X,y): return np.linalg.lstsq(X,y,rcond=None)[0]
res={}
for name in ['B','Dn','Fn','C','Bp']:
    S=P[name]; Ds=[d for d in S['D'] if d not in EXCL.get(name,[])]; idx={d:i for i,d in enumerate(S['D'])}; ii=[idx[d] for d in Ds]
    rec={x['D']:x for x in st[name]}; n=len(Ds)
    I=np.array([[float(abs(d)%p==0) for p in PR] for d in Ds])
    X=np.column_stack([np.ones(n),I]); Xi=np.column_stack([X,I[:,0]*I[:,1],I[:,0]*I[:,2],I[:,1]*I[:,2]])
    rng=np.random.default_rng(21); BI=[rng.integers(0,n,n) for _ in range(2000)]
    res[name]={}
    for stn in ('var','r'):
        y=np.array([rec[d][stn] for d in Ds]); g={v:np.array(S['pred'][v][stn])[ii] for v in ('law','law2','ignored')}
        c=ols(X,y-g['ignored']); cl=ols(X,g['law']-g['ignored']); cl2=ols(X,g['law2']-g['ignored'])
        ci=ols(Xi,y-g['ignored']); cli=ols(Xi,g['law']-g['ignored'])
        boot=[]; booti=[]
        for b in BI:
            boot.append((ols(X[b],(y-g['ignored'])[b]),ols(X[b],(g['law']-g['ignored'])[b]),ols(X[b],(g['law2']-g['ignored'])[b])))
            booti.append(ols(Xi[b],(y-g['law'])[b])[7:])
        res[name][stn]=dict(c=c,cl=cl,cl2=cl2,boot=boot,inter=ci[7:],inter_law=cli[7:],inter_se=np.array(booti).std(0),meanL=np.mean([rec[d]['Lmean'] for d in Ds]))
# pooled amplitudes with joint bootstrap; on p=2,3,5 and p=7,11,13; law and law2 pattern
def gam(c,l,w): return (w*c*l).sum()/(w*l*l).sum()
for stn in ('var','r'):
    for sel,lab in ((slice(1,4),'p=2,3,5'),(slice(4,7),'p=7,11,13')):
        for grp in (['B','Dn','Fn'],['C','Bp'],['B','Dn','Fn','C','Bp']):
            G=[];Gs=[];G2=[];G2s=[]; chiL=0; chiL2=0
            for nm in grp:
                R=res[nm][stn]; cb=np.array([x[0][sel] for x in R['boot']]); w=1/cb.std(0)**2
                gb=np.array([gam(x[0][sel],x[1][sel],w) for x in R['boot']]); g2b=np.array([gam(x[0][sel],x[2][sel],w) for x in R['boot']])
                G.append(gam(R['c'][sel],R['cl'][sel],w)); Gs.append(gb.std()); G2.append(gam(R['c'][sel],R['cl2'][sel],w)); G2s.append(g2b.std())
                # shape chi2 after amplitude fit, vs law and law2 (diag)
                d=np.array([x[0][sel]-x[1][sel] for x in R['boot']]); C=np.cov(d.T)
                for pat,key in ((R['cl'],'L'),(R['cl2'],'L2')):
                    a=gam(R['c'][sel],pat[sel],w); rr=R['c'][sel]-a*pat[sel]; ch=float(rr@np.linalg.solve(C,rr))
                    if key=='L': chiL+=ch
                    else: chiL2+=ch
            G,Gs,G2,G2s=map(np.array,(G,Gs,G2,G2s)); w=1/Gs**2; w2=1/G2s**2
            m=(w*G).sum()/w.sum(); s=1/math.sqrt(w.sum()); m2=(w2*G2).sum()/w2.sum(); s2=1/math.sqrt(w2.sum())
            print('%-3s %-9s %-16s amplitude on law %.3f +- %.3f (|1-g|/se %.1f; vs law2-expected 0.90: %.1f se) | on law2 pattern %.3f +- %.3f | shape chi2 (after amplitude) law %.1f law2 %.1f on %d dof'%(
                stn,lab,'+'.join(grp),m,s,abs(1-m)/s,abs(m-0.90)/s,m2,s2,chiL,chiL2,2*len(grp)))
    # gamma vs meanL slope (P+S, p=2,3,5)
    Ls=[];Gv=[];Gse=[]
    for nm in ['B','Dn','Fn','C','Bp']:
        R=res[nm][stn]; cb=np.array([x[0][1:4] for x in R['boot']]); w=1/cb.std(0)**2
        Gv.append(gam(R['c'][1:4],R['cl'][1:4],w)); Gse.append(np.std([gam(x[0][1:4],x[1][1:4],w) for x in R['boot']])); Ls.append(R['meanL'])
    Ls,Gv,Gse=map(np.array,(Ls,Gv,Gse)); W=np.diag(1/Gse**2); A=np.column_stack([np.ones(5),Ls-10.6])
    cov=np.linalg.inv(A.T@W@A); b=cov@A.T@W@Gv
    print('%-3s gamma vs mean L (5 sets): slope %+.3f +- %.3f per unit L; sets (L,gamma): %s'%(stn,b[1],math.sqrt(cov[1,1]),' '.join('(%.2f,%.3f)'%(l,g) for l,g in zip(Ls,Gv))))
# non-additivity: interaction terms 2x3, 2x5, 3x5 (data - law), pooled
for stn in ('var','r'):
    v=np.array([res[n][stn]['inter'] for n in ['B','Dn','Fn','C','Bp']]); vl=np.array([res[n][stn]['inter_law'] for n in ['B','Dn','Fn','C','Bp']]); se=np.array([res[n][stn]['inter_se'] for n in ['B','Dn','Fn','C','Bp']])
    # inter is coefficient of y-ign; compare data-law interaction: use bootstrap of y-law
    dl=v-vl; w=1/se**2; m=(w*dl).sum(0)/w.sum(0); s=1/np.sqrt(w.sum(0))
    print('%-3s interactions (2x3, 2x5, 3x5): law (mean over sets) %s ; data - law pooled %s +- %s'%(stn,' '.join('%+.5f'%x for x in vl.mean(0)),' '.join('%+.5f'%x for x in m),' '.join('%.5f'%x for x in s)))
