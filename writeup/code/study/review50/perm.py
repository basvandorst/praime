import numpy as np, json, math
st=json.load(open('base_stats.json'))
rows=[]
for si,name in enumerate(['B','Dn','Fn','C','Bp']):
    recs=st[name]; Ls=np.array([x['Lmean'] for x in recs]); med=np.median(Ls)
    for x in recs: rows.append(dict(set=si,D=x['D'],var=x['var'],r=x['r'],Lh=int(x['Lmean']>med)))
D=np.array([r['D'] for r in rows]); setid=np.array([r['set'] for r in rows]); Lh=np.array([r['Lh'] for r in rows])
I={p:(np.abs(D)%p==0).astype(int) for p in (2,3,5,7)}
rng=np.random.default_rng(2026)
def strat_diff(y,t,strata):
    num=0;den=0
    for s in np.unique(strata):
        m=strata==s; a=y[m&(t==1)]; b=y[m&(t==0)]
        if len(a)==0 or len(b)==0: continue
        w=1/(1/len(a)+1/len(b)); num+=w*(a.mean()-b.mean()); den+=w
    return num/den
def perm_test(stat,p,others,sets=None,NP=20000):
    sel=np.ones(len(D),bool) if sets is None else np.isin(setid,sets)
    y=np.array([r[stat] for r in rows])[sel]; t=I[p][sel]
    strata=setid[sel]*100+Lh[sel]*10
    for k,o in enumerate(others): strata=strata*2+I[o][sel]
    obs=strat_diff(y,t,strata)
    # permute t within strata
    idx=[np.where(strata==s)[0] for s in np.unique(strata)]
    null=np.empty(NP)
    for i in range(NP):
        tp=t.copy()
        for ix in idx: tp[ix]=rng.permutation(t[ix])
        null[i]=strat_diff(y,tp,strata)
    pval=(np.sum(np.abs(null)>=abs(obs))+1)/(NP+1)
    return obs,null.std(),pval
for stat in ('r','var'):
    for p,others in ((2,(3,5)),(3,(2,5)),(5,(2,3))):
        for sets,lab in ((None,'all five'),([0,1,2],'B,Dn,Fn'),([3,4],'C,Bp')):
            obs,sd,pv=perm_test(stat,p,others,sets,NP=5000 if lab!='all five' else 20000)
            print('%-3s p=%d %-9s  mean(p|D) - mean(p not|D), matched on %s, set, L-half: %+.5f   perm sd %.5f  z %+.1f  two-sided p %.1e'%(stat,p,lab,others,obs,sd,obs/sd,pv))
