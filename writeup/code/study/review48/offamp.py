# amplitude of the swapped (off-diagonal) term in the mean shift: mu - first = nuisance + g_off * off
import numpy as np, json
P=[2,3,5,7,11,13]
rows=[x for x in json.load(open('../stress_rows.json'))]
for band in ['all','lo','hi']:
    R=[x for x in rows if x['band']==band and x['q']!=1153]; n=len(R)
    cols=[np.ones(n)]
    for s in ['Fn','R120']: c=np.array([x['set']==s for x in R],float); cols.append(c-c.mean())
    c=np.array([x['lnL'] for x in R]); cols.append(c-c.mean()); ip=len(cols)
    Xn=np.array(cols).T
    mu=np.array([x['mu'] for x in R]); f1=np.array([x['first'] for x in R]); ff=np.array([x['full'] for x in R]); fo=ff-f1
    def est(i):
        X=np.column_stack([Xn[i],f1[i]-f1[i].mean(),fo[i]-fo[i].mean()])
        b=np.linalg.lstsq(X,mu[i],rcond=None)[0]
        X1=np.column_stack([Xn[i],fo[i]-fo[i].mean()])
        b1=np.linalg.lstsq(X1,(mu-f1)[i],rcond=None)[0]
        return np.array([b[-2],b[-1],b1[-1]])
    f=est(np.arange(n)); rng=np.random.default_rng(4)
    se=np.array([est(rng.integers(0,n,n)) for _ in range(1000)]).std(0,ddof=1)
    print('%-3s joint: first-order amp %.3f+-%.3f, swapped amp %.3f+-%.3f, corr; with first fixed at 1: swapped amp %.3f+-%.3f'%(band,f[0],se[0],f[1],se[1],f[2],se[2]))
