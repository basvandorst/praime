import json,glob,numpy as np
R={}
for f in glob.glob('out_*.json'): R.update(json.load(open(f)))
a=R['actual']; DS=('0.01','0.005','0.003')
EDGES=[0,50,100,200,300,500,1000,2000,5000,1e9]
for var in ('v1','v2'):
    S=[v for k,v in R.items() if k.startswith(var+'_')]
    print(var,'n sets',len(S),'mean count',np.mean([s['n'] for s in S]))
    for key in [d for d in DS]+['H5000_'+d for d in DS]+['H2000_'+d for d in DS]:
        g=np.array([s[key] for s in S]); print('  %-12s act %.5f gue %.5f +- %.5f  diff %.5f z %.2f below %d'%(key,a[key],g.mean(),g.std(ddof=1),a[key]-g.mean(),(a[key]-g.mean())/g.std(ddof=1),(g<a[key]).sum()))
    for d in DS:
        B=np.array([s['bins_'+d] for s in S]); ab=np.array(a['bins_'+d]); diff=ab-B.mean(0); se=B.std(0,ddof=1)
        print('  bins',d,' share of Q deficit by height band', ' '.join('%g-%g:%.2f(z%.1f)'%(EDGES[i],EDGES[i+1],diff[i]/diff.sum(),diff[i]/se[i]) for i in range(len(diff))))
