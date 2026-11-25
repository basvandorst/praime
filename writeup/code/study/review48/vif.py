import numpy as np, json
P=[2,3,5,7,11,13]
d=[r for r in json.load(open('../perL.json')) if r['q']!=1153]
n=len(d); cols={}; 
for s in ['Fn','R120']: cols['set_'+s]=np.array([r['set']==s for r in d],float)
cols['lnL']=np.array([r['lnL'] for r in d])
for p in P: cols['I%d'%p]=np.array([p in r['bad'] for r in d],float)
cols['nbad']=np.array([len(r['bad']) for r in d],float)
names=list(cols)[:-1]
M=np.column_stack([cols[k]-cols[k].mean() for k in names])
print('condition number of the centred design (scaled columns):', np.linalg.cond(M/M.std(0)))
for j,k in enumerate(names):
    others=np.delete(M,j,axis=1); b=np.linalg.lstsq(np.column_stack([np.ones(n),others]),M[:,j],rcond=None)[0]
    res=M[:,j]-np.column_stack([np.ones(n),others])@b; r2=1-res.var()/M[:,j].var()
    print('%-9s VIF %.2f  corr with lnL %+.3f'%(k,1/(1-r2),np.corrcoef(M[:,j],cols['lnL'])[0,1]))
# within-set correlation of lnL with the indicators
for s in ['R120','Dn','Fn']:
    m=np.array([r['set']==s for r in d])
    print(s, 'lnL range %.3f-%.3f'%(cols['lnL'][m].min(),cols['lnL'][m].max()), 'corr(lnL, I2) %.2f'%np.corrcoef(cols['lnL'][m],cols['I2'][m])[0,1])
