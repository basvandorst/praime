import numpy as np, json, sys, mpmath as mp
sys.path.insert(0,'.')
from gonek_k1 import poly_coeffs
z=np.load('../all.npz'); keys=sorted({k.rsplit('_',1)[0] for k in z.keys()})
j=json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/lprime/predictions_lprime.json'))
pj={}
for s in ['R120','Dn','Fn']:
    for r in j[s]: pj[(s,r['D'])]=r
rows=[]
for k in keys:
    s,D=k.split('_'); D=int(D); q=abs(D)
    y=z[k+'_y']; S1=float(np.sum(np.exp(2*y)))
    P=poly_coeffs(q)
    def S(T):
        x=mp.log(q*T/(2*mp.pi)); return T/(2*mp.pi)*sum(P[m]*x**m for m in range(5))
    pred=float(S(600)-S(100))
    pr=pj[(s,D)]['S1_recipe'] if (s,D) in pj else np.nan
    rows.append((s,D,S1,pred,pr,len(y)))
np.save('../k1_rows.npy',np.array([(r[1],r[2],r[3],r[4],r[5]) for r in rows]))
a=np.array([(r[2],r[3],r[4]) for r in rows])
lr=np.log(a[:,0]/a[:,1])
print('own pred vs project S1_recipe: max rel diff', np.nanmax(np.abs(a[:,1]/a[:,2]-1)))
keep=np.array([abs(r[1])!=1153 for r in rows])
print('pooled ratio (sum data/sum own pred), 547:', a[keep,0].sum()/a[keep,1].sum(), ' with 1153:', a[:,0].sum()/a[:,1].sum())
print('mean log ratio', lr[keep].mean(), 'sd per L-function', lr[keep].std(ddof=1), 'se', lr[keep].std(ddof=1)/np.sqrt(keep.sum()))
for s in ['R120','Dn','Fn']:
    m=np.array([r[0]==s for r in rows])&keep
    print(s, 'ratio', a[m,0].sum()/a[m,1].sum(), 'mean log', lr[m].mean(), '+-', lr[m].std(ddof=1)/np.sqrt(m.sum()))
print('1153:', [ (r[2]/r[3]) for r in rows if abs(r[1])==1153])
