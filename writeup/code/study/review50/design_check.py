import os, numpy as np, json
from sympy import primefactors
R={}; R.update(json.load(open('remap_B.json'))); R.update(json.load(open('remap_Dn_Fn_C_Bp_R120.json')))
st=json.load(open('base_stats.json')); P=json.load(open(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/spacing/predictions_spacing.json'))
PR=(2,3,5,7,11,13); EXCL={'Fn':[1153],'C':[197,-395]}
def ols(X,y): return np.linalg.lstsq(X,y,rcond=None)[0]
for stn,k in (('var',0),):
  for name in ['B','Dn','Fn','C','Bp']:
    S=P[name]; Ds=[d for d in S['D'] if d not in EXCL.get(name,[])]; idx={d:i for i,d in enumerate(S['D'])}; ii=[idx[d] for d in Ds]
    rec=[R[name][str(d)] for d in Ds]; n=len(Ds)
    y=np.array([r['v0'] for r in rec]); ign=np.array(S['pred']['ignored'][stn])[ii]; law=np.array(S['pred']['law'][stn])[ii]
    X0=np.array([[1.0]+[float(abs(d)%p==0) for p in PR] for d in Ds]); lp=np.array([sum(1/p for p in primefactors(abs(d)) if p>13) for d in Ds])
    X1=np.column_stack([X0,lp])
    out=[]
    for X in (X0,X1):
        c=ols(X,y-ign)[1:4]; l=ols(X,law-ign)[1:4]; a=[]
        for j,p in enumerate((2,3,5)):
            A=np.array([abs(d)%p==0 for d in Ds]); add=np.array([rec[i]['add%d'%p][k] if A[i] else 0 for i in range(n)])
            a.append(-ols(X,add)[j+1])
        out.append((c,l,np.array(a)))
    # plain group-mean ADD response
    gm=[-np.mean([rec[i]['add%d'%p][k] for i in range(n) if abs(Ds[i])%p==0]) for p in (2,3,5)]
    print(name,'data c: %s -> %s | law: %s -> %s | ADD: %s -> %s | ADD plain group mean %s | corr(I_p, largeprime load) %s'%(
      np.round(out[0][0],5),np.round(out[1][0],5),np.round(out[0][1],5),np.round(out[1][1],5),np.round(out[0][2],5),np.round(out[1][2],5),np.round(gm,5),
      np.round([np.corrcoef(X0[:,j],lp)[0,1] for j in (1,2,3)],2)))
