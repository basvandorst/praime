import json, numpy as np
from scipy.stats import spearmanr, pearsonr
from scipy.special import digamma
d=json.load(open('data.json'))
print('missing zeros:', [o['D'] for o in d if not o['z']])
D=np.array([o['D'] for o in d]); q=np.abs(D)
g1=np.array([o['z'][0] for o in d]); L1=np.array([o['L1'] for o in d]); dL=np.array([o['dL1'] for o in d])
g=g1*np.log(q)/(2*np.pi)
print('min zero', g1.min(), 'L1 min', L1.min())
for s,name in ((D<0,'neg'),(D>0,'pos')):
    print(name, s.sum(), 'spearman lnL1 vs g', spearmanr(np.log(L1[s]),g[s])[0], ' vs gamma1', spearmanr(np.log(L1[s]),g1[s])[0])
neg=D<0
order=np.argsort(g[neg]); Dn=D[neg]
r=list(Dn[order]).index(-163)+1
order2=np.argsort(g1[neg]); r2=list(Dn[order2]).index(-163)+1
print('rank -163 by g', r, 'by gamma1', r2, 'of', neg.sum(), 'gamma1', g1[D==-163], 'g', g[D==-163])
print('lowest 12 neg by g', [(int(x), round(float(y),4)) for x,y in zip(Dn[order][:12], g[neg][order][:12])])
print('chi-4 gamma1', g1[D==-4])
a=(D<0).astype(float)
S=2*dL/L1+np.log(q/np.pi)+digamma((1+a)/2)
lq=np.log(q)
for var in ('gamma','g'):
  for sgn,name in ((-1,'neg'),(1,'pos')):
    idx=np.where(np.sign(D)==sgn)[0]
    pred=np.zeros(len(idx)); base=np.zeros(len(idx))
    y = g1 if var=='gamma' else g
    for k,i in enumerate(idx):
        nb=idx[np.abs(lq[idx]-lq[i])<0.25]
        Sm=np.median(S[nb]); ym=np.median(y[nb])
        rhs=S[i]-Sm+2/(0.25+ym**2)
        pred[k]=np.sqrt(max(2/rhs-0.25,0)) if rhs>0 else np.inf
        base[k]=ym
    yt=y[idx]
    print(var,name,'spearman pred vs true',round(spearmanr(pred,yt)[0],4),' baseline(median only)',round(spearmanr(base,yt)[0],4),
          ' scaled: pred*ln q vs g', round(spearmanr(pred*lq[idx]/(2*np.pi) if var=='gamma' else pred, g[idx])[0],4), 'n_inf',np.isinf(pred).sum())
# within-window residual check: corr of (true - median) with (pred - median)
