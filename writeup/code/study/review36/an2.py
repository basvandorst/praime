import json, numpy as np
from scipy.stats import spearmanr, pearsonr
from scipy.special import digamma
d=json.load(open('data.json')); R=np.load('race.npy')
D=np.array([o['D'] for o in d]); q=np.abs(D); g1=np.array([o['z'][0] for o in d])
L1=np.array([o['L1'] for o in d]); dL=np.array([o['dL1'] for o in d]); a=(D<0)*1.0
S=2*dL/L1+np.log(q/np.pi)+digamma((1+a)/2); sh=R[1]
g=g1*np.log(q)/(2*np.pi)
print('share vs gamma1 spearman', spearmanr(sh,g1)[0], ' vs g', spearmanr(sh,g)[0])
print('share vs S', spearmanr(sh,S)[0], ' vs dL/L', spearmanr(sh,dL/L1)[0], ' vs lnL1', spearmanr(sh,np.log(L1))[0])
for lo,hi in ((0.05,0.1),(0.25,0.35),(0.35,0.5),(0.5,1),(1,2),(2,10)):
    s=(g1>=lo)&(g1<hi); print(lo,hi,s.sum(), 'mean share', sh[s].mean().round(3) if s.sum() else None, 'min',sh[s].min().round(3) if s.sum() else None)
# BH correction factor spread
import sympy
P=np.array(list(sympy.primerange(3,100001)))
chi=np.zeros((len(D),len(P)))
for j,p in enumerate(P):
    qr=np.full(p,-1.0); qr[(np.arange(1,p)**2)%p]=1; qr[0]=0; chi[:,j]=qr[D%p]
corr=np.sum(np.log((1-chi/(P-1))/(1-chi/P)),axis=1)
print('ln correction factor (odd p<1e5): sd',corr.std(),'range',corr.min(),corr.max(), ' sd lnL1', np.log(L1).std())
