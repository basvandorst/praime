import json, numpy as np, sympy
from scipy.stats import spearmanr
from scipy.special import digamma
d=json.load(open('data.json'))
D=np.array([o['D'] for o in d]); q=np.abs(D); g1=np.array([o['z'][0] for o in d]); lq=np.log(q)
L1=np.array([o['L1'] for o in d]); dL=np.array([o['dL1'] for o in d]); a=(D<0)*1.0
arch=np.log(q/np.pi)+digamma((1+a)/2)
def chi_of(n):
    return np.array([sympy.jacobi_symbol(0,1)*0 for _ in D]) # unused
pp=[]
for p in sympy.primerange(2,100001):
    pk=p
    while pk<=100000: pp.append((pk,p)); pk*=p
pp.sort()
from math import log
def kron(Dv,n): return np.array([int(sympy.ntheory.residue_ntheory.jacobi_symbol(int(x)%n,n)) if n%2 else 0 for x in Dv])
# compute chi(p) per prime then powers
P=sorted(set(p for _,p in pp))
chiP={}
m8=D%8
chiP[2]=np.where(D%2==0,0,np.where((m8==1)|(m8==7),1,-1)).astype(float)
for p in P[1:]:
    qr=np.full(p,-1.0); qr[(np.arange(1,p)**2)%p]=1; qr[0]=0; chiP[p]=qr[D%p]
def pred_from(S):
    out=np.zeros(len(D))
    for sgn in (-1,1):
        idx=np.where(np.sign(D)==sgn)[0]
        for i in idx:
            nb=idx[np.abs(lq[idx]-lq[i])<0.25]
            rhs=S[i]-np.median(S[nb])+2/(0.25+np.median(g1[nb])**2)
            out[i]=np.sqrt(max(2/rhs-0.25,0)) if rhs>0 else 1e9
    return out
Sfull=2*dL/L1+arch
acc=np.zeros(len(D)); cuts=[30,199,1000,10000,100000]; ci=0
for pk,p in pp:
    k=round(log(pk)/log(p)); acc+= -(chiP[p]**k)*log(p)/pk
    if ci<len(cuts) and (pk>=cuts[ci]):
        pass
res={}
for X in cuts:
    acc=np.zeros(len(D))
    for pk,p in pp:
        if pk>X: break
        k=round(log(pk)/log(p)); acc+= -(chiP[p]**k)*log(p)/pk
    S=2*acc+arch; pr=pred_from(S)
    r=[spearmanr(pr[np.sign(D)==s],g1[np.sign(D)==s])[0] for s in (-1,1)]
    # scaled comparison
    print('X',X,'spearman neg/pos',np.round(r,3), 'corr(Strunc,Sfull)',np.corrcoef(S,Sfull)[0,1].round(3))
pr=pred_from(Sfull); print('full', [round(spearmanr(pr[np.sign(D)==s],g1[np.sign(D)==s])[0],3) for s in (-1,1)])
# R^2 of linear regression of gamma1 (within-sign, with ln q) on chi(p), p<=199
for s in (-1,1):
    m=np.sign(D)==s
    X=np.column_stack([np.ones(m.sum()),lq[m]]+[chiP[p][m] for p in P if p<=199])
    b,*_=np.linalg.lstsq(X,g1[m],rcond=None); r2=1-np.var(g1[m]-X@b)/np.var(g1[m])
    X0=np.column_stack([np.ones(m.sum()),lq[m]]); b0,*_=np.linalg.lstsq(X0,g1[m],rcond=None); r20=1-np.var(g1[m]-X0@b0)/np.var(g1[m])
    print('sign',s,'R2 gamma1 ~ lnq + chi(p<=199):',round(r2,3),' lnq alone',round(r20,3))
