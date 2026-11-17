import numpy as np, math, sys
from lib import sieve
for (AM,MM) in ((2500,400),(10000,1000),(20000,3000)):
    m=np.arange(MM,dtype=np.int64); A=np.arange(1,AM,2,dtype=np.int64); V=m[None,:]**2+m[None,:]+A[:,None]; X=int(V.max()); s=sieve(X)
    for c in (0.4027,0.3,0.2):
        y0=int(round(X**c)); sv=np.ones(X+1,dtype=bool); sv[:2]=False
        for q in np.nonzero(s[:y0+1])[0]: sv[2*q::q]=False
        big=V>y0; P=(s[V]&big).sum(axis=1); S=(sv[V]&big).sum(axis=1); ok=S>=20; P=P[ok]; S=S[ok]; rho=P.sum()/S.sum(); res=P-rho*S; vc=(S*rho*(1-rho)).mean()
        print("m^2 + m + A, A odd < %d, m < %d (values to %d): sieve to %d (x^%.4f): rate of the coin %.3f; scatter of the prime counts around rate x survivors: variance %.2f; a coin would give %.2f; ratio %.2f"%(AM,MM,X,y0,c,rho,res.var(),vc,res.var()/vc))
