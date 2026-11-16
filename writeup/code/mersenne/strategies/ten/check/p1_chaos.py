import numpy as np, re
from lib import *
N=10**7; s=sieve(N); P=np.nonzero(s)[0]; P=P[P>=5]; g=np.diff(P)
sym=(g%6)//2
w=''.join(map(str,sym.tolist()))
print("primes",len(P),"gap classes 0,2,4 mod 6:",[int((sym==k).sum()) for k in range(3)])
print("forbidden '1 0* 1' (gap=2, zeros, gap=2):",len(re.findall(r'10*1',w)),"  '2 0* 2':",len(re.findall(r'20*2',w)))
def cgr(symbols,corners,size):
    x=np.empty(len(symbols)); y=np.empty(len(symbols)); cx=0.5; cy=0.5
    C=np.array(corners,float)
    # x_i = sum_{j<=i} C[s_j] 2^{-(i-j+1)} ; use last 24 symbols
    X=np.zeros(len(symbols)); Y=np.zeros(len(symbols))
    for k in range(24):
        sh=np.roll(symbols,k); 
        X+=C[sh,0]*2.0**-(k+1); Y+=C[sh,1]*2.0**-(k+1)
    X=X[24:]; Y=Y[24:]
    H,_,_=np.histogram2d(Y,X,bins=size,range=[[0,1],[0,1]]); return H
tri=[(0.5,0.02),(0.02,0.98),(0.98,0.98)]
H=cgr(sym,tri,700); png("chaos_gap6.png",seq(np.log1p(H),stops=HOT))
rng=np.random.default_rng(1); sh=sym.copy(); rng.shuffle(sh)
H2=cgr(sh,tri,700); png("chaos_gap6_shuffled.png",seq(np.log1p(H2),stops=HOT))
print("occupied pixels real",int((H>0).sum()),"shuffled",int((H2>0).sum()))
# four corners: p mod 12
cls={1:0,5:1,7:2,11:3}; c4=np.array([cls[int(x)] for x in (P%12)])
sq=[(0,0),(1,0),(0,1),(1,1)]
for depth in (2,3,4):
    Hc=cgr(c4,sq,2**depth); r=Hc/Hc.mean(); print("mod12 depth",depth,"cell ratio min %.3f max %.3f"%(r.min(),r.max()))
Hc=cgr(c4,sq,64); png("chaos_mod12.png",up(div(Hc/Hc.mean()-1,0.6),10))
c4s=c4.copy(); rng.shuffle(c4s); Hs=cgr(c4s,sq,64); print("shuffled depth6 ratio min %.3f max %.3f"%((Hs/Hs.mean()).min(),(Hs/Hs.mean()).max()))
png("chaos_mod12_shuffled.png",up(div(Hs/Hs.mean()-1,0.6),10))
# pair table
T=np.zeros((4,4),int); np.add.at(T,(c4[:-1],c4[1:]),1); print("pair counts rows=from 1,5,7,11:\n",T, "\nratio to uniform:\n",np.round(T/T.mean(),3))
