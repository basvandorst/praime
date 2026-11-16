import numpy as np
from lib import *
N=3*10**7; s=sieve(N); P=np.nonzero(s)[0]; print("primes to",N,len(P))
EG=0.5772156649015329
# ---------- rigidity: variance of theta in windows, real vs decoys
X0=10**7
def ratio_curve(ind,label):
    wt=np.log(np.arange(X0,N,dtype=np.float64)); wt*=ind[X0:N]; wt-=1.0
    c=np.empty(len(wt)+1); c[0]=0; np.cumsum(wt,out=c[1:]); L=len(wt); del wt; out=[]
    for k in range(1,19):
        H=2**k; st=max(8,H//4); idx=np.arange(0,L-H,st); d=c[idx+H]-c[idx]; lx=np.log(X0+idx+H/2)
        out.append((H,float(np.mean(d*d/(H*lx)))))
    print(label,"  ".join("%d:%.3f"%(H,r) for H,r in out[::2])); return out
real=ratio_curve(s,"real  ")
rng=np.random.default_rng(7)
def decoy(y):
    ind=np.zeros(N+1,dtype=bool); cand=np.ones(N+1,dtype=bool); cand[:2]=False; prod=1.0
    for p in P[P<=y]: cand[p::p]=False; prod*=1-1/p
    for a in range(X0,N,10**7):
        b=min(N,a+10**7); n=np.arange(a,b,dtype=np.float64); pr=(1/np.log(n))/prod
        ind[a:b]=cand[a:b]&(rng.random(b-a)<pr)
    return ind,prod
dc,_=decoy(1); cr=ratio_curve(dc,"cramer"); print("  decoy count",int(dc[X0:].sum()),"real",int(s[X0:].sum()))
d47,pr47=decoy(47); s47=ratio_curve(d47,"sieve47"); print("  decoy47 count",int(d47[X0:].sum()),"prod %.4f"%pr47)
lx=np.log(1.9e7)
print("MS prediction 1-(lnH+g+ln2pi)/lnx:","  ".join("%d:%.3f"%(2**k,1-(np.log(2**k)+EG+np.log(2*np.pi))/lx) for k in range(1,19,2)))
np.save("rigid.npy",np.array([[h for h,_ in real],[r for _,r in real],[r for _,r in cr],[r for _,r in s47]]))
del dc,d47
# ---------- race spectra (L-function zeros)
def race_spec(q,cls_pos,cls_neg,label,zeros):
    pm=P%q; step=np.zeros(N+1); 
    step[P[np.isin(pm,cls_pos)]]=1; step[P[np.isin(pm,cls_neg)]]-=1; walk=np.cumsum(step)
    u=np.linspace(np.log(200),np.log(N),2**14); x=np.exp(u); f=walk[np.floor(x).astype(int)]*u/np.sqrt(x)
    gam=np.linspace(0.3,30,1500); G=np.abs(np.exp(-1j*np.outer(gam,u))@((f-f.mean())*np.hanning(len(u))))
    pk=[gam[i] for i in range(1,len(gam)-1) if G[i]>G[i-1] and G[i]>G[i+1] and G[i]>0.3*G.max()]
    print(label,"mean of scaled race %.3f; spectrum peaks:"%f.mean(),np.round(pk,2),"\n    zeros:",zeros); return gam,G,u,f
g4=race_spec(4,[3],[1],"mod 4 :","6.0209 10.2437 12.9881 16.3426 18.2920 21.4506 23.2784 25.7288 28.3596")
g3=race_spec(3,[2],[1],"mod 3 :","8.0397 11.2492 15.7046 18.2620 20.4558 24.0594 26.5779 28.2182")
np.save("race_spec.npy",np.array([g4[0],g4[1],g3[1]]))
# ---------- race carpet mod 19
q=19; res=sorted({(a*a)%q for a in range(1,q)}); non=[a for a in range(1,q) if a not in res]; order=res+non
u=np.linspace(np.log(1000),np.log(N),900); x=np.floor(np.exp(u)).astype(int)
cnt=np.zeros((q,len(u))); pm=P%q
for a in range(1,q): cnt[a]=np.searchsorted(P[pm==a],x,side='right')
tot=cnt[1:].sum(axis=0); Ecar=np.array([(18*cnt[a]-tot)*u/np.sqrt(np.exp(u)) for a in order])
print("carpet mod 19: mean over residues rows %.3f (theory -1), non-residue rows %.3f (theory +1)"%(Ecar[:9].mean(),Ecar[9:].mean()))
print("  share of log-time residues ahead (sum of residue rows > 0): %.3f"%np.mean(Ecar[:9].sum(axis=0)>0))
img=div(Ecar,6.0); png("carpet19.png",np.repeat(img,22,axis=0))
