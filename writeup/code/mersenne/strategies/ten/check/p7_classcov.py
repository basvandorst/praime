import numpy as np
from lib import sieve
N=3*10**7; X0=10**7; s=sieve(N); n=np.arange(X0,N); ln=np.log(n.astype(np.float64))
def run(q,units,T):
    phi=len(units); cs=[]
    for a in units:
        w=np.where(s[X0:N]&(n%q==a),ln,0.0)-1.0/phi; c=np.empty(len(w)+1); c[0]=0; np.cumsum(w,out=c[1:]); cs.append(c)
    print("mod",q,"units",units)
    for k in (8,10,12,14,16):
        H=2**k; idx=np.arange(0,len(n)-H,H//4); D=np.array([c[idx+H]-c[idx] for c in cs]); D-=D.mean(axis=1,keepdims=True)
        C=(D@D.T)/D.shape[1]*phi*phi/H
        out=[]
        for r,t in T.items():
            vals=[C[i,j] for i,a in enumerate(units) for j,b in enumerate(units) if a!=b and (b*pow(a,-1,q))%q==r]
            out.append("r=%d: %.2f (pred %.2f)"%(r,np.mean(vals),t))
        Lm=np.log(2e7/H)-0.5772-np.log(2*np.pi)
        print("  H=%6d  off-diagonal, by ratio b/a:  %s   | diagonal %.2f (pred %.2f)  windows %d"%(H,"  ".join(out),np.mean(np.diag(C)),phi*Lm+sum(-t for t in T.values()),D.shape[1]))
l=np.log
run(12,[1,5,7,11],{5:l(4)-l(3)-l(12),7:-l(4)+l(3)-l(12),11:-l(4)-l(3)+l(12)})
run(8,[1,3,5,7],{3:-l(4)-l(8)+l(8),5:l(4)-l(8)-l(8),7:-l(4)+l(8)-l(8)})
run(5,[1,2,3,4],{2:-l(5),3:-l(5),4:-l(5)})
