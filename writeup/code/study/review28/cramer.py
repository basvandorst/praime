import numpy as np, math
import sys; rng=np.random.default_rng(int(sys.argv[1])); X=10**7
n=np.arange(3,int(1.11*X)+2); keep=rng.random(len(n))<1/np.log(n)
w=np.zeros(int(1.11*X)+2); w[n[keep]]=np.log(n[keep]); cum=np.cumsum(w)
x=np.arange(1,X)+rng.random(X-1)  # one random point per unit
for d in [1e-4,1e-5]:
    f=cum[np.floor(x*(1+d)).astype(int)]-cum[np.floor(x).astype(int)]-d*x
    R=np.dot(f,f)/(0.5*d*X*X); print(d,R,'lnX-1.5=',math.log(X)-1.5,'lnX-0.5=',math.log(X)-0.5)
