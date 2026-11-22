import numpy as np, math, sympy
from rq import *
# prime-power table Lambda(n)^2/n^3
N=10**5; L=np.zeros(N+1)
for p in sympy.primerange(2,N+1):
    q=p
    while q<=N: L[q]=math.log(p); q*=p
n=np.arange(N+1); w=np.zeros(N+1); w[2:]=L[2:]**2/n[2:]**3.0; lnn=np.zeros(N+1); lnn[2:]=np.log(n[2:])
print('sum Lambda^2/n^3 =', w.sum())
DS=(0.01,0.005,0.003); W=rows(Z,DS)
TH=np.log(Z/(2*math.pi))
S=np.array([w[(lnn>0)&(lnn<T)].sum() for T in TH])
I=0.25-np.exp(-2*TH)*(2*TH+1)/4
for d in DS:
    a=c2(Z,d); Q=4*W[d].sum()  # Q over both signs = 2*sum w (rows over positive) ... check
    Q=2*W[d].sum()
    dQ=np.sum(a*a*(S-I)/TH)   # Q_z - Q_GUE predicted, both signs
    rms=2/d*math.sqrt(2*Q)
    dQhi=np.sum((a*a*(S-I)/TH)[Z>500])
    print(d,'Q',Q,'rms',rms,'pred dRMS',rms*dQ/(2*Q),' share from zeros >500:',dQhi/dQ, ' >100:', np.sum((a*a*(S-I)/TH)[Z>100])/dQ)
