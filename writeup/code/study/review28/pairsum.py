import numpy as np, math
from sympy import primerange
N=2_000_000
C2=0.6601618158468695739
S=np.zeros(N+1); S[2::2]=2*C2
for p in primerange(3,N+1):
    S[p::p]*= (p-1)/(p-2)
g=0.5772156649015329; A=1-g-math.log(2*math.pi)
h=np.arange(N+1,dtype=float)
cS=np.cumsum(S); chS=np.cumsum(h*S)
for H in [1e3,1e4,1e5,1e6,2e6-0.5, 123456.7]:
    k=int(math.floor(H)); 
    if H==k: k-=0  # h<=H, term (H-h) zero at h=H anyway
    s=H*cS[k]-chS[k]   # sum_{1<=h<=H}(H-h)S(h)
    t=s-0.5*H*H+0.5*H*math.log(H)
    print(H, "sum(H-h)S - H^2/2 + H lnH/2 , /(H/2):", t/(H/2), " A=",A, " A+1=",A+1)
