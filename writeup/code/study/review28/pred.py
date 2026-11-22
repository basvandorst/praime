import numpy as np, math
from sympy import primerange
N=2_000_000
C2=0.6601618158468695739
S=np.zeros(N+1); S[2::2]=2*C2
for p in primerange(3,N+1): S[p::p]*=(p-1)/(p-2)
h=np.arange(N+1,dtype=float); cS=np.cumsum(S); chS=np.cumsum(h*S)
g=0.5772156649015329; A=1-g-math.log(2*math.pi)
def PS(H):  # sum_{1<=h<H}(H-h)S(h), vectorised
    H=np.asarray(H,float); out=0.5*H*H-0.5*H*np.log(np.maximum(H,1))+0.5*A*H
    m=H<N; k=np.floor(H[m]).astype(int); out[m]=H[m]*cS[k]-chS[k]; return out
def Rpred(X,d,model):
    # integrate V over x in [1,X] (log grid, trapezoid)
    u=np.linspace(0,math.log(X),400001); x=np.exp(u); H=d*x
    diag=(x+H)*np.log(x+H)-(x+H)-(x*np.log(x)-x)
    if model=='HL': V=diag+2*PS(H)-H*H
    else: V=diag-H   # Cramer Bernoulli: sum (ln n - 1)
    f=V*x; I=np.sum(0.5*(f[1:]+f[:-1])*np.diff(u))
    return I/(0.5*d*X*X)
for X in [1e8,1e9]:
    for d in [1e-1,1e-2,1e-3,1e-4,1e-5,1e-6]:
        print('X=%.0e d=%.0e  HL-pred R %.4f  ln(1/d)+B %.4f  diff %+.4f   Cramer %.3f (lnX-1.5=%.3f)'%(X,d,Rpred(X,d,'HL'),math.log(1/d)+A,Rpred(X,d,'HL')-math.log(1/d)-A,Rpred(X,d,'C'),math.log(X)-1.5))
