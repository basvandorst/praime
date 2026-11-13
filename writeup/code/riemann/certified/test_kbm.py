import math, numpy as np
from scipy.integrate import quad
from scipy.special import jv
from kbm import K, Khat, P, K0, majorant
from sympy import factorint
for (a,beta,m) in [(1.9,10.0,4),(1.9,14.0,3),(20.0,12.0,4)]:
    for w in [0.0, 1.0, 0.5*beta/a, 0.99*beta/a, 1.01*beta/a, 2*beta/a, 10*beta/a]:
        num = quad(lambda v: float(K(v,a,beta,m))*math.cos(w*v), -a, a, limit=500, epsabs=1e-13)[0]
        cf = float(Khat(w,a,beta,m))
        print(f'a={a} beta={beta} m={m} w={w:.4f}: closed {cf:.10e} quad {num:.10e} rel {abs(cf-num)/max(1e-300,abs(Khat(0,a,beta,m))):.1e}')
    x=1009.0
    num = quad(lambda v: float(K(v,a,beta,m))*math.exp(v/(2*x)), -a, a, limit=500)[0]
    print(f'  P({x}) closed {P(x,a,beta,m):.12e} quad {num:.12e}')
# majorant check: M(w) >= |Khat(w')| for w' >= w on a grid
a,beta,m=1.9,12.0,4
ws=np.linspace(1.0001*beta/a, 60*beta/a, 200000); kh=np.abs(Khat(ws,a,beta,m)); M=majorant(ws,a,beta,m)
sup=np.maximum.accumulate(kh[::-1])[::-1]
print('majorant valid on grid:', bool(np.all(M>=sup*(1-1e-12))), ' worst ratio sup/M', float(np.max(sup/M)))
# Landau bound sanity
y=np.linspace(0.01,500,200000); 
for nu in [3.5,4.5]: print('Landau', nu, float(np.max(np.abs(jv(nu,y))*y**(1/3))))
