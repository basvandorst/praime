# grid in r and Z(ir) = |zeta(1+ir)|^2 A(ir), A(a) = prod_p (1-p^{-1-a})(1-2/p+p^{-1-a})/(1-1/p)^2
import numpy as np, mpmath as mp, sympy, sys, time
mp.mp.dps=20
Pmax=int(sys.argv[1]); Rmax=float(sys.argv[2])
# nodes: r=u^2 on [0,0.5] (u in [0,sqrt .5]) with 40 GL nodes; then panels of width 0.1, 8 nodes
x,w=np.polynomial.legendre.leggauss(40); a=np.sqrt(0.5)
u=(x+1)/2*a; wu=w/2*a; r0=u**2; w0=wu*2*u
x8,w8=np.polynomial.legendre.leggauss(8)
edges=np.arange(0.5,Rmax+1e-9,0.1)
r1=((x8[None,:]+1)/2*0.1+edges[:-1,None]).ravel(); w1=np.tile(w8/2*0.1,len(edges)-1)
r=np.concatenate([r0,r1]); wr=np.concatenate([w0,w1])
t=time.time()
z2=np.array([abs(complex(mp.zeta(mp.mpc(1,ri))))**2 for ri in r])
print('zeta done',time.time()-t,len(r))
logA=np.zeros(len(r),complex)
ps=np.array(list(sympy.primerange(2,Pmax)),float)
for i in range(0,len(ps),500):
    p=ps[i:i+500][:,None]
    pa=p**-1*np.exp(-1j*r[None,:]*np.log(p))
    logA+=np.sum(np.log1p(-pa)+np.log((1-2/p)+pa)-2*np.log1p(-1/p),axis=0)
print('A done',time.time()-t)
np.savez(sys.argv[3],r=r,w=wr,Z=z2*np.exp(logA))
