# independent check of the sine-base closure cells at L=9.4: GUE bulk (tridiagonal) spectra + p's wave (random phase)
import numpy as np, math
from scipy.linalg import eigvalsh_tridiagonal
rng=np.random.default_rng(77); N=3000; L=9.4
def stats(x):
    s=np.diff(x); r=np.minimum(s[1:],s[:-1])/np.maximum(s[1:],s[:-1]); return np.var(s),r.mean()
acc={p:[] for p in (2,3,5,13)}; V0=[]
for it in range(250):
    d=rng.standard_normal(N)*math.sqrt(2); e=np.sqrt(rng.chisquare(2*np.arange(N-1,0,-1)))
    lam=eigvalsh_tridiagonal(d/math.sqrt(2),e/math.sqrt(2)); R=2*math.sqrt(np.mean(lam**2))
    F=0.5+(lam*np.sqrt(np.maximum(R**2-lam**2,0)))/(math.pi*R**2)+np.arcsin(np.clip(lam/R,-1,1))/math.pi
    y=N*F; y=y[N//2-700:N//2+700]
    v0,r0=stats(y); V0.append((v0,r0,np.diff(y).mean()))
    for p in acc:
        k=2*math.pi*math.log(p)/L; ph=rng.uniform(0,2*math.pi); x=y.copy()
        for _ in range(60): x=y-np.angle(1-p**-0.5*np.exp(1j*(k*x+ph)))/math.pi
        v,r=stats(x); acc[p].append((v0-v,r0-r))   # missing minus present
for p,a in acc.items():
    a=np.array(a); print('p=%2d  dVar %.5f +- %.5f   d<r~> %+.5f +- %.5f'%(p,a[:,0].mean(),a[:,0].std()/math.sqrt(len(a)),a[:,1].mean(),a[:,1].std()/math.sqrt(len(a))))
print('frozen law at L=9.4: dVar -0.01260 -0.01438 -0.01360 (13: -0.00872); d<r~> +0.00258 +0.00378 +0.00495 (13: +0.00574)')
print('base Var, r, mean s:',np.round(np.mean(V0,0),4))
