# control: adjacent-window correlation of Var(s) and <r~> for GUE bulk spectra (Dumitriu-Edelman tridiagonal, beta=2)
import numpy as np, math
from scipy.linalg import eigvalsh_tridiagonal
rng=np.random.default_rng(1)
N=4000; W=400; V=[];R=[]
for it in range(300):
    d=rng.standard_normal(N)*math.sqrt(2); e=np.sqrt(rng.chisquare(2*np.arange(N-1,0,-1)))
    lam=eigvalsh_tridiagonal(d/math.sqrt(2),e/math.sqrt(2)); Rr=math.sqrt(2*N)
    F=0.5+(lam*np.sqrt(np.maximum(Rr**2-lam**2,0)))/(math.pi*Rr**2)+np.arcsin(np.clip(lam/Rr,-1,1))/math.pi
    u=N*F; c=N//2; u=u[c-2*W:c+2*W+1]
    v=[];r=[]
    for j in range(4):
        s=np.diff(u[j*W:(j+1)*W+1]); v.append(np.var(s)); rr=np.minimum(s[1:],s[:-1])/np.maximum(s[1:],s[:-1]); r.append(rr.mean())
    V.append(v);R.append(r)
V=np.array(V);R=np.array(R)
print('GUE Var adjacent-window corr', np.round(np.corrcoef(V.T),2)[0:3,0:4]); print('GUE r~ ', np.round(np.corrcoef(R.T),2)[0:3,0:4])
