# zeta(1+i e), (zeta'/zeta)'(1+i e) on the eps grid, and A_zeta(i e) (Conrey-Snaith) for the deg-1 template
import mpmath as mp, numpy as np
mp.mp.dps=30
EPS=np.arange(1,1641)*0.005
Z0=[];ZD=[]
for e in EPS:
    s=mp.mpc(1,e); z=mp.zeta(s); z1=mp.zeta(s,derivative=1); z2=mp.zeta(s,derivative=2)
    Z0.append(complex(z)); ZD.append(complex(z2/z-(z1/z)**2))
np.savez('zeta_cache.npz',eps=EPS,Z0=np.array(Z0),ZD=np.array(ZD))
print('done', ZD[0], -1/EPS[0]**2)
