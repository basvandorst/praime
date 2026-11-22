import numpy as np, math, time
from rq import *
print('Nbar(g1), Nbar(gN)', Nbar(Z[0]), Nbar(Z[-1]), 'Ninv err', np.max(np.abs(Ninv(Nbar(Z))-Z)))
print('mean N(g_n)-n+1/2 :', np.mean(Nbar(Z)-np.arange(1,10001)+0.5))
t0=time.time(); r=rms(Z,(0.01,0.005,0.003)); print('actual', r, time.time()-t0)
u=gue_unf(1); s=np.diff(u[12000:24000]); print('gue mean spacing', s.mean(), 'var', s.var())
uz=Nbar(Z); sz=np.diff(uz); print('zeros unfolded spacing mean/var', sz.mean(), sz.var())
for lo_ in (0,):
  g=gue_set(1, Nbar(Z[0]), Nbar(Z[-1])); print(len(g), g[0], g[-1]); t0=time.time(); print('gue1', rms(g,(0.01,0.005,0.003)), time.time()-t0)
