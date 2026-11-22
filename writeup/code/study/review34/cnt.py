import sys; sys.path.insert(0,'/home/bas/praime/writeup/code/study')
import numpy as np, r100, r98, math
from scipy.linalg import eigh_tridiagonal
c=[];sp=[]
for s in range(100,108):
    g=r100.gue_stationary(s); c.append(len(g))
print('their stationary counts', c, np.mean(c), 'expected', float(np.interp(r98.Z[-1], r98.t, r98.Nmean))-float(np.interp(r98.Z[0], r98.t, r98.Nmean)))
# sample-R vs exact radius
n=30000; rs=[]
for s in range(100,106):
    rng=np.random.default_rng(s); dgl=rng.normal(0,math.sqrt(2),n); off=np.sqrt(rng.chisquare(2*np.arange(n-1,0,-1)))
    ev=eigh_tridiagonal(dgl,off,eigvals_only=True)/math.sqrt(2*n); rs.append(np.max(np.abs(ev)))
print('sample R / exact 2', np.array(rs)/2)
