import numpy as np, math, sys
sys.path.insert(0, '/home/bas/praime/writeup/code/study')
from gen import *
import r98
print('Nbar(g1)', N1, 'NT', NT, 'Ninv check', Ninv(np.array([N1, 100.5, NT])), Z[0], Z[-1])
print('Poisson r98 lower edge t at N=0:', float(Ninv(np.array([0.0]))[0]))
# compare r98.gue(0) with my reconstruction (sample radius)
g98 = r98.gue(0)
u, ev = gue_unfolded(0, radius='sample')
mid = 10000; us = u[mid:mid+10000] - u[mid]
mine = Ninv(N1 + us)
print('max |r98.gue(0) - mine|', np.max(np.abs(g98 - mine)))
# radius bias and spacing stats
for s in range(4):
    ue, ev = gue_unfolded(s, radius='exact'); usamp, _ = gue_unfolded(s, radius='sample')
    R = np.max(np.abs(ev))
    se = np.diff(ue[mid:mid+10000]); ss = np.diff(usamp[mid:mid+10000])
    print('seed', s, 'sample radius', R, 'mean spacing exactR', se.mean(), 'sampleR', ss.mean(), 'var', ss.var(), se.var()/se.mean()**2)
# zeros' normalised spacing variance
sp = np.diff(Nbar(Z)); print('zeros unfolded spacing mean', sp.mean(), 'var', sp.var())
# dense GUE check of spacing variance (central)
rng = np.random.default_rng(7); vs = []
for k in range(10):
    n = 1500; A = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n)); H = (A + A.conj().T)/2
    e = np.linalg.eigvalsh(H)/math.sqrt(n)   # radius 2
    x = np.clip(e/2, -1, 1); uu = n*(0.5 + (x*np.sqrt(1-x*x) + np.arcsin(x))/math.pi)
    d = np.diff(uu[500:1000]); vs.append(d.var()/d.mean()**2)
print('dense GUE central spacing var', np.mean(vs), np.std(vs)/math.sqrt(len(vs)))
