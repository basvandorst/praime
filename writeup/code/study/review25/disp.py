import numpy as np
T=np.load('twins.npy'); W=T[T<2*10**9]
BANDS = [(lo, lo + 24) for lo in range(6, 600, 30)] + [(lo, lo + 54) for lo in range(606, 1200, 60)]
rng=np.random.default_rng(1)
g=np.diff(W)
for k in (1,2,3):
    d=W[k:]-W[:-k]; 
    nb=400; idx=np.array_split(np.arange(len(d)),nb)
    ratios=[]
    for lo,hi in BANDS:
        ind=((d>=lo)&(d<=hi)&(d<=1200)).astype(float)
        c=np.array([ind[i].sum() for i in idx]); p=ind.mean()
        ratios.append(c.var(ddof=1)/(c.mean()*(1-p)))
    ratios=np.array(ratios)
    # covariance between adjacent bands' block counts is ignored; report the mean dispersion
    print(k,'mean dispersion %.3f (sd over bands %.3f)'%(ratios.mean(),ratios.std()))
# null check: shuffle gaps (independent) and compute dispersion -> should be ~1 for k=1 by construction
