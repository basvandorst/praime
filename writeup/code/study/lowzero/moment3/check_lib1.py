import os, sys, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import tmlib as T
# (1) pairs: per-prime series A_p vs closed form 1 - x1^2 x2^2
ml = T.multilinear(2); z=[ml.var(0), ml.var(1)]
for Ds in [(1001,-1003),(1005,1068),(-1048,1068)]:
    t0=time.time()
    la = T.jet_logA(ml, Ds, z, pmax=10**5)
    P = T.primes_upto(10**5).astype(float)
    c1 = T.kron_primes(Ds[0],10**5).astype(float); c2=T.kron_primes(Ds[1],10**5).astype(float)
    good = (c1*c2)**2
    lp = np.log(P)
    # exact: sum_p log(1 - good p^{-2} e^{-2 ln p (z1+z2)}) as multilinear jet
    u = ml.exp((z[0]+z[1])[None,:]*(-2*lp)[:,None]) * (good/P**2)[:,None]
    ex = ml.log(ml.const(1.0) - u).sum(axis=0)
    # add the same tail
    for nn in []: pass
    tail = ml.compose_taylor((z[0]+z[1])*2.0, [-T.tail_S(n,10**5)*((-1.0)**n)/math.factorial(n) for n in range(ml.deg+1)])
    print(Ds, 'series', la, '\n   closed', ex+tail, ' diff', np.abs(la-ex-tail).max(), '%.1fs'%(time.time()-t0))
# (2) single L: P vs (phi(q)/q)(L + 2 gamma + 2 sum log p/(p-1))
for D in [1001, -1003, 1068]:
    tr = T.Triple((D,))
    t = np.array([100.0, 500.0, 1000.0])
    q = abs(D); fac = T.prime_factors(q)
    phi = q*np.prod([1-1/p for p in fac])
    ref = phi/q*(np.log(q*t/(2*math.pi)) + 2*T.EULER_GAMMA + 2*sum(math.log(p)/(p-1) for p in fac))
    print(D, tr.P(t), ref, tr.diag(t))
