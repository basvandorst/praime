import sys, math, json
import numpy as np
from sympy import primerange

def tables(N):
    # F[n] = prod_{p|n, p>=5} r2(p) in logs; Gc[n] = prod_{p|n,p>=5} r1/r2^3
    lF = np.zeros(N+1); lG = np.zeros(N+1)
    for p in primerange(5, N+1):
        r1 = (1-1/p)/(1-3/p); r2 = (1-2/p)/(1-3/p)
        lF[p::p] += math.log(r2); lG[p::p] += math.log(r1) - 3*math.log(r2)
    return lF, lG

def f3(nu):
    return (1-nu/3)/(8/27)

def mi(lam, mult=24):
    G = int(mult*lam)//2*2
    lF, lG = tables(G)
    ev = np.arange(2, G+1, 2, dtype=np.int64)
    n = len(ev)
    m1 = np.zeros(n); m2 = np.zeros(n); tot = 0.0; Hj = 0.0
    for i in range(n):
        a = ev[i]
        b = ev[: (G - a)//2]
        if len(b) == 0: continue
        s = a + b
        # p = 3
        ra = a % 3; rb = b % 3
        nu = np.where((ra == 0) & (rb == 0), 1, np.where((ra != 0) & (rb == ra), 3, 2))
        f = f3(nu)
        g = np.gcd(a, b)
        lw = lF[a] + lF[b] + lF[s] + lG[g] - s/lam
        w = f*np.exp(lw)
        k = w > 0
        w = w[k]; bi = (b[k]//2 - 1)
        tot += w.sum(); Hj += float(np.sum(w*np.log(w)))
        m1[i] += w.sum(); np.add.at(m2, bi, w)
    H = lambda x: -float(np.sum((x[x>0]/tot)*np.log2(x[x>0]/tot)))
    Hjoint = -(Hj/tot - math.log(tot))/math.log(2)
    return H(m1) + H(m2) - Hjoint

def Iinf():
    s = 0.0
    for p in primerange(3, 10**6):
        q = 1-1/p
        fA = (1-2/p)/q**3; fB = (1-3/p)/q**3
        s += (p-1)/p**2*fA*math.log(q/(1-2/p))
        if p > 3: s += (p-1)*(p-2)/p**2*fB*math.log((1-3/p)*q/(1-2/p)**2)
    return s/math.log(2)

if __name__ == '__main__':
    I8 = Iinf(); print('I_inf bits', I8)
    for arg in sys.argv[1:]:
        lam, mult = map(float, arg.split(':'))
        I = mi(lam, mult)
        print(lam, mult, 'I', I, 'diff', I - I8, 'lam*diff', lam*(I - I8), flush=True)
