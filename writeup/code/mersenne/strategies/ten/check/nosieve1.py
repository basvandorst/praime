import numpy as np, math
from lib import sieve
import cypari2; pari=cypari2.Pari()
MEXP=[2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,136279841]
N=10**6; s=sieve(N); P=np.nonzero(s)[0]
def v2(x):
    k=0
    while x%2==0: x//=2; k+=1
    return k
print("== (a) primes in rows by the power of 2 in the class number h(-4p), p <= 1e6")
h=np.array([int(pari.qfbclassno(-4*int(p))) for p in P[1:]]); v=np.array([v2(int(x)) for x in h]); tot=len(v)
print(" rows 0..8 shares:"," ".join("%.4f"%np.mean(v==k) for k in range(9)),"  (halving: 0.5 0.25 0.125 ...)")
pp=P[1:]; print(" row 0 = p=3 mod 4:",bool(np.all((v==0)==(pp%4==3)))," row 1 = p=5 mod 8:",bool(np.all((v==1)==(pp%8==5)))," rows>=2 = p=1 mod 8:",bool(np.all((v>=2)==(pp%8==1))))
def is_x2_32y2(p):
    y=1
    while 32*y*y<p:
        r=p-32*y*y; x=math.isqrt(r)
        if x*x==r: return True
        y+=1
    return False
m8=pp%8==1; t=np.array([is_x2_32y2(int(p)) for p in pp[m8]]); print(" rows>=3 exactly the primes x^2+32y^2:",bool(np.all((v[m8]>=3)==t)), " count",int(t.sum()),"of",int(m8.sum()))
print(" largest power:",int(v.max()),"at p =",int(pp[np.argmax(v)]),"h =",int(h[np.argmax(v)]))
print(" composites n<=20000 (odd, squarefree): min v2(h(-4n)) by number of prime factors:")
from sympy import factorint
by={}
for n in range(3,20001,2):
    f=factorint(n)
    if any(e>1 for e in f.values()): continue
    by.setdefault(len(f),[]).append(v2(int(pari.qfbclassno(-4*n))))
for k in sorted(by): a=np.array(by[k]); print("   %d factors: min %d, mean %.2f, n=%d"%(k,a.min(),a.mean(),len(a)))
ev=[(p,v2(int(pari.qfbclassno(-4*p)))) for p in MEXP[1:]]
from collections import Counter
c=Counter(x for _,x in ev); print(" the 51 odd exponents by row:",dict(sorted(c.items()))," expected 25.5, 12.75, 6.4, 3.2, ...; deepest:",max(ev,key=lambda t:t[1]))
print("== (c) square/non-square walks: sum_{j<=n} (j/p)")
cnt=0; nonneg=0; ok=0; neg_first=None
for p in P[(P%4==3)&(P>3)&(P<20000)]:
    p=int(p); res=np.zeros(p,dtype=np.int8)-1; res[(np.arange(1,p,dtype=np.int64)**2)%p]=1; res[0]=0
    w=np.cumsum(res[1:]); mid=int(w[(p-1)//2-1]); hh=int(pari.qfbclassno(-p)); two=1 if p%8==7 else -1
    cnt+=1; ok+= (mid==(2-two)*hh); 
    if w.min()>=0: nonneg+=1
    elif neg_first is None: neg_first=p
print(" primes 3 mod 4 below 20000: %d; height at p/2 equals (2-(2/p)) h(-p) for %d; walk never below 0 for %d (%.1f%%); first prime whose walk dips below 0: %s"%(cnt,ok,nonneg,100*nonneg/cnt,neg_first))
print("== (d) Lehmer signs of the 51 odd exponents")
S=[1,1,-1,1,-1,-1,1,1,-1,-1,1,-1,-1,-1,1,-1,-1,1,-1,-1,1,-1,1,-1,1,-1,1,1,1,-1,1,-1,-1,1,1,1,1,1,1,-1,1,-1,-1,1,-1,1,-1,-1,1,-1,1]
E=MEXP[1:]; assert len(S)==len(E)
S=np.array(S); print(" plus %d, minus %d; runs %d (expected %.1f +- %.1f); longest run %d"%((S==1).sum(),(S==-1).sum(),1+int((np.diff(S)!=0).sum()),1+2*(S==1).sum()*(S==-1).sum()/len(S),math.sqrt(2*(S==1).sum()*(S==-1).sum()*(2*(S==1).sum()*(S==-1).sum()-len(S))/(len(S)**2*(len(S)-1))),max(len(list(g)) for _,g in __import__('itertools').groupby(S))))
for m in (8,3,5,12):
    print("  by p mod %d:"%m, {int(r):(int(((np.array(E)%m==r)&(S==1)).sum()),int(((np.array(E)%m==r)&(S==-1)).sum())) for r in sorted(set(int(e%m) for e in E))})
from scipy.stats import fisher_exact
print("  split p = 1 (mod 8) against p = 5 (mod 8), Fisher exact two-sided: p = %.3f"%fisher_exact([[7,13],[9,3]])[1])
