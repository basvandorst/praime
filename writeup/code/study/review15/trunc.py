# per-class exact Pi_a (inclusion-exclusion over admissible interior twin events) vs order-2 Bahadur, sampled classes
import math, sys, random, numpy as np
import law
from sympy import primerange
d=int(sys.argv[1]); nsamp=int(sys.argv[2]); random.seed(1)
lam=law.lam; SMALL=[2,3,5,7,11,13]; PP=30030
rho0=1/(lam*math.prod(1-1/q for q in SMALL))
QE=[q for q in primerange(17, d+5)]
GENq=np.array([q for q in law.PR if q > QE[-1]])
gtk={k: float(np.sum(np.log1p(-k/GENq)-k*np.log1p(-1/GENq))) for k in range(4, 200, 2)}
B=[0,2,d,d+2]; H=list(range(6,d,6))
def masks(pts): return [sum(1<<(x%q) for x in set(pts)) for q in QE]
MB=masks(B); MH=[masks([h,h+2]) for h in H]
lnq1=[math.log(1-1/q) for q in QE]
def lnS(m,k):
    s=0.0
    for q,mm,l1 in zip(QE,m,lnq1):
        nu=bin(mm).count('1')
        if nu>=q: return -math.inf
        s+=math.log(1-nu/q)-k*l1
    return s+gtk[k]
lB=lnS(MB,4)
ls=[lnS([a|b for a,b in zip(MB,MH[i])],6)-lB for i in range(len(H))]
p=[rho0**2*math.exp(x) for x in ls]; w=[x/(1-x) for x in p]
R={}
for i in range(len(H)):
    for j in range(i+1,len(H)):
        R[i,j]=math.expm1(lnS([a|b|c for a,b,c in zip(MB,MH[i],MH[j])],8)+lB-(ls[i]+lB)-(ls[j]+lB))
cls=[a for a in range(PP) if all(math.gcd(x,PP)==1 for x in (a,a+2,a+d,a+d+2))]
smp=random.sample(cls, min(nsamp,len(cls)))
num=0.0; den=0.0; den1=0.0
for a in smp:
    idx=[i for i,h in enumerate(H) if math.gcd(a+h,PP)==1 and math.gcd(a+h+2,PP)==1]
    L=sum(math.log1p(-p[i]) for i in idx); L1=L
    L+=sum(w[i]*w[j]*R[i,j] for x,i in enumerate(idx) for j in idx[x+1:])
    tot=0.0
    def rec(s, m, k, sign):
        global tot
        l=lnS(m,k)
        if l==-math.inf: return
        tot+=sign*math.exp(l-lB)*rho0**(k-4)
        for t in range(s,len(idx)):
            i=idx[t]; rec(t+1,[a|b for a,b in zip(m,MH[i])],k+2,-sign)
    rec(0,MB,4,1)
    num+=math.exp(L); den+=tot; den1+=math.exp(L1)
print('d=%d classes sampled %d of %d; order2/exact-1 = %.3e ; order1/exact-1 = %.3e'%(d,len(smp),len(cls),num/den-1,den1/den-1),flush=True)
