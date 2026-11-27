# Numerical spot check of Lemma 1 (projection of the high congruences to level n), exact fractions.
# For a system of congruences with distinct moduli m | 2^a*45, split low (v2(m) <= n) / high (v2(m) > n).
#   y(d, r mod 2^n d) = sum over high congruences (2^k d, rho) with rho = r (mod 2^n d) of 2^(n-k).
# Checked for every b mod 2^n*45 that lies in no low class:
#   (number of lifts of b mod 2^a*45 covered by a high class) <= 2^(a-n) * sum_d y(d, b mod 2^n d)
# and  sum_r y(d,r) <= 1 - 2^(n-a).   If the system covers, the left side is 2^(a-n), so sum_d y >= 1.
import random
from fractions import Fraction as F
DS=(1,3,5,9,15,45)
def v2(m):
    k=0
    while m%2==0: m//=2; k+=1
    return k
def check(system,a,n):
    L=2**a*45; Ln=2**n*45
    low=[(m,r) for m,r in system if v2(m)<=n]; high=[(m,r) for m,r in system if v2(m)>n]
    y={}
    for m,r in high:
        k=v2(m); d=m>>k
        y[(d,r%(2**n*d))]=y.get((d,r%(2**n*d)),0)+F(1,2**(k-n))
    for d in DS:
        assert sum(v for (dd,_),v in y.items() if dd==d) <= 1-F(1,2**(a-n))
    worst=None
    for b in range(Ln):
        if any(b%m==r%m for m,r in low): continue
        lifts=[b+Ln*j for j in range(2**(a-n))]
        cov=sum(1 for x in lifts if any(x%m==r%m for m,r in high))
        s=sum(y.get((d,b%(2**n*d)),0) for d in DS)
        assert cov <= 2**(a-n)*s, (system,b,cov,s)
        if cov==2**(a-n): assert s>=1
    return True
random.seed(1)
tests=0
for trial in range(400):
    a=random.choice([3,4,5]); n=random.choice([0,1,2])
    mods=[2**k*d for k in range(a+1) for d in DS if 2**k*d>=6]
    system=[(m,random.randrange(m)) for m in mods if random.random()<0.8]
    check(system,a,n); tests+=1
# a genuine covering system: 0 mod 2, 0 mod 3, 1 mod 4, 5 mod 6, 7 mod 12  (moduli | 2^2*45)
cov=[(2,0),(3,0),(4,1),(6,5),(12,7)]
assert all(any(x%m==r for m,r in cov) for x in range(180))
for n in (0,1): check(cov,2,n); tests+=1
print('Lemma 1 spot check passed on',tests,'systems (400 random with a in {3,4,5}, n in {0,1,2}; one real covering at n = 0, 1)')
