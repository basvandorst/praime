"""Beyond the list: numbers above 2^64 that fool the doubling wheel, built on purpose, on both kept wheels.
Family 1: Chernick's Carmichael numbers (6k+1)(12k+1)(18k+1) with all three factors prime.
Family 2: n = p q with q = k (p - 1) + 1 prime (k = 2..6) and 2^(n-1) = 1 (mod n).
Each n above 2^64 is put on the smallest Pascal ring it does not fit and on its two-cell wheel (plain Python, ref.py)."""
import random, time
from math import gcd
import gmpy2
from sympy import jacobi_symbol, primerange
from ref import ring_test
t0=time.time()
SP=list(primerange(3,200))
def misfit(n):
    for r in range(2,400):
        if gcd(n,r)==1 and (n*n)%r!=1: return r
    return 0
def duo(n):
    c=0
    for q in SP:
        j=jacobi_symbol(q,n)
        if j==0: return False
        if j==-1: c=q; break
    if not c: return False
    a,b=1,1
    for bit in bin(n)[3:]:
        a,b=(a*a+c*b*b)%n,(2*a*b)%n
        if bit=='1': a,b=(a+c*b)%n,(a+b)%n
    return a==1 and b==n-1
B=2**64
F1=[]
k=240000
while len(F1)<4000 and k<6000000:
    k+=1
    if gmpy2.is_prime(6*k+1) and gmpy2.is_prime(12*k+1) and gmpy2.is_prime(18*k+1):
        n=(6*k+1)*(12*k+1)*(18*k+1)
        if n>B: F1.append(n)
random.seed(64)
F2=[]
while len(F2)<6000:
    p=int(gmpy2.next_prime(random.randrange(3*10**9,10**12)))
    for kk in (2,3,4,5,6):
        q=kk*(p-1)+1
        if gmpy2.is_prime(q):
            n=p*q
            if n>B and pow(2,n-1,n)==1: F2.append(n)
for name,F in (("Chernick Carmichael numbers",F1),("p (k (p - 1) + 1), k = 2..6, base-2 pseudoprimes",F2)):
    assert all(pow(2,n,n)==2 and not gmpy2.is_prime(n) for n in F)
    rings={}; h1=0; h2=0
    for n in F:
        r=misfit(n); rings[r]=rings.get(r,0)+1
        if ring_test(n,r): h1+=1
        if duo(n): h2+=1
    print("%s: %d numbers from %.3g to %.3g; home on Pascal's wheel: %d; home on the two-cell wheel: %d; rings used: %s"%(name,len(F),min(F),max(F),h1,h2,sorted(rings.items())),flush=True)
print("time %.0f s"%(time.time()-t0))
