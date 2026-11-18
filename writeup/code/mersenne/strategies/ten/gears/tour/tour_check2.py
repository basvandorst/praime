"""The two-cell wheel in plain Python, with sympy's Jacobi symbol (the C program and the page use their own reciprocity rule).
Wheel n: step c = smallest odd prime with Jacobi(c, n) = -1; cells (a, b) = (1, 1); a turn is (a, b) -> (a + c b, a + b) mod n,
that is multiplication by 1 + sqrt(c); home after n turns means (1, -1).
Checks: 1. odd n from 3 to 5,000: home exactly for the primes   2. every odd composite below 2e6: none home
        3. a sample of primes of 8 to 64 bits: all home          4. the candidates below 1e10 and a sample up to 2^64: none home"""
import time
from math import isqrt
from sympy import jacobi_symbol, isprime, primerange
t0=time.time()
SP=list(primerange(3,200))
def step(n):
    for c in SP:
        if n==c: continue
        j=jacobi_symbol(c,n)
        if j==0: return -c
        if j==-1: return c
    return 0
def home(n):
    c=step(n)
    if c<=0: return False
    a,b=1,1
    for bit in bin(n)[3:]:
        a,b=(a*a+c*b*b)%n,(2*a*b)%n
        if bit=='1': a,b=(a+c*b)%n,(a+b)%n
    return a==1 and b==n-1
bad=[n for n in range(3,5001,2) if home(n)!=isprime(n)]
print("1. odd n from 3 to 5,000: mistakes %s; lamps lit %d (with 2 let through: %d)"%(bad,sum(1 for n in range(3,5001,2) if home(n)),1+sum(1 for n in range(3,5001,2) if home(n))))
N=2*10**6; sv=bytearray([1])*(N+1)
for i in range(2,isqrt(N)+1):
    if sv[i]: sv[i*i::i]=bytearray(len(sv[i*i::i]))
comp=[n for n in range(9,N,2) if not sv[n]]
print("2. odd composites below 2e6: %d tested, came home: %d; turned away because a step divides them: %d; squares or other numbers with no step: %d"%(len(comp),sum(1 for n in comp if home(n)),sum(1 for n in comp if step(n)<0),sum(1 for n in comp if step(n)==0)),flush=True)
P=[int(x) for x in open('primes_sample.txt')][::57]
print("3. primes of 8 to 64 bits: %d, came home: %d"%(len(P),sum(1 for p in P if home(p))))
L=[int(x) for x in open('psp2_to1e10.txt')]; S=[int(x) for x in open('sample64.txt')]
print("4. candidates below 1e10: %d, came home: %d, no step: %s; sample up to 2^64: %d, came home: %d"%(len(L),sum(1 for n in L if home(n)),[n for n in L if step(n)==0],len(S),sum(1 for n in S if home(n))))
use={}
for n in S: use[abs(step(n))]=use.get(abs(step(n)),0)+1
print("   steps used on the sample:",sorted(use.items()))
print("time %.0f s"%(time.time()-t0))
