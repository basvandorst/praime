"""Both kept wheels in plain Python on random odd numbers (seed 2026): the wheel must be home exactly for the primes.
300,000 numbers between 1e9 and 2e10 for the two-cell wheel, 20,000 of them also on Pascal's wheel; 20,000 numbers of 40 to 64 bits on both."""
import random, time
from math import gcd
from sympy import isprime
from tour_check2_lib import home
from ref import ring_test
t0=time.time(); random.seed(2026)
def misfit(n):
    for r in range(2,400):
        if gcd(n,r)==1 and (n*n)%r!=1: return r
    return 0
A=[random.randrange(10**9,2*10**10)|1 for _ in range(300000)]
bad=[n for n in A if home(n)!=isprime(n)]
print("two-cell wheel, 300,000 random odd numbers between 1e9 and 2e10: %d primes, wrong answers: %d"%(sum(1 for n in A if isprime(n)),len(bad)),flush=True)
bad=[n for n in A[:20000] if ring_test(n,misfit(n))!=isprime(n)]
print("Pascal's wheel on the first 20,000 of them: wrong answers: %d"%len(bad),flush=True)
B=[random.randrange(2**(b-1),2**b)|1 for b in range(40,65) for _ in range(800)]
print("20,000 random odd numbers of 40 to 64 bits: %d primes; wrong answers of the two-cell wheel: %d, of Pascal's wheel: %d"%(sum(1 for n in B if isprime(n)),sum(1 for n in B if home(n)!=isprime(n)),sum(1 for n in B if ring_test(n,misfit(n))!=isprime(n))))
print("time %.0f s"%(time.time()-t0))
