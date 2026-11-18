"""Pascal's wheel on a ring the number does not fit: plain-Python check of what the page shows.
1. every n from 2 to 5,000 (even numbers too): the wheel is home exactly when n is prime
2. the three results for 250 wheels (primes, Mersenne exponents, twin primes) and the ticks
3. the first number that cheats on a fixed ring (every odd composite, no prefilter)
4. even numbers that fool the doubling wheel
"""
import json
from math import gcd
from sympy import isprime
from ref import ring_test
def misfit(n):
    for r in range(2,400):
        if gcd(n,r)==1 and (n*n)%r!=1: return r
    return 0
def home(n):
    r=misfit(n); return r>0 and ring_test(n,r)
out={}
bad=[n for n in range(2,5001) if home(n)!=isprime(n)]
rings={}
for n in range(2,5001): rings[misfit(n)]=rings.get(misfit(n),0)+1
print("1. n = 2..5000: mistakes:",bad,"| lamps lit:",sum(1 for n in range(2,5001) if home(n)),"| ring sizes used:",sorted(rings.items()))
out['mistakes_to_5000']=len(bad); out['lit_to_5000']=sum(1 for n in range(2,5001) if home(n)); out['rings_to_5000']=sorted(rings.items())
MX={2,3,5,7,13,17,19,31,61,89,107,127}
lit=[]; best={'p':(0,0),'m':(0,0),'w':(0,0)}
for t in range(2,251):
    if home(t): lit.append(t)
    p=[k for k in lit if isprime(k)]; m=[k for k in p if k in MX]; w=[k for k in p if isprime(k-2) or isprime(k+2)]
    for key,a in (('p',p),('m',m),('w',w)):
        if len(a)>best[key][0]: best[key]=(len(a),t)
print("2. 250 wheels: primes %d (tick %d), Mersenne exponents %d (tick %d), twin primes %d (tick %d); composite lamps: %d"%(best['p']+best['m']+best['w']+(sum(1 for k in lit if not isprime(k)),)))
out['row250']=[best['p'][0],best['p'][1],best['m'][0],best['m'][1],best['w'][0],best['w'][1]]
first={}
for r in (5,7):
    for n in range(9,600000,2):
        if gcd(n,r)==1 and pow(2,n,n)==2 and not isprime(n) and ring_test(n,r): first[r]=(n,n%r); break
print("3. first odd composite that comes home on a fixed ring (n, n mod r):",first,"(base-2 condition used as a proved shortcut: x = 1 gives 2^n = 2)")
out['first_fixed']={str(k):v for k,v in first.items()}
ev=[n for n in range(4,3000000,2) if pow(2,n,n)==2]
print("4. even n below 3,000,000 with 2^n = 2 (mod n):",ev,"| home on their misfit ring:",[n for n in ev if home(n)])
out['even_psp']=ev; out['even_home']=[n for n in ev if home(n)]
json.dump(out,open('pascal_check.json','w'))
