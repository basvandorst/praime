import math
from sympy import primerange, totient, divisors
N=500000; P=list(primerange(2,N+1)); S=set(P)
g=[P[i+1]-P[i] for i in range(len(P)-1)]; mg=max(g); print("twin gaps %d; largest gap %d after %d; primes %d"%(g.count(2),mg,P[g.index(mg)],len(P)))
d=0; dc=[0,0,0,0]
for n in range(1,N+1):
    if n in S: dc[d]+=1
    d=(d+3 if ((n&-n)<<1)&n else d+1)&3
print("dragon: primes on pieces pointing right, up, left, down:",dc)
hi=0
for n in range(2,N+1):
    r=math.isqrt(n); a=next(x for x in range(r,0,-1) if n%x==0)
    if a**4>=n: hi+=1
print("fan: numbers whose largest divisor up to the root is at least the fourth root:",hi)
s=6*math.ceil(round(N**(1/3),9)/6); print("crystal: side %d, sheets that can hold primes %d; 7 levels of the flowsnake: %d steps"%(s,totient(s),7**7))
