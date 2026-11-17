import math, json, numpy as np, random
from sympy import isprime, nextprime, prevprime
ME=json.load(open('../round7/mexp.json')) if False else [2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,136279841]
REL={'p+2 prime':lambda p:isprime(p+2),'p-2 prime':lambda p:isprime(p-2),'twin (p-2 or p+2)':lambda p:isprime(p+2) or isprime(p-2),'p+4 or p-4 prime':lambda p:isprime(p+4) or isprime(p-4),'p+6 or p-6 prime':lambda p:isprime(p+6) or isprime(p-6),
 '2p+1 prime':lambda p:isprime(2*p+1),'2p-1 prime':lambda p:isprime(2*p-1),'(p+1)/2 prime':lambda p:p%2==1 and isprime((p+1)//2),'(p-1)/2 prime':lambda p:p%2==1 and isprime((p-1)//2),'4p+1 prime':lambda p:isprime(4*p+1),'4p-1 prime':lambda p:isprime(4*p-1),'6p+1 prime':lambda p:isprime(6*p+1),'6p-1 prime':lambda p:isprime(6*p-1),
 'p+2 and p+6 prime':lambda p:isprime(p+2) and isprime(p+6),'p^2+p+41-like: p+1 twice a prime or p-1 twice a prime':lambda p:(isprime((p+1)//2) or isprime((p-1)//2)),'p = 1 mod 4':lambda p:p%4==1,'p = 1 or 7 mod 8':lambda p:p%8 in (1,7)}
def survives(q,K=300):
    for k in range(1,K+1):
        f=2*k*q+1
        if f%8 in (1,7) and pow(2,q,f)==1: return False
    return True
odd=[p for p in ME if p>40]
res={r:[0,0.0,0.0,0.0,0.0] for r in REL}
for p in odd:
    nb=[]; q=p
    for _ in range(150): q=nextprime(q); nb.append(q)
    q=p
    for _ in range(150):
        if q<=43: break
        q=prevprime(q); nb.append(q)
    q=nb[149]
    while len(nb)<300: q=nextprime(q); nb.append(q)
    sv=[q for q in nb if survives(q)]
    for r,fn in REL.items():
        fa=sum(1 for q in nb if fn(q))/len(nb); fb=sum(1 for q in sv if fn(q))/max(1,len(sv))
        e=res[r]; e[0]+=1 if fn(p) else 0; e[1]+=fa; e[2]+=fa*(1-fa); e[3]+=fb; e[4]+=fb*(1-fb)
print("Mersenne exponents above 40: %d. For each property: how many exponents have it; expected from the 300 nearest primes; expected from those of them whose 2^q - 1 has no factor 2kq+1 with k <= 300"%len(odd))
for r,e in res.items():
    za=(e[0]-e[1])/math.sqrt(max(1e-9,e[2])); zb=(e[0]-e[3])/math.sqrt(max(1e-9,e[4]))
    print("   %-58s found %2d | nearby primes %5.1f (z %+.1f) | nearby survivors %5.1f (z %+.1f)"%(r,e[0],e[1],za,e[3],zb),flush=True)
