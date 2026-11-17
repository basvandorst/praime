"""Independent recheck of the results kept from bt1.py (other code: sympy, plain loops)."""
from sympy import isprime, totient, fibonacci
from math import gcd
best=(0,)
for N in range(2,251):
    for b in (1,2,3,4,5,6,8,10,12,15,20,24,30,40,60,120):
        for a in range(b):
            v=[k for k in range(1,N+1) if k%b==a]
            if v and all(isprime(x) for x in v) and len(v)>best[0]: best=(len(v),N,b,v)
print("dial: longest ray of primes only:",best,"; holds for N =",[N for N in range(2,251) if all(isprime(x) for x in range(7,N+1,30)) and len(range(7,N+1,30))==6][0],"to",[N for N in range(2,251) if all(isprime(x) for x in range(7,N+1,30)) and len(range(7,N+1,30))==6][-1])
tw=[]
for k in range(2,2001):
    t=1
    while fibonacci(t)%k: t+=1
    if t==2*k: tw.append(k)
    assert t<=2*k
print("Fibonacci wheels up to 2000: none later than 2k; exactly 2k for",tw)
print("accelerating: wheel 2^a first home at 2^(a+1)-1:",all(min(t for t in range(1,2**(a+2)) if (t*(t+1)//2)%(2**a)==0)==2**(a+1)-1 for a in range(1,11)),"; no other wheel up to 2000 later than k:",all(min(t for t in range(1,2*k+1) if (t*(t+1)//2)%k==0)<=k for k in range(2,2001) if k&(k-1)))
full=[]
for k in range(1,301):
    pos=0;t=0
    while True:
        t+=1
        if gcd(t,k)==1: pos=(pos+1)%k
        if t%k==0 and pos==0: break
    if t==k*k: full.append(k)
print("slipping: wheels up to 300 needing k rounds = numbers with gcd(k, phi(k)) = 1:",full==[k for k in range(1,301) if gcd(k,int(totient(k)))==1],"; count",len(full))
