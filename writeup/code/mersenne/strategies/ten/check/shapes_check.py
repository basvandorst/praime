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
print("== second batch of shapes")
J=N//6; both=one=none=0
for j in range(1,J+1):
    a=(6*j-1) in S; b=(6*j+1)<=N and (6*j+1) in S
    if a and b: both+=1
    elif a or b: one+=1
    else: none+=1
print("helix: rungs with both ends prime %d, one %d, none %d"%(both,one,none))
d=0; dc=[0,0,0]; bad=0
for n in range(1,N+1):
    m=n-1; ones=0
    while m: ones+=(m%3==1); m//=3
    if ones%3!=d: bad+=1
    if n in S: dc[d]+=1
    t=n
    while t%3==0: t//=3
    d=(d+(1 if t%3==1 else 2))%3
print("terdragon: primes by direction",dc,"; exceptions to 'direction = number of digits 1 of n-1 in base 3, mod 3':",bad)
from math import comb
w=[0,0,0]
for k in range(0,13,2): w[k%3]+=comb(12,k)*2**(12-k)
print("   even counts of the digit 1 among 12 digits, by count mod 3: %s (shares %s)"%(w,[round(x/sum(w),3) for x in w]))
fu=[0,1]+[0]*(2*N+3)
for i in range(1,N+2): fu[2*i]=fu[i]; fu[2*i+1]=fu[i]+fu[i+1]
print("fractions: largest a+b %d; odd n with fraction a/(a+b) <= 1/2: %d; distinct fractions %d of %d"%(max(fu[n]+fu[n+1] for n in range(1,N+1)),sum(1 for n in range(1,N+1,2) if 2*fu[n]<=fu[n]+fu[n+1]),len({(fu[n],fu[n+1]) for n in range(1,N+1)}),N))
f=[1,1]
while sum(x*x for x in f)<N: f.append(f[-1]+f[-2])
print("Fibonacci squares needed:",f,"; columns of width 144 that can hold primes:",totient(144))
