import math
from sympy import primerange, isprime
for N in (500000,1000000):
    y0=round(N**0.4027); R=math.isqrt(N); comp=0; rm1=0; rm3=0
    for q in primerange(y0+1,R+1):
        for c in primerange(q,N//q+1):
            m=q*c; comp+=1; hit1=False; hit3=False
            for fa,fb in ((1,1),(2,1),(3,1),(3,2),(4,1)):
                U=fa*q+fb*c; t=4*fa*fb*m; r=math.isqrt(t)
                if U-r<=2:
                    hit3=True
                    if (fa,fb)==(1,1): hit1=True
            rm1+=hit1; rm3+=hit3
    print("N=%d, sieve to %d: composite survivors %d; struck as products of two nearly equal primes %d (%.1f%%); with the ratios 2, 3, 3/2, 4 as well %d (%.1f%%)"%(N,y0,comp,rm1,100*rm1/comp,rm3,100*rm3/comp))
