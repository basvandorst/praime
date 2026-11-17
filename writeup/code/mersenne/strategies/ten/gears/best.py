import numpy as np, math, sys
exec(open('machines.py').read().split("RULES=[")[0])
TW=[k for k in range(3,NMAX+1) if isp[k] and (isp[k-2] or isp[k+2])]
def pat(f):
    return lambda k:(f(t)%k==0)
RULES=[('Wheels, one tooth per tick (1-6, 8-10, 18, 21)',lambda k:home('div',k)),('Dial (7), whole turns left out',lambda k:home('fast',k)&(t%120!=0)),('Doubling (11)',lambda k:home('double',k)),('Fibonacci (12)',lambda k:home('fib',k)),('Prime clock (13)',lambda k:home('pclock',k)&(PI>0)),('Accelerating (14)',lambda k:home('accel',k)),('Square jumps (15)',lambda k:home('square',k)),('Back and forth (16)',lambda k:home('flip',k)),('Slipping (17)',lambda k:home('slip',k)),
 ('Pattern: twin primes (20)',pat(lambda t:t*(t+2))),('Pattern: Sophie Germain (20)',pat(lambda t:t*(2*t+1))),('Pattern: quadruplets (20)',pat(lambda t:t*(t+2)*(t+6)*(t+8))),('Pattern: Euler (20)',pat(lambda t:t*t+t+41))]
print("250 gears, 100,000 ticks: most at home at once (tick)")
print("%-48s %-18s %-22s %-18s"%("machine","primes (of 53)","Mersenne exp. (of 12)","twin primes (of %d)"%len(TW)))
rows=[]
for name,fn in RULES:
    H=np.zeros((NMAX+1,T+1),dtype=bool)
    for k in range(2,NMAX+1): H[k]=fn(k)
    H[:,0]=False
    P=[k for k in range(2,NMAX+1) if isp[k]]
    a=H[P].sum(0); b=H[EXPS].sum(0); c=H[TW].sum(0)
    rows.append((name,int(a.max()),int(np.argmax(a)),int(b.max()),int(np.argmax(b)),int(c.max()),int(np.argmax(c))))
    print("%-48s %2d (%6d)        %2d (%6d)            %2d (%6d)"%rows[-1],flush=True)
from sympy import isprime
MXS={2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243}
bp=(0,0);be=(0,0);bt=(0,0)
for p in range(3,T+1):
    if not isp[p]: continue
    hm=[]
    for k in range(1,251):
        q=2*k*p+1
        if q%8 in (1,7) and pow(2,p,q)==1 and q!=2**p-1: hm.append(q)
    pr=[q for q in hm if isprime(q)]
    if len(pr)>bp[0]: bp=(len(pr),p)
    e=[q for q in pr if q in MXS]
    if len(e)>be[0]: be=(len(e),p)
    w=[q for q in pr if isprime(q-2) or isprime(q+2)]
    if len(w)>bt[0]: bt=(len(w),p)
rows.append(('Wheels that fit 2^p - 1 (22)',bp[0],bp[1],be[0],be[1],bt[0],bt[1])); print("%-48s %2d (%6d)        %2d (%6d)            %2d (%6d)"%rows[-1])
rows.append(('Lucas-Lehmer wheels (19)',1,1,1,1,1,1)); print("%-48s %2d (%6d)        %2d (%6d)            %2d (%6d)"%rows[-1])
import json
sc=lambda r:r[1]/53+r[3]/12+r[5]/len(TW)
rows.sort(key=sc,reverse=True)
print("\nranked by the three shares added up:")
for r in rows: print("   %.2f  %s"%(sc(r),r[0]))
json.dump({'rows':rows,'tw':len(TW)},open('best.json','w'))
