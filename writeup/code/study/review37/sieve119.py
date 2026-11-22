import numpy as np, math, sys, time, json
import sympy
P=np.array(list(sympy.primerange(3,int(math.isqrt(2*10**9))+2)),dtype=np.int64)
def counts(lo,hi,H):  # primes in [lo+iH, lo+(i+1)H), lo even, H even
    out=[]
    for a in range(lo,hi,10**7):
        n=10**7//2; s=np.ones(n,bool)  # odd numbers a+1+2j
        for p in P:
            st=(-(a+1))%p
            if st%2: st+=p
            st//=2
            s[st::p]=False
        idx=np.nonzero(s)[0]*2+a+1
        out.append(np.bincount((idx-a)//H,minlength=10**7//H))
    return np.concatenate(out)
Z=np.loadtxt('/home/bas/praime/writeup/data/riemann/zeros_10000_flint.txt')
def ei(w):
    s=np.ones_like(w);t=np.ones_like(w)
    for k in range(1,12): t=t*k/w; s=s+t
    return np.exp(w)/w*s
H=int(__import__("os").environ.get("HH","1000000")); res={}
for b in [int(a) for a in sys.argv[1:]]:
    t0=time.time(); Lb=10**9+b*10**8
    c=counts(Lb,Lb+10**8,H)
    x=Lb+np.arange(10**8//H)*H; e=H/np.log(x+H/2)
    rho=0.5+1j*Z; S=np.zeros(10**8//H)
    for i0 in range(0,len(Z),1000):
        r=rho[i0:i0+1000]; S+=-2*np.real(ei(np.outer(np.log(x+H),r))-ei(np.outer(np.log(x),r))).sum(axis=1)
    m=c-e
    print(b,'meas var/mean',round(np.var(m,ddof=1)/e.mean(),4),'corr(c-e, S)',round(np.corrcoef(m,S)[0,1],3),'var resid/mean',round(np.var(m-S,ddof=1)/e.mean(),4),'time',round(time.time()-t0,1),flush=True)
