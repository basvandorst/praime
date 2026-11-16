import numpy as np, math, time
from lib import sieve
from sympy import primerange
N=500000; s=sieve(1100000)
print("== F atlas recheck")
def h(D):
    c=0; a=1
    while 3*a*a<=-D:
        for b in range(-a+1,a+1):
            if (b*b-D)%(4*a): continue
            cc=(b*b-D)//(4*a)
            if cc<a or (b<0 and a==cc): continue
            if math.gcd(math.gcd(a,abs(b)),cc)==1: c+=1
        a+=1
    return c
R=400; A=np.arange(1,2500,2); m=np.arange(R); cnt=np.array([int(s[m*m+m+a].sum()) for a in A]); hs=np.array([h(1-4*a) if a>1 else 1 for a in A])
o=np.argsort(-cnt,kind='stable'); print(" richest:",[(int(A[i]),int(cnt[i]),int(hs[i])) for i in o[:6]]," poorest:",(int(A[o[-1]]),int(cnt[o[-1]])))
k=A>=3; x=np.log(np.sqrt(4*A[k]-1)/hs[k]); y=np.log(cnt[k]); print(" corr ln count vs ln(sqrt|D|/h): %.3f over %d columns"%(np.corrcoef(x,y)[0,1],k.sum()))
sp=list(primerange(3,3000))
def bh(a):
    D=1-4*a; C=2.0
    for p in sp:
        if D%p==0: nu=1
        else: nu=1+(1 if pow(D%p,(p-1)//2,p)==1 else -1)
        C*=(1-nu/p)/(1-1/p)
    v=m*m+m+a; v=v[v>2]; return C*np.sum(1/np.log(v))
B=np.array([bh(int(a)) for a in A]); z=(cnt-B)/np.sqrt(B)
print(" Bateman-Horn: total predicted %.0f, found %d; corr(count, BH) %.4f; z-scores mean %.2f sd %.2f"%(B.sum(),cnt.sum(),np.corrcoef(cnt,B)[0,1],z.mean(),z.std()))
rng=np.random.default_rng(5); cs=[np.corrcoef(x,np.log(np.maximum(1,rng.poisson(B[k]))))[0,1] for _ in range(20)]; print(" corr of ln(sqrt|D|/h) with Poisson draws around BH: %.3f (and with ln BH itself: %.3f)"%(np.mean(cs),np.corrcoef(x,np.log(B[k]))[0,1]))
print(" columns with first 40 values all prime:",[int(a) for a in A if s[m[:40]*m[:40]+m[:40]+a].all()], "; h=1:",[int(a) for a in A if a>1 and h(1-4*a)==1])
print("== G spectrum peak heights against 1/|rho| (primes to 500k)")
lam=np.zeros(N+1)
for p in primerange(2,N+1):
    q=p
    while q<=N: lam[q]=math.log(p); q*=p
psi=np.cumsum(lam); M=4096; u=np.linspace(math.log(50),math.log(N),M); xx=np.exp(u); f=(psi[np.minimum(N,np.floor(xx).astype(int))]-xx+math.log(2*math.pi))/np.sqrt(xx)
g=1/np.sqrt(xx); Am=np.vstack([np.ones(M),g]).T; c=np.linalg.lstsq(Am,f,rcond=None)[0]; f=(f-Am@c)*np.hanning(M)
Z=[14.134725,21.022040,25.010858,30.424876,32.935062,37.586178,40.918719,43.327073,48.005151,49.773832,52.970321,56.446248,59.347044]
amp=[abs(np.sum(f*np.exp(-1j*t*u))) for t in Z]
print(" height/first:"," ".join("%.2f"%(a/amp[0]) for a in amp)); print(" |rho1|/|rho|: "," ".join("%.2f"%(math.hypot(.5,Z[0])/math.hypot(.5,t)) for t in Z))
print("== E scales: the bottom rows against the zeta zeros")
import mpmath as mp
t0=time.time(); zs=[float(mp.zetazero(n).imag) for n in range(1,161)]; print(" 160 zeros in %.0f s, last %.2f"%(time.time()-t0,zs[-1]))
pi=np.cumsum(s[:N+1]); E=np.concatenate([[0,0],np.cumsum(1/np.log(np.arange(2,N+1)))])
def F(x):
    lx=np.log(x); out=np.zeros_like(x)
    for gmm in zs:
        rho=complex(0.5,gmm); out-=2*np.real(np.exp(rho*lx)/(rho*lx)*(1+1/(rho*lx)))
    return out
for H in (16384,8192,4096):
    x=np.arange(2*H,N-H,97); d=2*pi[x]-pi[x-H]-pi[x+H]-(2*E[x]-E[x-H]-E[x+H]); xf=x.astype(float); pr=2*F(xf)-F(xf-H)-F(xf+H)
    print(" H=%d: correlation of the row with the first 160 zeros %.3f; rms row %.2f, rms prediction %.2f, rms difference %.2f (needs zeros up to about %.0f)"%(H,np.corrcoef(d,pr)[0,1],d.std(),pr.std(),(d-pr).std(),2*math.pi*N/(2*H)*2))
