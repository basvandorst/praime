# Round 66 check: own listing of primes in [1e18, 1e18+2e6) by gmpy2, entropy, se, bootstrap, MM bias sim
import gmpy2, math, numpy as np, sys
L0=10**18; L1=L0+2*10**6
P=[]; x=gmpy2.next_prime(L0-1) if False else gmpy2.next_prime(L0)
while x < L1:
    P.append(int(x)); x=gmpy2.next_prime(x)
P_in=list(P); P.append(int(x))  # one beyond, as repo
rep=[int(s) for s in open('/home/bas/praime/writeup/data/study/r17_primes.txt').read().split()]
print('own primes in window',len(P_in),'with one beyond',len(P),'repo',len(rep),'equal',P==rep)
# spot-check with sympy (proof below 2^64)
import sympy, random
random.seed(1); print('sympy spot', all(sympy.isprime(q) for q in random.sample(P,200)))
g=np.diff(np.array(P,dtype=object)).astype(np.int64)
def H_mm(g):
    v,c=np.unique(g,return_counts=True); N=c.sum(); p=c/N
    return -np.sum(p*np.log2(p)), (len(c)-1)/(2*N*math.log(2)), len(c), N
h,mm,K,N=H_mm(g); print('gaps',N,'cells',K,'plugin',h,'MM term',mm,'H',h+mm)
v,c=np.unique(g,return_counts=True); p=c/N; l=-np.log2(p)
se=np.sqrt(np.sum(p*l**2)-np.sum(p*l)**2)/math.sqrt(N); print('iid se',se)
rng=np.random.default_rng(5)
for B in (1,50,500):
    nb=N//B; hs=[]
    for r in range(300):
        idx=rng.integers(0,nb,nb); gg=np.concatenate([g[i*B:(i+1)*B] for i in idx]) if B>1 else g[rng.integers(0,N,N)]
        a,b,_,_=H_mm(gg); hs.append(a+b)
    print('block bootstrap B=%d sd %.5f'%(B,np.std(hs)))
# MM-corrected bias by simulation from the data's own distribution (proxy for the truth)
hs=[]
for r in range(300):
    s=rng.choice(v,size=N,p=p); a,b,_,_=H_mm(s); hs.append(a+b)
print('MM-corrected bias (proxy truth=%.5f): mean %.5f'%(h,np.mean(hs)-h), 'plugin-only bias', np.mean(hs)-h-mm)
lam=math.log(P[0]+(P[-1]-P[0])/2); print('lambda',lam)
# law at the window's lambda by interpolation of r57's grid (quadratic in lambda through 34.5, 41.4, 50.0)
import json
r=json.load(open('/home/bas/praime/writeup/data/study/r57.json'))
xs=[34.5,41.4,50.0]; ys=[r['%.1f'%x]['H_law'] for x in xs]
co=np.polyfit(xs,ys,2); HL=np.polyval(co,lam); print('law at 41.4',ys[1],'interp at lam',HL,'slope',np.polyval(np.polyder(co),41.4))
Hd=h+mm; print('data-law(41.4) %.4f z %.2f ; data-law(lam) %.4f z %.2f'%(Hd-ys[1],(Hd-ys[1])/se,Hd-HL,(Hd-HL)/se))
K0=0.7627474; asym=(math.log(math.e*lam)-K0)/math.log(2); print('asym',asym,'data-asym',Hd-asym, 'law(lam)-asym', HL-asym)
