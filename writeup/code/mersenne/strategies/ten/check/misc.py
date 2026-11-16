import numpy as np, math
from lib import sieve, png
MEXP=[2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,136279841]
N=200000; s=sieve(N)
print("== continued fraction of sqrt(n): period length, n <=",N)
n=np.arange(2,N+1,dtype=np.int64); a0=np.floor(np.sqrt(n)).astype(np.int64); a0-= (a0*a0>n); a0+=((a0+1)**2<=n); nsq=a0*a0!=n
n=n[nsq]; a0=a0[nsq]; m=np.zeros(len(n),dtype=np.int64); d=np.ones(len(n),dtype=np.int64); a=a0.copy(); per=np.zeros(len(n),dtype=np.int64); act=np.ones(len(n),dtype=bool); step=0
while act.any():
    step+=1; i=np.nonzero(act)[0]; m[i]=d[i]*a[i]-m[i]; d[i]=(n[i]-m[i]*m[i])//d[i]; a[i]=(a0[i]+m[i])//d[i]
    done=i[a[i]==2*a0[i]]; per[done]=step; act[done]=False
isp=s[n]; r=per/np.sqrt(n)
print(" max period %d at n=%d; mean period/sqrt(n): primes %.3f composites %.3f"%(per.max(),n[per.argmax()],r[isp].mean(),r[~isp].mean()))
p1=isp&(n%4==1); p3=isp&(n%4==3)
print(" primes = 1 mod 4: odd period share %.4f ; primes = 3 mod 4: even period share %.4f"%((per[p1]%2==1).mean(),(per[p3]%2==0).mean()))
no3=np.ones(N+1,dtype=bool)
for q in np.nonzero(s)[0]:
    if q%4==3: no3[q::q]=False
sqf=np.ones(N+1,dtype=bool)
for q in range(2,int(N**.5)+1): sqf[q*q::q*q]=False
ok=no3[n]&sqf[n]&~isp; print(" squarefree composites with no prime factor = 3 mod 4: %d, odd period (x^2 - n y^2 = -1 solvable) share %.4f (Stevenhagen limit 0.5806)"%(ok.sum(),(per[ok]%2==1).mean()))
for lo,hi in ((2,2000),(2000,20000),(20000,200000)):
    k=ok&(n>=lo)&(n<hi); print("    n in [%d,%d): %.4f of %d"%(lo,hi,(per[k]%2==1).mean(),k.sum()))
H,W=500,1000; img=np.zeros((H,W,3))+np.array([14,24,40]); x=((n/N)*(W-1)).astype(int); y=np.clip((H-1)-(r/3.2*(H-1)),0,H-1).astype(int)
img[y[~isp],x[~isp]]=[38,70,120]; img[y[isp],x[isp]]=[160,200,240]; png("cf_period.png",img)
lr=np.log(r[isp]); h,e=np.histogram(lr,bins=30); print(" histogram of ln(period/sqrt p), primes:",h.tolist())
ex=[e_ for e_ in MEXP if 2<e_<=N]; print(" exponents: period/sqrt(p):",[round(float(r[n==e_][0]),2) for e_ in ex])
print("== Romanoff: odd n = p + 2^k")
M=10**6; s2=sieve(M); cnt=np.zeros(M+1,dtype=np.int16); k=1
while (1<<k)<M:
    cnt[(1<<k):]+=s2[:M+1-(1<<k)]; k+=1
odd=np.arange(3,M+1,2); c=cnt[odd]; print(" odd n <= 1e6: mean reps %.3f, share with none %.4f, Poisson with that mean would give %.4f; first holes:"%(c.mean(),(c==0).mean(),math.exp(-c.mean())),odd[c==0][:12].tolist())
for lo,hi in ((3,10**4),(10**4,10**5),(10**5,10**6)):
    kk=(odd>=lo)&(odd<hi); print("    [%d,%d): none %.4f mean %.3f"%(lo,hi,(c[kk]==0).mean(),c[kk].mean()))
print(" holes that are prime: %d of %d; holes among Mersenne exponents:"%(int(s2[odd[c==0]].sum()),int((c==0).sum())),[e_ for e_ in MEXP if e_<=M and e_%2 and cnt[e_]==0])
print("== exponents: order of 2")
from sympy import n_order, factorint
pr=0; idx=[]
for p in MEXP[1:]:
    o=n_order(2,p); idx.append((p-1)//o); pr+= (o==p-1)
print(" 2 is a primitive root for %d of 51 odd exponents (Artin 37.4%% -> 19.1); indices:"%pr,idx)
allp=np.nonzero(s2)[0]; 
from math import comb
exo=[e_ for e_ in MEXP if e_%2 and e_<=M]; holes=[e_ for e_ in exo if cnt[e_]==0]; pr_odd=odd[s2[odd]]; ph=float((cnt[pr_odd]==0).mean())
tail=sum(comb(len(exo),i)*ph**i*(1-ph)**(len(exo)-i) for i in range(len(holes),len(exo)+1))
print(" odd primes <= 1e6 that are not p + 2^k: share %.4f; odd exponents <= 1e6: %d, of them not p + 2^k: %d (expected %.1f), binomial tail %.3f"%(ph,len(exo),len(holes),ph*len(exo),tail))
print("== Collatz total stopping time, n <= 1e6 coprime to 6: primes against composites, size removed by a line in ln n")
st=np.zeros(M+1,dtype=np.int32)
for n_ in range(2,M+1):
    x=n_; c=0
    while x>=n_:
        x=x//2 if x%2==0 else 3*x+1; c+=1
    st[n_]=c+st[x]
nn=np.arange(M+1); k=(nn%2==1)&(nn%3!=0)&(nn>3)
for a_,b_ in ((5,10**4),(10**4,10**5),(10**5,10**6)):
    kk=k&(nn>=a_)&(nn<b_); x=np.log(nn[kk]); y=st[kk].astype(float); A=np.vstack([np.ones(len(x)),x]).T; co=np.linalg.lstsq(A,y,rcond=None)[0]; res=y-A@co; isp_=s2[nn[kk]]
    se=res.std()*math.sqrt(1/isp_.sum()+1/(~isp_).sum())
    print("  [%d,%d): raw difference %.3f; with size removed %.3f +- %.3f (primes %d)"%(a_,b_,y[isp_].mean()-y[~isp_].mean(),res[isp_].mean()-res[~isp_].mean(),se,isp_.sum()))
