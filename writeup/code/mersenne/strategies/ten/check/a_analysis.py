import numpy as np, math
from lib import sieve
from sympy import isprime, primerange
N=500000; s=sieve(1100000); P=np.nonzero(s[:N+1])[0]
MEXP=[2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,136279841]
print("== A chaoslane mod 12, transitions up to 500k (rows from 1,5,7,11), ratio to even share")
Q=P[P>=5]; c=np.array([{1:0,5:1,7:2,11:3}[int(x)] for x in Q%12]); T=np.zeros((4,4)); np.add.at(T,(c[:-1],c[1:]),1); print(np.round(T/T.mean(),2)); print(" repeats %.3f"%(np.trace(T)/T.sum()))
g=np.diff(Q); 
for a,b,gg in [(1,5,4),(5,7,2),(1,7,6),(1,1,12),(5,11,6),(11,1,2)]:
    m=(Q[:-1]%12==a)&(Q[1:]%12==b); print("  %d->%d: %d steps, %.0f%% of them by the gap %d"%(a,b,m.sum(),100*np.mean(g[m]==gg),gg))
print(" mirror pairs (a->b against -b->-a):", [(int(T[i,j]),int(T[3-j,3-i])) for i,j in [(0,1),(0,2),(1,3),(0,0),(1,1)]])
print("== A2 gap classes mod 6: transition counts rows 0,2,4")
gc=(g%6)//2; G3=np.zeros((3,3),int); np.add.at(G3,(gc[:-1],gc[1:]),1); print(G3)
print("== B returnmap symmetry (a,b) vs (b,a)")
a=g[:-1]; b=g[1:]; M=np.zeros((61,61)); m=(a<=60)&(b<=60); np.add.at(M,(a[m],b[m]),1)
iu=[(i,j) for i in range(2,61,2) for j in range(i+2,61,2) if M[i,j]+M[j,i]>=20]; chi=sum((M[i,j]-M[j,i])**2/(M[i,j]+M[j,i]) for i,j in iu); print(" pairs of cells %d, chi2 %.1f"%(len(iu),chi))
print(" top cells:",sorted([(int(M[i,j]),i,j) for i in range(2,61,2) for j in range(2,61,2)],reverse=True)[:8])
print(" equal gaps %d; of them 6|g: %d; share of primes with both gaps multiples of 6: %.3f"%((a==b).sum(),((a==b)&(a%6==0)).sum(),np.mean((a%6==0)&(b%6==0))))
E6=np.zeros((3,3)); np.add.at(E6,((a%6)//2,(b%6)//2),1); print(" class table (rows gap before 0,2,4):\n",E6.astype(int))
print("== C dejavu patterns from 11 (recheck) and Hardy-Littlewood")
pr=P.tolist(); i0=pr.index(11); cnt={}
for i in range(i0+1,len(pr)):
    L=0
    while i+L+1<len(pr) and i0+L+1<i and pr[i+L+1]-pr[i+L]==pr[i0+L+1]-pr[i0+L]: L+=1
    cnt[L]=cnt.get(L,0)+1
print(" rows:",dict(sorted(cnt.items())))
small=list(primerange(2,2000))
def hl(pat,X):
    k=len(pat); C=1.0
    for q in small:
        nu=len({x%q for x in pat}); C*=(1-nu/q)/(1-1/q)**k
    t=np.arange(20,X,dtype=float); return C*np.sum(1/np.log(t)**k)
off=[0,2,6,8,12,18,20,26]
for k in range(2,8):
    pat=off[:k]; real=sum(1 for p in pr if p>11 and all(s[p+o] for o in pat) and p+pat[-1]<=N)
    print("  %d primes spaced %s: found %d, Hardy-Littlewood %.1f"%(k,pat,real,hl(pat,N)))
print("== D goldfirst recheck and the sieve+coin law")
od=[p for p in pr if p>2]; J=np.zeros(N+1,int)
for n in range(6,N+1,2):
    j=0
    while not s[n-od[j]]: j+=1
    J[n]=j+1
ev=np.arange(6,N+1,2); m3=ev%3==0
print(" mean tries: multiples of 3 %.2f, others %.2f; others at once %.3f; record n=%d tries %d prime %d"%(J[ev[m3]].mean(),J[ev[~m3]].mean(),np.mean(J[ev[~m3]]==1),ev[np.argmax(J[ev])],J[ev].max(),od[J[ev].max()-1]))
n1=ev[ev%3==1]; print(" n = 1 mod 3: share whose partner is 1 mod 3: %.4f ; n = 2 mod 3: partner 2 mod 3: %.4f (exceptions are n - p = 3)"%(np.mean(np.array([od[J[n]-1]%3 for n in n1])==1),np.mean(np.array([od[J[n]-1]%3 for n in ev[ev%3==2]])==2)))
sp=[q for q in small if q<=47]; prod=np.prod([1-1/q for q in sp]); rng=np.random.default_rng(3)
samp=rng.choice(ev[ev>100000],4000,replace=False); surv=np.zeros(400); 
for n in samp:
    alive=1.0
    for j,p in enumerate(od[:400]):
        surv[j]+=alive; r=n-p
        ok=all(r%q for q in sp) 
        alive*=1-(min(1,1/(math.log(r)*prod)) if ok else 0)
surv/=len(samp); real=np.array([np.mean(J[samp]>j) for j in range(400)])
print(" P(more than j tries), model (sieve to 47 + coin) vs real: "+", ".join("j=%d: %.4f/%.4f"%(j,surv[j],real[j]) for j in (1,5,10,20,40,60)))
print(" mean tries model %.2f real %.2f"%(surv.sum(),J[samp].mean()))
print("== H survivors recheck")
for y,w in ((7,48),(197,300)):
    dead=np.zeros(N+1,bool); dead[:2]=True
    for q in primerange(2,y+1): dead[q::q]=True
    S=np.nonzero(~dead)[0]; ip=s[S]; pairs=int((ip[:-1]&ip[1:]).sum()); pe=sum(ip[i:i+2000].sum()**2/len(ip[i:i+2000]) for i in range(0,len(ip),2000))
    pr_=np.prod([1-1/q for q in primerange(2,y+1)]); print("  y=%d survivors %d primes %d pairs %d ratio %.3f; found/pred %.3f"%(y,len(S),ip.sum(),pairs,pairs/pe,ip.sum()/np.sum(1/(np.log(S)*pr_))))
print("== I exponents")
def rev(n,B):
    r=0
    while n: r=r*B+n%B; n//=B
    return r
for B in (10,2):
    K=B/sum(1 for x in range(B) if math.gcd(x,B)==1)
    for q in primerange(2,B+2):
        if (B*B-1)%q==0 and B%q: K*=q/(q-1)
    hit=[k for k in MEXP if isprime(rev(k,B))]; exp=sum(K/math.log(rev(k,B)) for k in MEXP if math.gcd(rev(k,B),B)==1 and rev(k,B)>3)
    print("  base %d: %d of 52 exponents have a prime reversal; plain expectation %.1f; they are %s"%(B,len(hit),exp,hit[:20]))
