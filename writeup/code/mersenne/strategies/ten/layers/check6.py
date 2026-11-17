"""Independent recomputation of the numbers quoted from add1.py, bits1.py, bits2.py, fold1.py and depth1.py: plain loops, no FFT, no Walsh transform, sympy/gmpy2 for primes and orders. Run "bits2.py 21" (writes W21.npy) and "depth2.py 200" (writes depth_rows_200.json) first."""
import math, random, json, numpy as np, gmpy2
from sympy import primerange, n_order, factorint
N=1_000_000
isp=bytearray(N+1)
for p in primerange(2,N+1): isp[p]=1
A=np.frombuffer(bytes(isp),dtype=np.uint8).astype(bool); PL=np.array([p for p in range(3,N+1) if isp[p]],dtype=np.int64)
def lowest(vals,ns,B):
    hit=tot=0
    for b in range(0,len(ns)-B+1,B):
        idx=list(range(b,b+B)); k=sum(isp[ns[i]] for i in idx); idx.sort(key=lambda i:vals[i]); hit+=sum(isp[ns[i]] for i in idx[:k]); tot+=k
    return hit,tot
ns=list(range(3,49999,2)); cum=np.concatenate([[0],np.cumsum(A[1:]&(np.arange(1,N+1)>2))])
v=[]
for n in ns:
    q=PL[PL<=N-2*n]; v.append(int(A[q+2*n].sum()))
h,t=lowest(v,ns,2048); print("pairs of primes 2n apart, primes to 1e6, odd n < 49999, blocks of 2048: lowest counts are prime in %d of %d = %.1f%%; share of primes among these n %.1f%%"%(h,t,100*h/t,100*sum(isp[n] for n in ns)/len(ns)))
N5=500_000; ns=list(range(3,10001,2)); v=[]
for n in ns:
    lim=N5-2*n; q=PL[PL<=lim]; v.append(int(A[q+2*n].sum())/max(1,len(q)))
h,t=lowest(v,ns,512); 
hit=tot=0
for b in range(0,len(ns),512):
    idx=list(range(b,min(len(ns),b+512))); k=sum(isp[ns[i]] for i in idx); idx.sort(key=lambda i:v[i]); hit+=sum(isp[ns[i]] for i in idx[:k]); tot+=k
print("page setting (primes to 500,000, odd n to 10,000, blocks of 512, last block short): %d of %d = %.1f%%"%(hit,tot,100*hit/tot))
ns=list(range(9,60001,2)); v=[]
for n in ns:
    q=PL[PL<=(n-3)//2]; v.append(int(A[n-2*q].sum()))
hit=tot=0
for b in range(0,len(ns),512):
    idx=list(range(b,min(len(ns),b+512))); k=sum(isp[ns[i]] for i in idx); idx.sort(key=lambda i:v[i]); hit+=sum(isp[ns[i]] for i in idx[:k]); tot+=k
print("n = p + 2q, page setting (odd n from 9 to 60,000, blocks of 512): %d of %d = %.1f%%; primes among these n %.1f%%"%(hit,tot,100*hit/tot,100*sum(isp[n] for n in ns)/len(ns)))
rnd=random.Random(3); blocks=sorted(rnd.sample(range(0,122),16)); hit=tot=0; ratios={3:[],5:[],7:[],11:[]}; pv=[]
for b in blocks:
    ns=[500001+2*(b*2048+i) for i in range(2048)]; v=[]
    for n in ns:
        q=PL[PL<=(n-3)//2]; v.append(int(A[n-2*q].sum()))
    k=sum(isp[n] for n in ns); idx=sorted(range(2048),key=lambda i:v[i]); hit+=sum(isp[ns[i]] for i in idx[:k]); tot+=k
    for i,n in enumerate(ns):
        if isp[n]: pv.append(v[i])
        else:
            for p in (3,5,7,11):
                if n%p==0:
                    ratios[p].append(v[i]); break
med=np.median(pv)
print("n = p + 2q between 500,000 and 1e6, 16 random blocks of 2048 odd numbers: %d of %d = %.1f%%; steps by smallest factor 3, 5, 7, 11: %s"%(hit,tot,100*hit/tot," ".join("%.3f"%(np.median(ratios[p])/med) for p in (3,5,7,11))))
K=21; P21=[p for p in range(1<<20,1<<21) if isp[p]] if (1<<21)<=N else list(primerange(1<<20,1<<21))
sh=[sum((p>>j)&1 for p in P21)/len(P21) for j in range(1,20)]; print("21-bit primes: %d; share of ones per layer 1..19: min %.4f max %.4f"%(len(P21),min(sh),max(sh)))
arr=np.array(P21,dtype=np.int64); best=0
for i in range(1,20):
    for j in range(i+1,20):
        c=np.bincount(2*((arr>>i)&1)+((arr>>j)&1),minlength=4)/len(arr); pa=[c[0]+c[1],c[2]+c[3]]; pb=[c[0]+c[2],c[1]+c[3]]
        best=max(best,sum(c[2*x+y]*math.log2(c[2*x+y]/(pa[x]*pb[y])) for x in (0,1) for y in (0,1)))
print("largest information between two layers: %.2e bits"%best)
M=1<<21; allp=np.array(list(primerange(3,M)),dtype=np.int64); odd=np.arange(1,M,2,dtype=np.int64)
def par(x,mask):
    y=x&mask; r=np.zeros(len(x),dtype=np.int64)
    for j in range(21): r^=(y>>j)&1
    return r
def bias(setarr,mask):
    mu=len(setarr)/len(odd); e_all=(par(odd,mask)==0).sum(); e=(par(setarr,mask)==0).sum(); o=len(setarr)-e
    return (e/e_all-o/(len(odd)-e_all))/(2*mu)
print("bias of the top layer alone: %.4f; of all layers 1..20: %.4f; of 110011001100100010000: %.4f"%(bias(allp,1<<20),bias(allp,(1<<21)-2),bias(allp,int('110011001100100010000',2))))
W=np.load('W21.npy'); top=np.argsort(-np.abs(W))[:1000]; rough=odd[(odd%3!=0)&(odd%5!=0)&(odd%7!=0)&(odd%11!=0)&(odd%13!=0)]
rs=rnd.sample(list(top),150); b1=[bias(allp,int(m)<<1) for m in rs]; b2=[bias(rough,int(m)<<1) for m in rs]
print("150 of the 1,000 strongest sets, biases recomputed by direct counting: correlation primes vs numbers free of 3, 5, 7, 11, 13: %.3f; largest difference from the transform %.1e"%(np.corrcoef(b1,b2)[0,1],max(abs(b1[i]-W[rs[i]]) for i in range(150))))
n21=np.arange((1<<20)+1,1<<21,2,dtype=np.int64); y=np.zeros(len(n21),bool); y[(arr-(1<<20)-1)//2]=True; mu=y.mean()
def Hb(p): return 0.0 if p<=0 or p>=1 else -p*math.log2(p)-(1-p)*math.log2(1-p)
for w in (2,4,8,12,5,7,13):
    Mw=(1<<w)-1; r=n21%Mw; tot=np.bincount(r,minlength=Mw); pr=np.bincount(r[y],minlength=Mw); h=sum(tot[i]/len(n21)*Hb(pr[i]/tot[i]) for i in range(Mw) if tot[i]); k=sum(1 for i in range(Mw) if tot[i])
    print("  fold %2d (mod %d): information %.4f bits"%(w,Mw,Hb(mu)-h-(k-1)/(2*len(n21)*math.log(2))))
print("full knowledge: %.4f bits"%Hb(mu))
Q=140000; fl=[]; full=0; cnt=0; ok=True
for q in primerange(3,Q+1):
    d=n_order(2,q); cnt+=1; full+=(d==q-1); ok&=(d>=math.log2(q+1)-1e-12)
    if 2**d-1==q: fl.append(q)
print("odd primes to %d: %d; depth >= log2(q+1) for all: %s; equality for %s; depth q-1 for %d (%.2f%%)"%(Q,cnt,ok,fl,full,100*full/cnt))
T=json.load(open('depth_rows_200.json')); bad=0; tot=0
for w in range(1,201):
    m=2**w-1
    for d in range(1,w+1):
        if w%d==0:
            for qs in T[str(d)]:
                q=int(qs)
                while m%q==0: m//=q
    if m!=1: bad+=1
    for qs in T[str(w)]:
        q=int(qs); tot+=1
        if not gmpy2.is_prime(q) or pow(2,w,q)!=1 or any(pow(2,w//r,q)==1 for r in factorint(w)): bad+=1
print("table to w = 200: %d primes; rows where 2^w - 1 is not fully split by the primes of the rows d | w, or entries that are not prime of order w: %d"%(tot,bad))
print("rows with a single prime equal to 2^w - 1:",[w for w in range(2,201) if len(T[str(w)])==1 and int(T[str(w)][0])==2**w-1])
