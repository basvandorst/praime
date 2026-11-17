"""Independent recomputation of the numbers quoted from lay2.py, two.py and link.py (other algorithms: coefficient lists, Rabin's irreducibility test, a fold instead of a division)."""
import json, math, itertools, numpy as np
from sympy import primerange, factorint, isprime
ME=json.load(open('mexp.json'))
def polymul(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if y: r[i+j]^=1
    return r
def bits(n): return [(n>>i)&1 for i in range(n.bit_length())]
for p,fs in ((7,[11,13]),(17,[313,471]),(31,[37,41,47,55,59,61]),(89,[2243,2359,2543,3121,3583,3785,3961,4091]),(127,[131,137,143,145,157,167,171,185,191,193,203,211,213,229,239,241,247,253])):
    r=[1]
    for f in fs: r=polymul(r,bits(f))
    print("p=%d: product without carries of %d patterns is %d ones: %s"%(p,len(fs),p,r==[1]*p))
def primroot2(p): return all(pow(2,(p-1)//r,p)!=1 for r in factorint(p-1))
od=[p for p in ME if p>2]; pr=[p for p in od if primroot2(p)]; print("odd exponents %d, with 2 a primitive root %d"%(len(od),len(pr)))
sh={1:[0,0],3:[0,0],5:[0,0],7:[0,0]}
for q in primerange(3,2_000_000):
    a=sh[q%8]; a[0]+=1
    if q%8 in (3,5) and primroot2(q): a[1]+=1
e=sum(sh[p%8][1]/sh[p%8][0] for p in od); v=sum((sh[p%8][1]/sh[p%8][0])*(1-sh[p%8][1]/sh[p%8][0]) for p in od); print("expected %.1f +- %.1f"%(e,math.sqrt(v)))
L=24; n=np.arange(1<<L,dtype=np.int64)
def strip_divides(w):
    g=n^(n<<1); acc=np.zeros(len(n),dtype=np.int64); M=(1<<w)-1
    while g.max()>0: acc^=g&M; g>>=w
    return acc==0
for w in (2,3,4,5,6,7,8):
    M=(1<<w)-1; y=strip_divides(w); x=(n%M==0); out=[]
    for LL in (16,20,24):
        T=1<<LL; a=int(x[:T].sum()); b=int(y[:T].sum()); c=int((x[:T]&y[:T]).sum()); out.append(c/(a*b/T))
    print("width %d (M=%d): with and without carries together, against chance: below 2^16 %.4f, 2^20 %.4f, 2^24 %.4f"%(w,M,*out))
N5=500000; y=strip_divides(7)[:N5+1]; x=(n[:N5+1]%127==0); x[0]=False; yy=y.copy(); yy[0]=False
print("page setting (1..500000, width 127): green among multiples %.1f%%, among the others %.1f%%"%(100*(x&yy).sum()/x.sum(),100*(yy&~x)[1:].sum()/((~x)[1:].sum())))
def pmulmod(a,b,f,df):
    r=0
    while b:
        if b&1: r^=a
        b>>=1; a<<=1
        if (a>>df)&1: a^=f
    return r
def pgcd(a,b):
    while b:
        while a.bit_length()>=b.bit_length() and a: a^=b<<(a.bit_length()-b.bit_length())
        a,b=b,a
    return a
def irreducible(f):
    d=f.bit_length()-1; x=2; t=2
    for i in range(1,d//2+1):
        t=pmulmod(t,t,f,d)
        if pgcd(f,t^x)!=1: return False
    return True
for k in (20,):
    a=b=c=m=0
    for q in range((1<<(k-1))+1,1<<k,2):
        if bin(q).count('1')%2==0: continue
        m+=1; P=isprime(q); I=irreducible(q); a+=P; b+=I; c+=(P and I)
    print("%d bits, odd with an odd number of ones: %d numbers, primes %d, patterns %d, both %d, expected %.1f, ratio %.3f"%(k,m,a,b,c,a*b/m,c/(a*b/m)))
full=[]; best=(0,0); s=0; t=0
for q in primerange(3,1<<20):
    b_=bin(q)[2:]; w=len(b_); k=sum(isprime(int(b_[i:]+b_[:i],2)) for i in range(w)); s+=k/w; t+=1
    if k==w: full.append(q)
    elif w==20 and k>best[0]: best=(k,q)
print("all rotations prime:",full,"; mean share %.4f; best of 20 bits not all ones: %d with %d"%(s/t,best[1],best[0]))
runs={}
for q in primerange(3,1<<20):
    r=len(list(itertools.groupby(bin(q)[2:]))); runs[r]=runs.get(r,0)+1
print("primes below 2^20 by runs:",sorted(runs.items()))
LIM=10**12; seen={}; dup=set(); cnt=0
w=3
while 2**w-1<LIM*2:
    b=2
    while True:
        v=(b**w-1)//(b-1)
        if v>=LIM: break
        cnt+=1
        if v in seen: dup.add((v,seen[v],(b,w)))
        else: seen[v]=(b,w)
        b+=1
    w+=1
print("all-ones numbers with 3 or more digits below 1e12: %d; twice: %s"%(cnt,sorted(dup)))
