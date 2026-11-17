import json, math, sys
from sympy import primerange, isprime, n_order, factorint
import numpy as np
ME=json.load(open('mexp.json'))
def clmul(a,b):
    r=0
    while b:
        if b&1: r^=a
        a<<=1; b>>=1
    return r
def pmod(a,m):
    dm=m.bit_length()
    while a.bit_length()>=dm: a^=m<<(a.bit_length()-dm)
    return a
def pdiv(a,m):
    q=0; dm=m.bit_length()
    while a.bit_length()>=dm:
        s=a.bit_length()-dm; q|=1<<s; a^=m<<s
    return q,a
print("1. WITHOUT CARRIES: 2^p - 1 as a product of patterns under carry-less multiplication (XOR instead of add)")
print("   theory: the strip of p ones splits into (p-1)/d irreducible patterns of d+1 bits, d = order of 2 mod p; it stays whole iff d = p-1")
for p in [q for q in ME if q<=127 and q>2]:
    d=n_order(2,p); R=(1<<p)-1; fs=[]
    if d<=14:
        rem=R
        for f in range((1<<d)+1,1<<(d+1),2):
            q,r=pdiv(rem,f)
            if r==0: fs.append(f); rem=q
            if rem==1: break
        prod=1
        for f in fs: prod=clmul(prod,f)
        assert prod==R and len(fs)==(p-1)//d,(p,fs)
    print("   p=%3d: order of 2 mod p = %3d -> %s"%(p,d,"stays whole (also prime without carries)" if d==p-1 else ("%d patterns of %d bits%s"%((p-1)//d,d+1,": "+" x ".join(map(str,fs)) if fs and len(fs)<=18 else ""))))
od=[p for p in ME if p>2]; pr=[p for p in od if n_order(2,p)==p-1]
print("   odd exponents: %d; with 2 a primitive root (strip whole without carries): %d = %s"%(len(od),len(pr),pr))
cls={}
for q in primerange(3,2_000_000):
    c=q%8; a=cls.setdefault(c,[0,0]); a[0]+=1; a[1]+=(n_order(2,q)==q-1)
exp_=sum(cls[p%8][1]/cls[p%8][0] for p in od); var=sum((cls[p%8][1]/cls[p%8][0])*(1-cls[p%8][1]/cls[p%8][0]) for p in od)
print("   expected for primes with the same remainders mod 8 (shares measured on the primes below 2e6: %s): %.1f +- %.1f"%({c:round(v[1]/v[0],3) for c,v in sorted(cls.items())},exp_,math.sqrt(var)))
print("2. NUMBERS BELOW 2^20: prime with carries against irreducible without carries")
K=20; N=1<<K; irr=bytearray([1])*N; irr[0]=irr[1]=0
for f in range(2,1<<(K//2+1)):
    if irr[f]:
        df=f.bit_length()-1
        for g in range(2,1<<(K-df)):
            irr[clmul(f,g)]=0
pri=bytearray(N)
for q in primerange(2,N): pri[q]=1
a=sum(pri); b=sum(irr); c=sum(1 for n in range(N) if pri[n] and irr[n])
oddpc=[n for n in range(3,N,2) if bin(n).count('1')%2==1]; po=sum(pri[n] for n in oddpc); io=sum(irr[n] for n in oddpc)
print("   primes %d, irreducible patterns %d, both %d; if independent among all numbers: %.0f; if independent among the odd numbers with an odd number of ones (%d; they hold %d of the primes and %d of the patterns): %.0f"%(a,b,c,a*b/N,len(oddpc),po,io,po*io/len(oddpc)))
print("3. ROTATIONS: share of the w rotations of a w-bit prime that are prime")
best={}; tot=0; s=0
for q in primerange(3,N):
    w=q.bit_length(); k=0
    for r in range(w):
        x=((q<<r)|(q>>(w-r)))&((1<<w)-1)
        if pri[x]: k+=1
    tot+=1; s+=k/w
    if q!=(1<<w)-1 and k/w>best.get(w,(0,0,0))[0]: best[w]=(k/w,q,k)
    if k==w: print("   all %d rotations prime: %d"%(w,q))
print("   mean share %.4f over %d primes; best that is not all ones, by length: %s"%(s/tot,tot," ".join("%d:%d/%d(%d)"%(w,v[2],w,v[1]) for w,v in sorted(best.items()) if w>=8)))
print("4. RUNS: primes below 2^20 by the number of runs of equal bits")
runs={}
for q in primerange(3,N):
    b_=bin(q)[2:]; r=1+sum(1 for i in range(1,len(b_)) if b_[i]!=b_[i-1]); runs[r]=runs.get(r,0)+1
print("   "+" ".join("%d:%d"%(r,c) for r,c in sorted(runs.items())))
print("5. ALL ONES IN TWO BASES: numbers 1 + b + ... + b^(w-1), w >= 3, below 1e13")
L=10**13; seen={}; dup={}
b=2
while 1+b+b*b<L:
    v=1+b+b*b; w=3
    while v<L:
        if v in seen: dup.setdefault(v,[seen[v]]).append((b,w))
        else: seen[v]=(b,w)
        v=v*b+1; w+=1
    b+=1
print("   values %d, bases up to %d; values that occur twice: %s"%(len(seen),b-1,{v:x for v,x in dup.items()}))
