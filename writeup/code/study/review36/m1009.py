import cypari2, math, cmath, json, sys, time
import numpy as np, sympy
pari=cypari2.Pari(); pari.allocatemem(600*10**6)
q=1009; G=pari('znstar(1009,1)'); g0=int(pari('znstar(1009,1).gen[1]'))
ind={}; v=1
for e in range(q-1): ind[v]=e; v=v*g0%q
lo,hi=int(sys.argv[1]),int(sys.argv[2]); rows=[]; t=time.time()
pk=[]
for p in sympy.primerange(2,1000):
    x=p
    while x<1000: pk.append((x,p)); x*=p
for k in range(lo,hi):
    L=pari.lfuncreate(pari('[znstar(1009,1),[%d]]'%k))
    T=4; z=[]
    while not z: z=[float(s) for s in pari.lfunzeros(L,pari('[%d,%d]'%(-T,T)))]; T*=2
    gmin=min(abs(s) for s in z)
    r=complex(pari.lfun(L,1,1))/complex(pari.lfun(L,1))
    tr={}
    for X in (30,199,999):
        s=0
        for n,p in pk:
            if n>X: continue
            if p==q: continue
            s+= -(cmath.exp(2j*math.pi*k*ind[n%q]/(q-1))).real*math.log(p)/n
        tr[X]=s
    rows.append(dict(k=k,a=k%2,gmin=gmin,ReLL=r.real,tr=tr))
json.dump(rows,open(f'm1009_{lo}_{hi}.json','w')); print(hi-lo, time.time()-t)
