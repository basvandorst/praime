import cypari2, math, random, sys, json, gmpy2
import numpy as np
from scipy.stats import mannwhitneyu
pari=cypari2.Pari(); pari.allocatemem(600*10**6)
seed=int(sys.argv[1]); NG=int(sys.argv[2]) if len(sys.argv)>2 else 15
M=3_000_000
s=np.ones(M+1,bool); s[:2]=False
for p in range(2,int(M**.5)+1):
    if s[p]: s[p*p::p]=False
def g1(d):
    Lf=pari.lfuncreate(d); T=3*2*math.pi/math.log(abs(d)); z=[]
    while not z:
        z=[float(t) for t in pari.lfunzeros(Lf,T) if float(t)>1e-9]; T*=2
    return min(z)*math.log(abs(d))/(2*math.pi)
Y=np.arange(-80,81)
def best(d):
    b=0;arg=None
    for A in range(1,301):
        for B in range(0,2*A):
            if (B*B-d)%(4*A): continue
            C=(B*B-d)//(4*A)
            v=np.abs(A*Y*Y+B*Y+C); assert v.max()<=M
            pr=s[v]
            # longest run of True
            r=0;m=0
            for t in pr:
                r=r+1 if t else 0
                if r>m:m=r
            if m>b:b=m;arg=(A,B,C)
    return b,arg
rng=random.Random(seed); ds=set()
while len(ds)<600:
    d=rng.randint(10**5,10**6)
    if d%8!=5 or int(pari.issquarefree(d))!=1: continue
    if all(gmpy2.kronecker(d,p)==-1 for p in (3,5,7,11)): ds.add(d)
ds=sorted(ds); gs={d:g1(d) for d in ds}
o=sorted(ds,key=lambda d:gs[d]); low=o[:NG]; high=o[-NG:]
L=[best(d) for d in low]; H=[best(d) for d in high]
lr=[x[0] for x in L]; hr=[x[0] for x in H]
p=mannwhitneyu(lr,hr,alternative='greater').pvalue
print('seed',seed,'LOW',sorted(lr),'med',np.median(lr),'HIGH',sorted(hr),'med',np.median(hr),'p %.2e'%p)
print('low g range %.3f-%.3f high %.3f-%.3f'%(gs[low[0]],gs[low[-1]],gs[high[0]],gs[high[-1]]))
mx=max(zip(lr,low)); print('max low',mx,gs[mx[1]],L[lr.index(mx[0])][1])
from scipy.stats import spearmanr
sub=rng.sample(ds,60); br=[best(d)[0] for d in sub]
print('spearman g vs best run over 60 random filtered d: %.3f p %.2e'%spearmanr([gs[d] for d in sub],br))
import sympy
def fnp(d):
    for p in sympy.primerange(13,2000):
        if gmpy2.kronecker(d,p)!=-1: return p
allD=sub+low+high; allR=br+lr+hr
F=[fnp(d) for d in allD]; Gg=[gs[d] for d in allD]
print('pooled n=%d: spearman run~g %.3f, run~first non-inert prime %.3f, g~fnp %.3f'%(len(allD),spearmanr(Gg,allR)[0],spearmanr(F,allR)[0],spearmanr(Gg,F)[0]))
print('sub only: run~fnp %.3f run~g %.3f'%(spearmanr([fnp(d) for d in sub],br)[0],spearmanr([gs[d] for d in sub],br)[0]))
print('LOW fnp',[fnp(d) for d in low]); print('HIGH fnp',[fnp(d) for d in high])
mid=o[NG:-NG]; R=[best(d)[0] for d in rng.sample(mid,15)]
print('MIDDLE 15 random from rest',sorted(R),'med',np.median(R),'p LOW>MID %.2e'%mannwhitneyu(lr,R,alternative='greater').pvalue)
json.dump(dict(seed=seed,low=list(zip(low,lr)),high=list(zip(high,hr)),g={str(d):gs[d] for d in ds}),open('r126_%d.json'%seed,'w'))
