import cypari2, numpy as np
from sympy import primerange
pari=cypari2.Pari(); pari.allocatemem(300*10**6)
Z={}
for D in [-4,8,-8,-3,12]:
    L=pari.lfuncreate(D)
    Z[D]=np.array([float(t) for t in pari.lfunzeros(L,50)])
    print(D,len(Z[D]),Z[D][:2])
P=np.array(list(primerange(2,10**7+1)),dtype=np.int64)
lx=np.log(P.astype(float))
W={}
for D,z in Z.items():
    acc=np.zeros(len(P))
    for c in range(0,len(z),8):
        gg=z[c:c+8]; acc+=(np.exp(1j*np.outer(lx,gg))/(0.5+1j*gg)).real.sum(1)
    W[D]=1+2*acc
def kron(D,a):
    return int(pari.kronecker(D,a))
for q,Ds in [(8,[-4,8,-8]),(12,[-4,-3,12])]:
    others=[a for a in range(2,q) if np.gcd(a,q)==1]
    sel=P>=23
    predlast=np.ones(len(P),bool); minmarg=np.inf
    for a in others:
        d=sum(W[D] for D in Ds if kron(D,a)==-1)
        predlast&=d>0
        minmarg=min(minmarg,d[sel].min())
    print(q,'pred share last in [23,1e7]',predlast[sel].mean(),'min margin',minmarg)
    r=P%q
    c={a:np.cumsum(r==a) for a in [1]+others}
    strict=np.ones(len(P),bool); weak=np.ones(len(P),bool)
    for a in others: strict&=c[a]>c[1]; weak&=c[a]>=c[1]
    print(q,'meas strict share',strict[sel].mean(),'weak',weak[sel].mean())
    bad=P[sel&~strict]; print(' non-strict primes >=23:',bad[:40], len(bad))
    if len(bad): print(' last non-strict',bad[-1], 'any not weak',P[sel&~weak][:10])
    for p in bad[:3].tolist()+bad[-3:].tolist():
        i=np.searchsorted(P,p); print('  ',p,{a:int(c[a][i]) for a in c})
