import numpy as np, math
K=21; N=1<<K
def sieve(n):
    s=np.ones(n+1,dtype=bool); s[:2]=False
    for i in range(2,int(n**0.5)+1):
        if s[i]: s[i*i::i]=False
    return s
s=sieve(N)[:N]; n=np.arange(N//2+1,N,2,dtype=np.int64); y=s[n]; mu=y.mean()
def H(p): return 0 if p<=0 or p>=1 else -p*math.log2(p)-(1-p)*math.log2(1-p)
def info(v,bins):
    c=np.bincount(v,minlength=bins).astype(float); c1=np.bincount(v,weights=y.astype(float),minlength=bins); ok=c>0
    p=c1[ok]/c[ok]; h=-(np.where(p>0,p*np.log2(np.where(p>0,p,1)),0)+np.where(p<1,(1-p)*np.log2(np.where(p<1,1-p,1)),0))
    return H(mu)-float((c[ok]/len(v)*h).sum())-(ok.sum()-1)/(2*len(v)*math.log(2))
def factors(m):
    f=[];d=2
    while d*d<=m:
        while m%d==0:
            if d not in f:f.append(d)
            m//=d
        d+=1
    if m>1 and m not in f:f.append(m)
    return f
print("odd numbers with %d bits: %d, primes %.4f; information about 'prime' in the folded number (bits, bias removed); a full knowledge is %.4f bits"%(K,len(n),mu,H(mu)))
print("fold in w layers = add the w-bit pieces of n = n mod (2^w - 1), without division")
for w in range(2,15):
    M=(1<<w)-1; v=np.zeros(len(n),dtype=np.int64); t=n.copy()
    while t.max()>0: v+=t&M; t>>=w
    while v.max()>M: v=(v&M)+(v>>w)
    v[v==M]=0
    fs=factors(M); pred=sum(-math.log2(1-1/q)*1 for q in fs if q>2)
    I=info(v,M+1)
    print("  w=%2d  2^w-1 = %-6d = %-14s information %.4f   (if only its prime factors matter: %.4f)"%(w,M,"·".join(map(str,fs)),I,H(mu)-( (lambda r: r*H(min(1,mu/r)))(np.prod([1-1/q for q in fs])) )))
