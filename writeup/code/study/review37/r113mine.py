import numpy as np, math, sympy
G=np.loadtxt('/home/bas/praime/writeup/data/riemann/lzeros/D-4.txt'); G=G[G<1500]; print('zeros',len(G))
rho=0.5+1j*G
pr=np.array(list(sympy.primerange(2,40000)))
def ei(w):
    s=np.ones_like(w);t=np.ones_like(w)
    for k in range(1,14): t=t*k/w; s=s+t
    return np.exp(w)/w*s
def Dx(x):
    x=np.atleast_1d(np.asarray(x,float))
    S=2*np.real(ei(np.outer(np.log(x),rho))).sum(axis=1)
    pis=np.searchsorted(pr,np.sqrt(x),side='right')
    return S+(pis-1)/2
def regions(lo,hi,n):
    xs=np.unique(np.round(np.geomspace(lo,hi,n)))
    v=Dx(xs); neg=v<0; out=[]
    i=0
    while i<len(xs):
        if neg[i]:
            j=i
            while j+1<len(xs) and neg[j+1]: j+=1
            # refine ends on integers
            a=int(xs[i-1]) if i>0 else int(xs[i]); b=int(xs[i])
            while b-a>1:
                m=(a+b)//2
                if Dx(m)[0]<0: b=m
                else: a=m
            st=b
            a=int(xs[j]); b=int(xs[j+1]) if j+1<len(xs) else int(xs[j])
            while b-a>1:
                m=(a+b)//2
                if Dx(m)[0]<0: a=m
                else: b=m
            out.append((st,a)); i=j+1
        else: i+=1
    return out
print(regions(20000,40000,4000))
print(regions(500000,800000,6000))
print(regions(100,1e9,60000))
