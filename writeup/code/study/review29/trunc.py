from rz import *
import json
prim={(0.1,1e8):0.9663146397385681,(0.1,1e9):1.2301149861229637,(0.1,1e10):1.2756640965621453,(0.01,1e8):2.9523712815587575,(0.01,1e9):3.179742624580468,(0.01,1e10):3.2544139014734137}
out=[]
for N in (2419,5000,10000):
    g=load(N); T=g[-1]
    for d in (0.1,0.01):
        te=tail_exact(T,d)
        for X in (1e8,1e9,1e10):
            rc,diag=R_closed(g,d,X)
            r=rc+te
            print('N %5d T %.0f d %.2f X %.0e  R %.5f (diag+tail %.5f off %+.5f) primes %.5f diff %+.5f'%(N,T,d,X,r,diag+te,rc-diag,prim[(d,X)],r-prim[(d,X)]),flush=True)
            out.append(dict(N=N,d=d,X=X,R=r,diag=diag+te,off=rc-diag,diff=r-prim[(d,X)]))
json.dump(out,open('trunc.json','w'),indent=1)
