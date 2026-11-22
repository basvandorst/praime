from rz import *
from tailchk import tail_num
prim={(3.268e9,0.1):0.8526746168469871,(3.268e9,0.01):3.110779664568085,(3.9e9,0.1):1.4809370680069456,(3.9e9,0.01):3.0299813372317757,(2e10,0.1):1.1261215316290516,(2e10,0.01):2.984755431005265}
for N in (2419,10000):
    g=load(N); T=g[-1]
    for d in (0.1,0.01):
        te=tail_num(T,d); ta=tail_avg(T,d)
        for X in (3.268e9,3.9e9,2e10):
            rc,diag=R_closed(g,d,X)
            p=prim[(X,d)]
            print('N %5d d %.2f X %.4g  zeros(avg tail) %.4f  zeros(exact tail) %.4f  primes %.4f  diff_exact %+.4f  off %+.4f'%(N,d,X,rc+ta,rc+te,p,p-(rc+te),rc-diag),flush=True)
