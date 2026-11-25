import cypari2, time, math
pari=cypari2.Pari(); pari.allocatemem(400*10**6)
pari.set_real_precision(30)
E=pari.ellinit([0,0,1,2,2]); 
t0=time.time()
Ls=pari.lfunsympow(E,2)
print('sympow', time.time()-t0)
t0=time.time(); Li=pari.lfuninit(Ls,[2,0.05,8.5],2); print('init', time.time()-t0)
t0=time.time()
for k in range(50):
    s=pari('2+%s*I'%(0.1*k)); v=pari.lfun(Li,s,0); d1=pari.lfun(Li,s,1); d2=pari.lfun(Li,s,2)
print('50x3 evals', time.time()-t0, v, d1, d2)
# compare L(3 in PARI = analytic 2) with Euler product
import sympy
s=3
prod=1.0
N=2699
for p in sympy.primerange(2,200000):
    if N%p==0:
        ap=int(pari.ellap(E,p)); fac=1/(1-ap*ap/p**s)  # arithmetic: (1 - a_p^2 p^{-s})?
    else:
        ap=int(pari.ellap(E,p)); a=ap/math.sqrt(p)
        # arithmetic sym2 factor: (1 - alpha^2 X)(1 - p X)(1 - beta^2 X) with alpha*beta=p, alpha+beta=ap
        X=p**(-s)
        fac=1/((1-p*X)*(1-(ap*ap-2*p)*X+p*p*X*X))
    prod*=fac
print('Euler', prod, 'PARI', pari.lfun(Ls,3))
print(pari.lfunparams(Ls) if hasattr(pari,'lfunparams') else '')
