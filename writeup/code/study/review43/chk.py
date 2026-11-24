import numpy as np, sympy, random
for f,x0 in [('pr_3e11.npy',3*10**11),('pr_6e11.npy',6*10**11),('pr_9e11.npy',9*10**11)]:
    pr=np.load(f)
    random.seed(1)
    # check random sample of listed numbers are prime, and random gaps contain no primes
    ok=all(sympy.isprime(int(pr[i])) for i in random.sample(range(len(pr)),300))
    gaps_ok=True
    for i in random.sample(range(len(pr)-1),100):
        a,b=int(pr[i]),int(pr[i+1])
        if sympy.nextprime(a)!=b: gaps_ok=False
    print(f, ok, gaps_ok, int(pr[0])==sympy.nextprime(x0-1))
