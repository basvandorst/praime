import mpmath as mp
from sympy import primerange
mp.mp.dps = 20
for chi, q in (([0,1,0,-1],4), ([0,1,-1],3)):
    d1 = mp.dirichlet(1,chi,1)/mp.dirichlet(1,chi); d0 = mp.dirichlet(0,chi,1)/mp.dirichlet(0,chi)
    terms = [(p**k, mp.mpf(chi[p%q])**k*mp.log(p)) for p in primerange(2,100) for k in range(1,8) if p**k <= 100]
    bad = []
    for i in range(0, 99*400+1):
        y = 1 + i/mp.mpf(400)
        s = sum(a*(y/n-1) for n,a in terms if n <= y)
        if abs(s - (d0 - d1*y)) > mp.mpf('0.156')*mp.sqrt(y): bad.append(float(y))
    print(q, 'in-words (0.156, T dropped) fails on', (bad[0], bad[-1]) if bad else None, len(bad))
    print(q, 'identity check d0 vs -d1 - ln(q/pi) + gamma + ln2:', mp.nstr(d0 - (-d1 - mp.log(q/mp.pi) + mp.euler + mp.log(2)), 5))
