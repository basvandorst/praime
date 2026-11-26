# Check the Guinand-Weil formula in the form used for the bootstrap, Gaussian test function.
from mpmath import mp, mpf, exp, log, sqrt, pi, euler, quad, cosh, zetazero, inf
from sympy import primerange
mp.dps = 30
s = mpf('0.45')
g = lambda u: exp(-u*u/(2*s*s))
h = lambda t: s*sqrt(2*pi)*exp(-s*s*t*t/2)
zs = 0
for k in range(1, 40):
    gam = zetazero(k).imag
    term = 2*h(gam); zs += term
    if term < mpf(10)**(-28): break
pole = quad(lambda u: 2*g(u)*cosh(u/2), [-inf, 0, inf])
arch = -(euler+log(pi))*g(0) + quad(lambda u: (2*exp(-2*u)*g(0) - 2*exp(-u/2)*g(u))/(1-exp(-2*u)), [0, 1, 5, inf])
pr = 0
for p in primerange(2, 3000):
    q = p
    while q < 10**6:
        pr += 2*log(p)/sqrt(mpf(q))*g(log(mpf(q))); q *= p
print('zeros ', zs)
print('rhs   ', pole + arch - pr)
print('pole', pole, 'arch', arch, 'primes', pr)
