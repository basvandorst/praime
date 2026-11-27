import sys, time, numpy as np, cypari2
from scipy.special import loggamma
from sympy import primerange
from loop9 import cth
pari = cypari2.Pari(); pari.allocatemem(600*10**6)
q, a, T = 4, 1, float(sys.argv[1])
chi = lambda n: 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)
g = np.array([float(z) for z in pari.lfunzeros(pari.lfuncreate(pari('Mod(3,4)')), T)])
theta = lambda t: t/2*np.log(q/np.pi) + np.imag(loggamma((0.5 + a + 1j*np.asarray(t, float))/2))
th = theta(g)
idx = np.arange(1, len(g) + 1.0); off = np.mean(idx - 0.5 - th/np.pi)
L = np.log(q*g/(2*np.pi)); sp = 2*np.pi/L
print('%d zeros of L(s, chi_-4) up to %.0f; mean of n - 1/2 - theta/pi = %.4f; rms with no primes %.4f spacings' % (len(g), T, off, np.std(idx - 0.5 - th/np.pi)))
for ucut in (0.5, 0.7, 0.9):
    nmax = int(np.exp(ucut*L.max())); A, F = [], []
    for p in primerange(2, nmax + 1):
        m = 1
        while p**m <= nmax:
            if chi(p): A.append(chi(p)**m/(np.pi*m*p**(m/2))); F.append(m*np.log(p))
            m += 1
    A = np.array(A); F = np.array(F)
    def Sp(y):
        U = F[None, :]/np.log(q*y/(2*np.pi))[:, None]
        C = np.where(U < ucut, cth(U.ravel()).reshape(U.shape), 0.0)
        return -((np.sin(np.outer(y, F))*C) @ A)
    y = g*0 + 2*np.pi*idx/np.log(idx + 10)
    for it in range(40):                                         # solve theta(y)/pi + Sp(y) = n - 1/2 (Newton on the smooth part, damped)
        r = theta(y)/np.pi + Sp(y) - (idx - 0.5)
        y = np.maximum(y - 0.7*r/(np.log(q*np.maximum(y, 3)/(2*np.pi))/(2*np.pi)), 0.5)
    e = (y - g)/sp; keep = slice(len(g)//10, None)                # drop the lowest tenth (L too small there)
    print('  cut %.1f (n <= %d, %d terms): rms error %.4f spacings' % (ucut, nmax, len(A), np.sqrt(np.mean(e[keep]**2))), flush=True)
