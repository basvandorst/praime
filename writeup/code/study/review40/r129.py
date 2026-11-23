import numpy as np, sympy
from scipy.integrate import quad
from scipy.special import sici

t = np.loadtxt('/home/bas/praime/writeup/data/riemann/zeros_40000_flint.txt')
lt = np.log(t/(2*np.pi))
sel = (lt >= 7.6) & (lt <= 8.6)
T = t[sel]
print('n zeros', sel.sum(), T[0], T[-1])

def Nsm(t):
    return t/(2*np.pi)*np.log(t/(2*np.pi*np.e)) + 7/8.

def sigma2(x, l, ngrid=400000, about='l'):
    # x sorted unfolded points; sample start y uniformly in [x0, x1-l]
    y = np.linspace(x[0], x[-1]-l, ngrid)
    n = np.searchsorted(x, y+l, side='right') - np.searchsorted(x, y, side='left')
    if about == 'l':
        return np.mean((n-l)**2), n
    return np.var(n), n

primes = np.array(list(sympy.primerange(2, 10**6)), dtype=float)

def berry(l, L, cut=1.0):
    s = 0.0
    lp = np.log(primes)
    k = 1
    while True:
        m = k*lp <= cut*L
        if not m.any(): break
        tau = k*lp[m]/L
        s += np.sum(np.sin(np.pi*l*tau)**2/(k**2*primes[m]**k))
        k += 1
    # int_cut^inf sin^2(pi l tau)/tau^2
    a = 2*np.pi*l*cut
    # int_c^inf sin^2(b tau)/tau^2 dtau, b=pi l : = (1/2)int (1-cos 2b tau)/tau^2
    b = np.pi*l
    c = cut
    si, ci = sici(2*b*c)
    # int_c^inf (1-cos(2bt))/(2t^2) dt = 1/(2c) - (1/2)[cos(2bc)/c - 2b*(pi/2 - Si(2bc))]
    I = 1/(2*c) - 0.5*(np.cos(2*b*c)/c - 2*b*(np.pi/2 - si))
    return 2/np.pi**2*(s+I)

def gue(l):
    f = lambda s: (l-s)*(np.sinc(s))**2
    return l - 2*quad(f, 0, l, limit=400)[0]

# check integral formula numerically
l=2.0
num = quad(lambda x: np.sin(np.pi*l*x)**2/x**2, 1, 2000, limit=20000)[0] + 0.5/2000
print('int check', num, (berry(l, 1e-9)*np.pi**2/2))

x = Nsm(T)
Lmean = np.mean(np.log(T/(2*np.pi)))
print('mean L over zeros', Lmean)
ls = [1, 2, 4.5, 6, 8]
for l in ls:
    s2l, n = sigma2(x, l)
    s2v, _ = sigma2(x, l, about='mean')
    # L-averaged Berry: average over unfolded position i.e. over zeros
    Lgrid = np.log(T[::200]/(2*np.pi))
    bav = np.mean([berry(l, LL) for LL in Lgrid])
    print(f'l={l}: data(about l)={s2l:.4f} var={s2v:.4f} meanN-l={n.mean()-l:+.4f} | Berry@meanL={berry(l,Lmean):.4f} Berry-avg={bav:.4f} '
          f'Berry@7.6={berry(l,7.6):.4f} @8.6={berry(l,8.6):.4f} | GUE={gue(l):.4f}')
