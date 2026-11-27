import re, sys, numpy as np, mpmath as mp
from sympy import primerange
mp.mp.dps = 50
def load(f):
    lines = open('/home/bas/praime/writeup/data/riemann/odlyzko_zeros%d.txt' % f).read().split('\n')
    return np.array([float(x) for x in lines if re.fullmatch(r'\s*\d+\.\d+\s*', x)])
OFF = {3: (mp.mpf('267653395647'), 10**12), 4: (mp.mpf('144176897509546973000'), 10**21), 5: (mp.mpf('1370919909931995300000'), 10**22)}
def S_mid(f):
    """S(t) at the midpoints between consecutive zeros: S = N(t) - theta(t)/pi - 1, N = index of the zero below."""
    T0, n0 = OFF[f]; g = load(f); mid = (g[:-1] + g[1:])/2
    # theta(T0 + x) = theta(T0) + theta'(T0) x + theta''(T0) x^2/2 + ... (x <= 3000 << T0): use mp for the expansion point
    th0 = mp.siegeltheta(T0); d1 = mp.log(T0/(2*mp.pi))/2; d2 = 1/(2*T0); d3 = -1/(2*T0**2)
    k0 = mp.floor(th0/mp.pi); frac0 = float(th0/mp.pi - k0)
    x = mid; th_over_pi = frac0 + (float(d1)*x + float(d2)*x*x/2 + float(d3)*x**3/6)/np.pi
    base = n0 + 1 - int(k0)                                # exact integers first (n0 is 10^21 or 10^22 in two files)
    S = (base + np.arange(len(mid))) - th_over_pi - 1
    return T0, mid, S
def P(T0, x, y):
    out = np.zeros(len(x)); var = 0.0
    for p in primerange(2, y + 1):
        lp = mp.log(p); ph0 = float(mp.fmod(T0*lp, 2*mp.pi)); lpf = float(lp); m = 1
        while p**m <= y:
            out -= np.sin(m*(ph0 + x*lpf))/(m*p**(m/2)); var += 1/(m*m*p**m); m += 1
    return out/np.pi, var/(2*np.pi**2)
if __name__ == '__main__':
    f = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    T0, x, S = S_mid(f)
    print('zeros%d: S at %d midpoints: mean %.4f var %.4f (min %.2f max %.2f)' % (f, len(S), S.mean(), S.var(), S.min(), S.max()))
    for y in (10, 100, 1000, 10**4, 10**5):
        Py, vth = P(T0, x, y); slope = np.cov(S, Py)[0, 1]/Py.var(); rem = (S - Py).var()
        print('  y = %6d: Var P_y data %.4f theory %.4f | slope %.3f | variance removed %.4f | remainder var %.4f | explained %.1f %%' % (y, Py.var(), vth, slope, S.var() - rem, rem, 100*(1 - rem/S.var())))
