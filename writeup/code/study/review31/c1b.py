import mpmath as mp, numpy as np, math
mp.mp.dps = 25
def closed(d):
    d = mp.mpf(d)
    return mp.log(1/d) - mp.euler - mp.log(2*mp.pi) + 2*mp.log1p(d) + ((2+d)/d)*mp.log1p(d/2)
def integ(f, d, a):
    """int_a^inf h(r) f(r) dr, split into smooth and oscillating parts"""
    d = mp.mpf(d); u = mp.log1p(d); s = mp.sqrt(1+d)
    sm = mp.quad(lambda r: (2+d)*f(r)/(mp.mpf(1)/4+r*r), [a, a+10, a+100, a+1000, mp.inf])
    osc = mp.quadosc(lambda r: -2*s*mp.cos(r*u)*f(r)/(mp.mpf(1)/4+r*r), [a, mp.inf], omega=u)
    return sm + osc
def weil(d):
    d = mp.mpf(d)
    g0 = 2*integ(lambda r: 1, d, 0)/(2*mp.pi)
    dig = 2*integ(lambda r: mp.re(mp.digamma(mp.mpf(1)/4 + 0.5j*r)), d, 0)/(2*mp.pi)
    return g0, (2*d*mp.log1p(d) - d*mp.log(mp.pi) + dig)/d   # using g0 = d exactly
Z = [mp.mpf(l) for l in open('/home/bas/praime/writeup/data/riemann/zeros_10000_flint.txt') if l.strip()]
Zf = np.array([float(z) for z in Z])
def zerosum(d):
    u = math.log1p(d); s = math.sqrt(1+d)
    S = 2*math.fsum(((2+d) - 2*s*np.cos(Zf*u))/(0.25+Zf*Zf))
    T = Z[-1]; n = len(Z); dd = mp.mpf(d)
    h = lambda t: (2+dd - 2*mp.sqrt(1+dd)*mp.cos(t*mp.log1p(dd)))/(mp.mpf(1)/4+t*t)
    edge = h(T)*(mp.siegeltheta(T)/mp.pi + 1 - (n - mp.mpf(1)/2))
    tail = integ(lambda t: (mp.log(t/(2*mp.pi))/2 + 1/(48*t*t))/mp.pi, d, T)
    return (S + 2*float(edge + tail))/d, len(Z)
for d in (0.1, 0.05, 0.01, 0.002, 0.3, 0.6, 0.9, 0.99):
    g0, Dw = weil(d); c = closed(d)
    line = 'd %.3f g0 %s  Weil %s closed %s diff %.2e' % (d, mp.nstr(g0, 12), mp.nstr(Dw, 14), mp.nstr(c, 14), float(Dw-c))
    zs, n = zerosum(d); line += '  zeros(%d) %.12f z-c %.2e' % (n, zs, zs-float(c))
    print(line, flush=True)
