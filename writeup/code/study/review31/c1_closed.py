"""Review 31, check 1: Weil pieces and closed form D(delta), against a zero sum with exact-theta tail."""
import mpmath as mp, numpy as np, math
mp.mp.dps = 30

def closed(d):
    d = mp.mpf(d)
    return mp.log(1/d) - mp.euler - mp.log(2*mp.pi) + 2*mp.log1p(d) + ((2+d)/d)*mp.log1p(d/2)

def hfun(r, d):
    u = mp.log1p(d)
    return (2 + d - 2*mp.sqrt(1+d)*mp.cos(r*u)) / (mp.mpf(1)/4 + r*r)

def weil(d):
    d = mp.mpf(d); u = mp.log1p(d)
    # h(i/2): limit
    hi = mp.limit(lambda e: hfun(mp.mpc(0, 0.5) + e, d), 0)
    g0 = mp.quad(lambda r: hfun(r, d), [0, mp.inf]) * 2 / (2*mp.pi)
    per = 2*mp.pi/u
    pts = [0] + [per*k for k in range(1, 400)]
    dig = 2*mp.quad(lambda r: hfun(r, d)*mp.re(mp.digamma(mp.mpf(1)/4 + 0.5j*r)), pts)
    # tail beyond pts[-1]: Re psi ~ ln(r/2), h ~ (2+d)/r^2 average
    R = pts[-1]
    dig += 2*mp.quad(lambda r: (2+d)/r**2*mp.log(r/2), [R, mp.inf])
    dig /= 2*mp.pi
    tot = 2*hi.real - g0*mp.log(mp.pi) + dig   # prime sum empty
    return hi, g0, tot/d

Z = np.array([float(l) for l in open('/home/bas/praime/writeup/data/riemann/zeros_10000_flint.txt') if l.strip()])
def zerosum(d):
    u = math.log1p(d); s = math.sqrt(1+d)
    w = lambda g: ((2+d) - 2*s*np.cos(g*u))/(0.25+g*g)
    S = 2*math.fsum(w(Z))
    T = mp.mpf(Z[-1]); n = len(Z)
    th = lambda t: mp.siegeltheta(t)
    # sum_{g>T} h = h(T)(theta(T)/pi + 1 - (n - 1/2)) + int_T^inf h theta'/pi
    dd = mp.mpf(d); per = 2*mp.pi/mp.log1p(dd)
    hh = lambda t: hfun(t, dd)
    edge = hh(T)*(th(T)/mp.pi + 1 - (n - mp.mpf(1)/2))
    # theta'(t) = 0.5 ln(t/2pi) + 1/(48 t^2)+...
    integ = mp.quad(lambda t: hh(t)*(mp.log(t/(2*mp.pi))/2 + 1/(48*t*t))/mp.pi, [T + k*per for k in range(0, 200)] + [mp.inf])
    return (S + 2*float(edge + integ))/d

for d in (0.1, 0.05, 0.01, 0.002, 0.3, 0.6, 0.9):
    hi, g0, Dw = weil(d)
    c = closed(d)
    line = 'd %.3f  h(i/2) %s vs d ln(1+d) %s  g0 %s  Weil %s  closed %s  diff %.2e' % (d, mp.nstr(hi.real, 12), mp.nstr(mp.mpf(d)*mp.log1p(d), 12), mp.nstr(g0, 12), mp.nstr(Dw, 14), mp.nstr(c, 14), float(Dw - c))
    if d <= 0.3:
        zs = zerosum(d); line += '  zeros %.12f  z-c %.2e' % (zs, zs - float(c))
    print(line, flush=True)
