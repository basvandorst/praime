import mpmath as mp, numpy as np, math
def closed(d):
    return math.log(1/d) - 0.5772156649015329 - math.log(2*math.pi) + 2*math.log1p(d) + ((2+d)/d)*math.log1p(d/2)
Zall = np.array([float(l) for l in open('/home/bas/praime/writeup/data/riemann/zeros_10000_flint.txt') if l.strip()])
xg, wg = np.polynomial.legendre.leggauss(48)
def gl(f, a, b, m):
    e = np.linspace(a, b, m+1); lo, hi = e[:-1], e[1:]
    t = (hi-lo)[:, None]/2*xg[None, :] + (hi+lo)[:, None]/2
    return float(np.sum(f(t)*wg[None, :]*(hi-lo)[:, None]/2))
for d in (0.3, 0.1, 0.05, 0.01, 0.002):
    u = math.log1p(d); s = math.sqrt(1+d); per = 2*math.pi/u
    h = lambda t: ((2+d) - 2*s*np.cos(t*u))/(0.25+t*t)
    f = lambda t: h(t)*(np.log(t/(2*math.pi))/2 + 1/(48*t*t))/math.pi
    out = []
    for n in (2000, 5000, 10000):
        Z = Zall[:n]; T = Z[-1]; thT = float(mp.siegeltheta(T))
        S = 2*math.fsum(h(Z)); edge = h(T)*(thT/math.pi + 1 - n)
        K = int(4e6/per); b = T + K*per
        tail = gl(f, T, b, 4*K)
        tail += (2+d)/(2*math.pi)*(math.log(b/(2*math.pi))+1)/b
        out.append('%d: %.2e' % (n, (S + 2*(edge+tail))/d - closed(d)))
    print('d', d, ' '.join(out), flush=True)
