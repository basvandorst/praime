import mpmath as mp, numpy as np, math
from scipy import integrate
mp.mp.dps = 20
def closed(d):
    return math.log(1/d) - 0.5772156649015329 - math.log(2*math.pi) + 2*math.log1p(d) + ((2+d)/d)*math.log1p(d/2)
Z = np.array([float(l) for l in open('/home/bas/praime/writeup/data/riemann/zeros_10000_flint.txt') if l.strip()])
n = len(Z); T = Z[-1]; thT = float(mp.siegeltheta(mp.mpf(T)))
for d in (0.3, 0.1, 0.05, 0.01, 0.002):
    u = math.log1p(d); s = math.sqrt(1+d)
    h = lambda t: ((2+d) - 2*s*np.cos(t*u))/(0.25+t*t)
    S = 2*math.fsum(h(Z))
    edge = h(T)*(thT/math.pi + 1 - (n - 0.5))
    f = lambda t: h(t)*(math.log(t/(2*math.pi))/2 + 1/(48*t*t))/math.pi
    per = 2*math.pi/u; a = T; tail = 0.0
    # integrate to a large bound in chunks of whole periods
    while a < 4e6:
        b = a + per*max(1, int(2000/per))
        tail += integrate.quad(f, a, b, limit=400, epsabs=0, epsrel=1e-13)[0]; a = b
    tail += (2+d)/(2*math.pi)*(math.log(a/(2*math.pi))+1)/a
    zs = (S + 2*(edge + tail))/d
    zs_noedge = (S + 2*tail)/d
    print('d %.3f closed %.10f zeros+tail %.10f diff %.1e  (without S(T) boundary term: diff %.1e)' % (d, closed(d), zs, zs-closed(d), zs_noedge-closed(d)), flush=True)
print('S(T) at last zero (midpoint conv):', n - 0.5 - thT/math.pi - 1)
