import mpmath as mp, numpy as np, math
from scipy import integrate
def closed(d):
    return math.log(1/d) - 0.5772156649015329 - math.log(2*math.pi) + 2*math.log1p(d) + ((2+d)/d)*math.log1p(d/2)
Zall = np.array([float(l) for l in open('/home/bas/praime/writeup/data/riemann/zeros_10000_flint.txt') if l.strip()])
print(len(Zall), Zall[-1])
def tailint(f, a, per, stop=4e6):
    t = 0.0
    while a < stop:
        b = a + per*max(1, int(2000/per)); t += integrate.quad(f, a, b, limit=400, epsabs=0, epsrel=1e-13)[0]; a = b
    return t, a
for d in (0.1, 0.01):
    u = math.log1p(d); s = math.sqrt(1+d)
    h = lambda t: ((2+d) - 2*s*np.cos(t*u))/(0.25+t*t)
    per = 2*math.pi/u
    for n in (2000, 5000, 8000, 10000):
        Z = Zall[:n]; T = Z[-1]; thT = float(mp.siegeltheta(T))
        S = 2*math.fsum(h(Z)); edge = h(T)*(thT/math.pi + 1 - (n - 0.5))
        f = lambda t: h(t)*(math.log(t/(2*math.pi))/2 + 1/(48*t*t))/math.pi
        tail, a = tailint(f, T, per)
        tail += (2+d)/(2*math.pi)*(math.log(a/(2*math.pi))+1)/a
        print(d, n, '%.3e' % ((S + 2*(edge+tail))/d - closed(d)), flush=True)
