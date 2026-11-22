import mpmath as mp
mp.mp.dps = 30
def integ(f, d, a=0):
    d = mp.mpf(d); u = mp.log1p(d); s = mp.sqrt(1+d)
    per = 2*mp.pi/u
    pts = [a] + [per*k/4 for k in range(1, 9)]
    sm = mp.quad(lambda r: (2+d)*f(r)/(mp.mpf(1)/4+r*r), pts + [per*10, per*100, mp.inf])
    osc = mp.quad(lambda r: -2*s*mp.cos(r*u)*f(r)/(mp.mpf(1)/4+r*r), pts) + mp.quadosc(lambda r: -2*s*mp.cos(r*u)*f(r)/(mp.mpf(1)/4+r*r), [pts[-1], mp.inf], omega=u)
    return sm + osc
for d in ('0.1', '0.05', '0.02', '0.01'):
    d = mp.mpf(d)
    D0 = 2/d*integ(lambda t: mp.log(t/(2*mp.pi))/(2*mp.pi), d)
    HL = mp.log(1/d) - mp.euler - mp.log(2*mp.pi) + (1+d)*mp.log1p(d)/d
    thp = 2/d*integ(lambda t: (mp.re(mp.digamma(mp.mpf(1)/4 + 0.5j*t)) - mp.log(mp.pi))/(2*mp.pi), d)
    D = mp.log(1/d) - mp.euler - mp.log(2*mp.pi) + 2*mp.log1p(d) + ((2+d)/d)*mp.log1p(d/2)
    print(mp.nstr(d, 3), 'D0-HL', mp.nstr(D0-HL, 10), ' /d^2', mp.nstr((D0-HL)/d**2, 8), ' (-1/144 = %.6f)' % (-1/144), ' D-(2ln(1+d)+theta part)', mp.nstr(D - 2*mp.log1p(d) - thp, 5))
