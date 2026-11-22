import mpmath as mp, json, math
from sympy import primerange, factorint
mp.mp.dps = 30
chi4 = [0,1,0,-1]; chi3=[0,1,-1]
def consts(chi): return mp.dirichlet(1,chi,1)/mp.dirichlet(1,chi), mp.dirichlet(0,chi,1)/mp.dirichlet(0,chi)
def T(y):
    y = mp.mpf(y); return mp.nsum(lambda k: y**(-(2*k+1))/((2*k+1)*(2*k+2)), [0, mp.inf])
def G(y, chi, q, d1):
    y = mp.mpf(y); s = mp.mpf(0)
    for p in primerange(2, int(y)+1):
        pk = p; k = 1
        while pk <= y:
            s += mp.mpf(chi[p % q])**k * mp.log(p) * (y/pk - 1); pk *= p; k += 1
    return s + y*d1
for chi, q, ys in ((chi4, 4, [1, 1.5, 2, 1511, 1511.5, 10**5+0.37]), (chi3, 3, [1, 2.5, 17509, 10**5+0.37])):
    d1, d0 = consts(chi); C = d1 - d0 + mp.log(2)
    for y in ys:
        g = G(y, chi, q, d1); D = g - d0 + T(y)
        print(q, y, 'G', mp.nstr(g, 12), 'ratio|G-d0+T|/sqrt', mp.nstr(abs(D)/mp.sqrt(y), 10), 'w/o T', mp.nstr(abs(g-d0)/mp.sqrt(y),10), 'C', mp.nstr(C, 8))
# zeros check chi4
zs = json.load(open('/home/bas/scratch/study-review33/zeros4.json'))
d1, d0 = consts(chi4)
for y in (1.5, 10, 50.5, 1511):
    y = mp.mpf(y)
    Zs = sum(2*mp.re(y**(mp.mpf(0.5)+1j*g))/(0.25+g*g) for g in zs)
    print('zeros formula y=%s' % y, mp.nstr(Zs + d0 - T(y), 10), 'direct', mp.nstr(G(y, chi4, 4, d1), 10), 'tail bound*sqrt y', mp.nstr(0.0093*mp.sqrt(y), 3))
# "in words" sentence: sum within 0.156 sqrt(y) of 0.78319 - 0.24561 y (T dropped)
bad = []
for i in range(0, 20001):
    y = 1 + i*0.005
    g = G(y, chi4, 4, d1) if y < 12 else None
    if y >= 12: break
    s = g - y*d1
    if abs(s - (d0 - y*d1)) > mp.mpf('0.155568')*mp.sqrt(y): bad.append(float(y))
print('in-words sentence fails for y in', (min(bad), max(bad)) if bad else None, 'count', len(bad))
# with T: full tube on [1, 12]
worst = 0
for i in range(0, 2201):
    y = 1 + i*0.005
    D = G(y, chi4, 4, d1) - d0 + T(y); worst = max(worst, abs(D)/mp.sqrt(y))
print('max |G-d0+T|/sqrt y on [1,12] chi4', mp.nstr(worst, 10))
