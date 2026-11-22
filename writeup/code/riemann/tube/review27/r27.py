import math, mpmath as mp, numpy as np
mp.mp.dps = 30
g = float(mp.euler); C = 2 + g - math.log(4*math.pi)
print("C", C, "check 2*sum 1/rho:", float(2*(1+mp.euler/2-mp.log(4*mp.pi)/2)))
# 1. table recompute with independent formula, and exclusion heights incl. peak-spacing correction
def ly(b,t):
    w = mp.mpc(b,t)*(1-mp.mpc(b,t)); return float(mp.log10(C*abs(w))/(b-0.5))
for t in (14.134725, 1e3, 1e6, 3e12, 3.000175332800e12):
    print(t, [round(ly(b,t),2) for b in (0.51,0.6,0.75,0.9,1.0)])
# rest bound: on-line rest <= (C - 4 Re 1/w) sqrt y, partners <= 2 y^(1-b)/|w| <= 2 sqrt y/|w|
for b in (0.51,0.75,1.0):
    for t in (14.13,100,1e3):
        w = complex(b,t)*(1-complex(b,t))
        print("rest<=C sqrt y?", b, t, C-4*(1/w).real+2/abs(w) <= C)
def bis(f,a=1.0,b=1e9):
    for _ in range(200):
        m=(a*b)**0.5
        (a,b)=(m,b) if f(m)<0 else (a,m)
    return a
for Y in (1e8,1e18):
    for b in (0.75,1.0):
        f = lambda t: ly(b,t)*math.log(10) + 2*math.pi/t - math.log(Y)
        t0 = bis(lambda t: ly(b,t)*math.log(10) - math.log(Y))
        t1 = bis(f)
        print("Y",Y,"beta",b,"height w/o peak corr",round(t0,2),"with",round(t1,2))
# zero-free region (Mossinghoff-Trudgian 2015): beta < 1 - 1/(5.558691 ln t)
t=3e12; bmax = 1-1/(5.558691*math.log(t)); print("bmax at 3e12", bmax, "log10 y*", ly(bmax,t))
# Schoenfeld-type analogue: amplitude 2x^b/|rho| vs sqrt(x) ln^2 x/(8 pi): threshold for b=1,0.75 at t=3e12
for b in (1.0,0.75):
    x = bis(lambda L: (b-0.5)*L - float(mp.log(t*L**2/(8*mp.pi))), 10, 1e5)   # L = ln x ; x^(b-1/2) = |rho| ln^2x /(8pi)
    print("Schoenfeld-type threshold b=",b," x ~ 10^", float(x/mp.log(10)))
