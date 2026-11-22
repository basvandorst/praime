import math, mpmath as mp
from scipy.integrate import quad
G=0.5772156649015329; B=1-G-math.log(2*math.pi)
# identity: int_0^inf ln t (1-cos at)/t^2 dt = (pi a/2)(1-gamma-ln a); use mpmath quad on [0,1] + quadosc tail, high dps
mp.mp.dps=25
for a in (0.5,1.0,3.0):
    f=lambda t: mp.log(t)*(1-mp.cos(a*t))/t**2
    v=mp.quad(f,[0,1])+mp.quadosc(f,[1,mp.inf],period=2*mp.pi/a)
    v0=mp.quad(lambda t:(1-mp.cos(a*t))/t**2,[0,1])+mp.quadosc(lambda t:(1-mp.cos(a*t))/t**2,[1,mp.inf],period=2*mp.pi/a)
    print('a',a,'ln-int',v,'formula',mp.pi*a/2*(1-mp.euler-mp.log(a)),' plain',v0,'pi a/2',mp.pi*a/2)
def D(d, lo=2*math.pi):
    u=math.log1p(d); s=math.sqrt(1+d)
    f=lambda t: math.log(t/(2*math.pi))/(2*math.pi)*((2+d)-2*s*math.cos(u*t))/(0.25+t*t)
    per=2*math.pi/u; tot=quad(f,lo,lo+per,limit=400)[0]; a=lo+per
    while a<3e7:
        b=a+per*max(1,int(a/per/20)); tot+=quad(f,a,b,limit=400)[0]; a=b
    tot+=(2+d)/(2*math.pi)*(math.log(a/(2*math.pi))+1)/a
    return 2/d*tot
for d in (1e-1,3e-2,1e-2,3e-3,1e-3,1e-4):
    v=D(d); print('d',d,'D-ln(1/d)',v-math.log(1/d),' minus B',v-math.log(1/d)-B, ' (from 0 instead of 2pi):',D(d,lo=1e-9)-math.log(1/d)-B)
