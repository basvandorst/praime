from rz import *
import mpmath as mp
from scipy.integrate import quad
def tail_num(T,d,Tmax=1e9):
    u=math.log1p(d); s=math.sqrt(1+d)
    f=lambda t: math.log(t/(2*math.pi))/(2*math.pi)*((2+d)-2*s*math.cos(u*t))/(0.25+t*t)
    # split into periods
    tot=0; a=T; per=2*math.pi/u
    while a<2e6:
        b=a+per*20; tot+=quad(f,a,b,limit=200)[0]; a=b
    # beyond: average
    tot+= (2+d)/(2*math.pi)*(math.log(a/(2*math.pi))+1)/a
    return 2/d*tot
g=load(10000)
for d in (0.1,0.01):
    c=coef(g,d); 
    for N,T in ((2419,g[2418]),(10000,g[-1])):
        diag=2*np.sum(np.abs(c[:N])**2)/d
        print(d,N,'diag',diag,'tail_avg',tail_avg(T,d),'tail_mp',tail_exact(T,d),'tail_quad',tail_num(T,d))
