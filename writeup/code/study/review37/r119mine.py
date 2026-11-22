import numpy as np, math, json
from scipy import integrate
Z=np.loadtxt('/home/bas/praime/writeup/data/riemann/zeros_10000_flint.txt'); T=Z[-1]
print(len(Z),T)
def ei_asym(w):  # Ei(w) for complex w, |w| large, without the i*pi constant (cancels in 2Re)
    s=np.ones_like(w); term=np.ones_like(w)
    for k in range(1,12):
        term=term*k/w; s=s+term
    return np.exp(w)/w*s
def tailvar(x,H):
    d=H/x;u=math.log1p(d);s=math.sqrt(1+d)
    f=lambda g:((2+d)-2*s*math.cos(g*u))/(0.25+g*g)*math.log(g/(2*math.pi))/(2*math.pi)
    # integrate to big cutoff then analytic remainder of the non-oscillating part
    G=T+200*2*math.pi/u
    I=integrate.quad(f,T,G,limit=5000)[0]
    I+= (2+d)*(math.log(G/(2*math.pi))+1)/(2*math.pi*G)
    return 2*x*I/math.log(x)**2
meas=json.load(open('/home/bas/praime/writeup/data/study/r10_var.json'))
for H in (10**6,10**5):
  tot=[];sh=[]
  for b in range(10):
    Lb=10**9+b*10**8
    x=Lb+np.arange(10**8//H)*H
    S=np.zeros(len(x))
    rho=0.5+1j*Z
    for i0 in range(0,len(Z),1000):
        r=rho[i0:i0+1000]
        S+=-2*np.real(ei_asym(np.outer(np.log(x+H),r))-ei_asym(np.outer(np.log(x),r))).sum(axis=1)
    e=H/np.log(x+H/2)
    vl=np.var(S,ddof=1); vt=np.mean([tailvar(float(xx),H) for xx in (Lb,Lb+5e7,Lb+1e8)])
    tot.append((vl+vt)/e.mean()); sh.append(vt/(vl+vt))
  print(H,'zeros',round(np.mean(tot),4),'tail share',round(np.mean(sh),3),'measured(repo)',round(np.mean(meas[str(H)]['meas']),4),'per-block sd of measured',round(np.std(meas[str(H)]['meas'],ddof=1)/math.sqrt(10),4))
  print(' per block zeros',np.round(tot,4)); print(' per block meas ',np.round(meas[str(H)]['meas'],4))
