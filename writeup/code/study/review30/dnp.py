# Independent D(delta) with numpy Gauss-Legendre per period; decomposition D - D0 = (7/8) f(0) + P3 + P4.
import numpy as np, math, json
from scipy.special import loggamma
Z=np.array([float(l) for l in open('/home/bas/praime/writeup/data/riemann/zeros_9999_flint.txt') if l.strip()])
T=Z[-1]; NT=len(Z); B=1-0.5772156649015329-math.log(2*math.pi)
xg,wg=np.polynomial.legendre.leggauss(60)
def gl(F, edges):
    tot=0.0
    for i in range(0,len(edges)-1,20000):
        lo=edges[i:i+20000+1][:-1]; hi=edges[i+1:i+20000+1]
        m=(lo+hi)/2; h=(hi-lo)/2
        t=m[:,None]+h[:,None]*xg[None,:]
        tot+=float(np.sum(h*(F(t)@wg)))
    return tot
def mk(d):
    u=math.log1p(d); s=math.sqrt(1+d)
    f=lambda t: 2*((1+d)+1-2*s*np.cos(t*u))/(d*(0.25+t*t))
    fp=lambda t: 2*(2*s*u*np.sin(t*u)*(0.25+t*t)-((2+d)-2*s*np.cos(t*u))*2*t)/(d*(0.25+t*t)**2)
    return u,f,fp
def theta(t): return np.imag(loggamma(0.25+0.5j*t))-0.5*t*math.log(math.pi)
def N0(t): return t/(2*math.pi)*np.log(t/(2*math.pi*math.e))
def r(t): return theta(t)/math.pi+0.125-N0(t)
def smooth_int(F,a,per,tmax):
    # [a, first period edge], then periods to tmax
    k0=math.ceil(a/per)+1
    edges=np.concatenate([[a],np.arange(k0,int(tmax/per)+1)*per])
    return gl(F,edges), edges[-1]
out=[]
for d in (0.1,0.05,0.02,0.01,0.005,0.002):
    u,f,fp=mk(d); per=2*math.pi/u
    tmax=max(2e7,4000*per)
    Fs=lambda t: f(t)*np.log(t/(2*math.pi))/(2*math.pi)
    # D0 from 0: fine subdivision of first period
    e0=np.concatenate([np.geomspace(1e-12,1,40),np.linspace(1,per,200)[1:]])
    D0=gl(Fs,np.concatenate([[0.0],e0]))
    rest,tm=smooth_int(Fs,per,per,tmax)
    tailan=lambda a:(2/d)*(2+d)/(2*math.pi)*(math.log(a/(2*math.pi))+1)/a
    D0+=rest+tailan(tm)
    hl=math.log(1/d)+B+(1+d)*math.log1p(d)/d-1
    zs=float(np.sum(f(Z)))
    tl,tm=smooth_int(Fs,T,per,tmax); tl+=tailan(tm)
    # r' correction small; S(T+) boundary
    ST=NT-theta(np.array([T]))[0]/math.pi-1
    D_exact=zs+tl-f(np.array([T]))[0]*ST
    T0=T+0.5*2*math.pi/math.log(T/(2*math.pi))
    tl2,tm=smooth_int(Fs,T0,per,tmax); tl2+=tailan(tm)
    D_half=zs+tl2
    f0=float(f(np.array([0.0]))[0]); fpl=2*s if False else 2*math.sqrt(1+d)*u*u/d
    # P3 = -int_0^inf f' r dt
    Fr=lambda t: -fp(t)*r(t)
    e=np.concatenate([np.linspace(0,1,50),np.linspace(1,per,400)[1:]])
    P3=gl(Fr,e); P3b,_=smooth_int(Fr,per,per,tmax); P3+=P3b
    P4=D_exact-D0-0.875*f0-P3
    row=dict(delta=d,D_exact=D_exact,D_halfspacing=D_half,D0=D0,HL=hl,D0mHL=D0-hl,S_Tplus=ST,fT=float(f(np.array([T]))[0]),
             ratio=(D_exact-hl)/d,ratio_half=(D_half-hl)/d,c2=((D_exact-hl)/d-1.75)/d,
             P78f0=0.875*f0,P_refined=0.875*fpl,P3=P3,P3_pred=d*d/144,P4=P4,
             resid_vs_refined=D_exact-hl-0.875*fpl,resid_vs_f0=D_exact-hl-0.875*f0)
    print(json.dumps(row),flush=True); out.append(row)
json.dump(out,open('dnp.json','w'),indent=1)
