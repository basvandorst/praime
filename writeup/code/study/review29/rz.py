# independent: R_zero via direct u-integration of S(u)^2 (not the closed-form double sum), and closed form with my own code
import numpy as np, math, sys, json
def load(n):
    if n<=2419: z=np.array([float(l) for l in open('/home/bas/praime/writeup/data/riemann/zeros_hp2419.txt') if l.strip()])
    else: z=np.loadtxt('zeros_10000.txt')
    return z[:n]
def coef(g,d):
    rho=0.5+1j*g
    return (np.exp(rho*math.log1p(d))-1)/rho
def tail_exact(T,d):
    # (1/d)*2*int_T^inf N'(t)|(1+d)^(1/2+it)-1|^2/(1/4+t^2) dt  with exact cos, via mpmath
    import mpmath as mp
    u=mp.log(1+d); s=mp.sqrt(1+d)
    f=lambda t: mp.log(t/(2*mp.pi))/(2*mp.pi)*((2+d)-2*s*mp.cos(u*t))/(mp.mpf(1)/4+t*t)
    return float(2/d*mp.quadosc(f,[T,mp.inf],omega=u))
def tail_avg(T,d): return (2+d)*(math.log(T/(2*math.pi))+1)/(math.pi*T*d)
def S_grid(g,c,u):
    out=np.zeros(len(u))
    for i in range(0,len(u),4000):
        uu=u[i:i+4000]
        out[i:i+4000]=2*np.real(np.exp(1j*np.outer(uu,g))@c)*np.exp(uu/2)  # S(x)/... we keep x^{1/2}
    return out
def R_direct(g,d,X,h=1e-4,span=14.0):
    c=coef(g,d); L=math.log(X)
    u=np.arange(L-span,L+h/2,h); u[-1]=L
    S=S_grid(g,c,u)                     # S(x)=sum c x^rho, x=e^u
    f=np.exp(u)*S**2                    # dx = e^u du ; integrand S(x)^2 dx
    w=np.ones(len(u)); w[0]=w[-1]=0.5   # simpson
    n=len(u)-1
    if n%2==0:
        w=np.ones(len(u)); w[1:-1:2]=4; w[2:-1:2]=2; w*=1/3
    I=h*np.sum(w*f)
    return 2*I/(d*X*X)
def R_closed(g,d,X):
    G=np.concatenate([g,-g]); c=coef(G,d); L=math.log(X); tot=0.0
    for i in range(0,len(G),400):
        D=G[i:i+400,None]-G[None,:]
        tot+=np.real(np.sum(c[i:i+400,None]*np.conj(c)[None,:]*(np.exp((2+1j*D)*L)-1)/(2+1j*D)))
    diag=np.sum(np.abs(c)**2)/d
    return 2*tot/(d*X*X), diag
if __name__=='__main__':
    g=load(2419); T=g[-1]
    for d in (0.1,0.01):
        print('d',d,'tail_avg',tail_avg(T,d),'tail_exact',tail_exact(T,d))
        for X in (1e8,1e9,1e10):
            rc,diag=R_closed(g,d,X); rd=R_direct(g,d,X)
            print(' X %.0e closed %.6f direct %.6f diag %.6f off %.6f  +tail_avg %.6f +tail_exact %.6f'%(X,rc,rd,diag,rc-diag,rc+tail_avg(T,d),rc+tail_exact(T,d)),flush=True)
