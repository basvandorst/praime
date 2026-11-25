# BK (2013) pair correlation at small r for primitive chi mod k, my own implementation from the published formula
# (bk/derivation.md section 0 & 3 for the formula statement). Predicts the small-gap rate P(s<s0) relative to the sine kernel.
import os, numpy as np, mpmath as mp, math, json
from sympy import primerange, primefactors
from numpy.polynomial.legendre import leggauss
mp.mp.dps=30
PR=np.array(list(primerange(2,200000)),dtype=float); lnp=np.log(PR)
EPS=np.linspace(1e-3,0.30,900)
# k-independent pieces on the eps grid
A=[]  # Re(zeta'/zeta)'(1+i e) + 1/e^2
Zm=[] # |zeta(1+ie)|^2  (complex factor times Phi_off handled below)
for e in EPS:
    s=mp.mpc(1,e); z=mp.zeta(s); z1=mp.zeta(s,derivative=1); z2=mp.zeta(s,derivative=2)
    A.append(float(mp.re(z2/z-(z1/z)**2))+1/e**2)
    Zm.append(complex(z))
A=np.array(A); Zc=np.array(Zm)
ZABS2=np.abs(Zc)**2
def phioff(e):
    pe=np.exp(1j*e*lnp); return np.prod(1-(pe-1)**2/(PR-1)**2)
PHI=np.array([phioff(e) for e in EPS])
ZREG=ZABS2*PHI-1/EPS**2
def d2lnphidiag(e):
    m=np.arange(2,40)[:,None]; return float(np.sum(2*lnp**2*np.sum((m-1)*np.cos(m*e*lnp)/PR**m,axis=0)))
D2PHI=np.array([d2lnphidiag(e) for e in EPS])
def Rpair(r,L,S):
    e=2*np.pi*r/L
    a=np.interp(e,EPS,A); zr=np.interp(e,EPS,ZREG.real)+1j*np.interp(e,EPS,ZREG.imag)
    d2=np.interp(e,EPS,D2PHI)
    psi_off=np.ones_like(e,dtype=complex); d2psi=np.zeros_like(e)
    for p in S:
        lp=math.log(p); m=np.arange(1,60)[:,None]
        d2psi+=2*lp**2*np.sum(np.cos(m*e*lp)/p**m,axis=0)
        psi_off*=1/(1+(p**(1j*e/2)-p**(-1j*e/2))**2/(p-p**(-1j*e)))
    C=(2/L**2)*(a+np.real(np.exp(2j*np.pi*r)*(zr*psi_off+(psi_off-1)/e**2)))-(1/L**2)*(d2+d2psi)
    return 1-np.sinc(r)**2+C
if __name__=='__main__':
    # sanity: C(r) small-r behaviour, and the leading expansion pi^2 r^2 [1/3 - 4 Lam/L^2 - 8Q/L^3]
    for S in ([],[2],[3,5]):
        L=10.0; r=np.array([0.02,0.05,0.1,0.2,0.3,0.5,1.0])
        Lam=1.573151+sum(math.log(p)**2/(p-1) for p in S); Q=2.315846-sum(math.log(p)**3/(p-1)**2 for p in S)
        print(S,'R/(sine R)',np.round(Rpair(r,L,S)/(1-np.sinc(r)**2),4),' leading', round(1-12*Lam/L**2-24*Q/L**3,4))
    # sine CDF and its R2 integral (for the higher-correlation correction)
    def E0(s,m=40):
        x,w=leggauss(m); x=(x+1)*s/2; w=w*s/2
        K=np.sinc(x[:,None]-x[None,:]); Am=np.sqrt(w)[:,None]*K*np.sqrt(w)[None,:]
        return np.linalg.det(np.eye(m)-Am)
    def cdf(s0,h=1e-5): return 1+(E0(s0+h)-E0(s0-h))/(2*h)
    S0=[0.05,0.1,0.15,0.2,0.3]
    kap={s0:cdf(s0)/np.trapezoid(1-np.sinc(np.linspace(0,s0,400))**2,np.linspace(0,s0,400)) for s0 in S0}
    print('sine CDF / int R2:',{k:round(v,4) for k,v in kap.items()})
    LZ=os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/data/riemann/lzeros/'
    SETS={'B':(['r134_B.npz'],100,600),'Dn':(['r136_Dn.npz'],100,600),'Fn':(['r138_Fn.npz'],100,600),
      'C':(['r134_C.npz'],600,1000),'Bp':(['r138_Bp.npz'],600,1000),'R120':(['r117_part0.npz','r117_part1.npz','r120_part0.npz','r120_part1.npz'],100,600)}
    st=json.load(open('base_stats.json'))
    xg,wg=leggauss(6)
    out={}
    for name,(files,lo,hi) in SETS.items():
        rows=[]
        for x in st[name]:
            D=x['D']; q=abs(D); S=primefactors(q)
            tt=(xg+1)/2*(hi-lo)+lo; Ls=np.log(q*tt/(2*np.pi)); ww=wg*Ls; ww=ww/ww.sum()
            pr={}
            for s0 in S0:
                rr=np.linspace(2e-3,s0,300)
                val=sum(w*np.trapezoid(Rpair(rr,L,S),rr) for w,L in zip(ww,Ls))
                pr[s0]=val*kap[s0]
            rows.append(dict(D=D,n=x['n'],pred=pr,obs={0.05:x['n05'],0.1:x['n10'],0.15:x['n15'],0.2:x['n20'],0.3:x['n30']}))
        out[name]=rows
        line=[]
        for s0 in S0:
            e=sum(r['n']*r['pred'][s0] for r in rows); o=sum(r['obs'][s0] for r in rows); es=sum(r['n']*cdf(s0) for r in rows)
            line.append('s<%.2f obs %5d BK %7.1f sine %7.1f obs/BK %.3f+-%.3f BK/sine %.3f'%(s0,o,e,es,o/e,math.sqrt(o)/e,e/es))
        print(name); print('   '+'\n   '.join(line))
    json.dump({k:[dict(D=r['D'],n=r['n'],pred={str(a):b for a,b in r['pred'].items()},obs={str(a):b for a,b in r['obs'].items()}) for r in v] for k,v in out.items()},open('bk_small.json','w'))
    # pooled over the five test sets
    for s0 in S0:
        e=sum(r['n']*r['pred'][s0] for n in ['B','Dn','Fn','C','Bp'] for r in out[n]); o=sum(r['obs'][s0] for n in ['B','Dn','Fn','C','Bp'] for r in out[n])
        print('POOLED five sets s<%.2f obs %d BK-pred %.1f ratio %.3f +- %.3f'%(s0,o,e,o/e,math.sqrt(o)/e))
