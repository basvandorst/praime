import os, numpy as np, math
from numpy.polynomial.legendre import leggauss
LZ=os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/data/riemann/lzeros/'
SETS={'B':('r134_B.npz',100,600),'Dn':('r136_Dn.npz',100,600),'Fn':('r138_Fn.npz',100,600),'C':('r134_C.npz',600,1000),'Bp':('r138_Bp.npz',600,1000)}
def unfold(t,q): return t/(2*math.pi)*np.log(q*t/(2*math.pi*math.e))
S=[];R=[]
for nm,(f,lo,hi) in SETS.items():
    z=np.load(LZ+f)
    for k in z.files:
        q=abs(int(k[1:])); t=z[k]; t=t[(t>=lo)&(t<=hi)]; s=np.diff(unfold(t,q)); S.append(s)
        R.append(np.minimum(s[1:],s[:-1])/np.maximum(s[1:],s[:-1]))
s=np.concatenate(S); r=np.concatenate(R)
# sine spacing density via E(0;s) second derivative
def E0(x,m=40):
    xg,w=leggauss(m); xg=(xg+1)*x/2; w=w*x/2
    K=np.sinc(xg[:,None]-xg[None,:]); A=np.sqrt(w)[:,None]*K*np.sqrt(w)[None,:]; return np.linalg.det(np.eye(m)-A)
grid=np.linspace(0,4,801); h=grid[1]-grid[0]
E=np.array([E0(x) if x>0 else 1.0 for x in grid]); cdf=1+np.gradient(E,h)   # F(s)=1+E'(s)
edges=[0,0.3,0.6,1.0,1.4,2.0,4.0]
print('pooled five sets: %d spacings, mean %.4f, Var %.5f (sine 0.17999)'%(len(s),s.mean(),s.var()))
tot=0
for a,b in zip(edges[:-1],edges[1:]):
    m=(s>=a)&(s<b); cd=np.sum((s[m]-1)**2)/len(s)
    # sine: int (x-1)^2 p(x) over [a,b] = via cdf on grid
    gi=(grid>=a)&(grid<=b); x=grid[gi]; p=np.gradient(cdf,h)[gi]; cs=np.trapezoid((x-1)**2*p,x)
    tot+=cd-cs
    print('  s in [%.1f,%.1f): data %.5f  sine %.5f  diff %+.5f'%(a,b,cd,cs,cd-cs))
print('  total diff %+.5f (Var data - sine, mean-1 convention)'%tot)
