import os, numpy as np, math
from scipy.special import loggamma
LZ=os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/data/riemann/lzeros/'
SETS={'B':(['r134_B.npz'],100,600),'Dn':(['r136_Dn.npz'],100,600),'Fn':(['r138_Fn.npz'],100,600),
      'C':(['r134_C.npz'],600,1000),'Bp':(['r138_Bp.npz'],600,1000)}
def theta(t,q,a):
    t=np.asarray(t,float); return np.imag(loggamma(0.25+a/2+0.5j*t))+t/2*math.log(q/math.pi)
def unfold(t,q): return t/(2*math.pi)*np.log(q*t/(2*math.pi*math.e))
for name,(files,lo,hi) in SETS.items():
    Z={}
    for f in files:
        z=np.load(LZ+f); Z.update({int(k[1:]):np.array(z[k]) for k in z.files})
    allg=[]
    for D,g in Z.items():
        q=abs(D); b=g[(g>=lo)&(g<=hi)]; s=np.diff(unfold(b,q))
        for i in np.argsort(s)[:3]:
            allg.append((s[i],D,b[i],b[i+1]))
    allg.sort()
    print(name, ' '.join('(%.4f D=%d t=%.4f)'%(x[0],x[1],x[2]) for x in allg[:8]))
# local count around the 3 smallest of Bp
