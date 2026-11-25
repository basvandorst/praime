# Independent recomputation: per-L-function statistics, completeness, small gaps.
import os, numpy as np, math, json
from scipy.special import loggamma
LZ=os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/data/riemann/lzeros/'
SETS={'B':(['r134_B.npz'],100,600),'Dn':(['r136_Dn.npz'],100,600),'Fn':(['r138_Fn.npz'],100,600),
      'C':(['r134_C.npz'],600,1000),'Bp':(['r138_Bp.npz'],600,1000),
      'R120':(['r117_part0.npz','r117_part1.npz','r120_part0.npz','r120_part1.npz'],100,600)}
def theta(t,q,a):
    t=np.asarray(t,float); return np.imag(loggamma(0.25+a/2+0.5j*t))+t/2*math.log(q/math.pi)
def unfold(t,q): return t/(2*math.pi)*np.log(q*t/(2*math.pi*math.e))
out={}
for name,(files,lo,hi) in SETS.items():
    Z={}
    for f in files:
        z=np.load(LZ+f); Z.update({int(k[1:]):np.array(z[k]) for k in z.files})
    Ds=sorted(Z.keys(),key=lambda d:(abs(d),d))
    rec=[]
    for D in Ds:
        g=Z[D]; q=abs(D); a=0 if D>0 else 1
        assert np.all(np.diff(g)>0)
        b=g[(g>=lo)&(g<=hi)]
        u=unfold(b,q); s=np.diff(u)
        r=np.minimum(s[1:],s[:-1])/np.maximum(s[1:],s[:-1])
        # count check on a fine grid: N(t) - theta/pi
        tt=np.linspace(lo,hi,2001)
        dev=np.searchsorted(g,tt)-theta(tt,q,a)/math.pi
        # tiny t-gaps anywhere in file
        dt=np.diff(g)
        L=np.log(q*np.array([lo,hi])/(2*math.pi))
        rec.append(dict(D=D,n=len(s),var=float(np.mean((s-s.mean())**2)),r=float(r.mean()),ms=float(s.mean()),
            m2=float(np.mean(s**2)),
            dev_max=float(np.abs(dev).max()),dev_mean=float(dev.mean()),
            oct=[float(x) for x in dev[:2000].reshape(8,250).mean(1)],
            maxgap=float(s.max()),mingap=float(s.min()),mintgap_file=float(dt.min()),
            n05=int((s<0.05).sum()),n10=int((s<0.1).sum()),n15=int((s<0.15).sum()),n20=int((s<0.2).sum()),n30=int((s<0.3).sum()),
            Lmean=float(np.mean(np.log(q*b/(2*math.pi)))),
            ))
    out[name]=rec
    v=np.array([x['var'] for x in rec]); rr=np.array([x['r'] for x in rec])
    print(name,len(rec),'meanVar %.5f mean r %.5f'%(v.mean(),rr.mean()),'n spacings mean %.1f'%np.mean([x['n'] for x in rec]),
          'max|dev| %.2f'%max(x['dev_max'] for x in rec),'max|oct| %.2f'%max(max(abs(o) for o in x['oct']) for x in rec),
          'mean dev %.3f'%np.mean([x['dev_mean'] for x in rec]),
          'maxgap %.2f'%max(x['maxgap'] for x in rec),'mingap %.4f'%min(x['mingap'] for x in rec),
          'min t-gap in file %.2e'%min(x['mintgap_file'] for x in rec))
json.dump(out,open('base_stats.json','w'))
