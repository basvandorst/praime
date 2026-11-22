# asymptotic (infinite ln X range) RMS of the off-diagonal part, from the zeros: groups ordered pairs by Delta
from rz import *
import json
out={}
for N in (2419,10000):
    g=load(N); G=np.concatenate([g,-g])
    for d in (0.1,0.03,0.01,0.003,0.001):
        c=coef(G,d); gen=0.0; selfc=0.0
        for i in range(0,len(G),500):
            D=G[i:i+500,None]-G[None,:]
            a2=np.abs(c[i:i+500,None])**2*np.abs(c[None,:])**2/(4+D**2)
            sc=np.isclose(G[i:i+500,None],-G[None,:])
            m=D>0
            gen+=a2[m & ~sc].sum(); selfc+=a2[m & sc].sum()
        ms=(2/d)**2*2*(2*gen+selfc)
        diag=2*np.sum(np.abs(c[:N])**2)/d
        out['%d_%g'%(N,d)]=dict(rms=float(np.sqrt(ms)),diag_trunc=float(diag))
        print(N,d,'rms off-diag',np.sqrt(ms),'diag(trunc)',diag,flush=True)
json.dump(out,open('swing.json','w'),indent=1)
