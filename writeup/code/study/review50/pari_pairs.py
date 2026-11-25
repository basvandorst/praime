import os, cypari2, numpy as np
pari=cypari2.Pari(); pari.allocatemem(200*10**6)
LZ=os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/data/riemann/lzeros/'
cases=[('r138_Bp.npz',-644,908.0174),('r134_B.npz',577,399.8012),('r138_Fn.npz',1265,314.3228),('r138_Bp.npz',677,652.8529)]
for f,D,t0 in cases:
    g=np.load(LZ+f)['D%d'%D]; i=np.argmin(abs(g-t0)); t1,t2=g[i],g[i+1]
    L=pari.lfuncreate(D)
    d=t2-t1
    pts=[t1-3*d,(t1+t2)/2,t2+3*d]
    vals=[float(pari.lfunhardy(L,pari(repr(float(x))))) for x in pts]
    # also check no extra zero: sign at midpoints with neighbours
    print(D,'pair',t1,t2,'gap_t %.5f'%d,'Z at t1-3d, mid, t2+3d:',['%.3e'%v for v in vals])
