import os, numpy as np
LZ=os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/data/riemann/lzeros/'
for f in ['r134_B.npz','r136_Dn.npz','r138_Fn.npz','r134_C.npz','r138_Bp.npz','r120_part0.npz','r120_part1.npz','r117_part0.npz','r117_part1.npz']:
    z=np.load(LZ+f)
    ks=list(z.keys())
    mins=[z[k].min() for k in ks]; maxs=[z[k].max() for k in ks]
    Ds=[int(k[1:]) for k in ks]
    srt=[np.all(np.diff(z[k])>0) for k in ks]
    print(f,len(ks),'min',min(mins),max(mins),'max',min(maxs),max(maxs),'|D| range',min(abs(d) for d in Ds),max(abs(d) for d in Ds),'sorted',all(srt), 'other keys',[k for k in ks if not k.startswith('D')][:5])
