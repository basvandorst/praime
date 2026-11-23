import cypari2, numpy as np, glob, sys
pari = cypari2.Pari(); pari.allocatemem(800*10**6)
pari.set_real_precision(25)
out = {}
for f in sorted(glob.glob('/home/bas/praime/writeup/data/riemann/lzeros/r128/D*.npy')):
    D = int(f.split('/D')[-1][:-4])
    g = np.sort(np.load(f))
    a, b = g[0]-0.01, g[-1]+0.01
    L = pari('lfuninit(lfuncreate(%d), [%f, %f])' % (D, 0, b+1))
    z = np.array([float(v) for v in pari.lfunzeros(L, [a, b])])
    np.save('pari_D%d.npy' % D, z)
    missing = [v for v in z if np.min(np.abs(g - v)) > 1e-5]
    extra = [v for v in g if np.min(np.abs(z - v)) > 1e-5]
    print(D, 'file', len(g), 'pari', len(z), 'missing', np.round(missing, 4), 'extra', np.round(extra, 4), flush=True)
