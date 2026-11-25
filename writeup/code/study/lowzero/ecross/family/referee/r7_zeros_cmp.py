import os, cypari2, numpy as np
pari = cypari2.Pari(); pari.allocatemem(300*10**6); pari.set_real_precision(38)
g = np.array([float(x) for x in pari('lfunzeros(lfunmf(mfinit([1,12],1),mfDelta()),420)')])
z = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'delta_zeros.npy'))
print(len(g), len(z), np.max(np.abs(g - z)))
