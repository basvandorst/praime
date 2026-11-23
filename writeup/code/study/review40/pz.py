import cypari2, numpy as np
pari = cypari2.Pari(); pari.allocatemem(400*10**6)
pari.set_real_precision(30)
for D, iv in [(-47, [(480.5, 484.5), (523.5, 527.5)]), (-59, [(341.5, 347.0), (501.5, 505.5)]), (-31,[(704.5,708.5)]), (29,[(447,451)])]:
    L = pari('lfuncreate(%d)' % D)
    g = np.sort(np.load('/home/bas/praime/writeup/data/riemann/lzeros/r128/D%d.npy' % D))
    for a, b in iv:
        z = [float(v) for v in pari.lfunzeros(L, [a, b])]
        have = g[(g > a) & (g < b)]
        print(D, (a, b), 'pari', np.round(z, 4), 'file', np.round(have, 4))
