"""Zeros of L(s, chi_D) up to height T (PARI lfunzeros), saved to writeup/data/riemann/lzeros/D<D>.txt.
usage: getzeros.py D T [D T ...]"""
import sys, os, time
import cypari2
pari = cypari2.Pari(); pari.allocatemem(900 * 10 ** 6)
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'data', 'riemann', 'lzeros')
a = sys.argv[1:]
for D, T in zip(a[::2], a[1::2]):
    t0 = time.time()
    z = [float(x) for x in pari.lfunzeros(pari.lfuncreate(int(D)), float(T))]
    open(os.path.join(OUT, 'D%s.txt' % D), 'w').write('\n'.join('%.12f' % x for x in z if x > 0) + '\n')
    print(D, len(z), 'zeros up to', T, round(time.time() - t0), 's', flush=True)
