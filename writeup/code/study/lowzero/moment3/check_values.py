"""check_values.py: independent check of stored |L(1/2+it, chi_D)| values (moment_values.py) at random grid points,
with L(s, chi_D) = q^{-s} sum_{a=1}^{q} chi_D(a) zeta(s, a/q) (mpmath Hurwitz zeta; kronecker symbol by PARI).
usage: check_values.py VALDIR [D ...] [--npts 2] [--seed 7]
Prints the largest relative difference; a stored file passes if every checked point agrees to 1e-10.
"""
import sys, os, math
import numpy as np, mpmath as mp
import cypari2

pari = cypari2.Pari()


def L_hurwitz(D, s):
    q = abs(D); tot = mp.mpc(0)
    for a in range(1, q + 1):
        c = int(pari.kronecker(D, a))
        if c:
            tot += c * mp.zeta(s, mp.mpf(a) / q)
    return tot * mp.mpf(q) ** (-s)


def main(argv):
    argv = list(argv)
    npts = int(argv[argv.index('--npts') + 1]) if '--npts' in argv else 2
    seed = int(argv[argv.index('--seed') + 1]) if '--seed' in argv else 7
    for flag in ('--npts', '--seed'):
        if flag in argv:
            i = argv.index(flag); del argv[i:i + 2]
    vd = argv[1]
    Ds = [int(x) for x in argv[2:]] or sorted(int(f[5:-4]) for f in os.listdir(vd) if f.startswith('absL_') and f.endswith('.npy'))
    mp.mp.dps = 20
    rng = np.random.default_rng(seed); worst = 0.0
    for D in Ds:
        v = np.load(os.path.join(vd, 'absL_%d.npy' % D))
        for k in rng.integers(0, len(v), npts):
            t = 100.0 + 0.05 * int(k)
            ref = float(abs(L_hurwitz(D, mp.mpc(0.5, t)))) if D != 1 else float(abs(mp.zeta(mp.mpc(0.5, t))))
            rel = abs(v[k] - ref) / max(ref, 1e-6); worst = max(worst, rel)
            print('D = %6d  t = %8.2f  stored %.12f  Hurwitz %.12f  rel %.1e %s' % (D, t, v[k], ref, rel, '' if rel < 1e-10 else 'FAIL'), flush=True)
    print('largest relative difference %.2e: %s' % (worst, 'PASS' if worst < 1e-10 else 'FAIL'))


if __name__ == '__main__':
    main(sys.argv)
