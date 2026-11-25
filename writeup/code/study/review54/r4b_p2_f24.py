import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+"/")
import math, numpy as np, cypari2, sympy
from r1_sums import theta
pari = cypari2.Pari(); pari.allocatemem(200 * 10 ** 6); pari.set_real_precision(40)
T2 = float(sys.argv[1])
pari('mf = mfinit([1,24],1)'); pari('F = mfeigenbasis(mf)[1]'); pari('v = mfcoefs(F, 500)'); pari('pol = mffields(mf)[1]'); pari('vt = variable(pol)'); pari('rts = polroots(pol)')
for i in (1, 2):
    pari('zz = real(rts[%d])' % i)
    co = [float(pari('subst(lift(v[%d]), vt, zz)' % (n + 1))) for n in range(0, 401)]
    lam = lambda p: co[p] / p ** 11.5
    mu = [11.5, 12.5]
    t = np.linspace(1e-6, T2, 400001); dt = t[1] - t[0]
    th2 = 2 * theta(t, 1.0, mu)
    P2 = 0j
    for p in sympy.primerange(2, 300):
        f = np.exp(1j * (th2 - t * 2 * np.log(p)))
        I = (f.sum() - 0.5 * (f[0] + f[-1])) * dt
        P2 += -(1 / (2 * np.pi)) * (lam(p) ** 2 - 2) * np.log(p) / p * I
    print('embedding %d: lambda_2 = %.4f ; prime-square term P2(T=%.0f) = %.2f%+.2fi' % (i, lam(2), T2, P2.real, P2.imag))
