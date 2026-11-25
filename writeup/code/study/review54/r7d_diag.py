import math, numpy as np, cypari2, mpmath, time
pari = cypari2.Pari(); pari.allocatemem(800 * 10 ** 6)
N = 13; PD = 100
pari.set_real_precision(PD); mpmath.mp.dps = PD
mf = pari('mfinit([13,2,Mod(4,13)],0)'); F = pari.mfeigenbasis(mf)[0]
co = pari.mfcoefs(F, 400)
print('exact a(1..6):', [str(co[n]) for n in range(1, 7)])
print('field:', pari.mffields(mf))
E = pari.mfembed(F, co)
A = [mpmath.mpc(mpmath.mpf(str(pari.real(E[n])).replace(' E', 'e')), mpmath.mpf(str(pari.imag(E[n])).replace(' E', 'e'))) for n in range(1, len(E))]
print('a(1..6) embedded:', [mpmath.nstr(x, 8) for x in A[:6]])
print('multiplicativity a(6) - a(2)a(3) =', mpmath.nstr(A[5] - A[1] * A[2], 5), ' a(4) - a(2)^2 + 2chi(2):', mpmath.nstr(A[3] - A[1] ** 2 + 2 * (A[1] ** 2 - A[3]) / 2, 5))
def fth(y, conj=False):
    return mpmath.fsum((mpmath.conj(A[k]) if conj else A[k]) * mpmath.exp(-2 * mpmath.pi * (k + 1) * y) for k in range(len(A)))
for yfac in ('0.6', '0.9', '1', '1.5'):
    y = mpmath.mpf(yfac) / mpmath.sqrt(N)
    eps = fth(1 / (N * y)) / (N * y * y * fth(y, True))
    print('y = %s/sqrt N: eps = %s' % (yfac, mpmath.nstr(eps, 40)))
