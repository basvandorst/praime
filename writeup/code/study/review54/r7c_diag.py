import math, numpy as np, cypari2, mpmath, time
pari = cypari2.Pari(); pari.allocatemem(800 * 10 ** 6)
N = 13
for PD in (60, 120):
    pari.set_real_precision(PD)
    mf = pari('mfinit([13,2,Mod(4,13)],0)'); F = pari.mfeigenbasis(mf)[0]
    E = pari.mfembed(F, pari.mfcoefs(F, 600))
    A = [mpmath.mpc(mpmath.mpf(str(pari.real(E[n])).replace(' E', 'e')), mpmath.mpf(str(pari.imag(E[n])).replace(' E', 'e'))) for n in range(1, len(E))]
    Ld = pari.lfunmf(mf, F)
    for t in (20.0, 30.0, 40.0, 45.0):
        t0 = time.time()
        lp = complex(pari.lfun(Ld, pari('1 + %s*I' % t)))
        tp = time.time() - t0
        print('PD=%d t=%.0f PARI(no init) %s (%.1fs)' % (PD, t, np.round(lp, 8), tp))
# own AFE with split parameter A at very high precision, check A-independence
def afe(t, A_split, dps, eps):
    with mpmath.workdps(dps):
        s = mpmath.mpc(1, t); x0 = 2 * mpmath.pi / mpmath.sqrt(N)
        tot = mpmath.mpc(0)
        for n in range(1, len(A) + 1):
            an = A[n - 1]
            tot += an * mpmath.power(N, s / 2) * mpmath.power(2 * mpmath.pi * n, -s) * mpmath.gammainc(s, n * x0 * A_split)
            tot += eps * mpmath.conj(an) * mpmath.power(N, (2 - s) / 2) * mpmath.power(2 * mpmath.pi * n, -(2 - s)) * mpmath.gammainc(2 - s, n * x0 / A_split) * mpmath.power(A_split, 0)
        return tot / (mpmath.power(N, s / 2) * mpmath.power(2 * mpmath.pi, -s) * mpmath.gamma(s))
mpmath.mp.dps = 120
yy = mpmath.mpf(9) / 10 / mpmath.sqrt(N)
th1 = mpmath.fsum(A[k] * mpmath.exp(-2 * mpmath.pi * (k + 1) / (N * yy)) for k in range(len(A)))
th2 = mpmath.fsum(mpmath.conj(A[k]) * mpmath.exp(-2 * mpmath.pi * (k + 1) * yy) for k in range(len(A)))
eps = th1 / (N * yy * yy * th2)
print('eps', mpmath.nstr(eps, 15))
for t in (20.0, 45.0):
    for dps in (60, 110):
        t0 = time.time()
        v = afe(t, 1, dps, eps)
        print('own AFE t=%.0f dps=%d: %s (%.1fs)' % (t, dps, mpmath.nstr(v, 12), time.time() - t0))
