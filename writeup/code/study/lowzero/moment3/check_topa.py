"""Compare the recipe's pair polynomial (tmlib, two-variable constant term) with Topacogullari's theorem
(arXiv 1909.11517, eq. (5576), q1 != q2): P(log t) = Res_{z=0} [ 2 (q1 q2)^z psi_z(q1)^2 psi_z(q2)^2 /
 ((psi_0(q2) q1^z psi_2z(q1) + psi_0(q1) q2^z psi_2z(q2)) psi_{1+2z}(q1 q2)) * zeta(1+z)^4 L(1+z)^2 / ((2pi)^z zeta(1+2z) zeta(2+2z)) t^z ]."""
import os, sys, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, mpmath as mp
import tmlib as T
pari = T.pari
mp.mp.dps = 30

def psi(z, q):
    out = mp.mpf(1)
    for p in T.prime_factors(q):
        out *= (1 - mp.mpf(p) ** (-1 - z))
    return out

def Limp(z, D1, D2):
    Ds = int(pari.coredisc(D1 * D2)); L = pari.lfuncreate(Ds)
    v = pari.lfun(L, pari('1+(%s)+I*(%s)' % (mp.nstr(z.real, 25), mp.nstr(z.imag, 25))))
    v = mp.mpc(float(pari.real(v)), float(pari.imag(v)))
    for p in T.prime_factors(D1 * D2):
        c = int(pari.kronecker(Ds, p))
        if c: v *= (1 - c * mp.mpf(p) ** (-1 - z))
    return v

def topa(D1, D2, ts, N=48, r=0.15):
    q1, q2 = abs(D1), abs(D2)
    out = np.zeros(len(ts))
    for k in range(N):
        z = r * mp.exp(2j * mp.pi * (k + 0.5) / N)
        F = (2 * (q1 * q2) ** z * psi(z, q1) ** 2 * psi(z, q2) ** 2
             / ((psi(0, q2) * mp.mpf(q1) ** z * psi(2 * z, q1) + psi(0, q1) * mp.mpf(q2) ** z * psi(2 * z, q2)) * psi(1 + 2 * z, q1 * q2))
             * mp.zeta(1 + z) ** 4 * Limp(z, D1, D2) ** 2 / ((2 * mp.pi) ** z * mp.zeta(1 + 2 * z) * mp.zeta(2 + 2 * z)))
        for i, t in enumerate(ts):
            out[i] += float(mp.re(F * mp.mpf(t) ** z * z)) / N      # (1/2 pi i) oint f dz = mean of f(z) z over the circle
    return out

ts = [100.0, 300.0, 1000.0]
for D1, D2 in [(1, 1001), (1, -1003), (1, 1068), (1001, -1003), (1005, 1068), (-1048, 1068), (1001, 1253), (-1195, 1005), (13, -1267)]:
    tr = T.Triple((D1, D2))
    mine = tr.P(np.array(ts)); th = topa(D1, D2, ts)
    print('%6d %6d  recipe %s  Topacogullari %s  rel.diff %s  lead h0 %.6f' % (D1, D2, np.round(mine, 6), np.round(th, 6), np.array2string((mine - th) / th, precision=2), tr.hm[0]))
