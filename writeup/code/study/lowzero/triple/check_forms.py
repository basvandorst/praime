"""Checks of the three-point function implementation (tplib.T3_parts):
 (1) the GUE limit of the six terms equals 2 S12 S23 S31;
 (2) R3 = 0 at a coincidence: T3(0, t -> 0, s) -> -2 T2(s), T2 = R2 - 1 of Bogomolny-Keating (bklib.C_of_r);
 (3) the ratios-conjecture form (Conrey-Snaith 2008, adapted) against Bogomolny-Keating 2013's independent form
     (their eqs. for R3^diag and R3^osc, random-matrix-universality route), both with p | k deleted, evaluated here with
     mpmath zeta, direct Euler products over p < 1e5 and no splines.
usage: check_forms.py -> check_forms_out.txt"""
import sys, math, os
import numpy as np, mpmath as mp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import tplib, bklib
mp.mp.dps = 20
PS = tplib.PrimeSums(10 ** 5)
P = tplib.primes_upto(10 ** 5).astype(float); LP = np.log(P)

def Lfun(s, S):          # L(s, chi_0) = zeta(s) prod_{p|k} (1 - p^-s)
    v = mp.zeta(s)
    for p in S: v *= (1 - mp.mpf(p) ** (-s))
    return complex(v)
def LdL(s, S):           # L'/L(s, chi_0)
    v = mp.zeta(s, derivative=1) / mp.zeta(s)
    for p in S: v += mp.log(p) / (mp.mpf(p) ** s - 1)
    return complex(v)

def T3_bk(e, L, S):
    """Bogomolny-Keating 2013 form, raw ordinates e = (e1, e2, e3); returns T3 (unfolded, connected)."""
    good = np.array([p not in S for p in P]); Pg = P[good]; lg = LP[good]
    kphi = np.prod([p / (p - 1.0) for p in S]) if S else 1.0; dbar = L / (2 * math.pi)
    def diag():
        tot = 0
        for (i, j, k) in ((0, 1, 2), (1, 0, 2), (2, 1, 0)):
            eij = e[i] - e[j]; eik = e[i] - e[k]
            tot += np.sum(lg ** 3 / ((Pg ** (1 - 1j * eij) - 1) * (Pg ** (1 - 1j * eik) - 1)))
        return -(2 * tot.real) / (2 * math.pi) ** 3
    def osc_term(i, j, k):
        e12 = e[i] - e[j]; e23 = e[j] - e[k]; e13 = e[i] - e[k]; e31 = -e13; e32 = -e23
        z2 = kphi ** 2 * abs(Lfun(mp.mpc(1, e12), S)) ** 2
        phi = np.prod(1 - (1 - Pg ** (1j * e12)) ** 2 / (Pg - 1) ** 2)
        # d/de3 ln|L(1+i e32)|^2 - ln|L(1+i e31)|^2 = 2 Re[i L'/L(1+i e32)] - 2 Re[i L'/L(1+i e31)]
        dln = 2 * (1j * LdL(mp.mpc(1, e32), S)).real - 2 * (1j * LdL(mp.mpc(1, e31), S)).real
        q = Pg; lq = lg
        ssum = np.sum(lq * ((q ** (1j * e12) - 1) / ((q ** (1 + 1j * e23) - 1) * (q ** (1 + 1j * e13) - 1))
                            + (q ** (1j * e12) - 1) / ((q ** (1 + 1j * e31) - 1) * (q ** (1 + 1j * e32) - 1))
                            + (1 - q ** (1j * e12)) ** 2 / (q - 2 + q ** (1j * e12)) * (1 / (q ** (1 - 1j * e31) - 1) + 1 / (q ** (1 - 1j * e23) - 1))))
        r = -np.exp(2j * math.pi * dbar * e12) / (2j * math.pi) ** 3 * z2 * phi * (dln - 1j * ssum)
        return 2 * r.real
    osc = osc_term(0, 1, 2) + osc_term(2, 0, 1) + osc_term(1, 2, 0)
    return (diag() + osc) / dbar ** 3

rng = np.random.default_rng(3)
print('(1) GUE limit: max |six GUE terms - 2 S S S| over 2000 random triples in [0, 8]^3:', end=' ')
u = np.zeros(2000); v = rng.uniform(0.01, 8, 2000); w = rng.uniform(0.01, 8, 2000)
d = tplib.T3_parts(u, v, w, 9.4, (), PS=PS)
print('%.2e' % np.abs(d['gue_terms'] - tplib.T3_gue(u, v, w)).max())

print('(2) coincidence: T3(0, t, s) at t = 1e-3, 1e-4 against -2 T2(s)')
for L, S in ((9.4, ()), (9.4, (2,)), (10.2, (3, 5)), (10.7, (2, 7, 13))):
    s = np.array([0.37, 1.0, 2.5, 4.4, 7.1])
    T2 = -np.sinc(s) ** 2 + bklib.C_of_r(tplib.G, L, S, r=s)
    for t in (1e-3, 1e-4):
        dd = tplib.T3_parts(np.zeros(5), np.full(5, t), s, L, S, PS=PS); tot = tplib.T3_gue(0 * s, 0 * s + t, s) + dd['Q'] + dd['P']
        print('   L=%.1f S=%-10s t=%.0e: T3 %s | -2T2 %s' % (L, S, t, ' '.join('%+.6f' % x for x in tot), ' '.join('%+.6f' % x for x in -2 * T2)))

print('(3) ratios form (ours) against Bogomolny-Keating 2013 form, T3 at random triples (unfolded)')
worst = 0
for L, S in ((9.4, ()), (9.4, (2,)), (10.2, (3, 5)), (10.7, (2, 7, 13))):
    for _ in range(4):
        r = np.sort(rng.uniform(0, 6, 3)); r = r - r[0]
        dd = tplib.T3_parts(np.array([r[0]]), np.array([r[1]]), np.array([r[2]]), L, S, PS=PS)
        ours = float(tplib.T3_gue(r[0], r[1], r[2]) + dd['Q'][0] + dd['P'][0])
        bk = T3_bk(2 * math.pi / L * r, L, S)
        worst = max(worst, abs(ours - bk))
        print('   L=%.1f S=%-10s r=(0, %.3f, %.3f): ours %+.9f  BK2013 %+.9f  GUE %+.9f  diff %.1e' % (L, S, r[1], r[2], ours, bk, tplib.T3_gue(*r), ours - bk))
print('   largest |ours - BK2013| = %.1e' % worst)
