"""Candidates (a), (b), (d): the finite-height corrections that the Bogomolny-Keating (BK) law leaves out, computed.
 (a1) the ratios-conjecture two-point function for L(s, chi), primitive chi mod k, written in its own (Conrey-Snaith)
      form with (L'/L)'(1+i eps, chi_0), B_k, the swap term (k/phi(k))^2 L(1+z,chi_0) L(1-z,chi_0) A_k(z), evaluated by
      direct Euler products (no splines), against bklib's C(r) (BK form);
 (a2) the exact off-diagonal phase 2[theta(t+eps/2) - theta(t-eps/2)] (true gamma factor, parity a) against eps L(t);
      the parity difference; the change of the p = 2 off-diagonal number-variance term and of gamma it causes;
 (b)  the density drift across a pair: d(t1) d(t2)/d(t)^2 - 1, and the pair-midpoint phase term eps^3/(24 t^2);
 (d)  the test's unfolding u = (t/2pi) ln(k t/(2 pi e)) against the exact smooth count theta(t)/pi: the non-constant part,
      the phase error it induces over 8 spacings, against the rms phase jitter a 7 % loss needs.
usage: theory.py -> theory_out.txt"""
import os, sys, math
import numpy as np
import mpmath as mp
from scipy.special import loggamma
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'bk'))
from common import *
import gridcache
mp.mp.dps = 25
G = gridcache.get()
PR = bklib.primes_upto(10 ** 6).astype(float); LP = np.log(PR)

def ratios_C(r, L, k):
    """R(r) - 1 + sinc^2(r) from the ratios-conjecture form, eps = 2 pi r/L, direct Euler products over p < 1e6."""
    e = 2 * math.pi * r / L; s = mp.mpc(1, e); bad = prime_factors(k) if k > 1 else []
    z = mp.zeta(s); z1 = mp.zeta(s, derivative=1); z2 = mp.zeta(s, derivative=2)
    dzz = z2 / z - (z1 / z) ** 2                                   # (zeta'/zeta)'(1+i eps)
    for p in bad:                                                  # (L'/L)'(s, chi_0) = (zeta'/zeta)' - sum (ln p)^2 p^-s/(1-p^-s)^2
        x = mp.mpf(p) ** (-s); dzz -= mp.log(p) ** 2 * x / (1 - x) ** 2
    good = np.array([p not in bad for p in PR])
    xs = np.exp(-(1 + 1j * e) * LP)                                # p^{-1-i eps}
    Bk = complex(((LP ** 2) * xs ** 2 / (1 - xs) ** 2)[good].sum())
    diag = (1 / (2 * math.pi ** 2)) * (complex(dzz).real - Bk.real)
    # swap term: z = -i eps pairs with e^{+i eps L};  A_p(z) = (1 - p^{-1-z})(1 - 2/p + p^{-1-z})/(1-1/p)^2, x = p^{-z} = p^{i eps}
    xp = np.exp(1j * e * LP)
    lnA = np.log((1 - xp / PR) * (1 - 2 / PR + xp / PR) / (1 - 1 / PR) ** 2)[good].sum()
    zz = complex(mp.zeta(mp.mpc(1, e)) * mp.zeta(mp.mpc(1, -e)))
    badf = 1.0 + 0j
    for p in bad:
        badf *= (p / (p - 1)) ** 2 * (1 - p ** (-1 + 1j * e)) * (1 - p ** (-1 - 1j * e))
    off = (1 / (2 * math.pi ** 2)) * (np.exp(1j * 2 * math.pi * r) * zz * badf * np.exp(lnA)).real
    total = diag + off                                             # = R2 - dbar^2 in raw units
    return total / (L / (2 * math.pi)) ** 2 + np.sinc(r) ** 2      # -> C(r) of bklib (R - 1 = -sinc^2 + C)

print('(a1) ratios-conjecture form (direct Euler products, mpmath zeta) against BK (bklib C(r), splines):')
mx = 0
for k in (1, 5, 8, 12, 105, 440, 664):
    B = tuple(prime_factors(k)) if k > 1 else ()
    for r in (0.25, 1.3, 3.7, 6.2):
        L = 10.2; a = ratios_C(r, L, k); b = float(bklib.C_of_r(G, L, B, r=np.array([r]))[0]); mx = max(mx, abs(a - b))
        if r in (1.3, 6.2): print('    k %4d r %.2f  ratios %+.8f  BK %+.8f  diff %.1e' % (k, r, a, b, a - b))
print('    largest |difference| over 7 conductors x 4 distances: %.1e' % mx)

def theta(t, k, a):
    t = np.asarray(t, float); return (t / 2) * np.log(k / np.pi) + np.imag(loggamma((0.5 + a + 1j * t) / 2))

print('\n(a2) exact off-diagonal phase 2[theta(t+eps/2) - theta(t-eps/2)] - eps ln(kt/2pi), k = 550 (radians):')
for t in (100.0, 300.0, 600.0):
    row = []
    for a in (0, 1):
        for e in (1.0, 3.0, 5.5):
            ph = 2 * (theta(t + e / 2, 550, a) - theta(t - e / 2, 550, a)) - e * math.log(550 * t / (2 * math.pi))
            row.append(float(ph))
    print('    t %5.0f  a=0: eps 1, 3, 5.5 -> %+.2e %+.2e %+.2e   a=1: %+.2e %+.2e %+.2e' % ((t,) + tuple(row)))
e = 5.5; t = 100.0
pa = [2 * (theta(t + e / 2, 550, a) - theta(t - e / 2, 550, a)) for a in (0, 1)]
print('    parity difference of the phase at t = 100, eps = 5.5: %.1e rad' % abs(pa[0] - pa[1]))
print('    predicted from Stirling: -eps/(24 t^2) - eps^3/(24 t^2) * ... ; at t=100, eps=5.5: eps^3/(24 t^2) = %.1e' % (e ** 3 / (24 * t * t)))
# effect on the p = 2 off-diagonal term and on gamma, at the worst height t = 100 (L = ln(550*100/2pi))
L = math.log(550 * 100 / (2 * math.pi)); r = bklib.R_NODES; ep = 2 * math.pi * r / L
dphi = 2 * (theta(t + ep / 2, 550, 0) - theta(t - ep / 2, 550, 0)) - ep * L
base = G('Z2', ep) * G('Phi', ep) * bklib.f_p(2, ep)
f0 = (2 / L ** 2) * (np.exp(1j * 2 * math.pi * r) * base).real
f1 = (2 / L ** 2) * (np.exp(1j * (2 * math.pi * r + dphi)) * base).real
s0 = bklib.cum_sigma(f0)[IDX]; s1 = bklib.cum_sigma(f1)[IDX]
print('    p = 2 off-diagonal Sigma2 term at t = 100 with the exact phase: max change %.1e (term itself up to %.1e); projection ratio %.6f' % (
    np.abs(s1 - s0).max(), np.abs(s0).max(), float((s1 * s0).sum() / (s0 * s0).sum())))

print('\n(b) density drift and pair midpoint (k = 550):')
for t in (100.0, 600.0):
    Lt = math.log(550 * t / (2 * math.pi)); e = 2 * math.pi * 8 / Lt
    d = lambda x: math.log(550 * x / (2 * math.pi)) / (2 * math.pi)
    print('    t %5.0f, r = 8 spacings (eps %.2f): d(t-eps/2) d(t+eps/2)/d(t)^2 - 1 = %+.1e;  eps^3/(24 t^2) = %.1e rad' % (
        t, e, d(t - e / 2) * d(t + e / 2) / d(t) ** 2 - 1, e ** 3 / (24 * t * t)))

print('\n(d) test unfolding u = (t/2pi) ln(kt/(2pi e)) against theta(t)/pi (exact smooth count), k = 550:')
ts = np.linspace(100, 600, 50001)
for a in (0, 1):
    du = theta(ts, 550, a) / np.pi - ts / (2 * np.pi) * np.log(550 * ts / (2 * np.pi * np.e))
    print('    a = %d: theta/pi - u = %+.6f at t=100, %+.6f at t=600 (constant %+.4f = %s1/8; the rest is 1/(48 pi t) = %.2e at 100)' % (
        a, du[0], du[-1], du[-1], '-' if a == 0 else '+', 1 / (48 * math.pi * 100)))
    # largest change of (theta/pi - u) across 8 spacings (raw width ~ 8 * 2pi/L)
    w = 8 * 2 * math.pi / math.log(550 * 100 / (2 * math.pi)); j = int(w / (ts[1] - ts[0]))
    print('           largest change across 8 spacings: %.1e spacings -> phase error 2pi x that = %.1e rad' % (np.abs(du[j:] - du[:-j]).max(), 2 * math.pi * np.abs(du[j:] - du[:-j]).max()))
print('    a 7 %% loss of a period-1 term by Gaussian phase jitter needs rms phase sqrt(2 x 0.07) = %.2f rad, i.e. %.3f spacings' % (
    math.sqrt(2 * 0.07), math.sqrt(2 * 0.07) / (2 * math.pi)))
