"""ps_check_lib.py: independent checks of ps_lib (sieve, enumerations, eigenvalues) against sympy and PARI."""
import sys, os, math
import numpy as np, sympy, cypari2
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ps_lib as L
pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)
fails = 0
def check(name, cond):
    global fails
    print(('PASS ' if cond else 'FAIL ') + name, flush=True); fails += (not cond)

# 1. sieve against sympy on three windows
for (Lo, U) in ((0, 200000), (10 ** 8, 10 ** 8 + 400000), (999000000, 999600000)):
    ok = L.sieve_odd(Lo, U); mine = set((Lo + 1 + 2 * np.nonzero(ok)[0]).tolist())
    ref = set(p for p in sympy.primerange(Lo + 1, U) if p % 2 == 1)
    check('sieve (%d, %d): %d primes' % (Lo, U, len(ref)), mine == ref)

# 2. all series on (0, 2e5) and (1e8, 1e8 + 2e5): lambda(p) against PARI coefficients / root counts
X = 200000
for (Lo, U) in ((0, X), (10 ** 8, 10 ** 8 + X)):
    ok = L.sieve_odd(Lo, U)
    pr = [p for p in sympy.primerange(max(Lo, 2), U)]
    # A: lambda(p) = a(p) / p^((k-1)/2) for the CM form; a(p) = 2 Re(pi^(k-1)), pi = a + bi any prime of norm p (exact integers)
    p, lams, ex = L.lam_series('A', Lo, U, ok, ['A4', 'A8', 'A12', 'A16', 'A20'])
    dct = {int(q): i for i, q in enumerate(p)}
    split = [q for q in pr if q % 4 == 1]
    check('A: every p = 1 mod 4 found once (%d)' % len(split), sorted(dct) == split and len(dct) == len(p))
    err = 0.0
    for q in split[:: max(1, len(split) // 300)]:
        a, b = [int(t) for t in pari.qfbsolve(pari.Qfb(1, 0, 1), q)]
        z = complex(a, b)
        for m in (4, 8, 12, 16, 20):
            ex2 = 2 * (complex(a, b) ** m).real / q ** (m / 2)   # (pi/|pi|)^m + conj
            err = max(err, abs(ex2 - lams['A%d' % m][dct[q]]))
    check('A: lambda = 2 Re((pi/|pi|)^m) (qfbsolve), max err %.1e' % err, err < 1e-9)
    # B
    p, lams, ex = L.lam_series('B', Lo, U, ok, ['B6', 'B12'])
    dct = {int(q): i for i, q in enumerate(p)}
    split = [q for q in pr if q % 3 == 1]
    check('B: every p = 1 mod 3 found once (%d)' % len(split), sorted(dct) == split and len(dct) == len(p))
    err = 0.0
    for q in split[:: max(1, len(split) // 300)]:
        a, b = [int(t) for t in pari.qfbsolve(pari.Qfb(1, 0, 3), q)]
        for m in (6, 12):
            ex2 = 2 * (complex(a, b * math.sqrt(3)) ** m).real / q ** (m / 2)
            err = max(err, abs(ex2 - lams['B%d' % m][dct[q]]))
    check('B: lambda = 2 Re((pi/|pi|)^m), max err %.1e' % err, err < 1e-9)
    # C23, C31: number of roots of the cubic mod p, minus 1
    for grp, pol in (('C23', [1, 0, -1, -1]), ('C31', [1, 0, 1, -1])):
        p, lams, ex = L.lam_series(grp, Lo, U, ok, [grp]); dct = {int(q): i for i, q in enumerate(p)}
        sample = [q for q in pr if q > 2][:: max(1, len(pr) // 400)]
        bad = 0
        for q in sample:
            nr = len(set(int(r) for r in pari.polrootsmod(pari.Pol(pol), q))) if True else 0
            want = nr - 1 if q != int(grp[1:]) else 1
            bad += abs(lams[grp][dct[q]] - want) > 1e-12
        check('%s: lambda = #roots - 1 on %d sampled primes' % (grp, len(sample)), bad == 0)
    # C47: against PARI's Hecke L-function coefficients (class group character [1], [2]) at primes
    p, lams, ex = L.lam_series('C47', Lo, U, ok, ['C47a', 'C47b']); dct = {int(q): i for i, q in enumerate(p)}
    if Lo == 0:
        for key, ch in (('C47a', '[1]'), ('C47b', '[2]')):
            an = pari('lfunan(lfuncreate([bnrinit(bnfinit(y^2+y+12,1),1), %s]), %d)' % (ch, 20000))
            bad = 0; n = 0
            for q in pr:
                if q > 20000: break
                if q == 2: continue
                n += 1; bad += abs(float(pari.real(an[q - 1])) - lams[key][dct[q]]) > 1e-9
            check('%s: lambda = PARI lfunan for p <= 2e4 (%d primes)' % (key, n), bad == 0)
    else:
        # class of p via qfbsolve on the three reduced forms
        bad = 0; n = 0
        for q in pr[:: max(1, len(pr) // 300)]:
            n += 1
            if pari.kronecker(-47, q) != 1: want = (0.0, 0.0)
            elif len(pari.qfbsolve(pari.Qfb(1, 1, 12), q)): want = (2.0, 2.0)
            elif len(pari.qfbsolve(pari.Qfb(2, 1, 6), q)): want = (L.C5a, L.C5b)
            else: want = (L.C5b, L.C5a)
            bad += abs(lams['C47a'][dct[q]] - want[0]) > 1e-12 or abs(lams['C47b'][dct[q]] - want[1]) > 1e-12
        check('C47: lambda from qfbsolve on %d primes near 1e8' % n, bad == 0)

# 3. lambda of A and B against the newforms' coefficients from PARI (p <= 2000)
for m, (N, chi) in ((4, (4, 'Mod(3,4)')), (8, (4, 'Mod(3,4)')), (6, (3, 'Mod(2,3)'))):
    k = m + 1; grp = 'A' if N == 4 else 'B'; key = '%s%d' % (grp, m)
    ok = L.sieve_odd(0, 4000); p, lams, ex = L.lam_series(grp, 0, 4000, ok, [key]); dct = {int(q): i for i, q in enumerate(p)}
    gen = '(2+I)' if N == 4 else '(2+sqrt(-3))'; p0 = 5 if N == 4 else 7; p1 = 3 if N == 4 else 5
    ap0 = pari('round(2*real(%s^%d))' % (gen, k - 1))
    r = pari('my(mf=mfinit([%d,%d,%s],0),B=mfeigenbasis(mf));select(F->mfcoefs(F,%d)[%d]==%s,B)[1]' % (N, k, chi, p0, p0 + 1, ap0))
    c = pari.mfcoefs(r, 2000); bad = 0
    for q in sympy.primerange(2, 2000):
        want = float(c[q]) / q ** ((k - 1) / 2)
        got = lams[key][dct[q]] if q in dct else L.special_primes(key).get(q, 0.0)
        bad += abs(want - got) > 1e-9
    check('%s: lambda(p) = a(p)/p^((k-1)/2) of the PARI newform for p < 2000 (incl. 2, ramified, inert)' % key, bad == 0)

# 4. Pythagorean directions and weights: sum rule c = (4/pi) prod (1 - 2cos(2 m theta_p)/(p-1)) equals the direct product (algebra)
dirs = L.pyth_directions(200)
w5 = [d for d in dirs if d[1] == 5]; w25 = [d for d in dirs if d[1] == 25]; w65 = [d for d in dirs if d[1] == 65]
check('directions: h=5 -> 2 rays weight -1/4; h=25 -> weight 0; h=65 -> 4 rays (2 product rays weight +1/48, 2 rays 13*(3,4)/5*(5,12) not primitive)',
      len(w5) == 2 and all(abs(d[2] + 0.25) < 1e-12 for d in w5) and all(d[2] == 0 for d in w25) and len(w65) == 4 and all(abs(d[2] - 1 / 48) < 1e-12 for d in w65))
# zone_F closed form against quadrature
for dl in (0.0, 0.7, -1.3):
    u = np.linspace(-40, 40, 400001); f = np.exp(-np.pi * (u + dl) ** 2 / L.BETA ** 2 - 1j * np.pi * u * u)
    q = np.trapezoid(f, u) if hasattr(np, 'trapezoid') else np.trapz(f, u)
    check('zone_F(%.1f) closed form vs quadrature: %.2e' % (dl, abs(q - L.zone_F(dl))), abs(q - L.zone_F(dl)) < 1e-6)
print('FAILS:', fails)
