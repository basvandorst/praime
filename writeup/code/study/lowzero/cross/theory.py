"""Theory numbers for the cross-correlation of zeros of two Dirichlet L-functions (derivation.md).
(T1) per-prime covariance A_p(w), B_p(w) on the raw grid w = 0.3 .. 4.8 and at fixed l for t = 100, 350, 600
(T2) the correlation coefficient against the BK number variance; spread of the full-law Cov over the 546 old pairs
(T3) the closed form (1/pi^2) ln|L(1,psi)/L(1+iw,psi)| - corr2 (PARI) against the prime sums (truncation)
(T4) the recipe's two-point function (1/2pi^2) Re D_psi(-i eps), integrated over the window square, against (T3)
(T5) the swapped (Gram-lattice) term: oddness for |D1| = |D2|, averaging factors for near-equal conductors, sizes
usage: theory.py  -> theory_out.txt (stdout)"""
import os, sys, json, math, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from crosslib import *
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', '..')); LOW = os.path.join(REPO, 'writeup/code/study/lowzero')
sys.path.insert(0, os.path.join(LOW, 'amp')); sys.path.insert(0, os.path.join(LOW, 'bk'))
import bklib, gridcache
import cypari2
pari = cypari2.Pari(); pari.allocatemem(500 * 10 ** 6)
np.set_printoptions(linewidth=200)
W = 0.3 * np.arange(1, 17)
PJ = json.load(open(os.path.join(LOW, 'triple', 'predictions_kappa3.json')))
OLD = PJ['R120']['D'] + PJ['B']['D'] + PJ['Dn']['D']; NEW = PJ['Fn']['D']
print('design: %d old real L-functions (R120 %d, B %d, Dn %d), %d new (Fn); heights 100-600' % (len(OLD), len(PJ['R120']['D']), len(PJ['B']['D']), len(PJ['Dn']['D']), len(NEW)))

# ---------------- T1
print('\n(T1) per-prime diagonal covariance, real characters: Cov_p = A_p(w) chi1(p)chi2(p) + B_p(w) 1[p not | k1k2]')
print('     raw window w (height units):   ' + ' '.join('%6.2f' % w for w in W))
for L in (9.08, 10.33, 10.87, 11.07):
    print('     l = wL/2pi at L=%5.2f:          ' % L + ' '.join('%6.2f' % (w * L / 2 / math.pi) for w in W))
for p in (2, 3, 5, 7, 11, 13):
    a, b = AB(p, W)
    print('  p=%-2d A_p ' % p + ' '.join('%+.4f' % x for x in a)); print('       B_p ' + ' '.join('%+.4f' % x for x in b))
print('  peak of the m = 1 term: (2/pi^2)/p at w = pi/ln p: ' + ', '.join('p=%d %.4f at w=%.2f' % (p, 2 / math.pi ** 2 / p, math.pi / math.log(p)) for p in (2, 3, 5, 7)))
print('\n  at fixed l (mean spacings of a conductor-550 L-function), heights t = 100 / 350 / 600:')
LS = 0.5 * np.arange(1, 17)
for t in (100, 350, 600):
    L = math.log(550 * t / 2 / math.pi); w = 2 * math.pi * LS / L
    print('   t=%d L=%.2f' % (t, L))
    for p in (2, 3, 5, 7):
        a, b = AB(p, w); print('     p=%d A ' % p + ' '.join('%+.4f' % x for x in a))

# ---------------- T2
print('\n(T2) correlation coefficient: Cov_p / Sigma2(l), Sigma2 from the BK law (bklib, frozen grid), k without small bad primes')
G = gridcache.get()
for L, B in ((10.33, ()), (10.33, (2,)), (11.07, ())):
    s2 = bklib.sigma2_full(G, L, B); lg = bklib.LGRID
    ell = W * L / 2 / math.pi; ok = ell <= 8.0
    V = np.interp(ell[ok], lg, s2)
    print('  L=%.2f bad primes %s: Sigma2 at w: ' % (L, B) + ' '.join('%.3f' % v for v in V))
    for p in (2, 3, 5, 7):
        a, _ = AB(p, W[ok]); print('     rho_%d = A_p/Sigma2: ' % p + ' '.join('%+.3f' % x for x in a / V))
# spread of the full law over the old pairs (P <= 2e5)
P = primes_upto(200000); X = kron_table(OLD, P)
t0 = time.time(); C = pair_cov_full(X, W, P); iu = np.triu_indices(len(OLD), 1)
print('  full-law Cov over the %d old pairs (primes <= 2e5, %.0fs): mean / sd / 1%% / 99%% / max|.|' % (len(iu[0]), time.time() - t0))
for k, w in enumerate(W):
    c = C[k][iu]; print('   w=%.1f  %+.4f %.4f %+.4f %+.4f %.4f' % (w, c.mean(), c.std(), np.percentile(c, 1), np.percentile(c, 99), np.abs(c).max()))
# decomposition of the spread: p<=7 vs the rest
X7 = X[:, :4]; P7 = P[:4]; C7 = pair_cov_full(X7, W, P7)
print('  sd over pairs of the part from p <= 7 / from p > 7:')
print('   ' + ' '.join('%.4f/%.4f' % (C7[k][iu].std(), (C[k] - C7[k])[iu].std()) for k in range(len(W))))

# ---------------- T3
print('\n(T3) closed form with PARI L-values against the prime sums')
def psi_L(D1, D2):
    d = int(pari.quaddisc(D1 * D2)); Lf = pari.lfuncreate(d)
    bad = [int(p) for p in pari.factor(abs(D1 * D2))[0]]; rem = [p for p in bad if d % p != 0]
    chi = lambda p: int(pari.kronecker(d, p)); return d, Lf, [(p, chi(p)) for p in rem]
def logabsL(Lf, rem, s):
    v = complex(pari.lfun(Lf, s)); out = math.log(abs(v))
    for p, c in rem: out += math.log(abs(1 - c * complex(p) ** (-complex(s))))
    return out
def corr2(D1, D2, w, PP):
    x = kc(D1, PP) * kc(D2, PP); s = 0.0
    for p, c in zip(PP, x):
        if c == 0: continue
        m = 2
        while float(p) ** m < 1e16:
            s += (1 / m - 1 / m ** 2) * (c ** m) * (float(p) ** -m) * (1 - math.cos(m * w * math.log(p))); m += 1
    return s
P6 = primes_upto(10 ** 6); PS = primes_upto(10 ** 4)
_KC = {}
def kc(D, PB):
    key = (D, len(PB))
    if key not in _KC: _KC[key] = kron_table([D], PB)[0]
    return _KC[key]
pairs = [(401, -403), (-995, 997), (929, -647), (-163, 1001), (8, -8 * 13)]
for D1, D2 in pairs:
    d, Lf, rem = psi_L(D1, D2); L1 = logabsL(Lf, rem, 1)
    for w in (0.6, 2.4, 4.5):
        cl = (L1 - 0.5 * (logabsL(Lf, rem, complex(1, w)) + logabsL(Lf, rem, complex(1, -w))) - corr2(D1, D2, w, PS)) / math.pi ** 2
        x2 = kron_table([D1, D2], P); s2 = pair_cov_full(x2, np.array([w]), P)[0, 0, 1]
        x6 = kron_table([D1, D2], P6); s6 = pair_cov_full(x6, np.array([w]), P6)[0, 0, 1]
        print('  D1=%5d D2=%5d (d*=%7d) w=%.1f: closed form %+.6f | primes<=2e5 %+.6f (diff %+.1e) | primes<=1e6 %+.6f (diff %+.1e)' % (D1, D2, d, w, cl, s2, s2 - cl, s6, s6 - cl))

# ---------------- T4
print('\n(T4) recipe two-point function integrated over the window square: (1/2pi^2) int_{-w}^{w} (w - |e|) Re D_psi(-i e) de')
def Dpsi(Lf, rem, d, D1, D2, s, PB):
    """D_psi at s = 1 + z: (L'/L)'(s, psi) - B_psi(z); psi = chi_d* without the primes dividing D1 D2"""
    v0 = complex(pari.lfun(Lf, s)); v1 = complex(pari.lfun(Lf, s, 1)); v2 = complex(pari.lfun(Lf, s, 2))
    out = v2 / v0 - (v1 / v0) ** 2; s = complex(s)
    for p, c in rem:
        u = c * p ** (-s); out += -c * math.log(p) ** 2 * p ** (-s) / (1 - u) ** 2
    x = (kc(D1, PB) * kc(D2, PB)).astype(float); lp = np.log(PB.astype(float)); u = x * PB.astype(float) ** (-s)
    out -= np.sum(lp ** 2 * u ** 2 / (1 - u) ** 2)
    return out
xg, wg = np.polynomial.legendre.leggauss(20)
for D1, D2 in pairs[:2]:
    d, Lf, rem = psi_L(D1, D2); L1 = logabsL(Lf, rem, 1)
    for w in (1.2, 4.5):
        e = w / 2 * (xg + 1); we = w / 2 * wg
        vals = np.array([Dpsi(Lf, rem, d, D1, D2, complex(1, -x), P6).real for x in e])
        rec = 2 * np.sum(we * (w - e) * vals) / (2 * math.pi ** 2)
        cl = (L1 - 0.5 * (logabsL(Lf, rem, complex(1, w)) + logabsL(Lf, rem, complex(1, -w))) - corr2(D1, D2, w, PS)) / math.pi ** 2
        print('  D1=%d D2=%d w=%.1f: recipe integral %+.6f, closed form %+.6f, diff %+.1e; Re D_psi(0) = %+.4f' % (D1, D2, w, rec, cl, rec - cl, Dpsi(Lf, rem, d, D1, D2, 1, P6).real))

# ---------------- T5
print('\n(T5) the swapped (Gram-lattice) term')
def H12(D1, D2, eps, Lf, rem, PB):
    """(k1/phi(k1))(k2/phi(k2)) L(1-i e, psi) L(1+i e, psi) A12(-i e) for real psi"""
    k1, k2 = abs(D1), abs(D2); f1 = np.prod([p / (p - 1) for p in set(int(q) for q in pari.factor(k1)[0])]); f2 = np.prod([p / (p - 1) for p in set(int(q) for q in pari.factor(k2)[0])])
    def Lpsi(s):
        v = complex(pari.lfun(Lf, s))
        for p, c in rem: v *= 1 - c * complex(p) ** (-complex(s))
        return v
    c1 = kc(D1, PB) != 0; c2 = kc(D2, PB) != 0; x = (kc(D1, PB) * kc(D2, PB)).astype(float)
    Pf = PB.astype(float); u = x * Pf ** (-1 + 1j * eps)
    A = np.prod((1 - u) * (1 - c1 / Pf - c2 / Pf + u) / ((1 - c1 / Pf) * (1 - c2 / Pf)))
    return f1 * f2 * Lpsi(complex(1, -eps)) * Lpsi(complex(1, eps)) * A
for D1, D2 in ((104, -104), (8 * 37, -8 * 37), (401, -403)):
    d, Lf, rem = psi_L(D1, D2)
    h = [H12(D1, D2, e, Lf, rem, P6) for e in (0.0, 1.0, -1.0, 2.5, -2.5)]
    print('  D1=%d D2=%d: H(0) = %.4f%+.4fi; H(1) = %.4f%+.4fi, H(-1) = %.4f%+.4fi; H(2.5) = %.4f%+.4fi, H(-2.5) = %.4f%+.4fi' % (D1, D2, *sum(([v.real, v.imag] for v in h), [])))
print('  -> H(-e) = conj H(e) (real psi); with |D1| = |D2| and opposite parity the phase is exp(+-i pi/2) exp(2 pi i r):')
print('     Re[+-i e^{i e L} H(e)] is odd in e, so its window-square integral is 0 for every w (equal windows).')
# averaging factors over the height range
tt = np.linspace(100, 600, 5001)
def favg(D1, D2):
    q1, q2 = abs(D1), abs(D2); a1 = 0 if D1 > 0 else 1; a2 = 0 if D2 > 0 else 1
    ph = 2 * (theta(tt, q2, a2) - theta(tt, q1, a1)); return np.abs(np.mean(np.exp(1j * ph)))
for name, A_, B_ in (('old x old', OLD, None), ('new x new', NEW, None), ('new x old', NEW, OLD)):
    F = []
    if B_ is None:
        for i in range(len(A_)):
            for j in range(i + 1, len(A_)):
                if abs(math.log(abs(A_[i]) / abs(A_[j]))) < 0.05: F.append((favg(A_[i], A_[j]), A_[i], A_[j]))
    else:
        for a in A_:
            for b in B_:
                if abs(math.log(abs(a) / abs(b))) < 0.05: F.append((favg(a, b), a, b))
    F = np.array([f[0] for f in F]); n_same = 0
    print('  %s: pairs with |ln(k1/k2)| < 0.05: %d; averaging factor |<exp 2i(theta2-theta1)>| > 0.1: %d, > 0.3: %d, > 0.6: %d (|k1| = |k2| pairs have 1 and give 0 by oddness)' % (name, len(F), (F > 0.1).sum(), (F > 0.3).sum(), (F > 0.6).sum()))
L = 10.3
for w in (0.3, 0.9, 1.5, 3.0, 4.5):
    print('  window factor (1/2pi^2)|int int e^{iL(y2-y1)}| = (1/2pi^2)(2 sin(Lw/2)/L)^2 at L=10.3, w=%.1f: %.5f' % (w, (2 * math.sin(L * w / 2) / L) ** 2 / (2 * math.pi ** 2)))
