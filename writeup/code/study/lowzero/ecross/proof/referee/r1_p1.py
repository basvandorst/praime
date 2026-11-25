"""Referee check 1 (own code): Proposition 1 (exact explicit formula) three ways, one curve, Gaussian-edged box.

  Z  = sum over the zero file of H(gamma) e^{2 i theta(gamma)}                      (theta by mpmath)
  R  = (1/2pi) int h(5/4+it) X^{-1}(5/4+it) (L'/L)(5/4+it) dt                       (L'/L by PARI lfun = approximate FE)
  E0 = -(1/2pi) int h(-1/4+it) X^{-1}(-1/4+it) (L'/L)(-1/4+it) dt                   (left edge, PARI lfun)
  P  = -(1/2pi) sum_n Lambda_E(n) n^{-1/2} int H(t) e^{i(2theta(t) - t ln n)} dt     (Lambda_E from ellan by the
        convolution recursion lambda(p^k) k ln p = sum_j Lambda_E(p^j) lambda(p^{k-j}); I_n by the trapezoid rule)
Residue theorem:  Z = R + E0.   Per-n shift (Prop. 1 step 4):  R = P.
usage: r1_p1.py CID T1 T2 Delta [kind]   kind = box (default) | hann  (hann: W = sin^2 on [T1,T2], convolved with g_Delta)
"""
import os, sys, json, math, cmath, time
import numpy as np, mpmath as mp, cypari2
from scipy.special import erfc

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), *(['..'] * 7)))
cid = sys.argv[1]; T1 = float(sys.argv[2]); T2 = float(sys.argv[3]); D = float(sys.argv[4])
kind = sys.argv[5] if len(sys.argv) > 5 else 'box'
C = {c['id']: c for c in json.load(open(REPO + '/writeup/code/study/lowzero/deg2/fresh_set.json'))['curves']}[cid]
N = C['N']
gam = np.sort(np.load(REPO + '/writeup/data/study/r144/fresh_zeros.npz')[cid])
mp.mp.dps = 25
pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6); pari.set_real_precision(38)
PI = math.pi
t_start = time.time()

def theta(t):
    return float(t / 2 * mp.log(N) - t * mp.log(2 * mp.pi) + mp.im(mp.loggamma(1 + 1j * mp.mpf(t))))
thetav = np.vectorize(theta)

def erfdiff(a, b):
    """erf(a) - erf(b) for complex arrays, stable in both tails."""
    a = np.asarray(a, complex); b = np.asarray(b, complex)
    right = (a.real + b.real) > 0
    return np.where(right, erfc(b) - erfc(a), erfc(-a) - erfc(-b))

def H(z):
    """H = W * g_D at complex z (entire)."""
    z = np.asarray(z, complex)
    if kind == 'box':
        return 0.5 * erfdiff((z - T1) / D, (z - T2) / D)
    # W(u) = sin^2(pi(u-T1)/L) on [T1,T2] = 1/2 - (e^{ik(u-T1)} + e^{-ik(u-T1)})/4
    L = T2 - T1; k = PI / L
    out = 0.25 * erfdiff((z - T1) / D, (z - T2) / D)
    for kk in (k, -k):
        # int_{T1}^{T2} e^{i kk (u-T1)} g_D(z-u) du = e^{i kk (z-T1)} e^{-kk^2 D^2/4} * 1/2 [erf((z-T1)/D + i kk D/2) - erf((z-T2)/D + i kk D/2)]
        c = 1j * kk * D / 2
        out = out - 0.25 * np.exp(1j * kk * (z - T1) - kk * kk * D * D / 4) * 0.5 * erfdiff((z - T1) / D + c, (z - T2) / D + c)
    return out

def Xinv(s):
    s = mp.mpc(s)
    return complex(mp.power(N, s - 0.5) * mp.power(2 * mp.pi, 1 - 2 * s) * mp.gamma(s + 0.5) / mp.gamma(1.5 - s))

lo, hi = T1 - 9 * D, T2 + 9 * D

# ---------------- Z
Zs = np.sum(H(gam) * np.exp(2j * thetav(gam)))
print('%s N=%d %s [%g,%g] Delta=%g; zeros in file %d (max %.2f); Z = %.12f %+.12fi  (%.0fs)' % (cid, N, kind, T1, T2, D, len(gam), gam.max(), Zs.real, Zs.imag, time.time() - t_start))

# ---------------- R and E0 by PARI lfun on the two vertical lines (trapezoid, step 0.05)
pari('E = ellinit(%s); Lf = lfuncreate(E); Li = lfuninit(Lf, [1, 0.8, %d]);' % (str(C['ainv']), int(hi + 5)))
def line(sig, step=0.05):
    tt = np.arange(lo, hi + step / 2, step)
    vals = np.empty(len(tt), complex)
    for j, t in enumerate(tt):
        sa = '%r + I*%r' % (float(sig + 0.5), float(t))
        LL = complex(pari('lfun(Li, %s, 1)/lfun(Li, %s)' % (sa, sa)))
        vals[j] = H(t - 1j * (sig - 0.5)) * Xinv(mp.mpc(sig, t)) * LL
    w_ = np.full(len(tt), step); w_[0] = w_[-1] = step / 2
    return np.sum(vals * w_) / (2 * PI), np.abs(vals[0]), np.abs(vals[-1])
R, r0, r1 = line(1.25)
E0n, e0, e1 = line(-0.25); E0 = -E0n
print('R (5/4 line) = %.12f %+.12fi | E0 (-1/4 line) = %.3e %+.3ei | endpoint |integrand| %.1e %.1e %.1e %.1e  (%.0fs)' % (
    R.real, R.imag, E0.real, E0.imag, r0, r1, e0, e1, time.time() - t_start))

# ---------------- P: prime powers with |ln n - 2theta'(t)| <= margin somewhere on [T1, T2]
def two_theta_prime(t):
    return float(mp.log(N) - 2 * mp.log(2 * mp.pi) + 2 * mp.re(mp.digamma(1 + 1j * mp.mpf(t))))
margin = float(sys.argv[6]) if len(sys.argv) > 6 else 3.5   # distance (in ln n) beyond the edge frequencies 2theta'(T1), 2theta'(T2)
xl = math.exp(two_theta_prime(T1) - margin); xh = math.exp(two_theta_prime(T2) + margin)
Xmax = int(xh) + 1
a = np.array(pari('Vecsmall(ellan(E, %d))' % Xmax), dtype=float)          # a[n-1] = a_n
lam = a / np.sqrt(np.arange(1, Xmax + 1, dtype=float))
primes = np.array(pari('Vecsmall(primes([2, %d]))' % Xmax), dtype=np.int64)
ns, Ls = [], []
for p in primes:
    lp = math.log(p)
    Lam = []                                                        # Lambda_E(p^j), j = 1..
    q = int(p); k = 1
    while q <= Xmax:
        # lambda(p^k) k ln p = sum_{j=1}^{k} Lambda_E(p^j) lambda(p^{k-j})
        val = lam[q - 1] * k * lp - sum(Lam[j - 1] * lam[int(p) ** (k - j) - 1] for j in range(1, k))
        Lam.append(val)
        if q >= xl:
            ns.append(q); Ls.append(val)
        q *= int(p); k += 1
    if p * p > Xmax:
        # remaining primes: only k = 1
        break
rest = primes[primes * primes > Xmax]
rest = rest[rest >= xl]
ns = np.concatenate([np.array(ns, float), rest.astype(float)])
Ls = np.concatenate([np.array(Ls), lam[rest - 1] * np.log(rest.astype(float))])
del a, lam
print('prime powers in [%.0f, %.0f]: %d  (%.0fs)' % (xl, xh, len(ns), time.time() - t_start))

ht = 1 / 16
tg = np.arange(lo, hi + ht / 2, ht)
F = H(tg) * np.exp(2j * thetav(tg)) * ht
F[0] *= 0.5; F[-1] *= 0.5
lnn = np.log(ns)
tot = 0j; tot_by_margin = {}
twoth = np.array([two_theta_prime(x) for x in (T1, T2)])
for s in range(0, len(ns), 1500):
    I = np.exp(-1j * np.outer(lnn[s:s + 1500], tg)) @ F
    contrib = Ls[s:s + 1500] * ns[s:s + 1500] ** -0.5 * I
    tot += contrib.sum()
    d = np.maximum(twoth[0] - lnn[s:s + 1500], lnn[s:s + 1500] - twoth[1])   # distance outside [2th'(T1), 2th'(T2)]
    for mg in (1.0, 2.0, 3.0, margin - 0.5):
        tot_by_margin[mg] = tot_by_margin.get(mg, 0j) + contrib[d <= mg].sum()
P = -tot / (2 * PI)
print('P = %.12f %+.12fi   (%.0fs)' % (P.real, P.imag, time.time() - t_start))
for mg in (1.0, 2.0, 3.0, margin - 0.5):
    Pm = -tot_by_margin[mg] / (2 * PI)
    print('   P restricted to n within %.1f of [2th\'(T1), 2th\'(T2)]: |P_m - P| = %.1e' % (mg, abs(Pm - P)))
print('|Z - P| = %.2e   |Z - (R + E0)| = %.2e   |R - P| = %.2e   |E0| = %.2e' % (abs(Zs - P), abs(Zs - R - E0), abs(R - P), abs(E0)))
ok = abs(Zs - P) < 1e-8 and abs(Zs - R - E0) < 1e-8 and abs(R - P) < 1e-8
print('r1', 'PASS' if ok else 'FAIL', '(%.0fs)' % (time.time() - t_start))
