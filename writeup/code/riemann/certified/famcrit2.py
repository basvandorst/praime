"""Family criterion for large X (see famcrit.py for the identity): every 3 <= n < X decided alone.

The class sum over the zeros, S(n) = sum_k chibar_k(n) F_k(ln n) with F_k(u) = sum_{gamma in chi_k} w(gamma) e^{i gamma u},
is evaluated by sampling every F_k on a grid in u (step du), one FFT over k per grid point, and degree-7 Lagrange
interpolation in u. Interpolation error, rigorously: |F_k^(8)| <= sum_{gamma in chi_k} |w| gamma^8, so the error of S is at
most 2 (sum over all zeros of |w| gamma^8)/8! du^8 max_{s in [0,1]} prod_{r=-3..4} |s - r|.
The zeros above T: counting bound of Bennett-Martin-O'Bryant-Rechnitzer; GRH to Platt's height; beyond, m integrations
by parts with the exact integral of |K^(m)| (kbderiv.py). zeta: its zeros up to T are added (here one, 14.1347).
Usage: famcrit2.py zeros.npz X c beta m cc eps out
"""
import sys, math, time
import numpy as np
from scipy.integrate import quad
from scipy.special import iv
from kbm import K, K0, Z as Zkb
from cert import tail_bound
from fambudget import J_L, H_platt
from kbderiv import Ij
import gmpy2

zf, X = sys.argv[1], int(sys.argv[2]); c, beta, m, cc, eps = map(float, sys.argv[3:8]); m = int(m); out = sys.argv[8]
D = np.load(zf); q = int(D['q']); T = float(D['T']); pw = D['pw']
assert abs(T - cc * beta / c) < 1e-9
assert (X - 1) * math.exp(c) < (X - 1) + q and 3 * math.exp(-c) > 1
ind = np.empty(q, dtype=np.int64); ind[pw] = np.arange(q - 1)
kz = D['k'].astype(np.int64); gz = D['gamma']; counts = D['count'][1:]
assert len(D['mism']) == 0, 'argument principle and zeros found disagree'
import os
ZZ = os.environ.get('ZETAZEROS', os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../data/riemann/zeros_hp2419.txt'))
zeta_zeros = np.array([float(l) for l in open(ZZ) if l.strip() and float(l) <= T])
k0 = K0(beta, m); C0 = math.sqrt(2 * math.pi) * beta ** m
wz = Zkb(c * gz, beta, m)
t0 = time.time()
print(f'q = {q}, T = {T:.6f}, {len(gz)} zeros with 0 < gamma <= T ({int(counts.sum())} with |gamma| <= T); zeta zeros below T: {len(zeta_zeros)}', flush=True)

# grid in u = ln n and interpolation weights
du = 0.02
ns = np.arange(3, X)
xg = (np.log(ns.astype(float)) - (math.log(3) - 4 * du)) / du
fl = np.floor(xg).astype(np.int64); sfrac = xg - fl
nodes = np.arange(-3, 5)
Wt = np.ones((len(ns), 8))
for a_, ra in enumerate(nodes):
    for rb in nodes:
        if rb != ra: Wt[:, a_] *= (sfrac - rb) / (ra - rb)
ss = np.linspace(0, 1, 100001); PImax = float(np.max(np.abs(np.prod([ss - r for r in nodes], axis=0))))
nu_grid = int(fl.max()) + 5
idn = ind[ns]
S = np.zeros(len(ns)); Fnorm = 0.0
for i in range(nu_grid):
    u = math.log(3) - 4 * du + i * du
    ph = gz * u
    F = np.bincount(kz, weights=wz * np.cos(ph), minlength=q - 1) + 1j * np.bincount(kz, weights=wz * np.sin(ph), minlength=q - 1)
    Fnorm = max(Fnorm, float(np.linalg.norm(F)))
    G = np.fft.fft(F)
    lo = np.searchsorted(fl, i - 4, side='left'); hi = np.searchsorted(fl, i + 3, side='right')    # i - 4 <= fl <= i + 3, i.e. node r = i - fl in -3 .. 4
    if hi > lo:
        r = i - fl[lo:hi]                             # node index in -3 .. 4
        S[lo:hi] += 2 * np.real(Wt[lo:hi, :][np.arange(hi - lo), r + 3] * G[idn[lo:hi]])
    if i % 100 == 0: print(f'grid point {i}/{nu_grid}  {time.time() - t0:.0f} s', flush=True)
print(f'zero sums for {len(ns)} numbers by FFT + interpolation: {time.time() - t0:.0f} s', flush=True)

# error budget, vectorised (all terms in units of (q-1) K0)
x = ns.astype(float); sq = np.sqrt(x); a = c * x; lnx = np.log(x)
nu = m + 0.5; qq = math.sqrt(beta ** 2 + (c / 2) ** 2)
Px = c * x * C0 * iv(nu, qq) / qq ** nu
khat0 = C0 * float(Zkb(0.0, beta, m))                 # Khat(0) / a
vals, mult = np.unique(counts, return_counts=True)
Jsum = sum(J_L(T, int(v), q, c, beta, m) * int(mu) for v, mu in zip(vals, mult))
tailL = c * sq * C0 * Jsum
tb1000 = tail_bound(1000.0, T, len(zeta_zeros), c * 1000.0, beta, m)
tailZ = 2 * sq * tb1000 / 1000.0
Iex = Ij(beta, m, m)
H = H_platt(q); j = m
SL = (j / math.pi) * H ** (1 - j) * (math.log(q * H) / (j - 1) + 1 / (j - 1) ** 2)
offL = (q - 2) * math.exp(c / 2) * Iex * c ** (1 - m) * x * SL
HZ = 3e12; SZ = (j / (2 * math.pi)) * HZ ** (1 - j) * (math.log(HZ) / (j - 1) + 1 / (j - 1) ** 2)
offZ = 2 * math.exp(c / 2) * Iex * c ** (1 - m) * x * SZ
Wabs = float(np.sum(np.abs(wz))); Wg8 = float(np.sum(np.abs(wz) * gz ** 8)); Nz = len(gz)
prec = (eps / sq) * (lnx * c * x * C0 * 2 * Wabs + 2 * Nz * c * c * x * khat0) * 1.01   # |d/dgamma Khat(gamma/x)| <= a Khat(0)/x
interp = c * sq * C0 * 2 * Wg8 / math.factorial(8) * du ** 8 * PImax
fft_err = c * sq * C0 * 2 * 8 * (3 * math.log2(q - 1) * 1.2e-16 * math.sqrt(q - 1) * Fnorm)
rnd = c * sq * C0 * 2 * Wabs * (T * lnx * 4.5e-16 + 3e-15) + fft_err + 1e-14 * Px
# trivial zeros: quadrature below 2000, a bound above (the terms are then left out of the estimate)
tmin = x * math.exp(-c)
Ee = np.zeros(len(ns)); Ez = np.zeros(len(ns))
trivb = math.exp(c / 2) * c * khat0 * (tmin ** 2 / (tmin ** 2 - 1) + 1 / (tmin ** 2 - 1))
for idx in range(min(len(ns), 2000 - 3)):
    xx = x[idx]; aa = a[idx]
    fz = lambda v: float(K(v, aa, beta, m)) * math.sqrt(xx * math.exp(v / xx)) / (xx * ((xx * math.exp(v / xx)) ** 3 - xx * math.exp(v / xx)))
    fe = lambda v: float(K(v, aa, beta, m)) * (xx * math.exp(v / xx)) ** 1.5 / (xx * ((xx * math.exp(v / xx)) ** 2 - 1))
    ez_, e1 = quad(fz, -aa, aa, limit=400, epsabs=0, epsrel=1e-12); ee_, e2 = quad(fe, -aa, aa, limit=400, epsabs=0, epsrel=1e-12)
    Ez[idx] = ez_ * math.sqrt(xx); Ee[idx] = ee_ * math.sqrt(xx); trivb[idx] = 10 * (e1 + e2) * math.sqrt(xx)
zsum = 2 * c * sq * C0 * np.sum(np.cos(np.outer(lnx, zeta_zeros)) * Zkb(c * zeta_zeros, beta, m), axis=1)
lam_hat = (Px - Ez + Ee - zsum - c * sq * C0 * S) / ((q - 1) * k0)
delta = (tailL + tailZ + offL + offZ + prec + interp + rnd + trivb) / ((q - 1) * k0)

truth = np.array([gmpy2.is_prime(int(v)) for v in ns])
lam_true = np.zeros(len(ns))
lam_true[truth] = lnx[truth]
for p in range(2, int(math.isqrt(X)) + 1):
    if gmpy2.is_prime(p):
        v = p * p
        while v < X:
            lam_true[v - 3] = math.log(p); v *= p
err = np.abs(lam_hat - lam_true)
decided = delta < 0.25 * lnx
prime = lam_hat > 0.75 * lnx
wrong = np.nonzero(decided & (prime != truth))[0]
assert np.all(err <= delta), ('bound violated at', ns[np.argmax(err - delta)])
iw = int(np.argmax(delta / lnx))
parts = {k_: float(v_[iw] / ((q - 1) * k0)) for k_, v_ in dict(tail_L=tailL, tail_zeta=tailZ, above_Platt=offL, above_3e12=offZ, precision=prec, interpolation=interp, rounding=rnd, trivial=trivb).items()}
lines = [f'family criterion mod {q}: c = {c}, beta = {beta}, m = {m}, T = {T:.6f}, eps = {eps}; Platt height {H:.1f}',
         f'zeros used: {int(counts.sum())} of {q - 2} L-functions with |gamma| <= T, and {len(zeta_zeros)} zero of zeta',
         f'numbers decided: {int(decided.sum())} of {len(ns)} (3 <= n < {X}); primes {int((decided & prime).sum())}, true {int(truth.sum())}; wrong {len(wrong)}',
         f'largest delta/(1/4 ln n) = {float(np.max(delta / (0.25 * lnx))):.3f} at n = {int(ns[iw])}: delta/K0 = {delta[iw]:.4f}; parts ' + ', '.join(f'{k_} {v_:.2e}' for k_, v_ in parts.items()),
         f'largest |Lambda_hat - Lambda| = {err.max():.3e} at n = {int(ns[np.argmax(err)])}; largest ratio to delta {float(np.max(err / delta)):.3e}; time {time.time() - t0:.0f} s']
# independent path: 300 random n summed directly over all zeros (no FFT, no interpolation)
rng = np.random.default_rng(11); pick = rng.choice(len(ns), 300, replace=False)
dmax = 0.0
for i_ in pick:
    nn = int(ns[i_]); ph = gz * math.log(nn) - 2 * math.pi * ((idn[i_] * kz) % (q - 1)) / (q - 1)
    Sd = 2 * float(np.sum(wz * np.cos(ph)))
    dmax = max(dmax, abs(Sd - S[i_]) * c * sq[i_] * C0 / ((q - 1) * k0))
lines.append(f'direct sum over all zeros for 300 random n: largest difference from the FFT + interpolation value {dmax:.2e} (in units of K0; interpolation bound there at most {float(np.max(interp)) / ((q - 1) * k0):.1e})')
print('\n'.join(lines))
open(out + '.txt', 'w').write('\n'.join(lines) + '\n')
np.save(out + '.npy', np.stack([ns, lam_hat, lam_true, delta]))
