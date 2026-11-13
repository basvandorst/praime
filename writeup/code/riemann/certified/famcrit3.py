"""Local prime criterion over the family mod q, for any prime q (also q < n): the test of the uncertainty law.

For n coprime to q take x = n, a = c n with e^c <= 1 + q/n and e^-c >= 1 - q/n, so the support n e^-c < m < n e^c holds no
other member of n's class. Orthogonality and the explicit formula give exactly
 (q-1) Lambda(n) K(0) = P(n) - E_zeta - Z_zeta(n) - corr_q(n) - sum_{chi != chi0} chibar(n) Z_chi(n) - E_even S_even(n) - E_odd S_odd(n)
(Z_.: zero sums; corr_q: the powers of q inside the support, which chi0 omits; S_even = (q-1)/2 [n = +-1] - 1,
S_odd = (q-1)/2 ([n = 1] - [n = -1]) mod q). Height T = cc beta / c. For several beta the script reports how many n are
certified (delta < (1/4) ln n) and how many are decided correctly (Lambda_hat > (3/4) ln n iff prime), and the zeros used.
Usage: famcrit3.py zeros.npz X c m cc eps betas(comma) out"""
import sys, math, time, json
import numpy as np
from scipy.integrate import quad
from scipy.special import iv
from kbm import K, K0, Z as Zkb
from cert import tail_bound
from fambudget import J_L, H_platt
from kbderiv import Ij
import gmpy2
zf, X = sys.argv[1], int(sys.argv[2]); c = float(sys.argv[3]); m = int(sys.argv[4]); cc = float(sys.argv[5]); eps = float(sys.argv[6])
betas = [float(b) for b in sys.argv[7].split(',')]; out = sys.argv[8]
D = np.load(zf); q = int(D['q']); Tl = float(D['T']); pw = D['pw']; kz = D['k'].astype(np.int64); gz = D['gamma']
assert ('mism' not in D.files) or len(D['mism']) == 0
N = q - 1; ind = np.empty(q, dtype=np.int64); ind[pw] = np.arange(N); kbar = (-np.arange(N)) % N
ns = np.array([n for n in range(3, X) if n % q != 0], dtype=np.int64)
assert np.all(ns * math.exp(c) < ns + q) and np.all(ns * math.exp(-c) > np.maximum(ns - q, 1)), 'support must isolate n in its class'
import os
zeta = np.array([float(l) for l in open(os.environ.get('ZETAZEROS', os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../data/riemann/zeros_hp2419.txt'))) if l.strip()])
x = ns.astype(float); lnx = np.log(x); sq = np.sqrt(x); idn = ind[ns % q]
r = ns % q; S_even = np.where((r == 1) | (r == q - 1), N / 2, 0.0) - 1.0; S_odd = np.where(r == 1, N / 2, 0.0) - np.where(r == q - 1, N / 2, 0.0)
truth = np.array([gmpy2.is_prime(int(v)) for v in ns])
lam_true = np.where(truth, lnx, 0.0)
for p in range(2, int(math.isqrt(X)) + 1):
    if gmpy2.is_prime(p):
        v = p * p
        while v < X:
            if v % q: lam_true[np.searchsorted(ns, v)] = math.log(p)
            v *= p
results = []
t0 = time.time()
for beta in betas:
    T = cc * beta / c
    assert T <= Tl + 1e-9, (T, Tl)
    k0 = K0(beta, m); C0 = math.sqrt(2 * math.pi) * beta ** m
    sel = gz <= T; g = gz[sel]; kk = kz[sel]; w = 2 * Zkb(c * g, beta, m)
    S = np.zeros(len(ns))
    for i in range(0, len(ns), 250):
        for jz in range(0, len(g), 40000):
            ph = np.outer(lnx[i:i + 250], g[jz:jz + 40000]) - (2 * math.pi / N) * (np.outer(idn[i:i + 250], kk[jz:jz + 40000]) % N)
            S[i:i + 250] += np.cos(ph) @ w[jz:jz + 40000]
    zz = zeta[zeta <= T]
    zpart = 2 * c * sq * C0 * (np.cos(np.outer(lnx, zz)) @ Zkb(c * zz, beta, m)) if len(zz) else np.zeros(len(ns))
    nu = m + 0.5; qq = math.sqrt(beta ** 2 + (c / 2) ** 2)
    Px = c * x * C0 * iv(nu, qq) / qq ** nu
    corr = np.zeros(len(ns))
    e = 1
    while q ** e < X * math.exp(c) + 1:
        mm = float(q ** e); inside = (x * math.exp(-c) < mm) & (mm < x * math.exp(c))
        if inside.any():
            corr[inside] += math.log(q) * np.sqrt(x[inside] / mm) * K(x[inside] * np.log(mm / x[inside]), c * x[inside], beta, m)
        e += 1
    # trivial zeros: quadrature for n < 300, a bound above (left out of the estimate, put into delta)
    khat0 = C0 * float(Zkb(0.0, beta, m)); tmin = x * math.exp(-c)
    # each of E_zeta, E_even, E_odd is at most (sqrt x / x) int K t^(3/2)/(t^2-1) dv <= e^(c/2) t_min^2/(t_min^2-1) c Khat(0)/a
    trivb = math.exp(c / 2) * c * khat0 * (tmin ** 2 / (tmin ** 2 - 1)) * (1 + np.abs(S_even) + np.abs(S_odd))
    Etriv = np.zeros(len(ns))
    for i in np.nonzero((x < 300) | (r == 1) | (r == q - 1))[0]:   # small n, and n = +-1 mod q where the even and odd sums are large
        xx = x[i]; aa = c * xx
        tt = lambda v: xx * math.exp(v / xx)
        Ez = quad(lambda v: float(K(v, aa, beta, m)) * math.sqrt(tt(v)) / (xx * (tt(v) ** 3 - tt(v))), -aa, aa, epsrel=1e-12, limit=200)
        Ee = quad(lambda v: float(K(v, aa, beta, m)) * tt(v) ** 1.5 / (xx * (tt(v) ** 2 - 1)), -aa, aa, epsrel=1e-12, limit=200)
        Eo = quad(lambda v: float(K(v, aa, beta, m)) * tt(v) ** 0.5 / (xx * (tt(v) ** 2 - 1)), -aa, aa, epsrel=1e-12, limit=200)
        Etriv[i] = math.sqrt(xx) * (Ez[0] + Ee[0] * S_even[i] + Eo[0] * S_odd[i])
        trivb[i] = 10 * math.sqrt(xx) * (Ez[1] + abs(S_even[i]) * Ee[1] + abs(S_odd[i]) * Eo[1])
    lam_hat = (Px - Etriv - zpart - corr - c * sq * C0 * S) / (N * k0)
    # error budget
    own = np.bincount(kk, minlength=N); cnt = own + own[kbar]
    vals, mult = np.unique(cnt[1:], return_counts=True)
    Jsum = sum(J_L(T, int(v_), q, c, beta, m) * int(mu) for v_, mu in zip(vals, mult))
    tailL = c * sq * C0 * Jsum
    tailZ = 2 * sq * tail_bound(1000.0, T, len(zz), c * 1000.0, beta, m) / 1000.0
    Iex = Ij(beta, m, m); H = H_platt(q); j = m
    offL = (q - 2) * math.exp(c / 2) * Iex * c ** (1 - m) * x * (j / math.pi) * H ** (1 - j) * (math.log(q * H) / (j - 1) + 1 / (j - 1) ** 2)
    HZ = 3e12; offZ = 2 * math.exp(c / 2) * Iex * c ** (1 - m) * x * (j / (2 * math.pi)) * HZ ** (1 - j) * (math.log(HZ) / (j - 1) + 1 / (j - 1) ** 2)
    Wabs = float(np.sum(np.abs(w)))
    prec = (eps / sq) * (lnx * c * x * C0 * Wabs + 2 * len(g) * c * c * x * khat0) * 1.01
    rnd = c * sq * C0 * Wabs * (T * lnx * 4.5e-16 + 3e-15) + 1e-14 * Px
    delta = (tailL + tailZ + offL + offZ + prec + rnd + trivb) / (N * k0)
    err = np.abs(lam_hat - lam_true)
    certified = delta < 0.25 * lnx
    correct = (lam_hat > 0.75 * lnx) == truth
    bound_ok = bool(np.all(err <= delta))
    zeros_used = int(cnt[1:].sum()) + 2 * len(zz)
    res = dict(q=q, c=c, beta=beta, T=T, Tc=T * c, zeros=zeros_used, n=len(ns), certified=int(certified.sum()), correct=int(correct.sum()),
               max_err_over_quarter=float(np.max(err / (0.25 * lnx))), max_delta_over_quarter=float(np.max(delta / (0.25 * lnx))), bound_holds=bound_ok,
               first_wrong=[int(v) for v in ns[~correct][:5]], secs=round(time.time() - t0))
    results.append(res)
    print(json.dumps(res), flush=True)
json.dump(results, open(out, 'w'))
