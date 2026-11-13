"""The localization law's prediction tested: the kernel that minimizes the out-of-band energy at fixed f(0) (random-phase
optimum) against the Kaiser-Bessel window, as a local prime test over the family mod q. Real decisions only (the optimal
kernel's spectrum decays too slowly for the proved bound). f on [-1, 1] (normalized u = v/a, a = c n):
  (q-1) Lambda(n) f(0) = P(n) - E(n) - Z_zeta(n) - corr_q(n) - c sqrt(n) sum_{chi != chi0, rho} chibar(n) n^{i gamma} fhat(c gamma)
with fhat(sigma) = int f(u) e^{i sigma u} du, P = a int f(u) e^{c u/2} du.
Usage: famcrit_opt.py zeros.npz X c s1,s2,... out"""
import sys, math, json, time
import numpy as np
import gmpy2
zf, X = sys.argv[1], int(sys.argv[2]); c = float(sys.argv[3]); svals = [float(v) for v in sys.argv[4].split(',')]; out = sys.argv[5]
D = np.load(zf); q = int(D['q']); Tl = float(D['T']); pw = D['pw']; kz = D['k'].astype(np.int64); gz = D['gamma']
N = q - 1; ind = np.empty(q, dtype=np.int64); ind[pw] = np.arange(N)
ns = np.array([n for n in range(3, X) if n % q], dtype=np.int64)
import os
if os.environ.get('NSAMPLE'):
    ns = np.sort(np.random.default_rng(int(os.environ.get('SEED', '1'))).choice(ns, int(os.environ['NSAMPLE']), replace=False))
x = ns.astype(float); lnx = np.log(x); sq = np.sqrt(x)
idn = ind[ns % q]; r = ns % q
S_even = np.where((r == 1) | (r == q - 1), N / 2, 0.0) - 1.0; S_odd = np.where(r == 1, N / 2, 0.0) - np.where(r == q - 1, N / 2, 0.0)
truth = np.array([gmpy2.is_prime(int(v)) for v in ns]); lam_true = np.where(truth, lnx, 0.0)
for p in range(2, int(math.isqrt(X)) + 1):
    if gmpy2.is_prime(p):
        v = p * p
        while v < X:
            i_ = np.searchsorted(ns, v)
            if i_ < len(ns) and ns[i_] == v: lam_true[i_] = math.log(p)
            v *= p
zeta = np.array([float(l) for l in open(os.environ.get('ZETAZEROS', os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../data/riemann/zeros_hp2419.txt'))) if l.strip()])
hu = 2 / 600; u = np.arange(-1 + hu / 2, 1, hu); Dm = u[:, None] - u[None, :]
def optimal(s):
    with np.errstate(invalid='ignore', divide='ignore'):
        B = np.where(Dm == 0, 2 * s, 2 * np.sin(s * Dm) / np.where(Dm == 0, 1, Dm))
    M = hu * (2 * math.pi * np.eye(len(u)) - hu * B); e0 = np.zeros(len(u)); i0 = len(u) // 2; e0[i0 - 1] = e0[i0] = 0.5
    f = np.linalg.solve(M, e0); return f / (e0 @ f)
def matched(s):
    # maximize (in-band signal)^2 / (out-of-band energy): f = M^-1 a, a_i = h sin(s u_i)/(pi u_i)
    with np.errstate(invalid='ignore', divide='ignore'):
        B = np.where(Dm == 0, 2 * s, 2 * np.sin(s * Dm) / np.where(Dm == 0, 1, Dm))
    M = hu * (2 * math.pi * np.eye(len(u)) - hu * B)
    a = hu * np.sin(s * u) / (math.pi * u)
    f = np.linalg.solve(M, a); return f / (a @ f)            # normalized: in-band signal S = 1
def inband(f, s): return float((hu * np.sin(s * u) / (math.pi * u)) @ f)
from scipy.special import iv
def kb(s, beta=None, m=4):
    beta = s / 1.2 if beta is None else beta
    t = np.clip(1 - u ** 2, 0, None); f = t ** (m / 2) * iv(m, beta * np.sqrt(t)); return f / iv(m, beta)
res = []
for s in svals:
    T = s / c; assert T <= Tl + 1e-9
    sel = gz <= T; g = gz[sel]; kk = kz[sel]; zz = zeta[zeta <= T]
    for name, f in (('matched filter', matched(s)), ('min-energy', optimal(s)), ('Kaiser-Bessel', kb(s))):
        S_in = inband(f, s)
        fhat = lambda sig: np.concatenate([hu * (np.cos(np.outer(sig[i:i + 20000], u)) @ f) for i in range(0, max(len(sig), 1), 20000)]) if len(sig) else np.zeros(0)   # f even: fhat real
        w = 2 * fhat(c * g)
        Ssum = np.zeros(len(ns))
        for i in range(0, len(ns), 250):
            for jz in range(0, len(g), 40000):
                ph = np.outer(lnx[i:i + 250], g[jz:jz + 40000]) - (2 * math.pi / N) * (np.outer(idn[i:i + 250], kk[jz:jz + 40000]) % N)
                Ssum[i:i + 250] += np.cos(ph) @ w[jz:jz + 40000]
        zpart = 2 * c * sq * (np.cos(np.outer(lnx, zz)) @ fhat(c * zz)) if len(zz) else 0.0
        a = c * x
        P = a * hu * np.sum(f[None, :] * np.exp(c * u[None, :] / 2), axis=1)
        tt = x[:, None] * np.exp(c * u[None, :])
        Ez = (sq / x) * a * hu * np.sum(f * np.sqrt(tt) / (tt ** 3 - tt), axis=1)
        Ee = (sq / x) * a * hu * np.sum(f * tt ** 1.5 / (tt ** 2 - 1), axis=1)
        Eo = (sq / x) * a * hu * np.sum(f * tt ** 0.5 / (tt ** 2 - 1), axis=1)
        corr = np.zeros(len(ns)); e = 1
        while q ** e < X * math.exp(c) + 1:
            mm = float(q ** e); uu = np.log(mm / x) / c; inside = np.abs(uu) < 1
            if inside.any(): corr[inside] += math.log(q) * np.sqrt(x[inside] / mm) * np.interp(uu[inside], u, f)
            e += 1
        lam_hat = (P - Ez - Ee * S_even - Eo * S_odd - zpart - corr - c * sq * Ssum) / (N * S_in)   # divide by the in-band part of the signal
        err = np.abs(lam_hat - lam_true)
        correct = (lam_hat > 0.75 * lnx) == truth
        rr = dict(q=q, s=s, T=T, kernel=name, inband=round(S_in, 4), zeros=int(2 * sel.sum()), correct=int(correct.sum()), n=len(ns),
                  max_err_over_margin=float(np.max(err / (0.25 * lnx))), rms_err_over_margin=float(np.sqrt(np.mean((err / (0.25 * lnx)) ** 2))))
        res.append(rr); print(json.dumps(rr), flush=True)
json.dump(res, open(out, 'w'))
