"""The family criterion: each n decided alone from the low zeros of all Dirichlet L-functions mod a prime q.

For 3 <= n < X with q prime, q > n(e^c - 1): take x = n, a = c n, K the Kaiser-Bessel window of order m. The support
n e^(-c) < t < n e^c lies in (1, n + q), so n is the only number = n (mod q) in it, and orthogonality gives exactly
  Lambda(n) K(0) = (1/(q-1)) { P(n) - E_zeta(n) + E_even(n) - Z_zeta(n) - S(n) },
  S(n) = sum over chi != chi0 of chibar(n) sum_{rho_chi} n^(i gamma) Khat(gamma/n)/sqrt(n),
with Z_zeta the zeta zero sum (no zeta zero lies below T here). The zeros with |gamma| <= T are summed; the rest is
bounded (counting bound of Bennett-Martin-O'Bryant-Rechnitzer; GRH to height 1e8/q by Platt; beyond that four
integrations by parts). n is prime iff Lambda(n) = ln n; every other n has Lambda(n) <= ln(n)/2, so the decision
'prime iff Lambda_hat > (3/4) ln n' is proved when delta(n) < (1/4) ln n.
Usage: famcrit.py zeros.npz X c beta m cc eps out.txt
"""
import sys, math, time
import numpy as np
from scipy.integrate import quad
from kbm import K, K0, Z as Zkb, P
from cert import tail_bound, offline_bound_j
from fambudget import J_L, offline_L, H_platt

zf, X = sys.argv[1], int(sys.argv[2]); c, beta, m, cc, eps = map(float, sys.argv[3:8]); m = int(m); out = sys.argv[8]
D = np.load(zf); q = int(D['q']); T = float(D['T']); pw = D['pw']
assert abs(T - cc * beta / c) < 1e-9, (T, cc * beta / c)
assert all(n * math.exp(c) < n + q for n in (3, X - 1)) and 3 * math.exp(-c) > 1
ind = np.empty(q, dtype=np.int64); ind[pw] = np.arange(q - 1)
kz = D['k'].astype(np.int64); gz = D['gamma']
counts = D['count'][1:]
k0 = K0(beta, m)
C0 = math.sqrt(2 * math.pi) * beta ** m
wz = 2 * Zkb(c * gz, beta, m)                        # each zero with 0 < gamma <= T and its mirror (chibar, -gamma)
print(f'q = {q}, T = {T:.6f}, {len(gz)} zeros with 0 < gamma <= T ({int(counts.sum())} with |gamma| <= T); zeta zeros below T: 0 (first at 14.1347)', flush=True)
assert T < 14.134725

# sums over the zeros for every n (chunked)
ns = np.arange(3, X)
S = np.zeros(len(ns)); Sabs = float(np.sum(np.abs(wz))); Sprec = float(np.sum(np.abs(wz)))
t0 = time.time()
for i in range(0, len(ns), 200):
    nn = ns[i:i + 200]; ln = np.log(nn.astype(float)); idn = ind[nn]
    acc = np.zeros(len(nn))
    for j in range(0, len(gz), 40000):
        g = gz[j:j + 40000]; kk = kz[j:j + 40000]
        ph = np.outer(ln, g) - (2 * math.pi / (q - 1)) * (np.outer(idn, kk) % (q - 1))
        acc += np.cos(ph) @ wz[j:j + 40000]
    S[i:i + 200] = acc
print(f'zero sums for {len(ns)} numbers: {time.time() - t0:.0f} s', flush=True)

# per-character tail bound, averaged over the characters (each has its own count N(T, chi))
vals, mult = np.unique(counts, return_counts=True)
Jmean = sum(J_L(T, int(v), q, c, beta, m) * int(mu) for v, mu in zip(vals, mult)) / (q - 1)

res = []; worst = 0.0; wrong = []; undecided = []
from sympy import isprime
for idx, n in enumerate(ns):
    x = float(n); a = c * x; sq = math.sqrt(x)
    def fz(v):
        t = x * math.exp(v / x); return float(K(v, a, beta, m)) * math.sqrt(t) / (x * (t ** 3 - t))
    def fe(v):
        t = x * math.exp(v / x); return float(K(v, a, beta, m)) * t ** 1.5 / (x * (t * t - 1))
    Ez, ez = quad(fz, -a, a, limit=400, epsabs=0, epsrel=1e-12); Ee, ee = quad(fe, -a, a, limit=400, epsabs=0, epsrel=1e-12)
    Ez *= sq; Ee *= sq; ez *= sq; ee *= sq
    Px = P(x, a, beta, m)
    SL = c * sq * C0 * S[idx]
    lam_hat = (Px - Ez + Ee - SL) / ((q - 1) * k0)
    # error budget (all in absolute units, then / ((q-1) K0) where it enters the class average)
    tailL = c * sq * C0 * Jmean * (q - 1)                       # sum over chi != chi0 of (1/sqrt x) tail
    tailZ = (2 / sq) * tail_bound(x, T, 0, a, beta, m)          # zeta: all its zeros lie above T
    offL = offline_L(x, a, beta, m, 4, q) * (q - 2)
    offZ = offline_bound_j(x, a, beta, m, 4)
    khat0 = a * C0 * float(Zkb(0.0, beta, m))
    prec = (1 / sq) * eps * (math.log(x) * c * x * C0 * Sprec + 2 * len(gz) * c * khat0) * 1.01
    rnd = c * sq * C0 * Sabs * ((T * math.log(x) + 2 * math.pi) * 4.5e-16 + 3e-15 + len(gz) * 1.2e-16) + 1e-14 * (Px + abs(Ez) + abs(Ee))
    triv = 10 * (ez + ee)
    delta = (tailL + tailZ + offL + offZ + prec + rnd + triv) / ((q - 1) * k0)
    L = math.log(x); truth = isprime(int(n))
    if delta < 0.25 * L:
        dec = lam_hat > 0.75 * L
        if dec != truth: wrong.append(int(n))
    else:
        undecided.append(int(n))
    lam_true = 0.0
    if truth: lam_true = L
    else:
        from sympy import factorint
        f = factorint(int(n)); lam_true = math.log(next(iter(f))) if len(f) == 1 else 0.0
    err = abs(lam_hat - lam_true)
    assert err <= delta, (n, lam_hat, lam_true, delta)
    worst = max(worst, err / delta)
    res.append((int(n), lam_hat, lam_true, delta, tailL / ((q - 1) * k0), offL / ((q - 1) * k0), prec / ((q - 1) * k0)))
R = np.array(res)
primes = [int(r[0]) for r in res if r[1] > 0.75 * math.log(r[0]) and r[3] < 0.25 * math.log(r[0])]
lines = [f'family criterion: q = {q}, c = {c}, beta = {beta}, m = {m}, T = {T:.6f} (no zeta zero), eps = {eps}',
         f'zeros used: {int(counts.sum())} (|gamma| <= T, {q - 2} L-functions); numbers decided: {len(ns) - len(undecided)} of {len(ns)} (3 <= n < {X}); undecided: {undecided[:10]}',
         f'primes found: {len(primes)}; equal to sympy: {primes == [p for p in range(3, X) if isprime(p)]}; wrong: {wrong[:10]}',
         f'largest delta/K0 = {R[:, 3].max():.4f} at n = {int(R[np.argmax(R[:, 3]), 0])} (threshold (1/4) ln n = {0.25 * math.log(R[np.argmax(R[:, 3]), 0]):.3f}); '
         f'parts there: tail {R[np.argmax(R[:, 3]), 4]:.3e}, above Platt height {R[np.argmax(R[:, 3]), 5]:.3e}, precision {R[np.argmax(R[:, 3]), 6]:.3e}',
         f'largest |Lambda_hat - Lambda| = {np.abs(R[:, 1] - R[:, 2]).max():.3e}; largest ratio to delta = {worst:.3e}',
         'samples: ' + '; '.join(f'n={int(r[0])}: Lambda_hat={r[1]:.6f} Lambda={r[2]:.6f} delta={r[3]:.4f}' for r in R[[0, 4, 1006, 7916, -1]])]
print('\n'.join(lines))
open(out, 'w').write('\n'.join(lines) + '\n')
np.save(out.replace('.txt', '.npy'), R)
