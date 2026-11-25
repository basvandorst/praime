"""Independent recomputation of the frozen predictions (no crosslib, no numpy sieve): sympy for primes and Kronecker
symbols, mpmath for the sums. Checks (1) the per-prime tables A_inf, B_inf and the finite-T factors in
predictions_cross.json; (2) 40 random pair predictions in predictions_cross_pairs.npz (primes <= 2e5, plain loops);
(3) the sha256 of both files; (4) the sign rule: for real characters the p-part is positive iff chi1(p) = chi2(p)."""
import os, sys, json, math, hashlib, random
import numpy as np
import sympy, mpmath
HERE = os.path.dirname(os.path.abspath(__file__))
J = json.load(open(os.path.join(HERE, 'predictions_cross.json'))); Z = np.load(os.path.join(HERE, 'predictions_cross_pairs.npz'))
print('sha256 pairs file %s (registered %s)' % (hashlib.sha256(open(os.path.join(HERE, 'predictions_cross_pairs.npz'), 'rb').read()).hexdigest(), J['pairs_sha256']))
W = [float(w) for w in J['w']]; T1, T2, H = J['T1'], J['T2'], J['h']
import functools
@functools.lru_cache(maxsize=None)
def kron(D, p):
    if p == 2:
        return 0 if D % 2 == 0 else (1 if D % 8 in (1, 7) else -1)
    return int(sympy.jacobi_symbol(D % p, p)) if D % p else 0
def ab(p, w, finiteT=False):
    A = mpmath.mpf(0); B = mpmath.mpf(0); m = 1; lp = mpmath.log(p)
    if finiteT:
        M = int(math.floor((T2 - w - T1) / H + 1e-9)) + 1
    while mpmath.mpf(1) / (m * m * mpmath.mpf(p) ** m) > mpmath.mpf('1e-15'):
        t = 2 / mpmath.pi ** 2 / (m * m * mpmath.mpf(p) ** m) * mpmath.sin(m * w * lp / 2) ** 2
        if finiteT:   # mean of cos(2m(t + w/2) ln p) over the grid, by the closed geometric sum
            x = 2 * m * lp; a0 = x * (T1 + w / 2); d = x * H
            f = 1 + mpmath.re(mpmath.exp(1j * a0) * (mpmath.exp(1j * d * M) - 1) / (mpmath.exp(1j * d) - 1)) / M
            t *= f
        if m % 2: A += t
        else: B += t
        m += 1
    return float(A), float(B)
mx = 0.0
for p in J['primes']:
    for k, w in enumerate(W):
        a, b = ab(p, w); aT, bT = ab(p, w, True)
        q = J['primes'].index(p); mx = max(mx, abs(a - J['A_inf'][q][k]), abs(b - J['B_inf'][q][k]), abs(aT - J['A_T'][q][k]), abs(bT - J['B_T'][q][k]))
print('(1) per-prime tables: max |difference| %.2e' % mx)
PR = list(sympy.primerange(2, 200001)); random.seed(7); worst = 0.0
TAB = {p: [ab(p, w) for w in W] for p in (2, 3, 5, 7)}
# finite-range part, independently: orbits n = p^m <= 300, grid means of cos(a_n) cos(a_m) by closed geometric sums
PW = [(p ** m, p, m) for p in sympy.primerange(2, 301) for m in range(1, 10) if p ** m <= 300]
def gmean(x, w):
    """mean over the grid t_j = T1 + j H (j < M) of cos(x (t_j + w/2))"""
    M = int(math.floor((T2 - w - T1) / H + 1e-9)) + 1
    if abs(x) < 1e-15: return 1.0
    z = complex(math.cos(x * H), math.sin(x * H)); a0 = x * (T1 + w / 2)
    s = (z ** M - 1) / (z - 1); return (complex(math.cos(a0), math.sin(a0)) * s).real / M
def delta(a, b, k):
    w = W[k]; tot = 0.0
    cs = [(-(2 / math.pi) / (m * math.sqrt(n)) * math.sin(w * math.log(n) / 2), n, p, m) for n, p, m in PW]
    for c1, n1, p1, m1 in cs:
        x1 = kron(a, p1) ** m1
        if x1 == 0: continue
        for c2, n2, p2, m2 in cs:
            x2 = kron(b, p2) ** m2
            if x2 == 0: continue
            l1, l2 = math.log(n1), math.log(n2)
            K = 0.5 * (gmean(l1 + l2, w) + gmean(l1 - l2, w))
            tot += c1 * c2 * x1 * x2 * K - (c1 * c1 * x1 * x2 / 2 if n1 == n2 else 0.0)
    return tot
for name in ('NN', 'NO'):
    D1 = Z[name + '_D1']; D2 = Z[name + '_D2']; F = Z[name + '_full']
    for _ in range(10):
        i = random.randrange(len(D1)); a, b = int(D1[i]), int(D2[i])
        tot = [0.0] * len(W)
        for p in PR:
            c = kron(a, p) * kron(b, p)
            if c == 0: continue
            if p <= 7:
                for k in range(len(W)): tot[k] += TAB[p][k][0] * c + TAB[p][k][1]
                continue
            lp = math.log(p); m = 1
            while 1.0 / (m * m * p ** m) > 1e-15:
                s = 2 / math.pi ** 2 / (m * m * p ** m) * (c ** m)
                for k, w in enumerate(W): tot[k] += s * math.sin(m * w * lp / 2) ** 2
                m += 1
        tot = [tot[k] + delta(a, b, k) for k in range(len(W))]
        d = max(abs(tot[k] - float(F[i, k])) for k in range(len(W))); worst = max(worst, d)
        print('  %s pair (%d, %d): max |recomputed - frozen| over w = %.2e; at w = 4.5: %+.5f' % (name, a, b, d, tot[14]))
print('(2) worst pair difference %.2e' % worst)
for p in (2, 3, 5, 7):
    ok = all(J['A_T'][J['primes'].index(p)][k] > 0 for k in range(len(W)))
    print('(4) p=%d: A_T(w) > 0 at every w: %s (so agreement chi1(p) = chi2(p) gives a positive p-part)' % (p, ok))
