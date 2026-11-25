"""Independent checks of deg2lib.py (separate formulas, no splines, different quadrature).
(a) O_p(x): power series of deg2lib against the closed form by partial fractions, and against the literal four-fold sum
    over (m, n, h, k) with a(p^m), mu(p^h) from the Euler factor (definition of the ratios recipe);
(b) Dcorr_p: closed form against direct sums over k of (ln p)^2 [W(k) - k(alpha^k + beta^k)^2] p^{-k(1+x)};
(c) PARI's L(s, Sym^2 E) at analytic s = 2 against the Euler product over p < 2*10^5 with the local factors of deg2lib;
    PARI's (L'/L)'(1 + i e) against a finite difference of PARI's log L;
(d) R(0) = 0 and the GUE limit: C(r) -> 0 as r -> 0;
(e) Sigma2 of one curve at one L recomputed with no spline (all functions evaluated at every r on a trapezoid grid of
    step 0.002, closed-form O_p) against deg2lib's Gauss-Legendre/spline value; and deg2lib with pmax = 10^5 against 3*10^5.
usage: validate.py"""
import sys, os, math, time, pickle
import numpy as np, mpmath as mp, cypari2
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deg2lib as D
pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6)
H = os.path.dirname(os.path.abspath(__file__)); out = []
def rep(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); out.append(s)


def O_closed(p, lam, x):
    th = math.acos(lam / 2); a = complex(math.cos(th), math.sin(th)); b = a.conjugate()
    y = np.exp(x * math.log(p)); Z = np.exp(-(1 + x) * math.log(p))
    A1 = (1 - 1 / y) * (a * y - b) / (y * (a - b)); A2 = (1 - 1 / y) * (b * y - a) / (y * (b - a))
    u = a * a * y * y * Z; v = y * y * Z; w = b * b * y * y * Z
    return 1 + A1 ** 2 * u / (1 - u) + 2 * A1 * A2 * v / (1 - v) + A2 ** 2 * w / (1 - w)


def O_literal(p, lam, x, J=40):
    a = [1.0, lam]
    for m in range(2, J + 1): a.append(lam * a[-1] - a[-2])
    mu = [1.0, -lam, 1.0] + [0.0] * (J - 2)
    tot = 0
    for m in range(J + 1):
        for h in range(3):
            j = m + h
            if j > J: continue
            for k in range(3):
                n = j - k
                if n < 0: continue
                tot = tot + a[m] * a[n] * mu[h] * mu[k] * p ** (-j) * np.exp(x * (m - k) * math.log(p))
    return tot


rep('(a) O_p: series vs closed form vs literal four-fold sum')
x = 1j * np.array([0.3, 1.7, 4.4, 7.9]); mx = 0
for p, ap in ((2, 0), (2, 1), (2, -2), (3, 3), (5, -4), (7, 1), (101, 13), (997, -40)):
    lam = ap / math.sqrt(p)
    th = math.acos(lam / 2)
    # deg2lib's series is inside local_good: recover O_p from its log ratio by adding back the denominators
    d, la, _, _, _ = D.local_good(p, lam, x)
    a2 = np.exp(2j * th); b2 = np.conj(a2); lp = math.log(p); zp = np.exp(-(1 + x) * lp); zm = np.exp(-(1 - x) * lp); z0 = 1 / p
    lnden = -(np.log(1 - zp) + np.log(1 - zm)) - (np.log(1 - a2 * zp) + np.log(1 - zp) + np.log(1 - b2 * zp)) \
            - (np.log(1 - a2 * zm) + np.log(1 - zm) + np.log(1 - b2 * zm)) + 2 * (np.log(1 - a2 * z0) + np.log(1 - z0) + np.log(1 - b2 * z0)).real
    Oser = np.exp(la + lnden + 2 * math.log1p(-1 / p))
    Oc = O_closed(p, lam, x); Ol = O_literal(p, lam, x, J=60 if p < 10 else 12)
    e1 = np.abs(Oser - Oc).max(); e2 = np.abs(Ol - Oc).max(); mx = max(mx, e1, e2)
    rep('  p=%d a_p=%d: |series - closed| %.1e, |literal - closed| %.1e' % (p, ap, e1, e2))
# multiplicative: closed form of deg2lib vs literal (a(p^m) = (e/sqrt p)^m, mu(p) = -e/sqrt p)
for p, e in ((2, 1), (3, -1), (11, 1)):
    al = e / math.sqrt(p); J = 50; a = [al ** m for m in range(J + 1)]; mu = [1.0, -al]
    tot = 0
    for m in range(J + 1):
        for h in range(2):
            j = m + h
            for k in range(2):
                n = j - k
                if 0 <= n <= J: tot = tot + a[m] * a[n] * mu[h] * mu[k] * p ** (-j) * np.exp(x * (m - k) * math.log(p))
    d, la, _, _, _ = D.local_mult(p, x); lp = math.log(p); zp = np.exp(-(1 + x) * lp); zm = np.exp(-(1 - x) * lp)
    lnden = -(np.log(1 - zp) + np.log(1 - zm)) - (np.log(1 - zp / p) + np.log(1 - zm / p)) + 2 * math.log1p(-1 / p ** 2)
    Om = np.exp(la + lnden + 2 * math.log1p(-1 / p)); rep('  mult p=%d: |deg2lib - literal| %.1e' % (p, np.abs(Om - tot).max()))

rep('(b) Dcorr_p closed form vs direct sums over k')
for p, ap in ((2, 0), (3, -3), (5, 2), (13, -7)):
    lam = ap / math.sqrt(p); th = math.acos(lam / 2); lp = math.log(p)
    d, _, j1, j2, j3 = D.local_good(p, lam, x)
    Wk = lambda k: (2 * math.cos(k * th)) ** 2
    direct = sum(lp ** 2 * (Wk(k) - k * Wk(k)) * np.exp(-k * (1 + x) * lp) for k in range(1, 80))
    J = sum(lp ** 2 * Wk(k) * np.exp(-k * (1 + x) * lp) for k in range(1, 80))
    rep('  p=%d a_p=%d: |Dcorr - direct| %.1e, |J1+J2+J3 - J| %.1e' % (p, ap, np.abs(d - direct).max(), np.abs(j1 + j2 + j3 - J).max()))

rep('(c) PARI Sym^2 against the Euler product at analytic s = 2, and (L\'/L)\' against finite differences')
for ainv in ([1, -1, 0, -17, 29], [0, 0, 1, -7, 6], [1, 1, 1, -15, 16]):
    E = pari.ellinit(ainv); N = int(pari.ellglobalred(E)[0]); bad = [int(q) for q in pari.factor(N)[0]]
    S = pari.lfunsympow(E, 2); v = float(pari.lfun(S, 3).real())
    P = D.primes_upto(200000); an = pari.ellan(E, 200000); lnp = 0.0
    for q in P:
        q = int(q); X = q ** -2.0
        if q in bad: lnp -= math.log(1 - X / q)
        else:
            lam = int(an[q - 1]) / math.sqrt(q); a2 = lam * lam - 2       # alpha^2 + beta^2
            lnp -= math.log((1 - X) * (1 - a2 * X + X * X))
    Li = pari.lfuninit(S, [2, 0.05, 3], 2); e = 1.3; h = 1e-3
    lnL = lambda t: complex(mp.log(complex(pari.lfun(Li, pari('2+%.10f*I' % t)))))
    fd = -(lnL(e + h) - 2 * lnL(e) + lnL(e - h)) / h ** 2       # d^2/de^2 log L(1+ie) = -(L'/L)'
    s = pari('2+%.10f*I' % e); v0 = complex(pari.lfun(Li, s)); v1 = complex(pari.lfun(Li, s, 1)); v2 = complex(pari.lfun(Li, s, 2))
    dd = v2 / v0 - (v1 / v0) ** 2
    rep('  N=%d: PARI L(Sym2, 2) = %.8f, Euler product (p < 2e5) %.8f, rel diff %.1e; (L\'/L)\'(1+1.3i): PARI %s, f.d. %s' % (
        N, v, math.exp(lnp), abs(v / math.exp(lnp) - 1), np.round(dd, 6), np.round(fd, 6)))

zf = os.path.join(H, 'zgrid.pkl'); Z = pickle.load(open(zf, 'rb'))
rep('(d) R(0) = 0: C(r) for r -> 0 (L = 12)')
c = D.Curve(pari, [0, 0, 1, -7, 6])
for r in (1e-4, 1e-3, 1e-2):
    rep('  N=%d r=%.0e C=%.2e' % (c.N, r, float(D.C_full(Z, c, 12.0, r=np.array([r]))[0])))

rep('(e) Sigma2 at L = 12.5, l = 1, 4, 8, 10: deg2lib (splines, Gauss-Legendre) vs direct (no spline, trapezoid step 0.002)')
L = 12.5; ls = (1.0, 4.0, 8.0, 10.0)
g1 = D.sigma2_gue() + D.cum_sigma(D.C_full(Z, c, L)); v1 = [g1[int(round(l / D.W)) - 1] for l in ls]
r = np.arange(0, 10.0 + 1e-9, 0.002); r[0] = 1e-5; e = 2 * math.pi * r / L
mp.mp.dps = 20
Lz = np.array([complex(mp.zeta(mp.mpc(1, t))) for t in e]); Lz1 = np.array([complex(mp.zeta(mp.mpc(1, t), derivative=1)) for t in e])
Lz2 = np.array([complex(mp.zeta(mp.mpc(1, t), derivative=2)) for t in e])
E = pari.ellinit([0, 0, 1, -7, 6]); S = pari.lfunsympow(E, 2); Li = pari.lfuninit(S, [2, 0.05, 6], 2)
LS0 = complex(pari.lfun(Li, 2)); vals = []
for t in e:
    s = pari('2+%.12f*I' % t); vals.append((complex(pari.lfun(Li, s)), complex(pari.lfun(Li, s, 1)), complex(pari.lfun(Li, s, 2))))
vals = np.array(vals)
J = Lz2 / Lz - (Lz1 / Lz) ** 2 + vals[:, 2] / vals[:, 0] - (vals[:, 1] / vals[:, 0]) ** 2       # (zeta'/zeta)' + (LS'/LS)'
lnA = np.zeros(len(e), complex); xx = 1j * e
P = D.primes_upto(10 ** 5); an = pari.ellan(E, 10 ** 5)
for q in P:
    q = int(q); lp = math.log(q)
    if q == 5077:
        d, la, _, _, _ = D.local_mult(q, xx); J += d; lnA += la; continue
    lam = int(an[q - 1]) / math.sqrt(q); th = math.acos(lam / 2); a2 = np.exp(2j * th); b2 = np.conj(a2)
    z = np.exp(-(1 + xx) * lp); zm = np.exp(-(1 - xx) * lp)
    J += -(lp ** 2) * ((a2 * z) ** 2 / (1 - a2 * z) ** 2 + 2 * z ** 2 / (1 - z) ** 2 + (b2 * z) ** 2 / (1 - b2 * z) ** 2)
    Oc = O_closed(q, lam, xx)
    lnA += np.log(Oc) - 2 * math.log1p(-1 / q) + np.log(1 - z) + np.log(1 - zm) + np.log(1 - a2 * z) + np.log(1 - z) + np.log(1 - b2 * z) \
           + np.log(1 - a2 * zm) + np.log(1 - zm) + np.log(1 - b2 * zm) - 2 * (np.log(1 - a2 / q) + np.log(1 - 1 / q) + np.log(1 - b2 / q)).real
K = Lz * np.conj(Lz) * np.abs(vals[:, 0]) ** 2 / abs(LS0) ** 2 * np.exp(lnA)
Rm1 = (2 / L ** 2) * (J + np.exp(-1j * e * L) * K).real
Rm1[0] = -1.0
v2 = []
for l in ls:
    m = r <= l + 1e-12; rr = r[m]; rr[0] = 0.0
    v2.append(l + 2 * np.trapezoid((l - rr) * Rm1[m], rr))
for l, a_, b_ in zip(ls, v1, v2): rep('  l=%.0f: deg2lib %.6f direct %.6f diff %.1e' % (l, a_, b_, a_ - b_))
c3 = D.Curve(pari, [0, 0, 1, -7, 6], pmax=3 * 10 ** 5)
g3 = D.sigma2_gue() + D.cum_sigma(D.C_full(Z, c3, L))
rep('  pmax 3e5 vs 1e5: ' + ' '.join('l=%.0f %.1e' % (l, g3[int(round(l / D.W)) - 1] - g1[int(round(l / D.W)) - 1]) for l in ls))
open(os.path.join(H, 'validate_out.txt'), 'w').write('\n'.join(out) + '\n')
