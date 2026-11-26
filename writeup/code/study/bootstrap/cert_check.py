"""Independent check of the certificate family ghat(t) = Xi(t)^2 B(t), B an even polynomial.

The certificates (rational coefficients of B) were produced by other code (theory/convex/certs/, Bessel closed form,
ball arithmetic).  Here nothing of that code is used: everything is computed from the zero side,

    g(u)  = (1/pi) int_0^oo Xi(t)^2 B(t) cos(tu) dt,
    W(g)  = 2 Xi(i/2)^2 B(i/2) + (1/pi) int_0^oo Xi(t)^2 B(t) Omega(t) dt,   Xi(i/2) = 1/2,
    Omega = Re psi(1/4 + it/2) - log pi,

and from the prime side W(g) = 2 sum_n Lambda(n) n^(-1/2) g(log n) (valid because ghat vanishes at every zero).

Parts:
 1. cost W(g)/g(0) of every class-A certificate (deg 4..48) and of the class-E optima, zero side and prime side;
 2. B has no real root (own root finder);
 3. g > 0 on the grid u = k/100, 0 <= k <= 380, for the degree-8, 24 and 48 certificates (sampled, not a proof;
    beyond u = 3.8 the values are below the working precision of 520 bits);
 4. the member J = 10 of the exact family of Johns, Part I (Proposition 4.5: P(0) = 1, deg 2J in t^2, double zeros of g
    at the log of the first J prime powers), rebuilt here and compared with Table 1 there.

Quadrature: Gauss-Legendre, 120 nodes on each panel [k, k+1], 0 <= k < 400, 520 bits.
Run: (ulimit -v 1500000; timeout 1500 nice python cert_check.py [cache.pkl])   about 10 minutes without cache.
"""
import sys, os, json, glob, ast, pickle
from fractions import Fraction
import numpy as np
import gmpy2
from gmpy2 import mpfr
import mpmath as mp

PREC = 520
mp.mp.prec = PREC
gmpy2.get_context().precision = PREC
T, NODES = 400, 120
HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.join(HERE, 'theory', 'convex', 'certs')
CACHE = sys.argv[1] if len(sys.argv) > 1 else None


def gl_nodes(n):
    xs, ws = [], []
    for i in range(1, n + 1):
        x = mp.cos(mp.pi * (i - mp.mpf(1) / 4) / (n + mp.mpf(1) / 2))
        for _ in range(60):
            p0, p1 = mp.mpf(1), x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < mp.mpf(2) ** (-PREC + 20):
                break
        p0, p1 = mp.mpf(1), x
        for k in range(2, n + 1):
            p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
        dp = n * (x * p1 - p0) / (x * x - 1)
        xs.append(x)
        ws.append(2 / ((1 - x * x) * dp * dp))
    return xs, ws


def Xi(t):
    s = mp.mpf(1) / 2 + 1j * t
    return (s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)).real


def Omega(t):
    return mp.digamma(mp.mpf(1) / 4 + 1j * t / 2).real - mp.log(mp.pi)


def G(x):
    return mpfr(mp.nstr(x, 160, strip_zeros=False))


def build():
    if CACHE and os.path.exists(CACHE):
        with open(CACHE, 'rb') as f:
            t, w, om = pickle.load(f)
        return (np.array([mpfr(x) for x in t], dtype=object), np.array([mpfr(x) for x in w], dtype=object),
                np.array([mpfr(x) for x in om], dtype=object))
    xs, ws = gl_nodes(NODES)
    t, w, om = [], [], []
    for k in range(T):
        for x, wt in zip(xs, ws):
            tt = k + (x + 1) / 2
            t.append(G(tt))
            w.append(G(wt / 2 * Xi(tt) ** 2 / mp.pi))     # weight * Xi^2 / pi
            om.append(G(Omega(tt)))
    if CACHE:
        with open(CACHE, 'wb') as f:
            pickle.dump(([str(x) for x in t], [str(x) for x in w], [str(x) for x in om]), f)
    return np.array(t, dtype=object), np.array(w, dtype=object), np.array(om, dtype=object)


t, w, om = build()
t2 = t * t
ZERO = mpfr(0)


def Bval(b):                       # B(t) = sum b_j t^(2j) on the nodes
    v = np.full(len(t), ZERO, dtype=object)
    for c in reversed(b):
        v = v * t2 + c
    return v


def tot(a):
    return gmpy2.fsum(list(a))


def zero_side(b):
    wb = w * Bval(b)
    g0 = tot(wb)
    arch = tot(wb * om)
    pole = sum(c * mpfr(-1) ** j / mpfr(4) ** j for j, c in enumerate(b)) / 2      # 2 Xi(i/2)^2 B(i/2)
    return g0, pole, arch, wb


def gat(wb, u):
    u = mpfr(u)
    return tot(wb * np.array([gmpy2.cos(x * u) for x in t], dtype=object))


def mangoldt(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0:
                m //= p
            return gmpy2.log(mpfr(p)) if m == 1 else ZERO
    return ZERO


def prime_side(wb, nmax=64):
    s = ZERO
    for n in range(2, nmax + 1):
        L = mangoldt(n)
        if L != 0:
            s += 2 * L / gmpy2.sqrt(mpfr(n)) * gat(wb, gmpy2.log(mpfr(n)))
    return s


def real_roots(b):                 # number of roots s >= 0 of q(s) = sum b_j s^j
    mp.mp.prec = 2000
    try:
        r = mp.polyroots([mp.mpf(str(c)) for c in reversed(b)], maxsteps=2000, extraprec=4000)
    finally:
        mp.mp.prec = PREC
    return sum(1 for z in r if abs(z.imag) < abs(z) * mp.mpf(10) ** (-30) and z.real >= 0), \
        min(abs(mp.arg(z)) for z in r)


def e(x, d=12):
    return mp.nstr(mp.mpf(str(x)), d)


print('quadrature: %d nodes, T = %d, %d bits' % (len(t), T, PREC))
g0, pole, arch, wb1 = zero_side([mpfr(1)])
print('B = 1: g(0) = %s  pole = %s  arch = %s  W/g(0) = %s' % (e(g0, 15), e(pole), e(arch, 15), e((pole + arch) / g0, 20)))
ps = prime_side(wb1)
print('B = 1: prime side W/g(0) = %s  (prime - zero)/W = %s' % (e(ps / g0, 20), e((ps - pole - arch) / ps, 3)))

print('\nclass A certificates (rational B from theory/convex/certs/certificate_A_d*.json, key b_rational, which includes the cushion eps in b_0)')
print('deg | W/g(0) zero side         | W/g(0) prime side        | rel. diff | certified ball centre     | roots s>=0 | min |arg root|')
keep = {}
for d in range(2, 25, 2):
    c = json.load(open(os.path.join(CERT, 'certificate_A_d%d.json' % d)))
    b = [mpfr(Fraction(x).numerator) / mpfr(Fraction(x).denominator) for x in c['b_rational']]
    g0, pole, arch, wb = zero_side(b)
    z = (pole + arch) / g0
    p = prime_side(wb) / g0
    nr, ma = real_roots(b)
    ball = c['cost_over_g0'].split(' ')[0].lstrip('[')
    print('%3d | %s | %s | %s | %s | %d | %s' % (2 * d, e(z, 20), e(p, 20), e(abs(z - p) / p, 3), ball[:24], nr, mp.nstr(ma, 4)))
    keep[2 * d] = (b, wb, g0)

print('\nclass E optima (floating B from theory/convex/certs/polish_E_d*.json), cost only')
print('deg | W/g(0) zero side         | W/g(0) prime side        | rel. diff | value in the file')
for d in range(2, 17, 2):          # degree 36 (d = 18) is left out: its stored coefficients are not an optimum
    c = json.load(open(os.path.join(CERT, 'polish_E_d%d.json' % d)))
    b = [mpfr(x) for x in (ast.literal_eval(c['b']) if isinstance(c['b'], str) else c['b'])]
    g0, pole, arch, wb = zero_side(b)
    z = (pole + arch) / g0
    p = prime_side(wb) / g0
    print('%3d | %s | %s | %s | %s' % (2 * d, e(z, 20), e(p, 20), e(abs(z - p) / abs(p), 3), c['cost'][:22]))

print('\nsampled positivity of g (class A certificates)')
c1 = np.array([gmpy2.cos(x / 100) for x in t], dtype=object)
for deg in (8, 24, 48):
    b, wb, g0 = keep[deg]
    cp, cc = np.array([mpfr(1)] * len(t), dtype=object), c1.copy()
    vals = [tot(wb)]
    for k in range(1, 601):
        vals.append(tot(wb * cc))
        cp, cc = cc, 2 * c1 * cc - cp
    chk = gat(wb, mpfr(3))                       # recurrence against a direct evaluation at u = 3
    KM = 380                                     # beyond u = 3.8 (e^u = 44.7) g/g(0) is below the working precision
    neg = sum(1 for v in vals[:KM + 1] if v <= 0)
    kmin = min(range(KM + 1), key=lambda k: vals[k])
    floor = max(abs(v) for v in vals[KM + 1:]) / g0
    print('deg %2d: non-positive values on u = k/100, k <= %d: %d; smallest g/g(0) = %s at u = %.2f; g(3)/g(0) = %s (recurrence) %s (direct); largest |g|/g(0) for 3.8 < u <= 6: %s'
          % (deg, KM, neg, e(vals[kmin] / g0, 4), kmin / 100, e(vals[300] / g0, 12), e(chk / g0, 12), e(floor, 3)))

print('\nexact family of Johns, Part I, member J = 10 (double zeros at log n, n = 2 3 4 5 7 8 9 11 13 16)')
PP = [2, 3, 4, 5, 7, 8, 9, 11, 13, 16]
J = len(PP)
TAU = mpfr(30)
s2 = t2 / (TAU * TAU)
pw = [np.array([mpfr(1)] * len(t), dtype=object)]
for k in range(2 * J):
    pw.append(pw[-1] * s2)
rows, rhs = [], []
for n in PP:
    u = gmpy2.log(mpfr(n))
    cu = np.array([gmpy2.cos(x * u) for x in t], dtype=object)
    su = np.array([gmpy2.sin(x * u) for x in t], dtype=object) * t
    wc, wsn = w * cu, w * su
    mk = [tot(wc * pw[k]) for k in range(2 * J + 1)]
    dk = [-tot(wsn * pw[k]) for k in range(2 * J + 1)]
    rows.append(mk[1:]); rhs.append(-mk[0])
    rows.append(dk[1:]); rhs.append(-dk[0])
A = mp.matrix([[mp.mpf(str(x)) for x in r] for r in rows])
for i in range(A.rows):                                    # row scaling
    m = max(abs(A[i, j]) for j in range(A.cols))
    for j in range(A.cols):
        A[i, j] /= m
    rhs[i] = mp.mpf(str(rhs[i])) / m
mp.mp.prec = 2 * PREC
sol = mp.lu_solve(A, mp.matrix(rhs))
mp.mp.prec = PREC
b = [mpfr(1)] + [G(sol[k]) / TAU ** (2 * (k + 1)) for k in range(2 * J)]
g0, pole, arch, wb = zero_side(b)
Wz = pole + arch
Wp = prime_side(wb, 80)
resid = max(abs(gat(wb, gmpy2.log(mpfr(n))) / g0) for n in PP)
print('integral of Xi^2 P over R = 2 pi g(0) = %s   (Johns: between 2.7884318 and 2.7884322)' % e(2 * mp.pi * mp.mpf(str(g0)), 10))
print('W zero side = %s   W prime side (n >= 17) = %s   rel. diff %s   (Johns Table 1, A(F_10) rounded up: 4.2125564e-19)'
      % (e(Wz, 10), e(Wp, 10), e(abs(Wz - Wp) / Wp, 3)))
print('W/(2 pi g(0)) = %s   (Johns kappa_10 rounded up: 1.5107258e-19);  W/g(0) = %s' % (e(Wp / g0 / (2 * gmpy2.const_pi()), 10), e(Wp / g0, 10)))
print('largest |g(log n)|/g(0) at the ten nodes = %s;  P has %d roots with t^2 >= 0' % (e(resid, 3), real_roots(b)[0]))
