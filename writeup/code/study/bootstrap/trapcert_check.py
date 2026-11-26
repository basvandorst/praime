"""Spot check, with separate code, of the twelve certified traps of degree 48 (theory/zeros/traps/trapcert_k*.json).

The certificates were produced with ball arithmetic and the Bessel closed form.  Here every function is evaluated from
the zero side (Gauss-Legendre quadrature of Xi(t)^2 B(t) cos(tu), 120 nodes per unit panel, 520 bits), as in
cert_check.py and trap_check.py.  This is a numerical check of what the certificates state, not a second proof.

For each node k the certificate gives: gt (rational coefficients, file theory/convex/certs/certificate_A_d24.json),
h (rational coefficients h_rational; the certificate uses h - DELTA * Phi*Phi), eta' (rational), and an interval
(a_k, b_k).  Its statement: w = eta' (h - DELTA g_0) - gt satisfies  B_w <= 0 on R,  W(w) > 0,  and  w <= 0 on
(0, oo) outside (a_k, b_k).  Checked here:
  1. W(h - DELTA g_0) from our prime sum against the ball W_h of the certificate;
  2. eta' W(h - DELTA g_0) - W(gt) > 0, against the certificate's margin;
  3. B_gt - eta' B_(h - DELTA g_0) has no root with t^2 >= 0 (own root finder);
  4. w(a_k) <= 0, w(b_k) <= 0, w(log n_k) > 0;
  5. w <= 0 at every point of the grid u = i/200, 0 < u <= 3.8, outside (a_k, b_k)   (a sample, not a proof).

Run: (ulimit -v 1500000; timeout 3000 nice python trapcert_check.py [cache.pkl])
"""
import sys, os, json, ast, pickle
from fractions import Fraction
import numpy as np
import gmpy2
from gmpy2 import mpfr
import mpmath as mp

PREC = 520
mp.mp.prec = PREC
gmpy2.get_context().precision = PREC
T, NODES, DH = 400, 120, 24
HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.join(HERE, 'theory', 'convex', 'certs')
TRAP = os.path.join(HERE, 'theory', 'zeros', 'traps')
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


def G(x):
    return mpfr(mp.nstr(x, 160, strip_zeros=False))


def build():
    if CACHE and os.path.exists(CACHE):
        with open(CACHE, 'rb') as f:
            t, w = pickle.load(f)
        return np.array([mpfr(x) for x in t], dtype=object), np.array([mpfr(x) for x in w], dtype=object)
    xs, ws = gl_nodes(NODES)
    t, w = [], []
    for k in range(T):
        for x, wt in zip(xs, ws):
            tt = k + (x + 1) / 2
            t.append(G(tt))
            w.append(G(wt / 2 * Xi(tt) ** 2 / mp.pi))
    if CACHE:
        with open(CACHE, 'wb') as f:
            pickle.dump(([str(x) for x in t], [str(x) for x in w]), f)
    return np.array(t, dtype=object), np.array(w, dtype=object)


t, w = build()
TAU = mpfr(30)
s2 = t * t / (TAU * TAU)
pw = [np.array([mpfr(1)] * len(t), dtype=object)]
for j in range(DH):
    pw.append(pw[-1] * s2)
wpw = [w * p for p in pw]
ZERO = mpfr(0)


def tot(a):
    return gmpy2.fsum(list(a))


def basis_at(u):
    u = mpfr(u)
    c = np.array([gmpy2.cos(x * u) for x in t], dtype=object)
    return [tot(wp * c) for wp in wpw]


def val(coef, bas):
    return sum(c * b for c, b in zip(coef, bas))


def mangoldt(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0:
                m //= p
            return gmpy2.log(mpfr(p)) if m == 1 else ZERO
    return ZERO


def fr(x):
    f = Fraction(x)
    return mpfr(f.numerator) / mpfr(f.denominator)


def e(x, d=6):
    return mp.nstr(mp.mpf(str(x)), d)


def ball(s):
    return mp.mpf(s.split(' ')[0].lstrip('['))


cg = json.load(open(os.path.join(CERT, 'certificate_A_d24.json')))
bg = [fr(x) * TAU ** (2 * j) for j, x in enumerate(cg['b_rational'])]
Wb = [ZERO] * (DH + 1)
for n in range(2, 48):
    L = mangoldt(n)
    if L != 0:
        bas = basis_at(gmpy2.log(mpfr(n)))
        f = 2 * L / gmpy2.sqrt(mpfr(n))
        Wb = [a + f * b for a, b in zip(Wb, bas)]
B0 = basis_at(0)
g0 = val(bg, B0)
kt = val(bg, Wb)
print('gt = cushioned certificate of degree 48: gt(0) = %s, W(gt) = %s' % (e(g0, 20), e(kt, 20)))
GRID = [mpfr(i) / 200 for i in range(1, 761)]
BG = [basis_at(u) for u in GRID]
gG = [val(bg, b) for b in BG]
print(' k |  n | W(h) ours - ball | eta\' W(h) - W(gt) | margin in file | roots s >= 0 | w(a_k)     | w(b_k)     | w(log n)   | grid points outside (a,b) with w > 0 | largest w/gt there | b_k - log n | a_k - log n')
for k in range(1, 13):
    c = json.load(open(os.path.join(TRAP, 'trapcert_k%d.json' % k)))
    hr = c['h_rational'] if not isinstance(c['h_rational'], str) else ast.literal_eval(c['h_rational'])
    ch = [fr(x) * TAU ** (2 * j) for j, x in enumerate(hr)]
    ch += [ZERO] * (DH + 1 - len(ch))
    delta = fr(c['DELTA'])
    eta = fr(c['eta_prime'])
    a, b = fr(c['a_k']), fr(c['b_k'])
    n = int(c['n_k'])
    Wraw = val(ch, Wb)
    Wg0 = Wb[0]
    # the certificate's h is h_rational - DELTA g_0, unless h_rational already contains it
    Wfile = ball(c['W_h'])
    d1 = abs(mp.mpf(str(Wraw - delta * Wg0)) - Wfile)
    d2 = abs(mp.mpf(str(Wraw)) - Wfile)
    if d1 <= d2:
        ch[0] = ch[0] - delta
        dW = d1
    else:
        dW = d2
    Wh = val(ch, Wb)
    cw = [eta * x - y for x, y in zip(ch, bg)]                 # w = eta' h - gt
    mp.mp.prec = 2500
    try:
        poly = [mp.mpf(str(-x)) for x in cw]                   # B_gt - eta' B_h, in s = t^2/TAU^2
        while abs(poly[-1]) == 0:
            poly.pop()
        r = mp.polyroots(list(reversed(poly)), maxsteps=4000, extraprec=8000)
        nreal = sum(1 for z in r if abs(z.imag) < abs(z) * mp.mpf(10) ** (-30) and z.real >= 0)
    finally:
        mp.mp.prec = PREC
    wa, wb_, wn = val(cw, basis_at(a)), val(cw, basis_at(b)), val(cw, basis_at(gmpy2.log(mpfr(n))))
    out = [(val(cw, bas), g) for u, bas, g in zip(GRID, BG, gG) if not (a < u < b) and u <= mpfr('3.8')]
    bad = sum(1 for v, g in out if v > 0)
    worst = max(v / g for v, g in out)
    ln = gmpy2.log(mpfr(n))
    print('%2d | %2d | %s | %s | %s | %d | %s | %s | %s | %d | %s | %s | %s'
          % (k, n, e(dW, 2), e(eta * Wh - kt, 8), c['ii_margin'].split(' ')[0].lstrip('[')[:14], nreal, e(wa, 3), e(wb_, 3), e(wn, 3), bad, e(worst, 4),
             e(b - ln, 6), e(a - ln, 6)), flush=True)
