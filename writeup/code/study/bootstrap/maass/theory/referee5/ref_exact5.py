"""Referee check of exact5.md (even character chi_5, Gamma_R(s), conductor 5).

Own code, no import of the author's scripts. mpmath at 50 digits.
Route: g_B(u) = (1/2pi) int Xi(t)^2 B(t) e^{-itu} dt by the trapezoid rule on the
whole line (step h = 2pi/20).  Exact aliasing identity:
    (h/2pi) sum_j ghat(jh) e^{-ijhu} = sum_k g_B(u + 2pi k/h),
so the error for |u| <= 7 is g_B at distance >= 13, i.e. below exp(-(2pi/5) e^13).
Xi(t) = Lambda(1/2+it), Lambda(s) = (5/pi)^{s/2} Gamma(s/2) L(s, chi_5), L from
mpmath.dirichlet (cross-checked against arb's acb.dirichlet_l through python-flint).
Run: (ulimit -v 1200000; timeout 600 nice .venv/bin/python -I ref_exact5.py)
"""
import json, sys, time
from fractions import Fraction
import mpmath as mp

mp.mp.dps = 50
HERE = '/home/bas/praime/writeup/code/study/bootstrap/maass/theory/'
NCHK = [0, 0]
LOG = []


def out(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


def check(name, ok, info=''):
    NCHK[0] += 1
    if not ok:
        NCHK[1] += 1
    out(('PASS ' if ok else 'FAIL ') + name + ('  ' + info if info else ''))


def nstr(x, n=12):
    return mp.nstr(x, n)


CHI5 = [0, 1, -1, -1, 1]
CHI8 = [0, 1, 0, -1, 0, -1, 0, 1]


def chi5(n):
    return CHI5[n % 5]


def Lam(s, chi, q):
    s = mp.mpc(s)
    return (mp.mpf(q) / mp.pi) ** (s / 2) * mp.gamma(s / 2) * mp.dirichlet(s, chi)


t0 = time.time()
# ---------------------------------------------------------------- (A) theta, Gauss sum, Mellin
out('== (A) theta function of the even character')


def theta0(x, N=260):
    return mp.fsum(chi5(n) * mp.exp(-mp.pi * n * n * x / 5) for n in range(1, N))


for x in ['0.07', '0.25', '0.769230769230769230769', '2.3']:
    x = mp.mpf(x)
    lhs, rhs = theta0(1 / x), mp.sqrt(x) * theta0(x)
    check('theta0(1/x) = x^(1/2) theta0(x) at x=%s' % nstr(x, 6), abs(lhs - rhs) < mp.mpf(10) ** -40 * abs(rhs),
          'rel %s' % nstr(abs(lhs - rhs) / abs(rhs), 3))
x = mp.mpf('0.25')
wrong = mp.power(x, 1.5) * theta0(x)
check('negative: x^(3/2) law (odd type) fails for chi_5', abs(theta0(1 / x) - wrong) > mp.mpf('0.01') * abs(wrong))
tau = mp.fsum(chi5(a) * mp.exp(2j * mp.pi * a / 5) for a in range(1, 5))
check('Gauss sum tau(chi_5) = sqrt5 (root number +1)', abs(tau - mp.sqrt(5)) < mp.mpf(10) ** -45, nstr(tau, 20))
check('chi_5 even: chi(-1)=chi(4)=+1, and sum_a chi(a)=0 (no n=0 term, no pole)',
      chi5(4) == 1 and sum(CHI5) == 0)
# Lambda real on the critical line (root number +1)
for t in ['3.3', '47.1']:
    v = Lam(mp.mpf(1) / 2 + 1j * mp.mpf(t), CHI5, 5)
    check('Im Lambda(1/2+it) = 0 at t=%s' % t, abs(v.imag) < mp.mpf(10) ** -40 * (abs(v.real) + mp.mpf(10) ** -30),
          'Re %s' % nstr(v.real, 8))
# Mellin integral at three points, incl. Re s < 0 (entire), with direct theta (no use of C1)


def mellin(s):
    f = lambda x: theta0(x) * x ** (s / 2 - 1)
    return mp.quad(f, [mp.mpf('0.005'), mp.mpf('0.05'), mp.mpf('0.3'), 1, 3, 10, 40, 120])


mp.mp.dps = 35
for s in [mp.mpc(2, 1), mp.mpc('0.5', 3), mp.mpc('-1.5', '0.7')]:
    a, b = mellin(s), Lam(s, CHI5, 5)
    check('Mellin int theta0 x^(s/2) dx/x = Lambda(s) at s=%s' % nstr(s, 4), abs(a - b) < mp.mpf(10) ** -22 * abs(b),
          'rel %s' % nstr(abs(a - b) / abs(b), 3))
mp.mp.dps = 50
lam0 = Lam(mp.mpf(10) ** -15, CHI5, 5)
check('Lambda(1e-15) finite, = Lambda(1) = 2 log golden ratio to O(1e-15) (pole of Gamma_R cancelled by L(0)=0)',
      abs(lam0.real - 2 * mp.log((1 + mp.sqrt(5)) / 2)) < mp.mpf(10) ** -13, nstr(abs(lam0.real - 2 * mp.log((1 + mp.sqrt(5)) / 2)), 3))
check('L(0, chi_5) = 0', abs(mp.dirichlet(0, CHI5)) < mp.mpf(10) ** -40)

# L cross-check with arb (independent implementation)
try:
    import flint
    flint.ctx.dps = 50
    c5 = flint.dirichlet_char(5, 4)
    for t in ['0.7', '37.3', '99.1', '149.9']:
        s = mp.mpf(1) / 2 + 1j * mp.mpf(t)
        a = mp.dirichlet(s, CHI5)
        b = flint.acb.dirichlet_l(flint.acb(flint.arb('0.5'), flint.arb(t)), c5)
        b = mp.mpc(mp.mpf(b.real.mid().str(50, radius=False)), mp.mpf(b.imag.mid().str(50, radius=False)))
        check('mpmath L(1/2+it, chi5) = arb L at t=%s' % t, abs(a - b) < mp.mpf(10) ** -40 * abs(a), nstr(abs(a - b), 3))
except Exception as e:
    out('flint not usable:', e)

# ---------------------------------------------------------------- trapezoid nodes
out('== nodes of the trapezoid rule')
H = 2 * mp.pi / 20
TMAX = 160
M = int(TMAX / H) + 1
TJ = [j * H for j in range(M + 1)]
X5 = [Lam(mp.mpf(1) / 2 + 1j * t, CHI5, 5).real ** 2 for t in TJ]
X8 = [Lam(mp.mpf(1) / 2 + 1j * t, CHI8, 8).real ** 2 for t in TJ]
out('M =', M, 'Xi5(T)^2 =', nstr(X5[-1], 3), 'time', round(time.time() - t0, 1))
W = [H / (2 * mp.pi) * (1 if j == 0 else 2) for j in range(M + 1)]


def Bval(b, t):
    t2 = t * t
    return mp.fsum(bk * t2 ** k for k, bk in enumerate(b))


COS = {}


def cosrow(u):
    key = mp.nstr(u, 45)
    if key not in COS:
        COS[key] = [mp.cos(t * u) for t in TJ]
    return COS[key]


def gfun(Xw, u):
    # Xw[j] = W_j Xi^2 B (t_j)
    return mp.fsum(a * c for a, c in zip(Xw, cosrow(u)))


def weights(X, b, extra=None):
    return [W[j] * X[j] * Bval(b, TJ[j]) * (1 if extra is None else extra(TJ[j])) for j in range(M + 1)]


def g0(Xw):
    return mp.fsum(Xw)


# GL nodes for the u-integral of A_0
GLd = mp.calculus.quadrature.GaussLegendre(mp.mp)
GLN = GLd.calc_nodes(5, mp.mp.prec)  # 48 nodes on [-1,1]


def glpanel(a, b):
    a, b = mp.mpf(a), mp.mpf(b)
    return [((b - a) / 2 * x + (b + a) / 2, (b - a) / 2 * w) for x, w in GLN]


UPAN = [(0, '0.25'), ('0.25', '0.5'), ('0.5', 1), (1, '1.5'), ('1.5', 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7)]
UNODES = [nw for a, b in UPAN for nw in glpanel(a, b)]
UMAX = mp.mpf(7)


def A_kappa(Xw, kappa):
    """A_kappa(g) = -(gamma+log pi) g(0) + int_0^inf [2e^{-2u} g(0) - 2 e^{-(1/2+kappa)u} g(u)]/(1-e^{-2u}) du."""
    G0 = g0(Xw)
    s = mp.fsum(w * (2 * mp.exp(-2 * u) * G0 - 2 * mp.exp(-(mp.mpf(1) / 2 + kappa) * u) * gfun(Xw, u)) / (1 - mp.exp(-2 * u))
                for u, w in UNODES)
    s += -G0 * mp.log(1 - mp.exp(-2 * UMAX))  # int_7^inf 2e^{-2u}/(1-e^{-2u}) du, g(u>7) negligible
    return -(mp.euler + mp.log(mp.pi)) * G0 + s


def vonmangoldt(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0:
                m //= p
            return mp.log(p) if m == 1 else 0
    return 0


NMAX = 1000
LAMB = {n: vonmangoldt(n) for n in range(2, NMAX + 1)}
PP = [n for n in range(2, NMAX + 1) if LAMB[n] != 0]


def primeside(Xw, chi):
    return 2 * mp.fsum(chi(n) * LAMB[n] / mp.sqrt(n) * gfun(Xw, mp.log(n)) for n in PP if chi(n) != 0)


def residual(Xw, kappa=0, chi=chi5, q=5):
    G0 = g0(Xw)
    return mp.log(q) * G0 + A_kappa(Xw, kappa) - primeside(Xw, chi), G0


# aliasing / truncation control: g(0) and g(log 2) with step 2pi/16 from a subset? use the identity sum g(u+2pi k/h)
out('== (B) Lemma 1 and the archimedean term')
B1 = [mp.mpf(1)]
Xw1 = weights(X5, B1)
G01 = g0(Xw1)
out('g(0) [B=1] =', nstr(G01, 30))
# step 2h (every second node): aliasing then at distance 10; compare g(0)
G0half = mp.fsum((2 * H / (2 * mp.pi)) * (1 if j == 0 else 2) * X5[j] for j in range(0, M + 1, 2))
check('trapezoid step h vs 2h agree for g(0) (aliasing g(10) tiny)', abs(G0half - G01) < mp.mpf(10) ** -40,
      nstr(abs(G0half - G01), 3))
check('truncation: Xi5(160)^2 B < 1e-90', X5[-1] < mp.mpf(10) ** -90)

R1, _ = residual(Xw1)
check('Lemma 1, B=1: (log5) g(0) + A_0(g) - 2 sum chi Lambda n^-1/2 g(log n) = 0',
      abs(R1) < mp.mpf(10) ** -35 * G01, 'residual %s, g(0)=%s' % (nstr(R1, 3), nstr(G01, 10)))

# Fourier-side A_0 for B=1: (1/2pi) int Xi^2 Omega, Omega = -log pi + Re psi(1/4+it/2); GL in t, own L evaluations
mp.mp.dps = 40
TPAN = [(0, '0.25'), ('0.25', '0.5'), ('0.5', 1)] + [(k, k + 1) for k in range(1, 20)] + [(k, k + 2) for k in range(20, 120, 2)]
GLN3 = GLd.calc_nodes(4, mp.mp.prec)  # 24 nodes


def glp3(a, b):
    a, b = mp.mpf(a), mp.mpf(b)
    return [((b - a) / 2 * x + (b + a) / 2, (b - a) / 2 * w) for x, w in GLN3]


TN = [nw for a, b in TPAN for nw in glp3(a, b)]
XI2T = [(t, w, Lam(mp.mpf(1) / 2 + 1j * t, CHI5, 5).real ** 2) for t, w in TN]
OM = [-mp.log(mp.pi) + mp.digamma(mp.mpf(1) / 4 + 0.5j * t).real for t, w, _ in XI2T]
A0F = mp.fsum(w * x * o for (t, w, x), o in zip(XI2T, OM)) / mp.pi
G0F = mp.fsum(w * x for (t, w, x) in XI2T) / mp.pi
mp.mp.dps = 50
A0U = A_kappa(Xw1, 0)
check('A_0 u-form (constants -(gamma+log pi), 2e^{-2u}g(0), kernel e^{-u/2}) = Fourier side (1/2pi) int ghat Omega',
      abs(A0U - A0F) < mp.mpf(10) ** -28 * abs(A0U), 'A0 = %s, diff %s' % (nstr(A0U, 20), nstr(abs(A0U - A0F), 3)))
check('g(0) Gauss-Legendre = g(0) trapezoid', abs(G0F - G01) < mp.mpf(10) ** -30, nstr(abs(G0F - G01), 3))

# Lemma 1 for other even B (any sign) -- B not >= 0 allowed in Lemma 1
for name, b in [('t^2', [0, 1]), ('1 - t^2/10', [1, mp.mpf(-1) / 10]), ('(t^2-44.2)^2', [mp.mpf('44.2') ** 2, -2 * mp.mpf('44.2'), 1]),
                ('t^6 - 3 t^2', [0, -3, 0, 1])]:
    Xw = weights(X5, b)
    R, G0 = residual(Xw)
    scale = mp.fsum(abs(a) for a in Xw)
    check('Lemma 1, B=%s' % name, abs(R) < mp.mpf(10) ** -35 * scale, 'residual %s, scale %s' % (nstr(R, 3), nstr(scale, 5)))

# negative control 1: wrong archimedean kernel (odd, kappa = 1) on the even kernel
Rk1, _ = residual(Xw1, kappa=1)
check('negative control: kernel e^{-3u/2} (kappa=1) instead of e^{-u/2} gives a nonzero residual', abs(Rk1) > mp.mpf('1e-3') * G01,
      nstr(Rk1, 10))

# negative control 2: ghat = Xi_{chi8}^2 inserted in chi_5's formula -> sum over zeros of L(s,chi5) of Xi8(gamma)^2
Xw8 = weights(X8, B1)
R8, G08 = residual(Xw8)
out('control chi8 kernel in chi5 formula: residual =', nstr(R8, 25), ' g8(0) =', nstr(G08, 10))
check('negative control: kernel Xi_{chi8}^2 does NOT give zero in chi_5 formula', abs(R8) > mp.mpf('1e-6') * G08)
# also: chi8 kernel in its own formula (log 8, chi8) is zero
R88, _ = residual(Xw8, chi=lambda n: CHI8[n % 8], q=8)
check('chi8 kernel in its own formula (log 8, chi8 primes) gives 0', abs(R88) < mp.mpf(10) ** -35 * G08, nstr(R88, 3))

# zeros of L(s, chi5) up to T=60 by sign changes, counted by the argument principle
mp.mp.dps = 25
TZ = 62
step = mp.mpf('0.02')
ts = [step * k for k in range(1, int(TZ / step) + 1)]
vals = [Lam(mp.mpf(1) / 2 + 1j * t, CHI5, 5).real for t in ts]
brackets = [(ts[k], ts[k + 1]) for k in range(len(ts) - 1) if vals[k] * vals[k + 1] < 0]
mp.mp.dps = 45
zeros = []
for a, b in brackets:
    z = mp.findroot(lambda t: Lam(mp.mpf(1) / 2 + 1j * t, CHI5, 5).real, (a, b), solver='anderson')
    zeros.append(z)
out('zeros of L(s,chi5) with 0<gamma<=62 found by sign change:', len(zeros), ' first', nstr(zeros[0], 15))


def Nplus(T):
    # (1/pi)[theta(T) + Delta arg L along 2 -> 2+iT -> 1/2+iT]
    th = T / 2 * mp.log(mp.mpf(5) / mp.pi) + mp.loggamma(mp.mpf(1) / 4 + 0.5j * T).imag
    # arg L on Re s = 2 is principal (|L-1| < 1); then track horizontally
    a = mp.arg(mp.dirichlet(2 + 1j * T, CHI5))
    prev = mp.dirichlet(2 + 1j * T, CHI5)
    nst = 600
    for k in range(1, nst + 1):
        s = 2 - mp.mpf(3) / 2 * k / nst + 1j * T
        cur = mp.dirichlet(s, CHI5)
        a += mp.arg(cur / prev)
        prev = cur
    return (th + a) / mp.pi


mp.mp.dps = 25
for T in [30, 62]:
    nT = Nplus(mp.mpf(T))
    cnt = sum(1 for z in zeros if z <= T)
    check('zero count up to %d: sign changes %d = argument principle %s' % (T, cnt, nstr(nT, 8)), abs(nT - cnt) < mp.mpf('1e-6'))
mp.mp.dps = 50
out('first zeros of L(s,chi5):', [nstr(zz, 12) for zz in zeros[:3]])
zs = 2 * mp.fsum(Lam(mp.mpf(1) / 2 + 1j * z, CHI8, 8).real ** 2 for z in zeros)
tailbound = Lam(mp.mpf(1) / 2 + 62j, CHI8, 8).real ** 2
check('control residual = sum over zeros of L(chi5) of Xi_8(gamma)^2 (both signs of gamma)',
      abs(zs - R8) < mp.mpf(10) ** -25, 'sum %s, diff %s' % (nstr(zs, 20), nstr(abs(zs - R8), 3)))


# B = 1 + t^2/7: the note's g(0) = 0.875 and Lemma 1
Xw7 = weights(X5, [1, mp.mpf(1) / 7])
R7, G07 = residual(Xw7)
check('B = 1 + t^2/7: g(0) = 0.87517073609894 (note: 0.875) and Lemma 1', abs(G07 - mp.mpf('0.87517073609894305955')) < mp.mpf(10) ** -19
      and abs(R7) < mp.mpf(10) ** -35, 'g(0) %s, residual %s' % (nstr(G07, 20), nstr(R7, 3)))

# (D)/(E) with a compactly supported test function: f(u) = (1-(u/a)^2)^k on [-a,a], true family chi_5 at q = 5.
# Q(f) computed in u-space (exact polynomial g = f*f~) must equal the zero side sum_rho |fhat(gamma)|^2 >= 0.
import sympy as sp
a_s, k_s = sp.Rational(6, 5), 8
v_s, w_s = sp.symbols('v w', real=True)
fpol = (1 - (v_s / a_s) ** 2) ** k_s
gpol = sp.expand(sp.integrate(sp.expand(fpol.subs(v_s, v_s + w_s) * fpol), (v_s, -a_s, a_s - w_s)))  # valid for 0 <= w <= 2a
gco = sp.Poly(gpol, w_s).all_coeffs()
gco = [mp.mpf(sp.Rational(c).p) / sp.Rational(c).q for c in gco]
gc = lambda u: mp.polyval(gco, u) if 0 <= u <= 2 * mp.mpf(6) / 5 else mp.mpf(0)
A = mp.mpf(6) / 5
G0c = gc(mp.mpf(0))
Ain = mp.fsum(w * (2 * mp.exp(-2 * u) * G0c - 2 * mp.exp(-u / 2) * gc(u)) / (1 - mp.exp(-2 * u))
              for p0, p1 in [(0, '0.3'), ('0.3', '0.8'), ('0.8', '1.6'), ('1.6', '2.4')] for u, w in glpanel(p0, p1))
A0c = -(mp.euler + mp.log(mp.pi)) * G0c + Ain - G0c * mp.log(1 - mp.exp(-4 * A))
Pc = 2 * mp.fsum(chi5(n) * LAMB[n] / mp.sqrt(n) * gc(mp.log(n)) for n in PP if n <= 12 and chi5(n) != 0)
Qc = mp.log(5) * G0c + A0c - Pc
fhat = lambda t: 2 * mp.quad(lambda u: (1 - (u / A) ** 2) ** 8 * mp.cos(t * u), [0, A / 2, A])
Zc = 2 * mp.fsum(fhat(z) ** 2 for z in zeros)
check('compact f=(1-(u/a)^2)^8, a=1.2: Q_{chi5,5}(f) = sum over zeros |fhat|^2 > 0 (sign and normalisation of Q)',
      Qc > 0 and abs(Qc - Zc) < mp.mpf(10) ** -12 * G0c, 'Q %s, zero side (gamma<62) %s, g(0) %s' % (nstr(Qc, 15), nstr(Zc, 15), nstr(G0c, 6)))
Qc49 = mp.log(mp.mpf('4.9')) * G0c + A0c - Pc
out('   same f at q = 4.9: Q =', nstr(Qc49, 10), '(4.9 lies below the B=1 bound 4.9058; this f is no witness, g_B approximants are)')


# Vacuity caveat: the trivial family mu_p = delta_0 (every p) is feasible at q >= exp(-Omega(0)), unconditionally,
# since A_0(g) = (1/2pi) int ghat Omega >= Omega(0) g(0) for ghat >= 0 (Re psi(x+iy) >= psi(x), x > 0).
Om0 = -mp.log(mp.pi) + mp.digamma(mp.mpf(1) / 4)
check('Omega(0) = -log pi + psi(1/4) = -log pi - gamma - pi/2 - 3 log 2', abs(Om0 - (-mp.log(mp.pi) - mp.euler - mp.pi / 2 - 3 * mp.log(2))) < mp.mpf(10) ** -45, nstr(Om0, 15))
Omg = [-mp.log(mp.pi) + mp.digamma(mp.mpf(1) / 4 + 0.5j * mp.mpf(k) / 20).real for k in range(0, 4001)]
check('Omega(t) increasing on a grid of [0, 200] (so Omega >= Omega(0))', all(Omg[k + 1] > Omg[k] for k in range(4000)))
qtriv = mp.exp(-Om0)
out('   trivial family delta_0 feasible at q >= exp(-Omega(0)) =', nstr(qtriv, 10), '; odd case exp(log pi - psi(3/4)) =',
    nstr(mp.exp(mp.log(mp.pi) - mp.digamma(mp.mpf(3) / 4)), 10))
check('A_0(g)/g(0) >= Omega(0) for B=1 and for the compact f', A0U / G01 >= Om0 and A0c / G0c >= Om0, '%s %s' % (nstr(A0U / G01, 8), nstr(A0c / G0c, 8)))
check('Q_{delta_0, q=215.4}(f) >= 0 for the compact f (delta_0 family: no prime sum)', mp.log(mp.mpf('215.4')) * G0c + A0c >= 0, nstr(mp.log(mp.mpf('215.4')) * G0c + A0c, 8))

# ---------------------------------------------------------------- (C) closed form and derivative recurrence
out('== (C) closed form of phi*phi')
mp.mp.dps = 50


def d(N):
    return sum(1 for k in range(1, N + 1) if N % k == 0)


DN = [0] + [d(N) for N in range(1, 400)]


def closed(w, j=0):
    # (d/dw)^j (phi*phi)(w) = 4 e^{w/2} sum chi(N) d(N) [p_j K0 + q_j K1](z), z = 2 pi N e^w / 5
    tot = 0
    for N in range(1, 400):
        if chi5(N) == 0:
            continue
        z = 2 * mp.pi * N * mp.exp(w) / 5
        if z > 400:
            break
        k0, k1 = mp.besselk(0, z), mp.besselk(1, z)
        if j == 0:
            v = k0
        elif j == 2:
            v = (mp.mpf(1) / 4 + z * z) * k0 - z * k1  # p_2 = 1/4 + z^2, q_2 = -z from the recurrence
        tot += chi5(N) * DN[N] * v
    return 4 * mp.exp(w / 2) * tot


for w in ['0.37', '-0.8', str(mp.log(7)), str(mp.log(13))]:
    w = mp.mpf(w)
    a, b = closed(w), gfun(Xw1, w)
    check('phi*phi closed form (4 e^{w/2} sum chi d(N) K0) = Fourier inversion at w=%s' % nstr(w, 5),
          abs(a - b) < mp.mpf(10) ** -38 * G01, 'diff %s, value %s' % (nstr(abs(a - b), 3), nstr(b, 8)))
Xw_t2 = weights(X5, [0, -1])  # FT of -t^2 Xi^2 = (phi*phi)''
for w in ['0.37', str(mp.log(3))]:
    w = mp.mpf(w)
    a, b = closed(w, 2), gfun(Xw_t2, w)
    check('(phi*phi)\'\' via p_2 = 1/4+z^2, q_2 = -z = Fourier of -t^2 Xi^2 at w=%s' % nstr(w, 5),
          abs(a - b) < mp.mpf(10) ** -36 * abs(b), 'diff %s' % nstr(abs(a - b), 3))
# recurrence itself: symbolic, against direct differentiation
import sympy as sp
z = sp.symbols('z', positive=True)
p, q = sp.Integer(1), sp.Integer(0)
f = sp.sqrt(z) * sp.besselk(0, z)
for j in range(1, 5):
    p, q = sp.expand(sp.Rational(1, 2) * p + z * sp.diff(p, z) - z * q), sp.expand(-sp.Rational(1, 2) * q + z * sp.diff(q, z) - z * p)
    f = z * sp.diff(f, z)
    diffexpr = f - sp.sqrt(z) * (p * sp.besselk(0, z) + q * sp.besselk(1, z))
    val = max(abs(complex(diffexpr.subs(z, zz).evalf(30))) for zz in [sp.Rational(3, 10), sp.Rational(17, 7), 9])
    check('recurrence D^%d[z^(1/2) K0] = z^(1/2)(p K0 + q K1)' % j, val < 1e-20, 'p=%s q=%s' % (p, q) if j <= 2 else '')

# ---------------------------------------------------------------- (D,E) defects and bounds for the certificates
out('== (D,E) Defect and the bound')
TABLE = {0: ('9.512e-3', '9.513e-3', '9.423e-2', 5), 2: ('1.593e-3', '1.594e-3', '1.591e-2', 7), 4: ('7.348e-5', '7.349e-5', '7.349e-4', 11),
         6: ('8.372e-6', '8.373e-6', '8.373e-5', 13), 8: ('2.206e-7', '2.207e-7', '2.207e-6', 19), 10: ('8.227e-9', '8.228e-9', '8.228e-8', 23),
         12: ('1.147e-10', '1.148e-10', '1.148e-9', 29), 14: ('1.918e-13', '1.919e-13', '1.919e-12', 37),
         16: ('2.763e-15', '2.764e-15', '2.764e-14', 41), 18: ('2.763e-15', '2.764e-15', '2.764e-14', 41)}
PRIMES = [p for p in range(2, NMAX + 1) if all(p % k for k in range(2, int(p ** 0.5) + 1))]
summary = {}
for deg in sorted(TABLE):
    c = json.load(open(HERE + 'cert5_d%d.json' % deg))
    b = [mp.mpf(Fraction(x).numerator) / Fraction(x).denominator for x in c['b_exact']]
    b[0] += mp.mpf(Fraction(c['eps_added_to_b0']).numerator) / Fraction(c['eps_added_to_b0']).denominator
    Xw = weights(X5, b)
    G0 = g0(Xw)
    shares = {}
    for p in PRIMES:
        coef = []
        m = 1
        while p ** m <= NMAX:
            coef.append(mp.log(p) * mp.power(p, -mp.mpf(m) / 2) * gfun(Xw, m * mp.log(p)))
            m += 1
        P = lambda x: mp.fsum(cm * x ** (k + 1) for k, cm in enumerate(coef))
        cand = [mp.mpf(-1), mp.mpf(1)]
        big = max(abs(cm) for cm in coef)
        ct = list(coef)
        while len(ct) > 1 and abs(ct[-1]) < mp.mpf(10) ** -70 * big:
            ct.pop()
        if len(ct) >= 2:
            dcoef = [(k + 1) * cm for k, cm in enumerate(ct)]  # P'(x) = sum (k+1) c_k x^k
            if len(dcoef) == 2:
                rts = [-dcoef[0] / dcoef[1]]
            else:
                rts = mp.polyroots(list(reversed(dcoef)), maxsteps=400, extraprec=400)
            cand += [mp.re(r) for r in rts if abs(mp.im(r)) < mp.mpf(10) ** -30 and abs(mp.re(r)) <= 1]
            grid = min(P(mp.mpf(k) / 200 - 1) for k in range(401))
        else:
            grid = None
        mn = min(P(x) for x in cand)
        if grid is not None:
            assert grid >= mn - mp.mpf(10) ** -45 * big, (p, grid, mn)
        shares[p] = P(mp.mpf(chi5(p))) - mn
    D = mp.fsum(shares.values())
    r = D / G0
    bound = 5 * mp.exp(-2 * r)
    lo, hi, gap, carrier = TABLE[deg]
    top = sorted(shares, key=lambda p: -shares[p])
    summary[deg] = (r, 5 - bound, top[:3])
    out('deg %2d: g(0)=%s  Defect/g(0)=%s  5-bound=%s  shares: %s' % (deg, nstr(G0, 12), nstr(r, 10), nstr(5 - bound, 10),
        ', '.join('%d:%s' % (p, nstr(shares[p] / G0, 3)) for p in top[:3])))
    check('deg %d: Defect/g(0) in table interval [%s, %s]' % (deg, lo, hi), mp.mpf(lo) <= r <= mp.mpf(hi), nstr(r, 10))
    check('deg %d: 5 - bound <= table value %s (table is a valid rounding up)' % (deg, gap), 5 - bound <= mp.mpf(gap)
          and 5 - bound > mp.mpf(gap) * mp.mpf('0.999'), nstr(5 - bound, 10))
    if deg > 0:
        check('deg %d: largest share outside p=5 at p=%d' % (deg, carrier), [p for p in top if p != 5][0] == carrier,
              str(top[:3]))
    else:
        check('deg 0: largest share p=5 then 11', top[0] == 5 and top[1] == 11, str(top[:3]))
    check('deg %d: Defect_2, Defect_3 negligible (< 1e-40 g(0))' % deg, shares[2] < mp.mpf(10) ** -40 * G0 and shares[3] < mp.mpf(10) ** -40 * G0,
          '%s %s' % (nstr(shares[2] / G0, 3), nstr(shares[3] / G0, 3)))
    g2, g3 = gfun(Xw, mp.log(2)), gfun(Xw, mp.log(3))
    check('deg %d: g(log 2), g(log 3) > 0' % deg, g2 > 0 and g3 > 0, '%s %s' % (nstr(g2 / G0, 5), nstr(g3 / G0, 5)))
    Dtop = len(b) - 1
    carr = [p for p in top if p != 5][0]
    out('   parity: carrier %d, chi5 = %d, (-1)^D = %d' % (carr, chi5(carr), (-1) ** Dtop))
    check('deg %d: carrier class chi5(p) = (-1)^D holds except at deg 18' % deg, (chi5(carr) == (-1) ** Dtop) == (deg != 18))
    if deg in (0, 8):
        Rz, _ = residual(Xw)
        check('deg %d: Lemma 1 residual' % deg, abs(Rz) < mp.mpf(10) ** -35 * G0, nstr(Rz, 3))
        check('deg %d: B > 0 on a grid t in [0,200]' % deg, min(Bval(b, mp.mpf(k) / 10) for k in range(2001)) > 0)
B1q = 5 * mp.exp(-2 * summary[0][0])
check('B=1: bound in [4.905776, 4.905777) (note: q >= 4.905776, 4.9057)', mp.mpf('4.905776') <= B1q < mp.mpf('4.905777'), nstr(B1q, 10))

out('time', round(time.time() - t0, 1))
out('%d checks, %d failed' % (NCHK[0], NCHK[1]))
open(HERE + 'referee5/ref_exact5_out.txt', 'w').write('\n'.join(LOG) + '\n')
