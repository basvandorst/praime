#!/usr/bin/env python3
"""check_echo.py: recompute and assert every number in echo_law.md.

Inputs: the zero table writeup/data/riemann/zeros_40000_flint.txt of the repo (read only), zeros of
L(s, chi_-4) and L(s, chi_-3) from PARI (lfunzeros), and the primes to 10^8 from the sieve below.
Nothing is imported from the repo. Run:
  (ulimit -v 1500000; timeout 600 nice ../../../../.venv/bin/python -I check_echo.py)
Every number printed with 'ok' is compared with the value written in the note (table REF below).
"""
import math, time, sys
import numpy as np, mpmath as mp, sympy as sp, cypari2
from scipy.special import sici
from scipy.optimize import minimize_scalar, brentq

T_START = time.time()
ZFILE = __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', '..', '..', 'data', 'riemann', 'zeros_40000_flint.txt')
NCHK = [0]; FAILED = []

def chk(name, ok, msg=''):
    NCHK[0] += 1
    if not ok: FAILED.append(name)
    print(('ok    ' if ok else 'FAIL  ') + name + ('   ' + msg if msg else ''))

def near(name, val, ref, tol):
    chk(name, abs(val - ref) <= tol, '%.12g  (note: %.12g, tol %.1g)' % (val, ref, tol))

def section(s):
    print('\n== %s  [%.0f s]' % (s, time.time() - T_START))

# ---------------------------------------------------------------------------------------------
section('0. constants')
mp.mp.dps = 40
EG = mp.euler; G1 = mp.stieltjes(1)
Cm = 2 + EG - mp.log(4 * mp.pi)
K2m = 6 + 2 * EG + 2 * EG**2 + 4 * G1 - mp.pi**2 / 4 - 2 * mp.log(4 * mp.pi)
C = float(Cm); K2 = float(K2m); EGf = float(EG); LN2PI = math.log(2 * math.pi)
near('C = 2 + gamma - ln 4pi', C, 0.0461914179322421, 1e-16)
near('gamma1 (first Stieltjes constant)', float(G1), -0.0728158454836767, 1e-15)
near('K2 = sum 1/|rho|^4 (closed form)', K2, 7.42012728749297e-5, 1e-18)
near('B = 1 - gamma - ln 2pi', float(1 - EG - mp.log(2 * mp.pi)), -1.41509273131, 1e-11)
zl = lambda s: mp.zeta(s, derivative=1) / mp.zeta(s)
near('zeta\'/zeta(0) = ln 2pi', float(zl(0) - mp.log(2 * mp.pi)), 0.0, 1e-30)
dzl0 = mp.diff(zl, 0)
near('(zeta\'/zeta)\'(0) = pi^2/12 - 2 gamma1 - gamma^2', float(dzl0 - (mp.pi**2 / 12 - 2 * G1 - EG**2)), 0.0, 1e-25)
near('sqrt(K2)/C (rms of the correlation coefficient)', math.sqrt(K2) / C, 0.1865, 5e-5)
near('3 - 3 K2/C^2 (kurtosis under LI)', 3 - 3 * K2 / C**2, 2.8957, 5e-5)
near('3 C^2 - 3 K2 (fourth moment under LI)', 3 * C**2 - 3 * K2, 0.0061783, 5e-8)
near('excess kurtosis -3 K2/C^2 under LI', -3 * K2 / C**2, -0.1043, 5e-5)

# Dirichlet characters, independently with mpmath (Hurwitz zeta, Stieltjes constants) and with PARI
def chi_val(q, m):
    m %= q
    if q == 4: return [0, 1, 0, -1][m]
    if q == 3: return [0, 1, -1][m]
    if q == 5: return [0, 1, -1, -1, 1][m]
def L_mp(q, s, der=0):
    f = lambda s: mp.mpf(q)**(-s) * sum(chi_val(q, r) * mp.zeta(s, mp.mpf(r) / q) for r in range(1, q))
    return mp.diff(f, s, der) if der else f(s)
def LpL1_mp(q):      # L'/L(1, chi) from the Laurent data of the Hurwitz zeta at s = 1
    g0 = [mp.stieltjes(0, mp.mpf(r) / q) for r in range(1, q)]; g1 = [mp.stieltjes(1, mp.mpf(r) / q) for r in range(1, q)]
    S0 = sum(chi_val(q, r) * g0[r - 1] for r in range(1, q)); S1 = sum(chi_val(q, r) * g1[r - 1] for r in range(1, q))
    L1 = S0 / q; dL1 = -mp.log(q) * S0 / q - S1 / q
    return dL1 / L1
pari = cypari2.Pari(); pari.allocatemem(120 * 10**6); pari.set_real_precision(38)
def LpL0_mp(q):     # L'/L(0, chi) from Lerch's formula zeta'(0, x) = ln Gamma(x) - (1/2) ln 2pi
    L0 = sum(chi_val(q, r) * (mp.mpf(1) / 2 - mp.mpf(r) / q) for r in range(1, q))
    dL0 = -mp.log(q) * L0 + sum(chi_val(q, r) * mp.loggamma(mp.mpf(r) / q) for r in range(1, q))
    return dL0 / L0
LP1, LP0 = {}, {}
for q, D in ((4, -4), (3, -3)):
    a1 = LpL1_mp(q); a0 = LpL0_mp(q)
    Lp = pari.lfuncreate(D)
    p1 = float(pari.lfun(Lp, 1, 1) / pari.lfun(Lp, 1)); p0 = float(pari.lfun(Lp, 0, 1) / pari.lfun(Lp, 0))
    chk('chi_-%d: L\'/L(1), L\'/L(0) agree mpmath vs PARI' % q, abs(float(a1) - p1) < 1e-14 and abs(float(a0) - p0) < 1e-14,
        '%.15f %.15f' % (float(a1), float(a0)))
    fe = a0 + a1 - (-mp.log(mp.mpf(q) / mp.pi) + EG + mp.log(2))
    near('chi_-%d: functional equation L\'/L(0) + L\'/L(1) = -ln(q/pi) + gamma + ln 2' % q, float(fe), 0.0, 1e-30)
    LP1[q] = float(a1); LP0[q] = float(a0)
near('chi_-4: L\'/L(1)', LP1[4], 0.2456095848, 1e-10)
near('chi_-4: L\'/L(0)', LP0[4], 0.7831887854, 1e-10)
near('chi_-3: L\'/L(1)', LP1[3], 0.3682816160, 1e-10)
near('chi_-3: L\'/L(0)', LP0[3], 0.9481988267, 1e-10)
CCHI = {q: LP1[q] - LP0[q] + math.log(2) for q in (4, 3)}
near('C_chi_-4 = L\'/L(1) - L\'/L(0) + ln 2', CCHI[4], 0.155568, 5e-7)
near('C_chi_-3', CCHI[3], 0.113230, 5e-7)
for q in (4, 3):
    near('C_chi_-%d = 2 L\'/L(1) + ln(q/pi) - gamma' % q, CCHI[q], 2 * LP1[q] + math.log(q / math.pi) - EGf, 1e-14)
# an even character, the Legendre symbol mod 5: B0 = lim (L'/L(s) - 1/s) at 0
tay = mp.taylor(lambda s: L_mp(5, s), 0, 2)
B0 = tay[2] / tay[1]
g0 = [mp.stieltjes(0, mp.mpf(r) / 5) for r in range(1, 5)]; g1 = [mp.stieltjes(1, mp.mpf(r) / 5) for r in range(1, 5)]
S0 = sum(chi_val(5, r) * g0[r - 1] for r in range(1, 5)); S1 = sum(chi_val(5, r) * g1[r - 1] for r in range(1, 5))
LP1_5 = (-mp.log(5) * S0 / 5 - S1 / 5) / (S0 / 5)
C5 = LP1_5 - B0 - mp.log(2)
near('(./5) even: B0 = b(chi)', float(B0), -0.0220247, 5e-8)
near('(./5) even: C_chi = L\'/L(1) - B0 - ln 2', float(C5), 0.156557, 5e-7)
near('(./5) even: B0 + L\'/L(1) = -ln(q/pi) + gamma + ln 2', float(B0 + LP1_5 + mp.log(5 / mp.pi) - EG - mp.log(2)), 0.0, 1e-25)
near('(./5) even: C_chi = 2 L\'/L(1) + ln(q/pi) + psi(1/2)', float(C5 - (2 * LP1_5 + mp.log(5 / mp.pi) + mp.digamma(0.5))), 0.0, 1e-25)

# ---------------------------------------------------------------------------------------------
section('1. closed forms and the tail identities')
def T_ser(a, K=400): return mp.nsum(lambda k: mp.mpf(a)**(-2 * k) / (2 * k * (2 * k + 1)), [1, mp.inf])
def To_ser(a): return mp.nsum(lambda k: mp.mpf(a)**(1 - 2 * k) / ((2 * k - 1) * 2 * k), [1, mp.inf])
def T_cl(a): a = mp.mpf(a); return 1 - mp.log(1 - a**-2) / 2 - a * mp.atanh(1 / a)
def To_cl(a): a = mp.mpf(a); return mp.atanh(1 / a) + a / 2 * mp.log(1 - a**-2)
mx = max(abs(T_cl(a) - T_ser(a)) for a in (1.3, 2, 7)); near('T(a): closed form = series (a = 1.3, 2, 7)', float(mx), 0, 1e-30)
mx = max(abs(To_cl(a) - To_ser(a)) for a in (1.3, 2, 7)); near('T_o(a): closed form = series', float(mx), 0, 1e-30)
near('T(1) = 1 - ln 2', float(T_ser(1) - (1 - mp.log(2))), 0, 1e-30)
near('T_o(1) = ln 2', float(To_ser(1) - mp.log(2)), 0, 1e-30)
mx = max(abs(-a / 2 * mp.quad(lambda t: mp.log(1 - t**-2) / t**2, [a, mp.inf]) - T_cl(a)) for a in (1.5, 3, 10))
near('-(a/2) int_a^oo ln(1 - t^-2) dt/t^2 = T(a)  (zeta: E0 tail = E tail - ln 2pi + T)', float(mx), 0, 1e-25)
mx = max(abs(a * mp.quad(lambda t: mp.atanh(1 / t) / t**2, [a, mp.inf]) - To_cl(a)) for a in (1.5, 3, 10))
near('a int_a^oo artanh(1/t) dt/t^2 = T_o(a)  (odd chi)', float(mx), 0, 1e-25)
a_, t_ = sp.symbols('a t', positive=True)
mx = max(abs(a * mp.quad(lambda t: mp.log(t) / t**2, [a, mp.inf]) - mp.log(a) - 1) for a in (1.5, 3, 10))
near('a int_a^oo ln t dt/t^2 = ln a + 1  (even chi)', float(mx), 0, 1e-25)
y_ = sp.symbols('y', positive=True)
Tsym = 1 - sp.log(1 - y_**-2) / 2 - y_ / 2 * sp.log((y_ + 1) / (y_ - 1))
chk('T\'(y) = 1/y - (1/2) ln((y+1)/(y-1))  (sympy; so dJ/da -> -oo as a -> 1+)',
    abs(float((sp.diff(Tsym, y_) - (1 / y_ - sp.log((y_ + 1) / (y_ - 1)) / 2)).subs(y_, sp.Rational(37, 10)))) < 1e-30)
# J2: the closed form of sum y^rho/(rho(1-rho))^2 (proof note C8)
def T2(a): return mp.nsum(lambda k: mp.mpf(a)**(-2 * k) / (4 * k**2 * (2 * k + 1)**2), [1, mp.inf])
near('T2(1) = pi^2/6 - 3 + 2 ln 2', float(T2(1) - (mp.pi**2 / 6 - 3 + 2 * mp.log(2))), 0, 1e-30)

# ---------------------------------------------------------------------------------------------
section('2. J(a) from 40,000 zeros of zeta against the closed form, 1 <= a <= 50')
gz = np.loadtxt(ZFILE); NZ = len(gz); TZ = float(gz[-1])
def isprime(p): return p > 1 and all(p % d for d in range(2, int(p**0.5) + 1))
PPL = sorted((p**k, p, k) for p in range(2, 2001) if isprime(p) for k in range(1, 12) if p**k <= 2000)
ppn = np.array([q for q, _, _ in PPL], float); ppb = np.array([p for _, p, _ in PPL]); ppk = np.array([k for _, _, k in PPL])
ppl = np.log(ppb.astype(float))
def chi_arr(q, m): return np.array([chi_val(q, int(v)) for v in m])
W = {'z': ppl, '4': ppl * chi_arr(4, ppb) ** ppk, '3': ppl * chi_arr(3, ppb) ** ppk}
CS1 = {k: np.concatenate(([0.0], np.cumsum(w / ppn))) for k, w in W.items()}
CS0 = {k: np.concatenate(([0.0], np.cumsum(w))) for k, w in W.items()}
def Pk(key, a):                       # sum_{n <= a} Lambda(n) chi(n) (a/n - 1)
    a = np.asarray(a, float); i = np.searchsorted(ppn, a, 'right'); return a * CS1[key][i] - CS0[key][i]
def Tf(a):
    a = np.atleast_1d(np.asarray(a, float)); out = np.full(a.shape, 1 - math.log(2)); m = a > 1; y = a[m]
    out[m] = 1 - 0.5 * np.log1p(-y**-2.0) - y * np.arctanh(1 / y); return out
def Tof(a):
    a = np.atleast_1d(np.asarray(a, float)); out = np.full(a.shape, math.log(2)); m = a > 1; y = a[m]
    out[m] = np.arctanh(1 / y) + 0.5 * y * np.log1p(-y**-2.0); return out
def Jc(key, a):                       # the closed forms J(a), J_chi(a) (odd chi)
    a = np.atleast_1d(np.asarray(a, float))
    if key == 'z': return Pk('z', a) - a * np.log(a) + (1 + EGf) * a - LN2PI + Tf(a)
    q = int(key); return Pk(key, a) + a * LP1[q] - LP0[q] + Tof(a)
near('J(1) = C (closed form)', Jc('z', 1.0)[0], C, 1e-15)
for q in (4, 3): near('J_chi_-%d(1) = C_chi (closed form)' % q, Jc(str(q), 1.0)[0], CCHI[q], 1e-15)
gfun = lambda t: (0.5 * mp.log(t / (2 * mp.pi)) - 1 / (48 * t**2)) / (mp.pi * (mp.mpf(1) / 4 + t**2))   # theta'(t)/(pi(1/4+t^2))
Tm = mp.mpf(TZ); gder = [complex(mp.diff(gfun, Tm, k)) for k in range(10)]
ST = NZ - 1 - float(mp.siegeltheta(Tm)) / math.pi
tail0 = float(mp.quad(gfun, [Tm, mp.inf]))
def Jz(a):                            # 2 sqrt(a) [sum_{n<=N} cos(g u)/(1/4+g^2) + smooth tail - f(T) S(T+)]
    out = []
    for ai in np.atleast_1d(a):
        u = math.log(ai); s = float(np.sum(np.cos(gz * u) / (0.25 + gz * gz)))
        if u == 0: I = tail0
        else:
            z = sum((-1)**k * gder[k] / (1j * u)**(k + 1) for k in range(10)); I = (-np.exp(1j * u * TZ) * z).real
        out.append(2 * math.sqrt(ai) * (s + I - math.cos(u * TZ) / (0.25 + TZ * TZ) * ST))
    return np.array(out)
def Ires(w):
    w = np.abs(w); si, ci = sici(w * TZ); return np.cos(w * TZ) / TZ - w * (np.pi / 2 - si)
A = np.round(np.arange(1.0, 50.0000001, 0.01), 10); A = A[(A == 1.0) | (A >= 1.01)]
d = Jz(A) - Jc('z', A)
res = np.zeros_like(A)
for q, l in zip(ppn[ppn <= 200], ppl[ppn <= 200]): res += (l / math.pi) * np.sqrt(A / q) * Ires(np.log(A) - math.log(q))
far = np.ones(len(A), bool)
for q in ppn[ppn <= 51]: far &= np.abs(np.log(A) - math.log(q)) > 0.004
print('     grid: %d values of a, T = %.3f, S(T+) = %.4f' % (len(A), TZ, ST))
near('grid size (step 0.01 on [1.01, 50], and a = 1)', len(A), 4901, 0)
near('height H of the last zero', TZ, 33190.012, 5e-4); near('S(H+)', ST, 0.4648, 5e-5)
near('zeros - closed at a = 1', d[0], -3.9e-14, 5e-15)
near('max |zeros - closed| with |ln a - ln n| > 0.004 for every prime power n', np.abs(d[far]).max(), 2.5e-7, 5e-9)
near('max |zeros - closed - resonance| over the whole grid', np.abs(d - res).max(), 5.0e-9, 1e-10)
mx = 0.0
for q, l in zip(ppn[ppn <= 50], ppl[ppn <= 50]):
    i = int(np.argmin(np.abs(A - q))); mx = max(mx, abs(d[i] * math.pi * TZ / l - 1))
chk('at each prime power n <= 50: zeros - closed = Lambda(n)/(pi T) to within 0.6%', mx < 0.006, 'max relative deviation %.4f' % mx)
near('Lambda(2)/(pi T)', math.log(2) / (math.pi * TZ), 6.648e-6, 5e-10)
rows = {}
for a0 in (1.5, 2.0, 2.5, 6.0, 10.0, 30.5, 50.0):
    i = int(np.argmin(np.abs(A - a0))); rows[a0] = (Jc('z', a0)[0], d[i], res[i], d[i] - res[i])
    print('     a = %5.1f  J closed %+.12f  zeros-closed %+.3e  resonance %+.3e  residual %+.1e' % ((a0,) + rows[a0]))
near('J(1.5)', rows[1.5][0], 0.0065637, 5e-8); near('J(2)', rows[2.0][0], -0.0245114, 5e-8)
near('J(2.5)', rows[2.5][0], 0.0157764, 5e-8); near('J(6)', rows[6.0][0], 0.0328967, 5e-8)
near('J(10)', rows[10.0][0], 0.0245927, 5e-8); near('J(30.5)', rows[30.5][0], -0.0614892, 5e-8); near('J(50)', rows[50.0][0], 0.0417043, 5e-8)
near('zeros - closed at a = 2', rows[2.0][1], 6.648e-6, 5e-10)
near('zeros - closed at a = 6', rows[6.0][1], -3.6e-9, 1e-10)
near('residual at a = 6 after the resonance term', rows[6.0][3], 5.7e-10, 1e-10)
# J2 (the tube echo kernel) against the zeros
def J2c(y):
    y = mp.mpf(y); P2 = sum(mp.log(p) * ((y / n + 1) * mp.log(y / n) - 2 * (y / n - 1)) for n, p, k in PPL if n <= y)
    return (y * (mp.log(y)**2 / 2 - (2 + EG) * mp.log(y) + 3 + 2 * EG + EG**2 + 2 * G1) - mp.log(2 * mp.pi) * (mp.log(y) + 2)
            + 2 * G1 + EG**2 - mp.pi**2 / 12 - T2(y) - P2)
mx = 0.0
for y in (1.0, 2.0, 5.0, 10.5, 40.0):
    zs = 2 * math.sqrt(y) * float(np.sum(np.cos(gz * math.log(y)) / (0.25 + gz * gz)**2)); mx = max(mx, abs(zs - float(J2c(y))))
near('J2(y) closed form vs 40,000 zeros at y = 1, 2, 5, 10.5, 40', mx, 2.6e-14, 5e-15)
near('J2(1) = K2', float(J2c(1)), K2, 1e-20)

# ---------------------------------------------------------------------------------------------
section('3. correlation coefficients r(a) = A(a)/(sqrt(a) A(1)), closed form')
CC = {'z': C, '4': CCHI[4], '3': CCHI[3]}
def r(key, a): a = np.atleast_1d(np.asarray(a, float)); return Jc(key, a) / (CC[key] * np.sqrt(a))
REF_R = {'z': [0.6540, -0.1369, 0.1160, -0.3752, -0.4540, -0.0457, -0.4533, 0.2907, -0.4418, -0.0445, 0.0285, 0.1684, -0.2515, -0.2101],
         '4': [0.8534, 0.2822, -0.2671, -0.1379, 0.4585, -0.1305, -0.5372, -0.0066, 0.4536, 0.2296, 0.0371, 0.1161, -0.4695, -0.2642],
         '3': [0.8112, 0.1345, -0.2300, 0.3122, -0.1021, -0.1850, 0.5034, -0.0591, -0.5185, -0.0414, 0.1301, 0.2879, 0.1209, 0.0258]}
AS = (1.01, 1.1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100, 1000)
for key in ('z', '4', '3'):
    vals = r(key, np.array(AS, float))
    print('     %s: ' % key + ' '.join('%g:%+.4f' % (a, v) for a, v in zip(AS, vals)))
    near('r table (%s), 14 values, max deviation' % key, float(np.max(np.abs(vals - np.array(REF_R[key])))), 0.0, 5.01e-5)
AG = np.exp(np.linspace(0, math.log(1000), 2_000_001))
EXT = {}
for key in ('z', '4', '3'):
    R = r(key, AG)
    def refine(i, sgn):
        lo, hi = AG[max(i - 2, 0)], AG[min(i + 2, len(AG) - 1)]
        f = lambda x: sgn * r(key, x)[0]
        cands = [lo, hi] + [q for q in ppn if lo <= q <= hi]
        cands.append(minimize_scalar(f, bounds=(lo, hi), method='bounded', options={'xatol': 1e-12}).x)
        b = min(cands, key=f); return float(b), float(r(key, b)[0])
    i1 = int(np.argmax(R < 0)); z0 = brentq(lambda x: r(key, x)[0], AG[i1 - 1], AG[i1])
    jmin = int(np.argmin(R)); after = AG > z0; jmax = int(np.argmax(np.where(after, R, -9)))
    EXT[key] = (z0, refine(jmin, 1), refine(jmax, -1), math.sqrt(float(np.mean(R**2))), float(np.mean(R)))
    print('     %s: first zero %.6f; min %s; largest value after the first zero %s; rms %.4f; mean %+.4f' % ((key,) + EXT[key]))
near('zeta: first zero of r', EXT['z'][0], 1.06585, 5e-6)
near('zeta: min of r on [1,1000] at a = 3', EXT['z'][1][0], 3.0, 0); near('zeta: min value', EXT['z'][1][1], -0.45400, 5e-6)
near('zeta: largest r after the first zero, location', EXT['z'][2][0], 208.115, 5e-4); near('zeta: its value', EXT['z'][2][1], 0.39341, 5e-6)
near('zeta: rms of r over ln a on [1,1000]', EXT['z'][3], 0.1854, 5e-5)
near('chi_-4: first zero of r', EXT['4'][0], 1.19747, 5e-6)
near('chi_-4: min at a = 5', EXT['4'][1][0], 5.0, 0); near('chi_-4: min value', EXT['4'][1][1], -0.53720, 5e-6)
near('chi_-4: largest r after the first zero at a = 491', EXT['4'][2][0], 491.0, 0); near('chi_-4: its value', EXT['4'][2][1], 0.60762, 5e-6)
near('chi_-4: rms of r on [1,1000]', EXT['4'][3], 0.2730, 5e-5)
near('chi_-3: first zero of r', EXT['3'][0], 1.13650, 5e-6)
near('chi_-3: min at a = 7', EXT['3'][1][0], 7.0, 0); near('chi_-3: min value', EXT['3'][1][1], -0.51852, 5e-6)
near('chi_-3: largest r after the first zero at a = 257', EXT['3'][2][0], 257.0, 0); near('chi_-3: its value', EXT['3'][2][1], 0.51542, 5e-6)
near('chi_-3: rms of r on [1,1000]', EXT['3'][3], 0.2304, 5e-5)
# the link to the variance theorem of the proof note: delta D(delta) = (2 + delta) A(1) - 2 A(1 + delta)
D1 = float(3 * mp.log(3) - mp.log(4 * mp.pi) - EG)
near('D(1) = ln(27/(4pi)) - gamma = (3 C - 2 J(2))/1', 3 * C - 2 * Jc('z', 2.0)[0], D1, 1e-14)
near('D(1) value', D1, 0.187596954133505, 1e-15)
# near a = 1: r(1 + d) = 1 - (d/(2C)) (ln(1/d) + B) + O(d^2 ln(1/d))
Bc = float(1 - EG - mp.log(2 * mp.pi))
for dd, ref in ((0.01, 0.6540), (0.001, 0.9406)):
    ex = r('z', 1 + dd)[0]; ap = 1 - dd / (2 * C) * (math.log(1 / dd) + Bc)
    near('r(1 + %g): the two-term expansion' % dd, ap, {0.01: 0.6547, 0.001: 0.9405}[dd], 5e-5)
    near('r(1 + %g) exact' % dd, ex, ref, 5e-5); chk('r(1 + %g): two-term expansion within 1e-3' % dd, abs(ex - ap) < 1e-3, 'expansion %.4f' % ap)
chk('T\'\'(y) = 1/(y^2 (y^2 - 1)) (sympy), so J\'\' = -1/a + 1/(a^2(a^2-1)) between prime powers',
    sp.simplify(sp.diff(Tsym, y_, 2) - 1 / (y_**2 * (y_**2 - 1))) == 0)
acv = brentq(lambda x: x**3 - x - 1, 1.1, 1.5); near('J convex on (1, a*), concave beyond: a*^3 = a* + 1', acv, 1.32472, 5e-6)
rm = minimize_scalar(lambda x: r('z', x)[0], bounds=(1.01, 1.32), method='bounded', options={'xatol': 1e-10})
near('zeta: interior local minimum of r on (1, a*): location', rm.x, 1.16482, 5e-5); near('zeta: interior local minimum of r: value', rm.fun, -0.20888, 5e-6)
ah = brentq(lambda x: r('z', x)[0] - 0.5, 1.0001, 1.05)
near('zeta: r(a) = 1/2 at a', ah, 1.01743, 5e-6)
# the race pi(x;4,3) - pi(x;4,1): covariance J_chi(a)/sqrt(a)
for a0, ref in ((2, -0.0214575), (3, 0.071325), (5, -0.083571)):
    near('race chi_-4: covariance J_chi(%d)/sqrt(%d)' % (a0, a0), Jc('4', a0)[0] / math.sqrt(a0), ref, 5e-7)

# ---------------------------------------------------------------------------------------------
section('4. zeros of L(s, chi_-4), L(s, chi_-3) from PARI to height 1000')
TCH = 1000.0
def theta_chi(q, t): t = mp.mpf(t); return t / 2 * mp.log(mp.mpf(q) / mp.pi) + mp.im(mp.loggamma((mp.mpf(3) / 2 + 1j * t) / 2))
CHZ = {}
for q, D in ((4, -4), (3, -3)):
    z = np.array([float(v) for v in pari.lfunzeros(pari.lfuncreate(D), TCH)]); CHZ[q] = z
    Sv = np.array([i + 1 - float(theta_chi(q, z[i])) / math.pi for i in range(len(z))])
    blocks = [Sv[(z > b) & (z <= b + 100)].mean() for b in range(0, 1000, 100)]
    chk('chi_-%d: %d zeros; mean of S(gamma_n+) in each block of height 100 rounds to 0.50' % (q, len(z)),
        all(abs(m - 0.5) < 0.005 for m in blocks), ' '.join('%.3f' % m for m in blocks))
    thp = lambda t: (mp.log(mp.mpf(q) / mp.pi) / 2 + mp.re(mp.digamma((mp.mpf(3) / 2 + 1j * t) / 2)) / 2) / mp.pi
    Tq = z[-1]; STq = len(z) - float(theta_chi(q, Tq)) / math.pi
    Cz = 2 * np.sum(1 / (0.25 + z * z)) + 2 * float(mp.quad(lambda t: thp(t) / (mp.mpf(1) / 4 + t**2), [Tq, mp.inf])) - 2 * STq / (0.25 + Tq * Tq)
    near('chi_-%d: C_chi from the zeros (sum 1/|rho|^2 + tail) = closed form' % q, Cz, CCHI[q], 2e-8)
    K2q = 2 * np.sum(1 / (0.25 + z * z)**2)
    print('     chi_-%d: K2_chi = sum 1/|rho|^4 = %.6e (zeros to 1000), rms of r: sqrt(K2)/C = %.4f' % (q, K2q, math.sqrt(K2q) / CCHI[q]))
    near('chi_-%d: sqrt(K2_chi)/C_chi' % q, math.sqrt(K2q) / CCHI[q], {4: 0.2754, 3: 0.2314}[q], 5e-5)
    near('chi_-%d: K2_chi from the zeros to 1000' % q, K2q, {4: 1.835740e-3, 3: 6.866039e-4}[q], 5e-10)
for q, D, ref in ((4, -4, 0.667691), (3, -3, 0.480868)):
    near('L(1/2, chi_-%d) (nonzero: no frequency 0, so the race has mean 1)' % q, float(pari.lfun(pari.lfuncreate(D), 0.5)), ref, 5e-7)
near('zeta: count of zeros below 1000', int(np.sum(gz <= 1000)), 649, 0)
near('chi_-4: zeros below 1000', len(CHZ[4]), 868, 0); near('chi_-3: zeros below 1000', len(CHZ[3]), 823, 0)
dmin = lambda u, v: float(np.min(np.abs(u[:, None] - v[None, :])))
gz1 = gz[gz <= 1000]
md1, md2 = dmin(gz1, CHZ[4]), dmin(CHZ[4], CHZ[3])
print('     closest ordinates below 1000: zeta vs chi_-4 %.2e, chi_-4 vs chi_-3 %.2e' % (md1, md2))
near('closest ordinates below 1000, zeta vs chi_-4', md1, 5.2e-4, 5e-6); near('closest ordinates below 1000, chi_-4 vs chi_-3', md2, 1.6e-4, 5e-6)

# ---------------------------------------------------------------------------------------------
section('5. primes to 10^8: windowed means against the law')
NMAX = 10**8
sv = np.ones(NMAX // 2, dtype=bool); sv[0] = False
for i in range(3, int(NMAX**0.5) + 1, 2):
    if sv[i // 2]: sv[i * i // 2::i] = False
PR = np.concatenate(([2], 2 * np.nonzero(sv)[0] + 1)).astype(np.int64); del sv
near('pi(10^8)', len(PR), 5761455, 0)
extra = []
for p in PR[PR <= 10**4].tolist():
    qq, k = p * p, 2
    while qq <= NMAX: extra.append((qq, p, k)); qq *= p; k += 1
extra = np.array(extra, dtype=np.int64)
nn = np.concatenate((PR, extra[:, 0])); bb = np.concatenate((PR, extra[:, 1])); kk = np.concatenate((np.ones(len(PR), np.int64), extra[:, 2]))
o = np.argsort(nn, kind='stable'); nn, bb, kk = nn[o], bb[o], kk[o]; del o, PR, extra
lam = np.log(bb.astype(float)); lnn = np.log(nn.astype(float))
c4 = np.where(bb % 2 == 0, 0, np.where(bb % 4 == 1, 1, -1)) ** kk
c3 = np.where(bb % 3 == 0, 0, np.where(bb % 3 == 1, 1, -1)) ** kk
CUM = {'z': np.cumsum(lam), '4': np.cumsum(lam * c4), '3': np.cumsum(lam * c3),
       'N': np.cumsum(np.where(kk == 1, -c4, 0)).astype(float), 'R': np.cumsum(np.where(kk == 1, -c4 * lam, 0.0))}
near('psi(10^8) - 10^8', CUM['z'][-1] - 1e8, -1757.2034, 5e-3)
X = 1e8; I1 = float(np.sum(lam / nn)) - CUM['z'][-1] / X - math.log(X)
del nn, bb, kk, lam, c4, c3
near('int_1^X E(t) dt/t^2 + 1 + gamma at X = 10^8', I1 + 1 + EGf, -1.1e-6, 5e-8)
near('its bound (C sqrt X + ln 2pi + 1)/X', (C * math.sqrt(X) + LN2PI + 1) / X, 4.6e-6, 5e-8)
Y0, Y1 = math.log(1e2), math.log(1e7); LW = Y1 - Y0; YC = Y0 + LW / 2
def fval(key, k, x):                 # normalised functions: E0/sqrt(x) etc.; k = index of the last prime power <= x (or -1)
    if key == '1': return np.ones_like(x)
    base = key[0]; cv = np.where(k >= 0, CUM[base][np.maximum(k, 0)], 0.0)
    if base == 'z': v = (cv - x + LN2PI + 0.5 * np.log1p(-x**-2.0)) / np.sqrt(x)
    elif base in '43': v = (cv + LP0[int(base)] - np.arctanh(1.0 / x)) / np.sqrt(x)
    elif base == 'N': v = np.log(x) * cv / np.sqrt(x)
    elif base == 'R': v = cv / np.sqrt(x)
    return v ** (int(key[1]) if len(key) > 1 else 1)
GLT, GLW = np.polynomial.legendre.leggauss(3)
def wmeans(a, pairs, chunk=1_000_000):
    """Hann-weighted mean over y = ln x in [Y0, Y1] of fA(a e^y) fB(e^y), exact up to 3-point Gauss-Legendre on each
    interval between consecutive jumps of the two step functions."""
    la = math.log(a)
    i0, i1 = np.searchsorted(lnn, Y0), np.searchsorted(lnn, Y1); j0, j1 = np.searchsorted(lnn, Y0 + la), np.searchsorted(lnn, Y1 + la)
    B = np.concatenate(([Y0, Y1], lnn[i0:i1], lnn[j0:j1] - la)); B.sort()
    out = np.zeros(len(pairs))
    for c0 in range(0, len(B) - 1, chunk):
        lo = B[c0:c0 + chunk + 1]; dd = np.diff(lo); ym = lo[:-1] + dd / 2
        k2 = np.searchsorted(lnn, ym, 'right') - 1; k1 = np.searchsorted(lnn, ym + la, 'right') - 1
        for gidx in range(3):
            y = ym + dd / 2 * GLT[gidx]; x = np.exp(y)
            wt = GLW[gidx] * dd / 2 * (1 - np.cos(2 * np.pi * (y - Y0) / LW)) / LW
            cache = {}
            for i, (ka, kb) in enumerate(pairs):
                if ('A', ka) not in cache: cache[('A', ka)] = fval(ka, k1, a * x)
                if ('B', kb) not in cache: cache[('B', kb)] = fval(kb, k2, x)
                out[i] += np.dot(wt, cache[('A', ka)] * cache[('B', kb)])
    return out
# zero side of the finite window: sum over pairs lambda != lambda' of beta conj(beta') e^{i lambda b} Kh(lambda - lambda')
def Kh(th):
    z = th * LW / 2; small = np.abs(z) < 1e-8; big = np.abs(np.abs(z) - np.pi) < 1e-6
    zz = np.where(small | big, 1.0, z); out = np.sin(zz) / zz * np.pi**2 / (np.pi**2 - zz**2)
    out = np.where(small, 1.0, np.where(big, 0.5, out)); return np.exp(1j * th * YC) * out
def fr(pos): l = np.concatenate((pos, -pos)); return l, 1.0 / (0.5 + 1j * l)
def vdense(l1, b1, l2, b2, same):
    M = Kh(l1[:, None] - l2[None, :]) * np.conj(b2)[None, :]
    if same: np.fill_diagonal(M, 0)
    return M.sum(axis=1)
def vband(l, b, Dw=60.0):
    o = np.argsort(l); ls, bs = l[o], b[o]; v = np.zeros(len(ls), complex); k = 1
    while k < len(ls):
        m = (ls[k:] - ls[:-k]) < Dw
        if not m.any(): break
        i = np.nonzero(m)[0]; j = i + k; Kij = Kh(ls[i] - ls[j])
        np.add.at(v, i, np.conj(bs[j]) * Kij); np.add.at(v, j, np.conj(bs[i]) * np.conj(Kij)); k += 1
    vv = np.zeros(len(ls), complex); vv[o] = v; return vv
def offd(l, b, v, a): return float((math.sqrt(a) * np.sum(b * np.exp(1j * l * math.log(a)) * v)).real)
LZ, BZ = fr(gz); VZ_all = vband(LZ, BZ); VZ_120 = vband(LZ, BZ, 120.0)
near('zeta window term: band 60 vs band 120, max over a = 1, 2, 10', max(abs(offd(LZ, BZ, VZ_all, a) - offd(LZ, BZ, VZ_120, a)) for a in (1.0, 2.0, 10.0)), 0.0, 1e-6)
del VZ_120
LZ1, BZ1 = fr(gz1); VZ1 = vdense(LZ1, BZ1, LZ1, BZ1, True)
L4, B4 = fr(CHZ[4]); V44 = vdense(L4, B4, L4, B4, True)
L3, B3 = fr(CHZ[3]); V33 = vdense(L3, B3, L3, B3, True)
VZ4 = vdense(LZ1, BZ1, L4, B4, False); V43 = vdense(L4, B4, L3, B3, False)
RATIOS = (1.0, 1.5, 2.0, 2.9, 3.0, 3.1, 4.9, 5.0, 5.1, 7.0, 10.0)
WM = {}
for a in RATIOS:
    pairs = [('z', 'z'), ('4', '4'), ('3', '3')]
    if a in (1.0, 2.0, 3.0): pairs += [('z', '4'), ('4', '3')]
    if a in (1.0, 2.0, 3.0, 5.0): pairs += [('N', 'N'), ('N', '1'), ('1', 'N'), ('R', 'R'), ('R', '1'), ('1', 'R')]
    if a == 1.0: pairs += [('z', 'z2'), ('z', 'z3'), ('4', '42'), ('4', '43')]
    vals = wmeans(a, pairs); WM[a] = dict(zip(pairs, vals))
print('     window: ln x in [ln 1e2, ln 1e7], length %.4f, Hann weight; zeta pairs with %d zeros (band 60), chi pairs with zeros to 1000' % (LW, len(gz)))
near('window length', LW, 11.5129, 5e-5)
TAB = {}
for a in RATIOS:
    sa = math.sqrt(a)
    pz = WM[a][('z', 'z')] * sa; p4 = WM[a][('4', '4')] * sa; p3 = WM[a][('3', '3')] * sa
    jz, j4, j3 = Jc('z', a)[0], Jc('4', a)[0], Jc('3', a)[0]
    oz, oz1, o4, o3 = offd(LZ, BZ, VZ_all, a), offd(LZ1, BZ1, VZ1, a), offd(L4, B4, V44, a), offd(L3, B3, V33, a)
    TAB[a] = (pz, jz, oz, pz - jz - oz, oz1, p4, j4, o4, p4 - j4 - o4, p3, j3, o3, p3 - j3 - o3)
    print('     a=%4.1f zeta: primes %+.7f law %+.7f window %+.2e resid %+.1e (T=1000 window %+.2e) | chi4: %+.7f %+.7f %+.2e %+.1e | chi3: %+.7f %+.7f %+.2e %+.1e'
          % ((a,) + TAB[a]))
mz = max(abs(TAB[a][3]) for a in RATIOS); m4 = max(abs(TAB[a][8]) for a in RATIOS); m3 = max(abs(TAB[a][12]) for a in RATIOS)
mz1 = max(abs(TAB[a][2] - TAB[a][4]) for a in RATIOS)
md = max(max(abs(TAB[a][0] - TAB[a][1]), abs(TAB[a][5] - TAB[a][6]), abs(TAB[a][9] - TAB[a][10])) for a in RATIOS)
near('max |primes - law| before the window correction, 11 ratios, 3 functions', md, 7.6e-4, 5e-6)
near('zeta: max |primes - law - window| (40,000 zeros)', mz, 3.8e-6, 1e-7)
near('zeta: max change of the window term from 1000 to 40,000 zeros', mz1, 2.8e-5, 1e-6)
near('chi_-4: max |primes - law - window| (zeros to 1000)', m4, 4.8e-5, 1e-6)
near('chi_-3: max |primes - law - window| (zeros to 1000)', m3, 4.7e-5, 1e-6)
for a0, refs in ((1.0, (0.046168, 0.046191, 0.155546, 0.155568, 0.113120, 0.113230)), (2.0, (-0.024498, -0.024511, -0.030498, -0.030346, 0.049976, 0.049988)),
                 (3.0, (-0.036279, -0.036323, 0.123648, 0.123539, -0.019783, -0.020029)), (10.0, (0.024159, 0.024593, 0.057319, 0.057107, 0.102915, 0.103070))):
    vals = (TAB[a0][0], TAB[a0][1], TAB[a0][5], TAB[a0][6], TAB[a0][9], TAB[a0][10])
    print('     a = %g:  zeta %+.6f (law %+.6f)  chi_-4 %+.6f (law %+.6f)  chi_-3 %+.6f (law %+.6f)' % ((a0,) + vals))
    near('table row a = %g (primes and law, 3 functions), max deviation' % a0, max(abs(u - v) for u, v in zip(vals, refs)), 0.0, 5.01e-7)
# corners: second differences with h = 0.1 at a = 3 and 5
JUMP = {('z', 3.0): math.log(3) / 3, ('z', 5.0): math.log(5) / 5, ('4', 3.0): -math.log(3) / 3, ('4', 5.0): math.log(5) / 5,
        ('3', 3.0): 0.0, ('3', 5.0): -math.log(5) / 5}
CORN = {}
for key, ix in (('z', 0), ('4', 5), ('3', 9)):
    for n0 in (3.0, 5.0):
        sd = lambda j: TAB[n0 + 0.1][j] - 2 * TAB[n0][j] + TAB[n0 - 0.1][j]
        dp, dl, dw = sd(ix), sd(ix + 1), sd(ix + 2)
        CORN[(key, n0)] = (dp, dl, dl + dw, dl - 0.1 * JUMP[(key, n0)])
        print('     corner %s at %g: second difference primes %+.5f, law %+.5f, law + window %+.5f, law without the corner %+.5f' % ((key, n0) + CORN[(key, n0)]))
REFC = {('z', 3.0): (0.03340, 0.03343, -0.00319), ('z', 5.0): (0.02956, 0.03021, -0.00198), ('4', 3.0): (-0.03635, -0.03620, 0.00042),
        ('4', 5.0): (0.03054, 0.03227, 0.00008), ('3', 3.0): (-0.00002, 0.00042, 0.00042), ('3', 5.0): (-0.03131, -0.03211, 0.00008)}
for kn, (rp, rl, rn) in REFC.items():
    dp, dl, dlw, dn = CORN[kn]
    near('second difference at %g (%s): primes' % (kn[1], kn[0]), dp, rp, 5e-6)
    near('second difference at %g (%s): law' % (kn[1], kn[0]), dl, rl, 5e-6)
    near('second difference at %g (%s): law without its corner' % (kn[1], kn[0]), dn, rn, 5e-6)
    chk('second difference at %g (%s): |primes - (law + window)| < 1e-4' % (kn[1], kn[0]), abs(dp - dlw) < 1e-4, '%.1e' % (dp - dlw))
# cross-correlations
for a in (1.0, 2.0, 3.0):
    sa = math.sqrt(a)
    xz4 = WM[a][('z', '4')] * sa; x43 = WM[a][('4', '3')] * sa
    pz4 = offd(LZ1, BZ1, VZ4, a); p43 = offd(L4, B4, V43, a)
    print('     cross a=%g: zeta x chi4 primes %+.6f, window prediction %+.6f, diff %+.1e;  chi4 x chi3 primes %+.6f, prediction %+.6f, diff %+.1e'
          % (a, xz4, pz4, xz4 - pz4, x43, p43, x43 - p43))
    TAB[('x', a)] = (xz4, pz4, x43, p43)
REFX = {1.0: (-0.001063, -0.001029, 0.005724, 0.005705), 2.0: (-0.007410, -0.007409, 0.005138, 0.005096), 3.0: (0.003865, 0.003892, 0.005583, 0.005543)}
for a in (1.0, 2.0, 3.0):
    near('cross a = %g: primes and window predictions (4 numbers), max deviation' % a, max(abs(u - v) for u, v in zip(TAB[('x', a)], REFX[a])), 0.0, 5.01e-7)
mxc = max(max(abs(TAB[('x', a)][0] - TAB[('x', a)][1]), abs(TAB[('x', a)][2] - TAB[('x', a)][3])) for a in (1.0, 2.0, 3.0))
near('cross: max |primes - window prediction| (a = 1, 2, 3)', mxc, 4.3e-5, 1e-6)
# the race
for a in (1.0, 2.0, 3.0, 5.0):
    out = []
    for key in ('N', 'R'):
        m2 = WM[a][(key, key)]; ma = WM[a][(key, '1')]; mb = WM[a][('1', key)]
        out += [ma, mb, m2 - ma * mb]
    TAB[('race', a)] = out
    print('     race a=%g: N: mean(ax) %.4f mean(x) %.4f cov %+.5f | theta-race: mean %.4f %.4f cov %+.5f | law cov J_chi(a)/sqrt(a) = %+.5f'
          % ((a,) + tuple(out) + (Jc('4', a)[0] / math.sqrt(a),)))
for a, refN, refR, refL in ((1.0, 0.15696, 0.15997, 0.15557), (2.0, -0.02010, -0.01761, -0.02146), (3.0, 0.07114, 0.07317, 0.07133), (5.0, -0.08349, -0.08168, -0.08357)):
    near('race a = %g: window covariance of the normalised pi-race' % a, TAB[('race', a)][2], refN, 5e-6)
    near('race a = %g: window covariance of the theta-race' % a, TAB[('race', a)][5], refR, 5e-6)
    near('race a = %g: law J_chi(a)/sqrt(a)' % a, Jc('4', a)[0] / math.sqrt(a), refL, 5e-6)
near('race: window mean of the normalised pi-race', TAB[('race', 1.0)][1], 1.2233, 5e-5)
near('race: window mean of the theta-race', TAB[('race', 1.0)][4], 0.9214, 5e-5)
near('race: max |window covariance - law| over N, R_theta and a = 1, 2, 3, 5', max(max(abs(TAB[('race', a)][2] - Jc('4', a)[0] / math.sqrt(a)), abs(TAB[('race', a)][5] - Jc('4', a)[0] / math.sqrt(a))) for a in (1.0, 2.0, 3.0, 5.0)), 0.0044, 5e-5)
near('theta(10)/10 - 1 (the prime-square term at x = 10^2)', math.log(2 * 3 * 5 * 7) / 10 - 1, -0.4653, 5e-5)

# ---------------------------------------------------------------------------------------------
section('6. third and fourth moments (a = 1) against the LI prediction')
m2z = WM[1.0][('z', 'z')]; m3z = WM[1.0][('z', 'z2')]; m4z = WM[1.0][('z', 'z3')]
m24 = WM[1.0][('4', '4')]; m34 = WM[1.0][('4', '42')]; m44 = WM[1.0][('4', '43')]
del CUM, lnn
def Qgrid(pos, yg, chunk=100):
    Q = np.zeros(len(yg))
    for i in range(0, len(pos), chunk):
        p = pos[i:i + chunk]; Q += 2 * np.real((np.exp(1j * np.outer(yg, p)) / (0.5 + 1j * p)).sum(axis=1))
    return Q
def moments(pos, ng):
    yg = np.linspace(Y0, Y1, ng); hg = yg[1] - yg[0]; wg = (1 - np.cos(2 * np.pi * (yg - Y0) / LW)) / LW
    Q = Qgrid(pos, yg); wint = lambda f: float(np.sum(wg * f) * hg)
    PT = 2 * np.sum(1 / (0.25 + pos**2)); K2T = 2 * np.sum(1 / (0.25 + pos**2)**2)
    return -wint(Q**3), wint(Q**4) - 3 * PT**2 + 3 * K2T, K2T
gz4 = gz[gz <= 4000]
m3a, dTa, _ = moments(gz1, 60001); m3b, dTb, _ = moments(gz[gz <= 2000], 60001); m3c, dTc, _ = moments(gz4, 120001)
pred4 = 3 * C**2 - 3 * K2 + dTc
print('     zeta: M2 %.7f  M3 primes %+.4e (window value from zeros to 1000, 2000, 4000: %+.4e %+.4e %+.4e)' % (m2z, m3z, m3a, m3b, m3c))
print('     zeta: M4 primes %.7f; LI + window from zeros to 1000, 2000, 4000: %.7f %.7f %.7f; window part %+.3e; without the -3K2 term %.7f; kurtosis %.4f'
      % (m4z, 3 * C**2 - 3 * K2 + dTa, 3 * C**2 - 3 * K2 + dTb, pred4, dTc, pred4 + 3 * K2, m4z / m2z**2))
near('zeta: number of zeros to 4000', len(gz4), 3474, 0)
near('zeta: M3 primes', m3z, -1.689e-4, 5e-8); near('zeta: M3 window value (zeros to 4000)', m3c, -1.668e-4, 5e-8)
near('zeta: M3 window value (zeros to 2000)', m3b, -1.703e-4, 5e-8); near('zeta: M3 window value (zeros to 1000)', m3a, -1.751e-4, 5e-8)
near('zeta: M4 primes', m4z, 0.0059981, 5e-8)
near('zeta: M4 = 3C^2 - 3K2 + window part (zeros to 4000)', pred4, 0.0059975, 5e-8)
near('zeta: M4 window part (zeros to 4000)', dTc, -1.808e-4, 5e-8)
near('zeta: M4 with zeros to 2000', 3 * C**2 - 3 * K2 + dTb, 0.0059949, 5e-8)
near('zeta: M4 prediction without the -3K2 term', pred4 + 3 * K2, 0.0062201, 5e-8)
near('zeta: M4 without the -3K2 term minus primes', pred4 + 3 * K2 - m4z, 2.2e-4, 5e-6)
near('zeta: kurtosis M4/M2^2 in the window', m4z / m2z**2, 2.8140, 5e-5)
m3x, dTx, K24T = moments(CHZ[4], 60001)
p4x = 3 * CCHI[4]**2 - 3 * K24T + dTx
print('     chi_-4: M2 %.6f  M3 %+.4e (zeros to 1000 %+.4e)  M4 %.6f (LI + window %.6f)  kurtosis %.4f (LI limit %.4f)'
      % (m24, m34, m3x, m44, p4x, m44 / m24**2, 3 - 3 * K24T / CCHI[4]**2))
near('chi_-4: M3 primes', m34, -9.452e-4, 5e-8); near('chi_-4: M3 window value (zeros to 1000)', m3x, -8.870e-4, 5e-8)
near('chi_-4: M4 primes', m44, 0.070489, 5e-7); near('chi_-4: M4 = 3C^2 - 3K2 + window part (zeros to 1000)', p4x, 0.070487, 5e-7)
near('chi_-4: kurtosis in the window', m44 / m24**2, 2.9134, 5e-5)
near('chi_-4: LI kurtosis 3 - 3 K2_chi/C_chi^2', 3 - 3 * K24T / CCHI[4]**2, 2.7724, 5e-5)

# ---------------------------------------------------------------------------------------------
print('\n%d checks, %d failed  [%.0f s]' % (NCHK[0], len(FAILED), time.time() - T_START))
if FAILED: print('failed:', FAILED); sys.exit(1)
