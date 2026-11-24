"""check_proof.py: recomputes and asserts every number quoted in D_delta_proof.md (same folder).

usage:  python3 -I check_proof.py REPO_ROOT        (REPO_ROOT holds writeup/data/riemann/zeros_40000_flint.txt)

Self-contained: numpy, mpmath, sympy, cypari2 only; nothing is imported from the repo or from older scripts.
Every check prints one line, "ok" or "FAIL", and the run ends with "N checks, F failed".
A quoted number passes when the recomputed value rounds to it (half a unit in its last digit).
Sections of this script (the note's section numbers differ; each check line names what it checks):
  1  closed forms and constants              5  zeta: K2, sigma2, Dt against the zeros
  2  exact algebra and series (sympy)        6  D0 - HL_exact (two quadratures)
  3  the Weil route by quadrature            7  characters: zeros from PARI, closed forms, resonance
  4  zeta: closed form against 40,000 zeros  8  primes against zeros, end to end
"""
import sys, math, time
import numpy as np
import mpmath as mp
import sympy as sp

T_START = time.time()
ROOT = sys.argv[1] if len(sys.argv) > 1 else '.'
ZFILE = ROOT.rstrip('/') + '/writeup/data/riemann/zeros_40000_flint.txt'
mp.mp.dps = 30
EG = mp.euler
PI = mp.pi
G1 = mp.stieltjes(1)
CC = 2 + EG - mp.log(4 * PI)              # C = sum_rho 1/(rho(1-rho))
BB = 1 - EG - mp.log(2 * PI)              # B = 1 - gamma - ln 2pi

# ---------------------------------------------------------------- check helpers
NCHK = [0, 0]


def check(label, ok, detail=''):
    NCHK[0] += 1
    if not ok:
        NCHK[1] += 1
    print('%s  %s%s' % ('ok  ' if ok else 'FAIL', label, (': ' + detail) if detail else ''), flush=True)


def ulp_of(s):
    s = s.strip().lstrip('+')
    if 'e' in s:
        mant, ex = s.split('e')
        dec = len(mant.split('.')[1]) if '.' in mant else 0
        return mp.mpf(10) ** (int(ex) - dec)
    dec = len(s.split('.')[1]) if '.' in s else 0
    return mp.mpf(10) ** (-dec)


def quoted(label, value, s):
    """value must round to the decimal string s (fixed point or scientific)."""
    v = mp.mpf(value) if not isinstance(value, mp.mpf) else value
    target = mp.mpf(s.strip().lstrip('+'))
    ok = abs(v - target) <= ulp_of(s) / 2 * (1 + mp.mpf('1e-12'))
    check(label, ok, 'computed %s, quoted %s' % (mp.nstr(v, 17), s))


def small(label, value, bound):
    v = abs(float(value))
    check(label, v <= bound, '|%.3e| <= %.1e' % (float(value), bound))


# ---------------------------------------------------------------- arithmetic helpers
def lambda_table(M):
    lam = [mp.mpf(0)] * (M + 1)
    sieve = bytearray([1]) * (M + 1)
    for p in range(2, M + 1):
        if sieve[p]:
            for k in range(p * p, M + 1, p):
                sieve[k] = 0
            lp = mp.log(p)
            v = p
            while v <= M:
                lam[v] = lp
                v *= p
    return lam


LAM = lambda_table(1100)


def Psum(y, chi=None):
    """P(y) = sum_{n <= y} Lambda(n) Re chi(n) (y/n - 1)."""
    y = mp.mpf(y)
    tot = mp.mpf(0)
    for n in range(2, int(mp.floor(y)) + 1):
        if LAM[n] != 0:
            c = 1 if chi is None else mp.re(chi(n))
            tot += LAM[n] * c * (y / n - 1)
    return tot


def D_free(d):
    """prime-free part of the closed form."""
    d = mp.mpf(d)
    return mp.log(1 / d) - EG - mp.log(2 * PI) + 2 * mp.log1p(d) + ((2 + d) / d) * mp.log1p(d / 2)


def D_closed(d):
    d = mp.mpf(d)
    return D_free(d) - (2 / d) * Psum(1 + d)


def HL_exact(d):
    d = mp.mpf(d)
    return mp.log(1 / d) - EG - mp.log(2 * PI) + (1 + d) * mp.log1p(d) / d


def Tfun(y):
    """T(y) = sum_k y^-2k/(2k(2k+1)) (trivial zeros), closed form for y > 1, 1 - ln 2 at y = 1."""
    y = mp.mpf(y)
    if y == 1:
        return 1 - mp.log(2)
    return 1 - mp.log(1 - y ** -2) / 2 - (y / 2) * mp.log((y + 1) / (y - 1))


def Gfun(y):
    """fact 2048: G(y) = sum_{n<=y} Lambda(n)(y/n - 1) - y ln y + (1 + gamma) y."""
    y = mp.mpf(y)
    return Psum(y) - y * mp.log(y) + (1 + EG) * y


def Dchi_closed(d, q, a, chi, pole=False, parity=None):
    d = mp.mpf(d)
    a = a if parity is None else parity
    P0 = ((2 + d) / d) * mp.log1p(d / 2)
    P1 = (2 * (1 + d) / d) * mp.log1p(d) - P0
    v = mp.log(q / (2 * PI * d)) - EG + (P0 if a == 0 else P1) - (2 / d) * Psum(1 + d, chi)
    return v + (2 * mp.log1p(d) if pole else 0)


def h_mp(t, d):
    d = mp.mpf(d)
    u = mp.log1p(d)
    s = mp.sqrt(1 + d)
    return ((s - 1) ** 2 + 4 * s * mp.sin(u * t / 2) ** 2) / (mp.mpf(1) / 4 + t * t)


def h_np(t, d):
    u = math.log1p(d)
    s = math.sqrt(1 + d)
    return ((s - 1) ** 2 + 4 * s * np.sin(u * t / 2) ** 2) / (0.25 + t * t)


def osc_integral(f, u, a, K=60, nibp=8):
    """int_a^inf f(t) cos(u t) dt for smooth f decaying like (ln t)/t^2: tanh-sinh quadrature on half periods
    up to R (u R >= 200), then the integration-by-parts series -Re e^{iuR} sum_k (-1)^k f^(k)(R)/(iu)^(k+1)."""
    a = mp.mpf(a)
    hp = PI / u
    K = max(K, int(mp.ceil((200 - u * a) / PI)) + 1)
    R = a + K * hp
    pts = [a]
    if a < 1:
        x = mp.mpf('0.25')
        while x < min(hp, R):
            pts.append(a + x)
            x *= 2
    k0 = 1
    while a + k0 * hp <= R:
        if a + k0 * hp > pts[-1]:
            pts.append(a + k0 * hp)
        k0 += 1
    g = lambda t: f(t) * mp.cos(u * t)
    tot = mp.fsum(mp.quad(g, [pts[i], pts[i + 1]]) for i in range(len(pts) - 1))
    R = pts[-1]
    e = mp.exp(1j * u * R)
    ibp = -e * mp.fsum((-1) ** k * mp.diff(f, R, k) / (1j * u) ** (k + 1) for k in range(nibp))
    return tot + mp.re(ibp)


def split_integral(d, weight, a):
    """int_a^inf h_delta(t) weight(t) dt with h = (A - 2 s cos(ut)) q(t), q = 1/(1/4 + t^2)."""
    d = mp.mpf(d)
    u = mp.log1p(d)
    s = mp.sqrt(1 + d)
    A = 2 + d
    f = lambda t: weight(t) / (mp.mpf(1) / 4 + t * t)
    a = mp.mpf(a)
    br = [a] + ([mp.mpf(x) for x in ('0.25', '0.5', '1', '2', '4', '8', '16', '64')] if a < 1 else []) \
        + [2 * a + 1, 8 * a + 8, mp.inf]
    return A * mp.quad(f, br) - 2 * s * osc_integral(f, u, a)


def theta_prime(t):
    """theta'(t) for zeta: (1/2) Re psi(1/4 + it/2) - (1/2) ln pi."""
    return mp.re(mp.digamma(mp.mpf(1) / 4 + 0.5j * t)) / 2 - mp.log(PI) / 2


def resonance_I(w, T):
    """I(T) = int_T^inf cos(w t)/t^2 dt = cos(wT)/T - |w| (pi/2 - Si(|w| T))."""
    w = abs(mp.mpf(w))
    T = mp.mpf(T)
    if w == 0:
        return 1 / T
    return mp.cos(w * T) / T - w * (PI / 2 - mp.si(w * T))


def sec(title):
    print('', flush=True)
    print('== %s  (%.0f s)' % (title, time.time() - T_START), flush=True)


# ================================================================ 1. closed forms and constants
sec('1. closed forms and constants')
quoted('C = 2 + gamma - ln 4pi', CC, '0.0461914179322421')
quoted('C (short)', CC, '0.0461914')
quoted('B = 1 - gamma - ln 2pi', BB, '-1.41509273131')
quoted('B (short)', BB, '-1.41509')
quoted('gamma_1', G1, '-0.0728158454836767')
D1 = D_closed(1)
quoted('D(1) = ln(27/(4pi)) - gamma', mp.log(mp.mpf(27) / (4 * PI)) - EG, '0.18759695413350542')
small('D(1): closed form equals ln(27/(4pi)) - gamma', D1 - (mp.log(mp.mpf(27) / (4 * PI)) - EG), 1e-28)
TABLE_D = {'0.001': '5.49491150668117', '0.01': '3.21247396010240', '0.1': '1.10270616884989',
           '0.5': '0.204702422036445', '1': '0.187596954133505', '1.5': '0.0867447359919386',
           '2': '0.128705436305314', '2.999': '0.0797441945246622', '3': '0.0798007375292860',
           '3.001': '0.0797419206506647', '5': '0.0515093215827293', '6': '0.0795872835034556',
           '10': '0.0663954730327453'}
for dstr, val in TABLE_D.items():
    quoted('D(%s)' % dstr, D_closed(dstr), val)
quoted('D(0.1) short', D_closed('0.1'), '1.10271')
quoted('D(0.5) short', D_closed('0.5'), '0.20470')
quoted('D(1) short', D_closed('1'), '0.18760')
# delta = 2: the prime-free formula against the full one
quoted('delta = 2: prime-free formula', D_free(2), '0.4753')
quoted('delta = 2: full closed form', D_closed(2), '0.1287')
# prime-power term values
for dstr, val in (('1.5', '-0.231'), ('2', '-0.347'), ('3', '-0.706'), ('5', '-1.261'), ('10', '-2.162')):
    quoted('prime term -(2/delta) P(1+delta) at %s' % dstr, -(2 / mp.mpf(dstr)) * Psum(1 + mp.mpf(dstr)), val)
# 1 < delta <= 2: prime term = -((delta-1)/delta) ln 2 ; 2 < delta <= 3 adds -(2/delta) ln 3 (delta-2)/3
for dstr in ('1.3', '1.7', '2'):
    d = mp.mpf(dstr)
    small('1 < delta <= 2 prime term = -((delta-1)/delta) ln 2 at %s' % dstr,
          D_closed(d) - (D_free(d) - ((d - 1) / d) * mp.log(2)), 1e-28)
for dstr in ('2.5', '3'):
    d = mp.mpf(dstr)
    small('2 < delta <= 3 prime term at %s' % dstr,
          D_closed(d) - (D_free(d) - ((d - 1) / d) * mp.log(2) - (2 / d) * mp.log(3) * (d - 2) / 3), 1e-28)
# fact 2048 form: delta D = (2+delta) C - 2 [G(1+delta) - ln 2pi + T(1+delta)]
for dstr in ('0.3', '1', '2', '3.7', '10', '57.2'):
    d = mp.mpf(dstr)
    rhs = ((2 + d) * CC - 2 * (Gfun(1 + d) - mp.log(2 * PI) + Tfun(1 + d))) / d
    small('fact-2048 form equals the closed form at delta = %s' % dstr, rhs - D_closed(d), 1e-26)
# shape on (0, 1]
phi = lambda d: d ** 2 - (1 + d) * mp.log1p(d / 2)
d0 = mp.findroot(phi, mp.mpf('0.74'))
quoted('delta_0 (zero of phi)', d0, '0.7403218')
quoted('min D on (0,1] = D(delta_0)', D_closed(d0), '0.1594331')
Dp = lambda d: 2 / (1 + d) - (2 / d ** 2) * mp.log1p(d / 2)
small("D'(0.3) formula against numerical derivative", mp.diff(D_closed, mp.mpf('0.3')) - Dp(mp.mpf('0.3')), 1e-20)
quoted("D'(1-)", Dp(mp.mpf(1)), '0.18907')
quoted("D'(1-) = 1 - 2 ln(3/2)", 1 - 2 * mp.log(mp.mpf(3) / 2), '0.18907')
Dp_plus = mp.diff(D_closed, mp.mpf(1), direction=1)
quoted("D'(1+) (one-sided numerical derivative)", Dp_plus, '-0.50408')
small("D'(1+) - D'(1-) = -ln 2", Dp_plus - Dp(mp.mpf(1)) + mp.log(2), 1e-15)
# slope jump -2 Lambda(n)/(n(n-1)) at delta = n - 1 for n = 3, 4, 5
for n in (3, 4, 5):
    dl = mp.diff(D_closed, mp.mpf(n - 1), direction=-1)
    dr = mp.diff(D_closed, mp.mpf(n - 1), direction=1)
    small('slope jump at delta = %d equals -2 Lambda(%d)/(%d*%d)' % (n - 1, n, n, n - 1),
          (dr - dl) + 2 * LAM[n] / (n * (n - 1)), 1e-15)
# grid on (0, 60], step 0.05: minimum; RH bounds on (0, 100], step 0.05
grid = [mp.mpf(k) / 20 for k in range(1, 2001)]
vals = [D_closed(x) for x in grid]
i60 = 1200
m60 = min(vals[:i60])
quoted('min of D on the grid (0, 60] step 0.05', m60, '0.04320')
quoted('... attained at delta', grid[vals[:i60].index(m60)], '33.5')
check('the grid minimum on (0, 60] lies below C', m60 < CC)
lo_gap = min(v - CC * (mp.sqrt(1 + x) - 1) ** 2 / x for x, v in zip(grid, vals))
hi_gap = min(CC * (mp.sqrt(1 + x) + 1) ** 2 / x - v for x, v in zip(grid, vals))
check('C(sqrt(1+d)-1)^2 <= d D(d) <= C(sqrt(1+d)+1)^2 on the grid (0, 100] step 0.05 (2000 points)',
      lo_gap > 0 and hi_gap > 0, 'min margins %s, %s' % (mp.nstr(lo_gap, 4), mp.nstr(hi_gap, 4)))
check('(1 + |l|)^2 <= 8(1/4 + l^2): 7l^2 - 2l + 1 has discriminant 4 - 28 < 0', 4 - 28 < 0)

# ================================================================ 2. exact algebra and series (sympy)
sec('2. exact algebra and series (sympy)')
dS = sp.symbols('delta', positive=True)
gS, g1S, PS, lq = sp.symbols('gamma gamma1 P lq', real=True)   # P = P(1 + delta), lq = ln(q/pi)
aS, bS = sp.symbols('a b', real=True)                          # L'/L(1, chi) = a + i b
ld, l1, l2d, l2, lpi = sp.symbols('ld l1 l2d l2 lpi', real=True)  # ln delta, ln(1+delta), ln(2+delta), ln 2, ln pi
# closed forms of T and T_o at y = 1 + delta, written in the five logarithms:
#   ln(1 - y^-2) = ln delta + ln(2+delta) - 2 ln(1+delta),   ln((y+1)/(y-1)) = ln(2+delta) - ln delta
yS = 1 + dS
T_y = 1 - (ld + l2d - 2 * l1) / 2 - (yS / 2) * (l2d - ld)
To_y = (l2d - ld) / 2 + (yS / 2) * (ld + l2d - 2 * l1)
T_1, To_1 = 1 - l2, l2
# zeta: J(y) = P(y) - y ln y + (1+gamma) y - ln 2pi + T(y), P(1) = 0
J1z = (1 + gS) - (l2 + lpi) + T_1
Jyz = PS - yS * l1 + (1 + gS) * yS - (l2 + lpi) + T_y
target_z = dS * (-ld - gS - l2 - lpi + 2 * l1) + (2 + dS) * (l2d - l2) - 2 * PS
check('Riesz route, zeta: (2+d) J(1) - 2 J(1+d) = d D(d) identically (sympy, symbolic)',
      sp.expand((2 + dS) * J1z - 2 * Jyz - target_z) == 0)
# odd chi: J(y) = P(y) + y L1 - L0 + T_o(y),  L0 = -ln(q/pi) + gamma + ln 2 - conj(L1)
L1 = aS + sp.I * bS
L0 = -lq + gS + l2 - sp.conjugate(L1)
J1o = L1 - L0 + To_1
ReJyo = PS + yS * aS - sp.re(L0) + To_y
target_o = dS * (lq - l2 - ld - gS) + 2 * (1 + dS) * l1 - (2 + dS) * (l2d - l2) - 2 * PS
check('Riesz route, odd chi: J(1) is real', sp.simplify(sp.im(J1o)) == 0)
check('Riesz route, odd chi: (2+d) J(1) - 2 Re J(1+d) = d D_chi(d) identically, L\'/L(1) drops out',
      sp.expand((2 + dS) * sp.re(J1o) - 2 * ReJyo - target_o) == 0)
# even chi: J(y) = P(y) + y L1 - ln y - 1 - B0 + T(y),  B0 = -ln(q/pi) + gamma + ln 2 - conj(L1)
B0 = -lq + gS + l2 - sp.conjugate(L1)
J1e = L1 - 1 - B0 + T_1
ReJye = PS + yS * aS - l1 - 1 - sp.re(B0) + T_y
target_e = dS * (lq - l2 - ld - gS) + (2 + dS) * (l2d - l2) - 2 * PS
check('Riesz route, even chi: J(1) is real', sp.simplify(sp.im(J1e)) == 0)
check('Riesz route, even chi: (2+d) J(1) - 2 Re J(1+d) = d D_chi(d) identically',
      sp.expand((2 + dS) * sp.re(J1e) - 2 * ReJye - target_e) == 0)
check('J(1) for odd chi equals 2 Re L\'/L(1) + ln(q/pi) - gamma',
      sp.expand(sp.re(J1o) - (2 * aS + lq - gS)) == 0)
# the closed forms of T and T_o against their defining series
wS = sp.symbols('w', positive=True)
Tcl = 1 - sp.log(1 - wS ** 2) / 2 - sp.atanh(wS) / wS
Tocl = sp.atanh(wS) + sp.log(1 - wS ** 2) / (2 * wS)
serT = sp.series(Tcl, wS, 0, 14).removeO()
serTo = sp.series(Tocl, wS, 0, 14).removeO()
check('T(y) = 1 - (1/2) ln(1 - w^2) - artanh(w)/w, w = 1/y: coefficients of w^2 .. w^12 are 1/(2k(2k+1))',
      all(serT.coeff(wS, 2 * k) == sp.Rational(1, 2 * k * (2 * k + 1)) for k in range(1, 7)))
check('T_o(y) = artanh(w) + ln(1 - w^2)/(2w): coefficients of w .. w^13 are 1/((2k-1)2k)',
      all(serTo.coeff(wS, 2 * k - 1) == sp.Rational(1, (2 * k - 1) * 2 * k) for k in range(1, 8)))
y3 = mp.mpf(3)
small('T(3) closed form against the series', Tfun(3) - mp.nsum(lambda k: y3 ** (-2 * k) / (2 * k * (2 * k + 1)), [1, mp.inf]), 1e-28)
small('T(1) = 1 - ln 2 (series)', mp.nsum(lambda k: 1 / (2 * k * (2 * k + 1)), [1, mp.inf]) - (1 - mp.log(2)), 1e-25)
small('T_o(1) = ln 2 (series)', mp.nsum(lambda k: 1 / ((2 * k - 1) * 2 * k), [1, mp.inf]) - mp.log(2), 1e-25)
# the Riesz kernels: residues of z^s/(s(s-1)) and z^s/(s(s-1))^2
sS, zS = sp.symbols('s z', positive=True)
r1 = sp.residue(zS ** sS / (sS * (sS - 1)), sS, 1) + sp.residue(zS ** sS / (sS * (sS - 1)), sS, 0)
check('kernel 1/(s(s-1)): residues at 1 and 0 give z - 1', sp.simplify(r1 - (zS - 1)) == 0)
r2 = sp.residue(zS ** sS / (sS * (sS - 1)) ** 2, sS, 1) + sp.residue(zS ** sS / (sS * (sS - 1)) ** 2, sS, 0)
check('kernel 1/(s(s-1))^2: residues give (z+1) ln z - 2(z-1)', sp.simplify(sp.expand_log(r2 - ((zS + 1) * sp.log(zS) - 2 * (zS - 1)), force=True)) == 0)
# residue at s = 1 of (-zeta'/zeta)(s) y^s/(s(s-1)) with -zeta'/zeta = 1/(s-1) - gamma + (2 gamma1 + gamma^2)(s-1)
yy = sp.symbols('y', positive=True)
mz = 1 / (sS - 1) - gS + (2 * g1S + gS ** 2) * (sS - 1)
res1 = sp.residue(sp.series(mz * yy ** sS / (sS * (sS - 1)), sS, 1, 3).removeO(), sS, 1)
check('residue at s = 1 (kernel 1/(s(s-1))): y ln y - (1 + gamma) y', sp.simplify(res1 - (yy * sp.log(yy) - (1 + gS) * yy)) == 0)
res2 = sp.residue(sp.series(mz * yy ** sS / (sS * (sS - 1)) ** 2, sS, 1, 4).removeO(), sS, 1)
check('residue at s = 1 (kernel 1/(s(s-1))^2): y[ln^2 y/2 - (2+gamma) ln y + 3 + 2gamma + gamma^2 + 2gamma1]',
      sp.simplify(sp.expand(res2 - yy * (sp.log(yy) ** 2 / 2 - (2 + gS) * sp.log(yy) + 3 + 2 * gS + gS ** 2 + 2 * g1S))) == 0)
# Laurent data used above, numerically
eps = mp.mpf('1e-6')
lz = lambda s: mp.zeta(s, derivative=1) / mp.zeta(s)
small("zeta'/zeta(s) = -1/(s-1) + gamma - (2 gamma1 + gamma^2)(s-1) + O((s-1)^2) at s = 1 + 1e-6",
      lz(1 + eps) - (-1 / eps + EG - (2 * G1 + EG ** 2) * eps), 1e-10)
small("zeta'/zeta(0) = ln 2pi", lz(0) - mp.log(2 * PI), 1e-25)
small("(zeta'/zeta)'(0) = pi^2/12 - 2 gamma1 - gamma^2", mp.diff(lz, 0) - (PI ** 2 / 12 - 2 * G1 - EG ** 2), 1e-20)
# derivative of D and the convexity of phi
DS = -sp.log(dS) - gS - l2 - lpi + 2 * sp.log(1 + dS) + ((2 + dS) / dS) * sp.log(1 + dS / 2)
check("D'(d) = 2/(1+d) - (2/d^2) ln(1+d/2) = 2 phi/(d^2 (1+d)), phi = d^2 - (1+d) ln(1+d/2)",
      sp.simplify(sp.diff(DS, dS) - 2 * (dS ** 2 - (1 + dS) * sp.log(1 + dS / 2)) / (dS ** 2 * (1 + dS))) == 0)
phiS = dS ** 2 - (1 + dS) * sp.log(1 + dS / 2)
check("phi'' = 2 - 1/(2+d) - 1/(2+d)^2, phi(0) = 0, phi'(0) = -1/2",
      sp.simplify(sp.diff(phiS, dS, 2) - (2 - 1 / (2 + dS) - 1 / (2 + dS) ** 2)) == 0
      and phiS.subs(dS, 0) == 0 and sp.diff(phiS, dS).subs(dS, 0) == sp.Rational(-1, 2))
# series
Dz = -sp.log(dS) - gS - sp.log(2 * sp.pi) + 2 * sp.log(1 + dS) + ((2 + dS) / dS) * sp.log(1 + dS / 2)
Bs = 1 - gS - sp.log(2 * sp.pi)
HLs = -sp.log(dS) - gS - sp.log(2 * sp.pi) + (1 + dS) * sp.log(1 + dS) / dS


def coeffs(expr, n):
    ser = sp.expand(sp.series(expr, dS, 0, n + 1).removeO())
    return [sp.nsimplify(ser.coeff(dS, k)) for k in range(0, n + 1)]


c = coeffs(sp.expand_log(Dz - Bs + sp.log(dS), force=True), 4)
check('D - ln(1/d) - B = (9/4)d - (25/24)d^2 + (65/96)d^3 - (161/320)d^4 + ...',
      c == [0, sp.Rational(9, 4), sp.Rational(-25, 24), sp.Rational(65, 96), sp.Rational(-161, 320)], str(c))
c = coeffs(sp.expand_log(Dz - HLs, force=True), 5)
check('D - HL_exact = (7/4)d - (7/8)d^2 + (19/32)d^3 - (29/64)d^4 + (353/960)d^5 - ...',
      c == [0, sp.Rational(7, 4), sp.Rational(-7, 8), sp.Rational(19, 32), sp.Rational(-29, 64), sp.Rational(353, 960)], str(c))
c = coeffs(((2 + dS) / dS) * sp.log(1 + dS / 2) - (1 + dS) * sp.log(1 + dS) / dS, 3)
check('Gamma part beyond HL_exact: ((2+d)/d) ln(1+d/2) - (1+d) ln(1+d)/d = -d/4 + d^2/8 - 7d^3/96 + ...',
      c == [0, sp.Rational(-1, 4), sp.Rational(1, 8), sp.Rational(-7, 96)], str(c))
c = coeffs(2 * sp.log(1 + dS), 2)
check('pole part: 2 ln(1+d) = 2d - d^2 + ...', c == [0, 2, -1], str(c))
c = coeffs((8 / dS) * (sp.sqrt(1 + dS) - 1) ** 2, 3)
check('(2/d) h(0) = (8/d)(sqrt(1+d) - 1)^2 = 2d - d^2 + (5/8)d^3 - ...', c == [0, 2, -1, sp.Rational(5, 8)], str(c))
c = coeffs(2 * sp.log(1 + dS) - (8 / dS) * (sp.sqrt(1 + dS) - 1) ** 2, 4)
check('S part: 2 ln(1+d) - (8/d)(sqrt(1+d) - 1)^2 = d^3/24 - d^4/16 + ...', c == [0, 0, 0, sp.Rational(1, 24), sp.Rational(-1, 16)], str(c))
c = coeffs((2 / dS) * ((1 + dS) * sp.log(1 + dS) - (2 + dS) * sp.log(1 + dS / 2)), 3)
check('odd minus even: (2/d)[(1+d) ln(1+d) - (2+d) ln(1+d/2)] = d/2 - d^2/4 + 7d^3/48 - ...',
      c == [0, sp.Rational(1, 2), sp.Rational(-1, 4), sp.Rational(7, 48)], str(c))
Pe = ((2 + dS) / dS) * sp.log(1 + dS / 2)
Po = (2 * (1 + dS) / dS) * sp.log(1 + dS) - Pe
c = coeffs(Pe - 1, 2)          # D_even - ln(q/d) - B = P_0 - 1
check('D_chi even - ln(q/d) - B = d/4 - d^2/24 + ...', c == [0, sp.Rational(1, 4), sp.Rational(-1, 24)], str(c))
c = coeffs(Po - 1, 2)
check('D_chi odd - ln(q/d) - B = 3d/4 - 7d^2/24 + ...', c == [0, sp.Rational(3, 4), sp.Rational(-7, 24)], str(c))
check('9/4 = 2 (pole) + 1/4 (even Gamma part)', sp.Rational(9, 4) == 2 + sp.Rational(1, 4))
# D0 - HL_exact: D0 = -ln 2pi + (1/d) int_0^u (d e^-x - g(x))/x dx + E1(u), g = (1+d) e^(-x/2) - e^(x/2), u = ln(1+d)
uS, xS = sp.symbols('u x', positive=True)
KS = 10
dl_u = sp.exp(uS) - 1
integrand = sp.series((dl_u * sp.exp(-xS) - sp.exp(uS) * sp.exp(-xS / 2) + sp.exp(xS / 2)) / xS, xS, 0, KS).removeO()
I1 = sp.integrate(sp.expand(integrand), (xS, 0, uS))
E1reg = sum(sp.Integer(-1) ** (k + 1) * uS ** k / (k * sp.factorial(k)) for k in range(1, KS + 1))  # E1(u) = -gamma - ln u + E1reg
logratio = sp.series(sp.log(dl_u / uS), uS, 0, KS).removeO()                                       # ln((e^u - 1)/u)
diff_u = sp.series(I1 / dl_u + E1reg + logratio - sp.exp(uS) * uS / dl_u, uS, 0, 7).removeO()       # D0 - HL_exact in u
c = coeffs(diff_u.subs(uS, sp.log(1 + dS)), 4)
check('D0 - HL_exact = -d^2/144 + d^3/144 - (703/115200) d^4 + ...',
      c == [0, 0, sp.Rational(-1, 144), sp.Rational(1, 144), sp.Rational(-703, 115200)], str(c))
D0_series = lambda d: -d ** 2 / 144 + d ** 3 / 144 - mp.mpf(703) / 115200 * d ** 4
# the epsilon part of C3: -(2/d) int h' eps = [Gamma part beyond HL_exact] + (1/d)(sqrt(1+d)-1)^2 - (D0 - HL_exact)
eps_part = sp.series(((2 + dS) / dS) * sp.log(1 + dS / 2) - (1 + dS) * sp.log(1 + dS) / dS + (sp.sqrt(1 + dS) - 1) ** 2 / dS, dS, 0, 5).removeO() \
    - diff_u.subs(uS, sp.log(1 + dS))
c = coeffs(eps_part, 3)
check('epsilon part of C3: -(2/d) int_0^inf h\' eps dt = d^2/144 - d^3/576 + ...', c == [0, 0, sp.Rational(1, 144), sp.Rational(-1, 576)], str(c))
check('(7/8)(2/d) h(0) = (7/d)(sqrt(1+d) - 1)^2 = (7/4) d - (7/8) d^2 + ...', coeffs(7 * (sp.sqrt(1 + dS) - 1) ** 2 / dS, 2) == [0, sp.Rational(7, 4), sp.Rational(-7, 8)])
# K2 = 2 sigma2 + 2C and the closed forms
sig2S = 1 + gS ** 2 + 2 * g1S - sp.pi ** 2 / 8
K2S = 6 + 2 * gS + 2 * gS ** 2 + 4 * g1S - sp.pi ** 2 / 4 - 2 * sp.log(4 * sp.pi)
check('K2 = 2 sigma2 + 2C (symbolic)', sp.simplify(sp.expand_log(K2S - 2 * sig2S - 2 * (2 + gS - sp.log(4 * sp.pi)), force=True)) == 0)
K2 = 6 + 2 * EG + 2 * EG ** 2 + 4 * G1 - PI ** 2 / 4 - 2 * mp.log(4 * PI)
SIG2 = 1 + EG ** 2 + 2 * G1 - PI ** 2 / 8
quoted('K2', K2, '7.42012728749297e-5')
quoted('sigma2', SIG2, '-0.0461543172958')
quoted('sqrt(K2)', mp.sqrt(K2), '0.00861')
quoted('sqrt(K2)/C', mp.sqrt(K2) / CC, '0.1865')
T21 = mp.nsum(lambda k: 1 / (4 * k ** 2 * (2 * k + 1) ** 2), [1, mp.inf])
small('T2(1) = pi^2/6 - 3 + 2 ln 2', T21 - (PI ** 2 / 6 - 3 + 2 * mp.log(2)), 1e-25)


def T2fun(y):
    y = mp.mpf(y)
    return mp.nsum(lambda k: y ** (-2 * k) / (4 * k * k * (2 * k + 1) ** 2), [1, mp.inf])


def P2sum(y):
    y = mp.mpf(y)
    return mp.fsum(LAM[n] * ((y / n + 1) * mp.log(y / n) - 2 * (y / n - 1)) for n in range(2, int(mp.floor(y)) + 1) if LAM[n] != 0)


def J2(y):
    y = mp.mpf(y)
    ly = mp.log(y)
    return (y * (ly ** 2 / 2 - (2 + EG) * ly + 3 + 2 * EG + EG ** 2 + 2 * G1) - mp.log(2 * PI) * (ly + 2)
            + 2 * G1 + EG ** 2 - PI ** 2 / 12 - T2fun(y) - P2sum(y))


small('J2(1) = K2', J2(1) - K2, 1e-25)

# ================================================================ 3. the Weil route by quadrature
sec('3. the Weil route: formula (4.3) and the archimedean integral by quadrature')
mp.mp.dps = 25


def g_weil(x, d):
    d = mp.mpf(d)
    x = abs(mp.mpf(x))
    return (1 + d) * mp.exp(-x / 2) - mp.exp(x / 2) if x <= mp.log1p(d) else mp.mpf(0)


def gauss_43(d, alpha):
    """right side of (4.3): -gamma g(0) + 2 int_0^inf [g(0) e^-2x - e^(-2 alpha x) g(x)]/(1 - e^-2x) dx."""
    d = mp.mpf(d)
    u = mp.log1p(d)
    f = lambda x: (d * mp.exp(-2 * x) - mp.exp(-2 * alpha * x) * g_weil(x, d)) / (-mp.expm1(-2 * x))
    return -EG * d + 2 * mp.quad(f, [0, u / 2, u, u + 1, mp.inf])


def arch_closed(d, a):
    d = mp.mpf(d)
    base = -EG * d - d * mp.log(2 * d)
    return base + ((2 + d) * mp.log1p(d / 2) if a == 0 else 2 * (1 + d) * mp.log1p(d) - (2 + d) * mp.log1p(d / 2))


worst43 = 0
for dstr in ('0.01', '0.1', '0.5', '1', '1.5', '3', '10'):
    for a, alpha in ((0, mp.mpf(1) / 4), (1, mp.mpf(3) / 4)):
        worst43 = max(worst43, abs(gauss_43(dstr, alpha) - arch_closed(dstr, a)))
small('(4.3) as an x-integral against its closed form, alpha = 1/4 and 3/4, seven delta in [0.01, 10]', worst43, 1e-20)
small('g(0) = delta and g(ln(1+delta)) = 0 at delta = 0.7', abs(g_weil(0, '0.7') - mp.mpf('0.7')) + abs((1 + mp.mpf('0.7')) * mp.exp(-mp.log1p(mp.mpf('0.7')) / 2) - mp.exp(mp.log1p(mp.mpf('0.7')) / 2)), 1e-24)
# h = Fourier transform of g, at three r
for r in ('0', '2.5', '17'):
    d = mp.mpf('0.6')
    u = mp.log1p(d)
    ft = 2 * mp.quad(lambda x: g_weil(x, d) * mp.cos(mp.mpf(r) * x), [0, u])
    small('h(r) = int g(x) e^{irx} dx at delta = 0.6, r = %s' % r, ft - h_mp(mp.mpf(r), d), 1e-22)
worstW = 0
for dstr in ('0.01', '0.1', '1', '3'):
    d = mp.mpf(dstr)
    arch = split_integral(d, lambda r: mp.re(mp.digamma(mp.mpf(1) / 4 + 0.5j * r)), 0) / PI
    weil = 2 * d * mp.log1p(d) - d * mp.log(PI) + arch - 2 * Psum(1 + d)
    worstW = max(worstW, abs(weil / d - D_closed(d)))
small('Weil formula for zeta with the digamma integral done by quadrature in r, delta = 0.01, 0.1, 1, 3', worstW, 1e-17)
worstWo = 0
for dstr in ('0.1', '1'):
    d = mp.mpf(dstr)
    arch = split_integral(d, lambda r: mp.re(mp.digamma(mp.mpf(3) / 4 + 0.5j * r)), 0) / PI
    worstWo = max(worstWo, abs(arch - arch_closed(d, 1)))
small('odd parity: digamma integral at alpha = 3/4 by quadrature in r, delta = 0.1, 1', worstWo, 1e-17)

# ================================================================ 4. zeta: closed form against 40,000 zeros
sec('4. zeta: the closed form against 40,000 zeros')
Z = np.loadtxt(ZFILE)
check('zero file: 40000 increasing ordinates', len(Z) == 40000 and bool(np.all(np.diff(Z) > 0)))
worstz = 0
for n in (1, 2, 1000, 20000, 40000):
    zz = mp.zetazero(n)
    worstz = max(worstz, abs(float(zz.imag) - Z[n - 1]))
small("mpmath zetazero(n), n = 1, 2, 1000, 20000, 40000, against the file", worstz, 1e-12)
th = np.array([float(mp.siegeltheta(mp.mpf(z))) for z in Z[49::50]])
Sv = np.arange(50, len(Z) + 1, 50) - 1 - th / math.pi
quoted('S(gamma_n+) at every 50th zero: mean', Sv.mean(), '0.50')
quoted('... mean over the first half', Sv[:len(Sv) // 2].mean(), '0.51')
quoted('... mean over the second half', Sv[len(Sv) // 2:].mean(), '0.49')
TZ = mp.mpf(Z[-1])
quoted('T = gamma_40000', TZ, '33190.012')


def zero_side(d, N):
    d = mp.mpf(d)
    df = float(d)
    T = mp.mpf(Z[N - 1])
    main = math.fsum(h_np(Z[:N], df).tolist())
    smooth = split_integral(d, lambda t: theta_prime(t) / PI, T)
    S_T = N - 1 - mp.siegeltheta(T) / PI
    Dz = (2 / d) * (mp.mpf(main) + smooth - h_mp(T, d) * S_T)
    u = mp.log1p(d)
    s = mp.sqrt(1 + d)
    A = 2 + d
    # Trudgian: |int_t1^t2 S| <= 0.059 ln t2 + 2.067 (t2 > t1 > 168 pi); one more integration by parts, |h''| bound
    a_, b_ = mp.mpf('0.059'), mp.mpf('2.067')
    I = lambda k: T ** (1 - k) * ((a_ * mp.log(T) + b_) / (k - 1) + a_ / (k - 1) ** 2)
    bound1 = (2 / d) * (2 * s * u ** 2 * I(2) + 8 * s * u * I(3) + 6 * (A + 2 * s) * I(4))
    # Rosser: |N(t) - F(t)| < 0.137 ln t + 0.443 ln ln t + 1.588, |S| <= that + 1/(40 t); |h'| <= 2su/t^2 + 2(A+2s)/t^3
    Bt = lambda t: mp.mpf('0.137') * mp.log(t) + mp.mpf('0.443') * mp.log(mp.log(t)) + mp.mpf('1.588') + 1 / (40 * t)
    bound2 = (2 / d) * mp.quad(lambda t: (2 * s * u / t ** 2 + 2 * (A + 2 * s) / t ** 3) * Bt(t), [T, 2 * T, 10 * T, mp.inf])
    return Dz, bound1, bound2, S_T


def nearest_pp(y):
    n = int(mp.nint(y))
    return n if n < len(LAM) and LAM[n] != 0 else None


ZROWS = [  # delta, N, quoted zeros - closed, quoted bound S1, quoted bound S, resonance n, quoted prediction, quoted residual
    ('0.001', 40000, '+1.2e-9', '3.5e-7', '5.4e-4', None, None, None),
    ('0.01', 40000, '-5.1e-10', '3.3e-6', '5.1e-4', None, None, None),
    ('0.1', 40000, '+2.3e-10', '3.1e-5', '5.1e-4', None, None, None),
    ('0.5', 40000, '-1.6e-9', '1.3e-4', '5.1e-4', None, None, None),
    ('1', 10000, '-4.47e-5', None, None, 2, '-4.47e-5', None),
    ('1', 20000, '-2.45e-5', None, None, 2, '-2.45e-5', None),
    ('1', 40000, '-1.330e-5', '2.2e-4', '5.0e-4', 2, '-1.330e-5', '-1.6e-10'),
    ('1.5', 40000, '-8.1e-10', '2.9e-4', '4.9e-4', None, None, None),
    ('2', 40000, '-1.054e-5', '3.5e-4', '4.8e-4', 3, '-1.054e-5', '-6.3e-10'),
    ('2.999', 40000, '+4.98e-7', '4.2e-4', '4.7e-4', 4, '+4.97e-7', '+1.6e-9'),
    ('3', 40000, '-4.43e-6', '4.2e-4', '4.7e-4', 4, '-4.43e-6', '-1.9e-9'),
    ('3.001', 40000, '+4.97e-7', '4.2e-4', '4.7e-4', 4, '+4.97e-7', '+6.3e-11'),
    ('5', 40000, '+1.4e-9', '5.2e-4', '4.5e-4', None, None, None),
    ('6', 40000, '-6.22e-6', '5.5e-4', '4.4e-4', 7, '-6.22e-6', '+2.2e-10'),
    ('10', 40000, '-4.60e-6', '6.3e-4', '4.1e-4', 11, '-4.60e-6', '+5.8e-10'),
]
for dstr, N, q_diff, q_b1, q_b2, npp, q_pred, q_res in ZROWS:
    d = mp.mpf(dstr)
    Dz, b1, b2, S_T = zero_side(d, N)
    diff = Dz - D_closed(d)
    tag = 'delta %s, N %d' % (dstr, N)
    quoted('zeros - closed, ' + tag, diff, q_diff)
    if q_b1:
        quoted('bound from Trudgian S1, ' + tag, b1, q_b1)
        quoted('bound from Rosser S, ' + tag, b2, q_b2)
        check('|zeros - closed| inside both bounds, ' + tag, abs(diff) <= b1 and abs(diff) <= b2)
    if npp:
        w = mp.log1p(d) - mp.log(npp)
        T = mp.mpf(Z[N - 1])
        pred = -(2 * LAM[npp] / (PI * d)) * resonance_I(w, T)
        quoted('predicted resonance -(2 Lambda(n)/(pi delta)) int_T^inf cos(wt)/t^2, n = %d, %s' % (npp, tag), pred, q_pred)
        if q_res:
            quoted('residual after the resonance term, ' + tag, diff - pred, q_res)
        if dstr == '2.999':
            quoted('|w| T at delta = 2.999, n = 4', abs(w) * T, '8.3')
        if w == 0:
            quoted('(zeros - closed) * T / (2 Lambda(n)/(pi delta)) at exact resonance, ' + tag, diff * T / (2 * LAM[npp] / (PI * d)), '-1.00')
    else:
        small('non-resonant row: |zeros - closed| <= 2e-9, ' + tag, diff, 2e-9)

# ================================================================ 5. zeta: K2, sigma2, C and Dt against the zeros
sec('5. zeta: sums over zeros with the second Riesz kernel')
ST40 = len(Z) - 1 - mp.siegeltheta(TZ) / PI
thp_np_tail = lambda f: mp.quad(lambda t: f(t) * theta_prime(t) / PI, [TZ, 2 * TZ, 10 * TZ, mp.inf])
for name, f, f_np, exact, tol in (
        ('C = sum 1/|rho|^2', lambda t: 2 / (mp.mpf(1) / 4 + t * t), lambda t: 2 / (0.25 + t * t), CC, 1e-12),
        ('sigma2 = sum rho^-2', lambda t: 2 * (mp.mpf(1) / 4 - t * t) / (mp.mpf(1) / 4 + t * t) ** 2,
         lambda t: 2 * (0.25 - t * t) / (0.25 + t * t) ** 2, SIG2, 1e-12),
        ('K2 = sum 1/|rho|^4', lambda t: 2 / (mp.mpf(1) / 4 + t * t) ** 2, lambda t: 2 / (0.25 + t * t) ** 2, K2, 1e-18)):
    zs = mp.mpf(math.fsum(f_np(Z).tolist())) + thp_np_tail(f) - f(TZ) * ST40
    small('%s from 40,000 zeros + smooth tail + boundary term, against the closed form' % name, zs - exact, tol)
Dt_worst = 0
for dstr in ('0.01', '0.1', '0.5', '1', '1.5', '3', '10'):
    d = mp.mpf(dstr)
    df = float(d)
    closed = ((2 + d) * K2 - 2 * J2(1 + d)) / d
    zsum = math.fsum((h_np(Z, df) / (0.25 + Z * Z)).tolist())
    tail = split_integral(d, lambda t: theta_prime(t) / PI / (mp.mpf(1) / 4 + t * t), TZ)
    Dtz = (2 / d) * (mp.mpf(zsum) + tail - h_mp(TZ, d) / (mp.mpf(1) / 4 + TZ * TZ) * ST40)
    Dt_worst = max(Dt_worst, abs(Dtz - closed))
    if dstr == '1':
        quoted('Dt zeros - closed at delta = 1', Dtz - closed, '-4.0e-15')
        quoted('resonance -2 Lambda(2)/(3 pi delta T^3) at delta = 1', -2 * mp.log(2) / (3 * PI * TZ ** 3), '-4.0e-15')
    if dstr == '0.01':
        quoted('Dt(0.01)/0.01', closed / d, '0.0393')
small('Dt(delta) = [(2+delta) K2 - 2 J2(1+delta)]/delta against the zeros, delta = 0.01 ... 10 (7 values)', Dt_worst, 5e-15)

# ================================================================ 6. D0 - HL_exact
sec('6. the smooth-density diagonal D0 against HL_exact')
for dstr, q_t, q_s3 in (('0.1', '-6.306e-5', '-6.311e-5'), ('0.05', '-1.653e-5', '-1.653e-5'), ('0.02', '-2.723e-6', '-2.723e-6')):
    d = mp.mpf(dstr)
    u = mp.log1p(d)
    D0_t = (2 / d) * split_integral(d, lambda t: mp.log(t / (2 * PI)) / (2 * PI), 0)
    D0_x = -mp.log(2 * PI) + mp.quad(lambda x: (d * mp.exp(-x) - g_weil(x, d)) / x, [0, u]) / d + mp.e1(u)
    small('D0 by the t-integral (period by period) and by the x-integral agree, delta = %s' % dstr, D0_t - D0_x, 1e-15)
    quoted('D0 - HL_exact at delta = %s' % dstr, D0_x - HL_exact(d), q_t)
    quoted('-d^2/144 + d^3/144 - 703 d^4/115200 at delta = %s' % dstr, D0_series(d), q_s3)
    quoted('-d^2/144 at delta = %s' % dstr, -d ** 2 / 144, {'0.1': '-6.944e-5', '0.05': '-1.736e-5', '0.02': '-2.778e-6'}[dstr])

# ================================================================ 7. characters
sec('7. characters: zeros from PARI lfunzeros, closed forms, resonance')
import cypari2
pari = cypari2.Pari()
pari.allocatemem(300 * 10 ** 6, silent=True)
pari.set_real_precision(30)
TMAXC = 1000
W3 = mp.exp(2j * PI / 3)


def p2mp(z):
    """PARI number -> mpmath mpc (via strings, full precision)."""
    rs = lambda x: mp.mpf(str(x).replace(' E', 'e').replace('E', 'e'))
    return mp.mpc(rs(pari.real(z)), rs(pari.imag(z)))


def p2frac(z):
    from fractions import Fraction
    fr = Fraction(str(z))
    return mp.mpf(fr.numerator) / fr.denominator
CHARS = [  # name, q, a, PARI L(chi), PARI L(conj chi) or None (real), chi(n) for n mod q
    ('chi4 (odd, mod 4)', 4, 1, 'lfuncreate(-4)', None, {0: 0, 1: 1, 2: 0, 3: -1}),
    ('chi3 (odd, mod 3)', 3, 1, 'lfuncreate(-3)', None, {0: 0, 1: 1, 2: -1}),
    ('chi5 = (./5) (even)', 5, 0, 'lfuncreate(5)', None, {0: 0, 1: 1, 2: -1, 3: -1, 4: 1}),
    ('chi5c (odd, mod 5, chi(2) = i)', 5, 1, 'lfuncreate([znstar(5,1),[1]])', 'lfuncreate([znstar(5,1),[3]])',
     {0: 0, 1: 1, 2: 1j, 3: -1j, 4: -1}),
    ('chi7c (even, mod 7, order 3, chi(3) = e^(2 pi i/3))', 7, 0, 'lfuncreate([znstar(7,1),[2]])', 'lfuncreate([znstar(7,1),[4]])',
     {0: 0, 1: 1, 3: W3, 2: W3 ** 2, 6: 1, 4: W3, 5: W3 ** 2}),
]
# character tables and primitivity, checked with PARI
for name, q, a, cmd, cmdbar, tab in CHARS[3:]:
    Gs = 'znstar(%d,1)' % q
    chi_s = cmd.split(',', 1)[1].rstrip(')').rstrip(']').split('[')[-1]
    vals_ok = all(abs(mp.exp(2j * PI * (p2frac(pari('chareval(%s,[%s],%d)' % (Gs, chi_s, n))))) - tab[n % q]) < 1e-20
                  for n in range(1, q) if math.gcd(n, q) == 1)
    cond = int(pari('zncharconductor(%s,[%s])' % (Gs, chi_s)))
    odd = int(pari('zncharisodd(%s,[%s])' % (Gs, chi_s)))
    check('%s: table matches PARI chareval, conductor %d (primitive), parity %d' % (name, cond, odd), vals_ok and cond == q and odd == a)
REF = {  # closed values D_chi(delta) quoted in the note, delta = 0.01, 0.1, 0.5, 1, 1.5, 2.5, 3.5
    'chi4 (odd, mod 4)': ['4.5838428044', '2.3460172309', '0.9814217024', '0.5273950277', '0.3142687894', '0.2494544590', '0.3102969029'],
    'chi3 (odd, mod 3)': ['4.2961607320', '2.0583351585', '0.6937396300', '0.2397129553', '0.2576357772', '0.2311790564', '0.1543202212'],
    'chi5 = (./5) (even)': ['4.8020112108', '2.5215237217', '1.0032101183', '0.4107405054', '0.3256993051', '0.3000987853', '0.2907246295'],
    'chi5c (odd, mod 5, chi(2) = i)': ['4.8069863557', '2.5691607823', '1.2045652538', '0.7505385790', '0.5374123407', '0.3261163718', '0.2690617417'],
    'chi7c (even, mod 7, order 3, chi(3) = e^(2 pi i/3))': ['5.1384834475', '2.8579959583', '1.3396823549', '0.7472127421', '0.5466470116', '0.3553860485', '0.2969654583'],
}
DCHI = ('0.01', '0.1', '0.5', '1', '1.5', '2.5', '3.5')
QDIFF = {}   # computed zeros - closed per character and delta
QDIFF_QUOTED = {  # as in the table of the note
    'chi4 (odd, mod 4)': ['+1.7e-7', '-1.5e-7', '-1.8e-7', '-3.5e-7', '+4.8e-7', '-6.1e-7', '-2.9e-6'],
    'chi3 (odd, mod 3)': ['-5.8e-7', '+4.9e-7', '+1.1e-6', '+4.42e-4', '-1.1e-6', '+2.0e-6', '+1.7e-6'],
    'chi5 = (./5) (even)': ['-5.2e-7', '+4.6e-7', '+7.6e-7', '+4.42e-4', '-8.2e-7', '+1.6e-6', '-9.7e-7'],
    'chi5c (odd, mod 5, chi(2) = i)': ['-4.1e-7', '+4.0e-7', '-6.1e-8', '+3.6e-7', '-2.5e-7', '-6.2e-7', '+1.7e-6'],
    'chi7c (even, mod 7, order 3, chi(3) = e^(2 pi i/3))': ['-3.8e-7', '+3.3e-7', '+5.4e-7', '+2.21e-4', '-4.1e-7', '+2.0e-7', '+2.0e-6'],
}
QRES_QUOTED = {'chi3 (odd, mod 3)': '+9.7e-7', 'chi5 = (./5) (even)': '+8.1e-7', 'chi7c (even, mod 7, order 3, chi(3) = e^(2 pi i/3))': '+3.3e-7'}


def theta_chi(t, q, a):
    return (t / 2) * mp.log(q / PI) + mp.im(mp.loggamma((mp.mpf(1) / 2 + a + 1j * t) / 2))


def thetap_chi(t, q, a):
    return mp.log(q / PI) / 2 + mp.re(mp.digamma((mp.mpf(1) / 2 + a + 1j * t) / 2)) / 2


def arg_half(L):
    """arg L(1/2) by continuous variation along the real axis from sigma = 30 (where arg ~ 0) down to 1/2."""
    sig = [mp.mpf(30) - k * (mp.mpf(30) - mp.mpf(1) / 2) / 240 for k in range(241)]
    vals = [p2mp(pari.lfun(L, pari(mp.nstr(s, 25)))) for s in sig]
    ang = mp.arg(vals[0])
    for v0, v1 in zip(vals[:-1], vals[1:]):
        step = mp.arg(v1 / v0)
        assert abs(step) < 0.5
        ang += step
    return ang, vals[-1]


def zeros_side(L, q, a, argh):
    """ordinates in (0, TMAXC] of L, blockwise; blocks whose mean S(gamma_n+) drops below 0.3 are recomputed
    with a finer subdivision (PARI can miss close pairs).  S(t) = N(t) - theta(t)/pi + arg L(1/2)/pi."""
    Li = pari.lfuninit(L, [TMAXC + 5])
    blocks = {b: [float(z) for z in pari.lfunzeros(Li, [b, b + 100])] for b in range(0, TMAXC, 100)}
    divz = {b: 8 for b in blocks}
    def collect():
        zs = np.array(sorted(set(round(z, 10) for b in blocks for z in blocks[b])))
        S = np.arange(1, len(zs) + 1) - np.array([float(theta_chi(mp.mpf(z), q, a)) for z in zs]) / math.pi + float(argh) / math.pi
        return zs, S
    zs, S = collect()
    for it in range(6):
        bad = [b for b in sorted(blocks) if np.mean(S[(zs > b) & (zs <= b + 100)]) < 0.3]
        if not bad:
            break
        for b in (bad[0] - 100, bad[0]):
            if b in blocks:
                divz[b] *= 4
                blocks[b] = [float(z) for z in pari.lfunzeros(Li, [b, b + 100], divz[b])]
        zs, S = collect()
    redone = sum(1 for b in divz if divz[b] > 8)
    return zs, S, redone


def tail_chi(d, T, q, a):
    return split_integral(d, lambda t: thetap_chi(t, q, a) / PI, T)


worst_nonres = {}
for name, q, a, cmd, cmdbar, tab in CHARS:
    t0 = time.time()
    chi = lambda n, tab=tab, q=q: tab[n % q]
    sides = []
    for c in ([cmd] if cmdbar is None else [cmd, cmdbar]):
        L = pari(c)
        argh, Lh = arg_half(L)
        zs, S, redone = zeros_side(L, q, a, argh)
        sides.append((zs, S, argh, redone))
    real = cmdbar is None
    check('%s: blocks recomputed with a finer subdivision: %d (the note: two for chi5c, none for the others)'
          % (name, sum(s_[3] for s_ in sides)), sum(s_[3] for s_ in sides) == (2 if name.startswith('chi5c') else 0))
    check('%s: %s zeros up to height %d (%d blocks of 100 recomputed with finer subdivision; %.0f s)'
          % (name, ' + '.join(str(len(s[0])) for s in sides) + ('' if not real else ' (x2, symmetric)'), TMAXC,
             sum(s[3] for s in sides), time.time() - t0), True)
    for k, (zs, S, argh, redone) in enumerate(sides):
        blocks_mean = [np.mean(S[(zs > b) & (zs <= b + 100)]) for b in range(0, TMAXC, 100)]
        check('%s side %d: arg L(1/2)/pi = %+.3f; mean S(gamma_n+) = %.3f, every block of height 100 in [0.3, 0.7]'
              % (name, k + 1, float(argh) / math.pi, S.mean(), ), all(0.3 <= m <= 0.7 for m in blocks_mean) and abs(S.mean() - 0.5) < 0.05)
        RvM = theta_chi(mp.mpf(TMAXC), q, a) / PI - argh / PI
        check('%s side %d: count %d against theta(T)/pi - arg L(1/2)/pi = %.2f (|S(T)| < 1.5)' % (name, k + 1, len(zs), float(RvM)),
              abs(len(zs) - RvM) < 1.5)
    if not real:
        quoted('%s: arg L(1/2, chi)/pi' % name, sides[0][2] / PI, {5: '0.088', 7: '-0.073'}[q])
    # J_chi(1) = sum_rho 1/(rho(1-rho)) against 2 Re L'/L(1) + ln(q/pi) - gamma (odd) or ... - gamma - 2 ln 2 (even)
    L = pari(cmd)
    L1v = p2mp(pari.lfun(L, 1, 1)) / p2mp(pari.lfun(L, 1))
    Cchi = 2 * mp.re(L1v) + mp.log(q / PI) - EG - (2 * mp.log(2) if a == 0 else 0)
    Jz = 0
    for zs, S, argh, redone in (sides if not real else sides * 2):
        Tl = mp.mpf(zs[-1])
        Sl = len(zs) - theta_chi(Tl, q, a) / PI + argh / PI
        Jz += math.fsum((1 / (0.25 + zs * zs)).tolist()) + mp.quad(lambda t: thetap_chi(t, q, a) / PI / (mp.mpf(1) / 4 + t * t), [Tl, 2 * Tl, 10 * Tl, mp.inf]) \
            - Sl / (mp.mpf(1) / 4 + Tl * Tl)
    small('%s: sum over zeros of 1/|rho|^2 against the residue formula for J_chi(1) = %s' % (name, mp.nstr(Cchi, 8)), Jz - Cchi, 1e-9)
    # functional equation at s = 0
    if a == 1:
        lhs = p2mp(pari.lfun(L, 0, 1)) / p2mp(pari.lfun(L, 0))
    else:
        lhs = p2mp(pari.lfun(L, 0, 2)) / (2 * p2mp(pari.lfun(L, 0, 1)))       # B0 = L''(0)/(2 L'(0))
    small("%s: %s + conj L'/L(1) = -ln(q/pi) + gamma + ln 2" % (name, "L'/L(0)" if a == 1 else 'B0'),
          abs(lhs + mp.conj(L1v) - (-mp.log(q / PI) + EG + mp.log(2))), 1e-15)
    for j, dstr in enumerate(DCHI):
        d = mp.mpf(dstr)
        df = float(d)
        Dc = Dchi_closed(d, q, a, chi)
        if REF[name][j]:
            quoted('%s: D_chi(%s) closed' % (name, dstr), Dc, REF[name][j])
        tot = mp.mpf(0)
        for zs, S, argh, redone in (sides if not real else sides * 2):
            Tl = mp.mpf(zs[-1])
            Sl = len(zs) - theta_chi(Tl, q, a) / PI + argh / PI
            tot += mp.mpf(math.fsum(h_np(zs, df).tolist())) + tail_chi(d, Tl, q, a) - h_mp(Tl, d) * Sl
        Dz = tot / d
        diff = Dz - Dc
        n1 = nearest_pp(1 + d)
        line = '%s, delta %s: zeros - closed' % (name, dstr)
        if n1 and mp.almosteq(1 + d, n1, 1e-20) and mp.re(chi(n1)) != 0:
            pred = -(mp.re(chi(n1)) * LAM[n1] / (PI * d)) * sum(resonance_I(0, mp.mpf(s_[0][-1])) for s_ in (sides if not real else sides * 2))
            QDIFF[(name, dstr)] = (diff, pred)
            quoted(line + ': residual after the resonance term (Re chi(2) Lambda(2)/(pi delta))(1/T+ + 1/T-)', diff - pred, QRES_QUOTED[name])
        else:
            QDIFF[(name, dstr)] = (diff, None)
            worst_nonres[name] = max(worst_nonres.get(name, 0), abs(diff))
        quoted(line, diff, QDIFF_QUOTED[name][j])
        # the alternatives: pole term and wrong parity
        if dstr in ('0.01', '3.5'):
            pole_shift = Dchi_closed(d, q, a, chi, pole=True) - Dc
            par_shift = Dchi_closed(d, q, a, chi, parity=1 - a) - Dc
            quoted('%s, delta %s: the pole term 2 ln(1+delta) would move D_chi by' % (name, dstr), pole_shift, {'0.01': '0.020', '3.5': '3.01'}[dstr])
            quoted('%s, delta %s: the wrong parity would move D_chi by (absolute value)' % (name, dstr), abs(par_shift), {'0.01': '0.0050', '3.5': '0.69'}[dstr])
    small('%s: all non-resonant rows, |zeros - closed|' % name, worst_nonres.get(name, 0), 3e-6)
print('   zeros - closed per character and delta (for the table):')
for name, *_ in CHARS:
    print('   %-50s ' % name + '  '.join('%+.1e' % float(QDIFF[(name, dd)][0]) for dd in DCHI), flush=True)
c4 = lambda n: {0: 0, 1: 1, 2: 0, 3: -1}[n % 4]
c3 = lambda n: {0: 0, 1: 1, 2: -1}[n % 3]
for dstr, q4, q3 in (('0.1', '2.346', '2.058'), ('0.01', '4.584', '4.296'), ('0.1', '2.3460172', '2.0583352'), ('0.01', '4.5838428', '4.2961607')):
    quoted('repo value D_chi4(%s)' % dstr, Dchi_closed(dstr, 4, 1, c4), q4)
    quoted('repo value D_chi3(%s)' % dstr, Dchi_closed(dstr, 3, 1, c3), q3)
    small('D_chi4 - D_chi3 = ln(4/3) at %s' % dstr, Dchi_closed(dstr, 4, 1, c4) - Dchi_closed(dstr, 3, 1, c3) - mp.log(mp.mpf(4) / 3), 1e-25)
quoted('ln(4/3)', mp.log(mp.mpf(4) / 3), '0.287682072452')
for dstr in ('0.3', '0.9', '1'):
    d = mp.mpf(dstr)
    small('even chi, delta <= 1: D_chi = D_zeta - 2 ln(1+delta) + ln q at %s (q = 5)' % dstr,
          Dchi_closed(d, 5, 0, lambda n: {0: 0, 1: 1, 2: -1, 3: -1, 4: 1}[n % 5]) - (D_closed(d) - 2 * mp.log1p(d) + mp.log(5)), 1e-25)
for dstr in ('1.3', '1.9'):
    d = mp.mpf(dstr)
    small('1 < delta <= 2: chi3 gains +((delta-1)/delta) ln 2 over its prime-free form at %s' % dstr,
          Dchi_closed(d, 3, 1, c3) - (Dchi_closed(d, 3, 1, lambda n: 0) + ((d - 1) / d) * mp.log(2)), 1e-25)
    small('1 < delta <= 2: chi4 gains nothing at %s' % dstr, Dchi_closed(d, 4, 1, c4) - Dchi_closed(d, 4, 1, lambda n: 0), 1e-25)

# ================================================================ 8. primes against zeros, end to end
sec('8. primes against zeros: R(X, delta) by exact integration of the step function')
XMAX = 10 ** 6
NS = 2 * XMAX + 10
sv = np.ones(NS + 1, bool)
sv[:2] = False
for p in range(2, int(NS ** 0.5) + 1):
    if sv[p]:
        sv[p * p::p] = False
primes = np.nonzero(sv)[0]
pp, wt = [], []
for p in primes.tolist():
    v = p
    while v <= NS:
        pp.append(v)
        wt.append(math.log(p))
        v *= p
order = np.argsort(pp)
pp = np.array(pp, float)[order]
wt = np.array(wt)[order]
cw = np.concatenate([[0.0], np.cumsum(wt)])
psi = lambda x: cw[np.searchsorted(pp, x, side='right')]
small('psi(100) = 94.0453...', psi(np.array([100.0]))[0] - 94.04531, 1e-4)
Z1 = Z[:1000]
GG = np.concatenate([Z1, -Z1])
RHO = 0.5 + 1j * GG
thm = GG[:, None] - GG[None, :]
WIN = {'0.1': ('1.10865', '1.11008', '1.10271', '0.69', '1.71', '+0.007'),
       '0.5': ('0.20541', '0.20561', '0.20470', '0.13', '0.29', '+0.001'),
       '1': ('0.18819', '0.18835', '0.18760', '0.12', '0.27', '+0.001')}
pt_diffs = []
for dstr, (q_p, q_z, q_D, q_lo, q_hi, q_res) in WIN.items():
    d = float(dstr)
    Bk = np.unique(np.concatenate([[1.0, XMAX], pp[pp <= XMAX], pp[pp / (1 + d) >= 1] / (1 + d)]))
    Bk = Bk[(Bk >= 1) & (Bk <= XMAX)]
    lo, hi = Bk[:-1], Bk[1:]
    mid = (lo + hi) / 2
    K = psi((1 + d) * mid) - psi(mid)
    seg = ((d * hi - K) ** 3 - (d * lo - K) ** 3) / (3 * d)          # int of (d x - K)^2 on each piece
    cum = np.concatenate([[0.0], np.cumsum(seg)])

    def I_at(X):
        j = np.searchsorted(Bk, X, side='right') - 1
        Kj = psi((1 + d) * (Bk[j] + X) / 2) - psi((Bk[j] + X) / 2)
        return cum[j] + ((d * X - Kj) ** 3 - (d * Bk[j] - Kj) ** 3) / (3 * d)

    R = lambda X: I_at(X) / (0.5 * d * X * X)
    c = ((1 + d) ** RHO - 1) / RHO
    Dd = float(D_closed(dstr))
    diag1000 = float(np.sum(np.abs(c) ** 2)) / d
    M = c[:, None] * np.conj(c)[None, :]
    for X in (1e4, 3e5, 1e6):
        lx = math.log(X)
        Rz = (2 / (d * X * X)) * np.real(np.sum(M * (np.exp((2 + 1j * thm) * lx) - 1) / (2 + 1j * thm))) + (Dd - diag1000) * (1 - 1 / X ** 2)
        pt_diffs.append(abs(R(np.array(X)) - Rz))
    U0, U1 = math.log(1e2), math.log(1e6)
    Lw = U1 - U0
    us = U0 + (np.arange(20000) + 0.5) * Lw / 20000
    Rp = R(np.exp(us))
    with np.errstate(divide='ignore', invalid='ignore'):
        ker = np.where(thm == 0, 1.0 + 0j, (np.exp(1j * thm * U1) - np.exp(1j * thm * U0)) / (1j * thm * Lw))
    kexp = (math.exp(-2 * U0) - math.exp(-2 * U1)) / (2 * Lw)
    Rz_avg = (2 / d) * np.real(np.sum(M * (ker - kexp) / (2 + 1j * thm))) + (Dd - diag1000) * (1 - kexp)
    quoted('delta %s: window mean of R from the primes, ln X in [ln 1e2, ln 1e6]' % dstr, Rp.mean(), q_p)
    quoted('delta %s: same window from 1000 zeros + diagonal tail' % dstr, Rz_avg, q_z)
    quoted('delta %s: D(delta)' % dstr, Dd, q_D)
    quoted('delta %s: smallest R in the window' % dstr, Rp.min(), q_lo)
    quoted('delta %s: largest R in the window' % dstr, Rp.max(), q_hi)
    quoted('delta %s: off-diagonal window residue (zeros window mean - D)' % dstr, Rz_avg - Dd, q_res)
    del M, ker
quoted('pointwise |R primes - R zeros| at X = 1e4, 3e5, 1e6: smallest', min(pt_diffs), '1.3e-4')
quoted('... largest', max(pt_diffs), '2.3e-3')

# ---------------------------------------------------------------- summary
print('', flush=True)
print('%d checks, %d failed  (%.0f s)' % (NCHK[0], NCHK[1], time.time() - T_START), flush=True)
sys.exit(1 if NCHK[1] else 0)
