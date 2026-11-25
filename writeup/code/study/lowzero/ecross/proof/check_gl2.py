"""check_gl2.py: independent checks of every identity in gl2_onepoint_proof.md.

Own code (does not import gl2lib): theta by mpmath, derivatives of theta by the Stirling series, I_n by direct
exponentials on a uniform grid (Simpson), coefficients by PARI.  Asserts:
  A1  X(1/2+it)^{-1} = e^{2i theta(t)};  L(1/2+it) = w e^{-2i theta(t)} L(1/2-it) (PARI lfun, centre 1);  Z = w conj Z
  A2  Stirling forms of theta, theta', theta'' (errors O(t^-3), O(t^-4), O(t^-5))
  A3  stationary phase for one n: quadrature vs leading term vs second order; the simplified constant
  A4  Prop. 1 (exact explicit formula), Gaussian-edged box, two curves:  |zeros - primes| < 1e-6
  A5  Voronoi identity sum lambda(n) n^-1/2 I_n = w sum lambda(m) m^-1/2 Hhat(-ln m):  < 1e-6
  A6  Theorem 2 (Hann [20,400]), 5 curves: |Z - P| < 0.3 (leading), < 0.1 (second order)
  A7  twist r = 2, 1/2 (Hann [20,400]), 2 curves: |Z - P| < 0.3
  A8  the line factors from the Voronoi heuristic equal the recipe's rho(p), rho(1/p)   (sympy, exact)
  A9  constants: partial summation (2 sqrt2/3), prime-side density <-> zero-side law (quadrature)
  A10 sharp Voronoi partial sums: ratio to (2sqrt2/3) e^{-3pi i/4} w (x^3/N)^{1/4} within 3 % at x = 4000 N
Run: (ulimit -v 1500000; timeout 600 nice python -I check_gl2.py)"""
import os, math, json, cmath, time, sys
PARTS = sys.argv[1] if len(sys.argv) > 1 else 'all'   # 'a' = A1-A5, A8-A10; 'b' = A6-A7; 'all'
import numpy as np, mpmath as mp, cypari2, sympy as sp
from scipy.special import erf

t_start = time.time()
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', '..', '..'))
ZF = np.load(REPO + '/writeup/data/study/r144/fresh_zeros.npz')
CV = {c['id']: c for c in json.load(open(REPO + '/writeup/code/study/lowzero/deg2/fresh_set.json'))['curves']}
pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)
PI = math.pi
ok = []
def check(name, cond, msg):
    print(('PASS ' if cond else 'FAIL ') + name + ': ' + msg); ok.append(bool(cond))

# ------------------------------------------------------------------ theta (mpmath) and Stirling-series derivatives
def theta_mp(N, t):
    return float(t / 2 * mp.log(N) - t * mp.log(2 * mp.pi) + mp.im(mp.loggamma(1 + 1j * t)))
theta_v = np.vectorize(theta_mp)
def d2theta(t):   # 2 theta''
    return 2 / t - 1 / (3 * t ** 3) - 1 / (15 * t ** 5) - 1 / (21 * t ** 7)
def d3theta(t):   # 2 theta'''
    return -2 / t ** 2 + 1 / t ** 4 + 1 / (3 * t ** 6)
def d4theta(t):
    return 4 / t ** 3 - 4 / t ** 5 - 2 / t ** 7
def Lfreq(N, t):  # 2 theta'
    return math.log(N) - 2 * math.log(2 * PI) + 2 * np.log(t) + 1 / (6 * t ** 2) + 1 / (60 * t ** 4) + 1 / (126 * t ** 6)

if PARTS in ('a', 'all'):
    # ------------------------------------------------------------------ A1
    E = CV['N307_0_0_1_m8_m9']; N = E['N']; w = E['w']
    pari('E1 = ellinit(%s); L1 = lfuncreate(E1);' % str(E['ainv']))
    errs = []
    for t in (5.3, 47.1, 213.7):
        s = mp.mpf(0.5) + 1j * mp.mpf(t)
        Xinv = mp.power(N, s - 0.5) * mp.power(2 * mp.pi, 1 - 2 * s) * mp.gamma(s + 0.5) / mp.gamma(1.5 - s)
        errs.append(abs(complex(Xinv) - cmath.exp(2j * theta_mp(N, t))))
    fe = []
    for t in (3.7, 11.2, 25.0):
        Lp = complex(pari('lfun(L1, 1 + I*%r)' % t)); Lm = complex(pari('lfun(L1, 1 - I*%r)' % t))
        fe.append(abs(Lp - w * cmath.exp(-2j * theta_mp(N, t)) * Lm) / abs(Lp))
        Zt = cmath.exp(1j * theta_mp(N, t)) * Lp
        fe.append(abs(Zt - w * Zt.conjugate()) / abs(Zt))
    check('A1', max(errs) < 1e-10 and max(fe) < 1e-8, 'X^-1 vs e^{2i theta}: %.1e; FE with root number and Z = w conj(Z): %.1e' % (max(errs), max(fe)))

    # ------------------------------------------------------------------ A2
    r0, r1, r2 = [], [], []
    for t in (20.0, 50.0, 100.0, 400.0):
        th2 = 2 * theta_mp(N, t)
        r0.append(abs(th2 - (t * math.log(N * t * t / (4 * PI * PI)) - 2 * t + PI / 2 - 1 / (6 * t))) * t ** 3)
        d1 = float(mp.log(N) - 2 * mp.log(2 * mp.pi) + 2 * mp.re(mp.digamma(1 + 1j * t)))      # 2 theta'
        r1.append(abs(d1 - math.log(N * t * t / (4 * PI * PI)) - 1 / (6 * t * t)) * t ** 4)
        d2 = float(-2 * mp.im(mp.psi(1, 1 + 1j * t)))                                           # 2 theta''
        r2.append(abs(d2 - 2 / t + 1 / (3 * t ** 3)) * t ** 5)
    check('A2', max(r0) < 0.2 and max(r1) < 0.2 and max(r2) < 1.0,
          'max t^3|2th - stirling| = %.3f, t^4|2th\' - ...| = %.3f, t^5|2th\'\' - 2/t + 1/(3t^3)| = %.3f' % (max(r0), max(r1), max(r2)))

# ------------------------------------------------------------------ weights
def hann(t, a, b, d=0):
    L = b - a; k = PI / L; inside = (t > a) & (t < b)
    if d == 0: return np.where(inside, np.sin(k * (t - a)) ** 2, 0.0)
    if d == 1: return np.where(inside, k * np.sin(2 * k * (t - a)), 0.0)
    return np.where(inside, 2 * k * k * np.cos(2 * k * (t - a)), 0.0)
def gbox(t, a, b, D):
    return 0.5 * (erf((t - a) / D) - erf((t - b) / D))

# ------------------------------------------------------------------ stationary phase per n (own implementation)
def tstat(N, lnm):
    t = 2 * PI * np.sqrt(np.exp(lnm) / N)
    for _ in range(5):
        t = t - (Lfreq(N, t) - lnm) / d2theta(t)
    return t
def sp_terms(N, n, coef, r=1.0, lo=0, hi=1e9):
    """stationary-phase pieces per n: t_n, base_n = -(1/2pi) coef n^-1/2 sqrt(2pi/f'') e^{i(phi(t_n)+pi/4)}, f'', f''', f''''."""
    lnm = np.log(n * r); t0 = tstat(N, lnm)
    m = (t0 > lo) & (t0 < hi); n, coef, t0, lnm = n[m], coef[m], t0[m], lnm[m]
    f2, f3, f4 = d2theta(t0), d3theta(t0), d4theta(t0)
    ph = 2 * theta_v(N, t0) - t0 * lnm
    base = -(1 / (2 * PI)) * coef * n ** -0.5 * np.sqrt(2 * PI / f2) * np.exp(1j * (ph + PI / 4))
    return t0, base, f2, f3, f4
def apply_weight(terms, Hf, H1=None, H2=None, second=False):
    t0, base, f2, f3, f4 = terms
    g = Hf(t0); amp = g
    if second:
        g1, g2 = H1(t0), H2(t0)
        amp = g + (1j / (2 * f2)) * (g2 - g1 * f3 / f2 - g * f4 / (4 * f2) + 5 * g * f3 ** 2 / (12 * f2 ** 2))
    return np.sum(base * amp)

def In_quad(N, Hf, lo, hi, lnn, h=0.02):
    """Simpson rule, direct exponentials (no recursion)."""
    t = np.arange(lo, hi + h / 2, h)
    if len(t) % 2 == 0: t = np.append(t, t[-1] + h)
    wS = np.ones(len(t)); wS[1:-1:2] = 4; wS[2:-1:2] = 2; wS *= h / 3
    F = Hf(t) * np.exp(2j * theta_v(N, t)) * wS
    out = np.empty(len(lnn), complex)
    for s in range(0, len(lnn), 800):
        out[s:s + 800] = np.exp(-1j * np.outer(lnn[s:s + 800], t)) @ F
    return out

def prime_powers(ainv, bad, X):
    pari('Ec = ellinit(%s); Pc = primes([2, %d]);' % (str(ainv), X))
    ap = np.array(pari('Vecsmall(apply(p -> ellap(Ec, p), Pc))'), float)
    pr = np.array(pari('Vecsmall(Pc)'), float)
    lam = ap / np.sqrt(pr); isb = np.isin(pr.astype(np.int64), bad)
    th = np.arccos(np.clip(lam / 2, -1, 1))
    ns, Ls = [], []; k = 1
    while True:
        m = pr ** k <= X
        if not m.any(): break
        ns.append(pr[m] ** k); Ls.append(np.where(isb[m], lam[m] ** k, 2 * np.cos(k * th[m])) * np.log(pr[m])); k += 1
    n = np.concatenate(ns); L = np.concatenate(Ls); o = np.argsort(n)
    return n[o], L[o]

if PARTS in ('a', 'all'):
    # ------------------------------------------------------------------ A3: one n, Hann [20,400], resonant n with t_n ~ 250
    n0 = round(N * 250.0 ** 2 / (4 * PI ** 2)) | 1
    Hf = lambda t: hann(t, 20, 400)
    Iq = In_quad(N, Hf, 20, 400, np.array([math.log(n0)]), h=0.005)[0]
    t0 = tstat(N, np.array([math.log(n0)]))[0]
    f2, f3, f4 = d2theta(t0), d3theta(t0), d4theta(t0)
    lead = math.sqrt(2 * PI / f2) * cmath.exp(1j * (2 * theta_mp(N, t0) - t0 * math.log(n0) + PI / 4)) * Hf(np.array([t0]))[0]
    g, g1, g2 = hann(np.array([t0]), 20, 400)[0], hann(np.array([t0]), 20, 400, 1)[0], hann(np.array([t0]), 20, 400, 2)[0]
    sec = lead / g * (g + (1j / (2 * f2)) * (g2 - g1 * f3 / f2 - g * f4 / (4 * f2) + 5 * g * f3 ** 2 / (12 * f2 ** 2)))
    simple = -(cmath.exp(3j * PI / 4) / math.sqrt(2)) * N ** -0.25 * n0 ** -0.25 * cmath.exp(-4j * PI * math.sqrt(n0 / N)) * Hf(np.array([2 * PI * math.sqrt(n0 / N)]))[0]
    exact_term = -(1 / (2 * PI)) * n0 ** -0.5 * Iq
    check('A3', abs(Iq - lead) / abs(Iq) < 3 / t0 and abs(Iq - sec) / abs(Iq) < abs(Iq - lead) / abs(Iq) and abs(simple - exact_term) / abs(exact_term) < 3 / t0,
          'n=%d t_n=%.2f: |I-lead|/|I| = %.2e, |I-2nd|/|I| = %.2e, simplified constant vs exact term: %.2e (1/t = %.2e)' % (
              n0, t0, abs(Iq - lead) / abs(Iq), abs(Iq - sec) / abs(Iq), abs(simple - exact_term) / abs(exact_term), 1 / t0))

    # ------------------------------------------------------------------ A4, A5: exact identities (Gaussian-edged box)
    for cid, a, b, D in (('N307_0_0_1_m8_m9', 60.0, 100.0, 4.0), ('N611_0_0_1_m1_1', 150.0, 200.0, 5.0)):
        C = CV[cid]; N = C['N']; w = C['w']
        lo, hi = a - 6 * D, b + 6 * D
        Hf = lambda t: gbox(t, a, b, D)
        gam = np.sort(ZF[cid]); Zs = np.sum(Hf(gam) * np.exp(2j * theta_v(N, gam)))
        xlo = N * lo ** 2 / (4 * PI ** 2) * math.exp(-1.6); xhi = N * hi ** 2 / (4 * PI ** 2) * math.exp(1.6)
        n, L = prime_powers(C['ainv'], C['bad'], int(xhi)); m = n >= xlo
        P = -(1 / (2 * PI)) * np.sum(L[m] * n[m] ** -0.5 * In_quad(N, Hf, lo, hi, np.log(n[m]), h=0.04))
        check('A4', abs(Zs - P) < 1e-6, '%s box[%g,%g] Delta=%g: zeros %.8f%+.8fi, primes %.8f%+.8fi, |diff| %.1e (%d prime powers)' % (
            cid, a, b, D, Zs.real, Zs.imag, P.real, P.imag, abs(Zs - P), m.sum()))
        if cid != 'N307_0_0_1_m8_m9':
            continue                                      # Voronoi quadrature only on the small case (cost)
        a_n = np.array(pari('Vecsmall(ellan(Ec, %d))' % int(xhi)), float); nn = np.arange(1, len(a_n) + 1, dtype=float)
        lam = a_n / np.sqrt(nn); mm = nn >= xlo
        left = np.sum(lam[mm] * nn[mm] ** -0.5 * In_quad(N, Hf, lo, hi, np.log(nn[mm]), h=0.04))
        k = np.arange(1, 60, dtype=float); u = np.log(k)
        Hh = np.where(k == 1, b - a, (np.exp(1j * b * u) - np.exp(1j * a * u)) / (1j * np.where(k == 1, 1, u)) * np.exp(-D * D * u * u / 4))
        right = w * np.sum(lam[:59] * k ** -0.5 * Hh)
        check('A5', abs(left - right) < 1e-6, 'Voronoi identity: %.8f%+.8fi vs w*(dual) %.8f%+.8fi (m = 1 term %g), |diff| %.1e' % (
            left.real, left.imag, right.real, right.imag, w * (b - a), abs(left - right)))

if PARTS in ('b', 'all'):
    # ------------------------------------------------------------------ A6, A7: Theorem 2 with the Hann taper, twists
    pred = np.load(REPO + '/writeup/code/study/lowzero/ecross/pred_fresh.npz')
    cmap = dict(zip([str(s) for s in pred['ids']], pred['c']))
    rows = []
    for cid in ('N307_0_0_1_m8_m9', 'N335_0_0_1_m2_2', 'N373_0_1_1_m2_m2', 'N611_0_0_1_m1_1', 'N870_1_0_0_m5_m3'):
        C = CV[cid]; N = C['N']; w = C['w']
        rmax = 2.0 if cid in ('N307_0_0_1_m8_m9', 'N611_0_0_1_m1_1') else 1.0
        n, L = prime_powers(C['ainv'], C['bad'], int(rmax * N * 402.0 ** 2 / (4 * PI ** 2)))
        gam = np.sort(ZF[cid])
        H0 = lambda t: hann(t, 20, 400); H1 = lambda t: hann(t, 20, 400, 1); H2 = lambda t: hann(t, 20, 400, 2)
        Zs = np.sum(H0(gam) * np.exp(2j * theta_v(N, gam)))
        terms = sp_terms(N, n, L, lo=19.5, hi=400.5)
        P1 = apply_weight(terms, H0); P2 = apply_weight(terms, H0, H1, H2, second=True)
        law = -w * cmap[cid] * 190 / (2 * PI)
        check('A6', abs(Zs - P1) < 0.3 and abs(Zs - P2) < 0.1, '%s Hann[20,400]: zeros %.3f%+.3fi, primes %.3f%+.3fi (|d| %.3f), 2nd order |d| %.3f; law -w c intH/2pi = %.2f' % (
            cid, Zs.real, Zs.imag, P1.real, P1.imag, abs(Zs - P1), abs(Zs - P2), law))
        rows.append((cid, Zs, P1, P2, law))
        if rmax == 2.0:
            for r in (2.0, 0.5):
                Zr = np.sum(H0(gam) * np.exp(2j * theta_v(N, gam) - 1j * gam * math.log(r)))
                Pr = apply_weight(sp_terms(N, n, L, r=r, lo=19.5, hi=400.5), H0)
                check('A7', abs(Zr - Pr) < 0.3, '%s twist r=%g: zeros %.3f%+.3fi, primes %.3f%+.3fi, |d| %.3f' % (cid, r, Zr.real, Zr.imag, Pr.real, Pr.imag, abs(Zr - Pr)))

if PARTS in ('a', 'all'):
    # ------------------------------------------------------------------ A8: line factors (exact algebra)
    lp, X = sp.symbols('lambda_p X')
    l2, l3 = lp ** 2 - 1, lp ** 3 - 2 * lp                          # lambda(p^2), lambda(p^3) at a good prime
    mu = {0: 1, 1: -lp, 2: 1}
    num_p = sum(mu[j] * [1, lp, l2, l3][j + 1] * X ** j for j in range(3))      # sum_j mu(p^j) lambda(p^{j+1}) X^j
    rec_p = (1 - X) * (lp - l3 * X)                                               # recipe numerator of rho(p)
    num_1p = sum([mu[j + 1] * [1, lp, l2][j] * X ** j for j in range(2)])         # sum_j mu(p^{j+1}) lambda(p^j) X^j
    rec_1p = -lp * (1 - X)
    Fp = sum(mu[j] * [1, lp, l2][j] * X ** j for j in range(3))
    check('A8', sp.expand(num_p - rec_p) == 0 and sp.expand(num_1p - rec_1p) == 0 and sp.expand(Fp - (1 - X) * (1 - l2 * X)) == 0,
          'sum mu(p^j)lambda(p^{j+1})X^j = (1-X)(lambda_p - lambda(p^3)X); sum mu(p^{j+1})lambda(p^j)X^j = -lambda_p(1-X); F_p = (1-X)(1-lambda(p^2)X)')

    # ------------------------------------------------------------------ A9: constants
    x_ = sp.symbols('x', positive=True); n_ = sp.symbols('n', positive=True)
    dens = sp.Rational(1, 1) / sp.sqrt(2) * n_ ** sp.Rational(-1, 4)           # |mean density| of Lambda_E e(..) per unit c w N^{-1/4}
    S_int = sp.integrate(dens, (n_, 0, x_))
    c1 = sp.simplify(S_int / x_ ** sp.Rational(3, 4) - 2 * sp.sqrt(2) / 3) == 0
    # prime-side smoothed law -> zero side: -(e^{3pi i/4}/sqrt2) N^{-1/4} int n^{-1/4} mu(n) H(t_n) dn = -(w c/2pi) int H dt
    Nn, wc = 500.0, 1.7
    tt = np.linspace(20, 400, 200001); nn = Nn * tt ** 2 / (4 * PI ** 2); dn = Nn * tt / (2 * PI ** 2)
    mu_t = (wc / math.sqrt(2)) * cmath.exp(-3j * PI / 4) * (nn * Nn) ** -0.25          # mean density of Lambda_E(n) e(-2 sqrt(n/N))
    val = -(cmath.exp(3j * PI / 4) / math.sqrt(2)) * Nn ** -0.25 * np.trapezoid(nn ** -0.25 * mu_t * hann(tt, 20, 400) * dn, tt)
    target = -wc / (2 * PI) * 190
    check('A9', c1 and abs(val - target) < 1e-6 * abs(target), 'int_0^x n^{-1/4}/sqrt2 dn = (2sqrt2/3) x^{3/4}: %s; density -> zero-side law: %.6f%+.6fi vs %.6f' % (c1, val.real, val.imag, target))

    # ------------------------------------------------------------------ A10: sharp Voronoi partial sums
    for cid in ('N307_0_0_1_m8_m9', 'N870_1_0_0_m5_m3'):
        C = CV[cid]; N = C['N']; w = C['w']; Xs = 4000 * N
        pari('Ev = ellinit(%s);' % str(C['ainv']))
        a_n = np.array(pari('Vecsmall(ellan(Ev, %d))' % Xs), float); nn = np.arange(1, Xs + 1, dtype=float)
        S = np.sum(a_n / np.sqrt(nn) * np.exp(-4j * PI * np.sqrt(nn / N)))
        main = (2 * math.sqrt(2) / 3) * cmath.exp(-3j * PI / 4) * w * (Xs ** 3 / N) ** 0.25
        check('A10', abs(S / main - 1) < 0.03, '%s x=%d: sum lambda(n) e(-2sqrt(n/N)) = %.1f%+.1fi, main %.1f%+.1fi, |S/main - 1| = %.4f, |S-main|/sqrt(x) = %.3f' % (
            cid, Xs, S.real, S.imag, main.real, main.imag, abs(S / main - 1), abs(S - main) / math.sqrt(Xs)))
        del a_n, nn

print('%d/%d checks passed (%.0fs)' % (sum(ok), len(ok), time.time() - t_start))
assert all(ok)
